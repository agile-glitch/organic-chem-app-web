/* RCRoute — whole routes built from single steps, for the RetroChimera tab (js/retrochimera.js supplies the steps).

   A best-first search over partial routes, in the spirit of Retro* (Chen et al., ICML 2020) and AiZynthFinder:
     state     the steps chosen so far and the molecules still to be made ("open"); a state with nothing open is a route
     expand    the first open molecule gets its candidate steps (RetroChimera's suggestions merged with the app's taught
               reactions; the caller labels each with its evidence and says which reactants can be bought)
     order     cost so far + H × (molecules still open): the cheapest-looking partial route is continued first
     cost      per step: ½·(−ln p) for RetroChimera's probability p (capped), plus a penalty by evidence:
                 taught 0 · recorded 0 · confirmed by the forward check 0.3 · RetroChimera alone 2.0,
               and 0.3 more when only one of its two sub-models proposed the step
     stop      a molecule stops being open when it can be bought; a branch stops at the depth limit; the search stops
               at the route count, the expansion budget, the time limit, or when the caller aborts
     counted   only routes with no 'model' step count as routes; a finished route with one goes to `unconfirmed`
               (at most as many again) and the search carries on looking for supported ones
   The costs are a ranking device, set by hand, not a measured probability that a route works: the page says so.
   Molecules to expand are sent in small batches (one GPU call each). A molecule already expanded is not asked again.
   Nothing here touches the page: RCRoute.search(target, deps, opts, onProgress, signal) → Promise<result>, where
     deps.expandMany(smilesList) → Promise<[options[] per molecule]>, option = {reactants: [{smi, key, buyable}], p,
       agree, tier: 'taught'|'recorded'|'confirmed'|'model', …anything else the caller wants back}
     deps.key(smiles) → a canonical key (for cycles and repeats)
   result = {routes: [{steps: [{product: {smi, key, depth}, option}], open: [], cost}], unconfirmed: [same],
             partial: {…} | null (only when neither list has anything),
             stats: {expansions, molecules, ms, stopped}} */
window.RCRoute = (() => {
  'use strict';
  const DEFAULTS = { maxDepth: 4, maxSteps: 8, maxExpansions: 40, maxMs: 240000, branch: 6, routes: 3, batch: 4, H: 1.5, perLastStep: 2 };
  const TIER = { taught: 0, recorded: 0, confirmed: 0.3, model: 2.0 };
  function stepCost(o) {
    const base = o.p == null ? 0.3 : Math.min(3.5, -0.5 * Math.log(Math.max(o.p, 1e-3)));
    return base + (TIER[o.tier] ?? 2) + (o.p != null && o.agree === false && o.tier !== 'taught' ? 0.3 : 0);
  }
  const sig = s => s.steps.map(x => x.product.key + '<' + x.option.reactants.map(r => r.key).sort().join('+')).sort().join('|');

  async function search(target, deps, opts = {}, onProgress = () => {}, signal = null) {
    const O = Object.assign({}, DEFAULTS, opts);
    const t0 = Date.now(), cache = new Map();
    let expansions = 0, stopped = 'searched everything within the limits';
    const f = s => s.cost + O.H * s.open.length;
    const tk = deps.key(target);
    let queue = [{ open: [{ smi: target, key: tk, depth: 0, anc: [] }], steps: [], cost: 0 }];
    const routes = [], unconfirmed = [], seen = new Set(), lastStep = new Map();
    let partial = null;
    const keepPartial = s => {
      if (!s.steps.length) return;
      if (!partial || s.open.length < partial.open.length || (s.open.length === partial.open.length && s.cost < partial.cost)) partial = s;
    };
    const budgetLeft = () => expansions < O.maxExpansions && Date.now() - t0 < O.maxMs;
    const progress = () => onProgress({ expansions, molecules: cache.size, routes: routes.length, unconfirmed: unconfirmed.length, ms: Date.now() - t0, queue: queue.length });

    while (queue.length && routes.length < O.routes) {
      if (signal && signal.aborted) { stopped = 'stopped by you'; break; }
      queue.sort((a, b) => f(a) - f(b));
      /* expand, in one batch, the first open molecule of the best few states that has not been expanded yet */
      if (budgetLeft()) {
        const need = [];
        for (const s of queue.slice(0, 16)) {
          const m = s.open[0];
          if (!m || m.depth >= O.maxDepth || s.steps.length >= O.maxSteps || cache.has(m.key) || need.some(x => x.key === m.key)) continue;
          need.push(m);
          if (need.length >= Math.min(O.batch, O.maxExpansions - expansions)) break;
        }
        if (need.length) {
          let res;
          try { res = await deps.expandMany(need.map(m => m.smi)); } catch (e) { stopped = 'RetroChimera stopped answering: ' + (e.message || e); break; }
          need.forEach((m, i) => cache.set(m.key, (res[i] || []).map(o => Object.assign(o, { cost: stepCost(o) })).sort((a, b) => a.cost - b.cost)));
          expansions += need.length;
          progress();
          if (signal && signal.aborted) { stopped = 'stopped by you'; break; }
        }
      }
      const s = queue.shift();
      if (!s.open.length) {                                              // a complete route
        const k = sig(s), last = s.steps[0] && s.steps[0].option.reactants.map(r => r.key).sort().join('+');
        if (seen.has(k) || (lastStep.get(last) || 0) >= O.perLastStep) continue;   // no repeats; some variety in the final step
        seen.add(k);
        if (s.steps.some(x => x.option.tier === 'model')) { if (unconfirmed.length < O.routes) unconfirmed.push(s); continue; }
        lastStep.set(last, (lastStep.get(last) || 0) + 1);
        routes.push(s); progress();
        continue;
      }
      const m = s.open[0];
      if (m.depth >= O.maxDepth || s.steps.length >= O.maxSteps || !cache.has(m.key)) { keepPartial(s); continue; }
      for (const o of cache.get(m.key).slice(0, O.branch)) {
        if (o.reactants.some(r => r.key === m.key || m.anc.includes(r.key))) continue;      // no cycles
        const open = s.open.slice(1);
        for (const r of o.reactants) {
          if (r.buyable || open.some(x => x.key === r.key)) continue;
          open.push({ smi: r.smi, key: r.key, depth: m.depth + 1, anc: m.anc.concat(m.key) });
        }
        queue.push({ open, steps: s.steps.concat({ product: m, option: o }), cost: s.cost + o.cost });
      }
      if (queue.length > 4000) { queue.sort((a, b) => f(a) - f(b)); queue.length = 2000; }
      if (!budgetLeft() && !queue.some(q => !q.open.length || cache.has(q.open[0].key))) break;
    }
    if (routes.length >= O.routes) stopped = `found ${O.routes} routes`;
    else if (stopped === 'searched everything within the limits') {
      if (expansions >= O.maxExpansions) stopped = `used its budget of ${O.maxExpansions} molecules`;
      else if (Date.now() - t0 >= O.maxMs) stopped = `reached the time limit (${Math.round(O.maxMs / 60000)} min)`;
    }
    for (const q of queue) keepPartial(q);
    routes.sort((a, b) => a.cost - b.cost); unconfirmed.sort((a, b) => a.cost - b.cost);
    return { routes, unconfirmed, partial: routes.length || unconfirmed.length ? null : partial, stats: { expansions, molecules: cache.size, ms: Date.now() - t0, stopped } };
  }
  return { search, stepCost, DEFAULTS, TIER };
})();
