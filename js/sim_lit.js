/* Simulate → Literature: a published experiment played out on the 3D bench.

   The experiments are in data/lab_experiments.js and every quantity, time, temperature and yield in them is the
   published source's (each names its source and links it). This file only plays them: it walks the steps, sets what
   is in each vessel, turns heaters on and off, raises and lowers baths, pours from one vessel into another, moves
   the camera to what the step is about, and shows the words. Nothing is predicted here.

   The words keep their origins apart, as the data does:
     the step text, quantities, times and temperatures   the source (paraphrased)
     "The source reports"                                 an observation or result the source states
     "Worked out"                                         arithmetic from the source's numbers
     "Why"                                                this app's explanation, not the source's
     "Drawn, not in the source"                           what the picture shows that the source does not say
   and every experiment carries a yield badge: measured (a checked, published yield), expected (stated as an
   expectation), example (one teaching run), or none (the source gives no yield).

   The bench is js/lab3d.js (three.js, loaded on first use like the Molecule tab's 3D view). Click a liquid to see what
   is in it; hover anything for its name and amounts. */
(() => {
  'use strict';
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const THREE_SRC = 'vendor/three/three.js';
  const STEP_MS = 3600;                                   // how long Play dwells on a step
  const SPEEDS = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 2];            // playback speeds on offer
  const S = { speed: 1, view: null, loading: null, exp: null, at: -1, playing: false, timer: 0, pinned: null, hoverEvt: null };
  const HEATERS = new Set(['burner', 'hotplate', 'stirplate', 'mantle', 'oilbath', 'waterbath', 'sandbath', 'steambath', 'alblock', 'overheadstirrer', 'microwave', 'meltemp', 'uvlamp']);
  const GAUGES = new Set(['thermometer', 'meltemp']);            // pieces with a red column that step.reading sets
  const YIELD_BADGE = {
    measured: ['measured', 'A yield the source reports from its own run of the procedure.'],
    expected: ['expected', 'The source states this as what to expect, not as a measured result.'],
    example: ['one example run', 'A single teaching run, not a published, checked yield.'],
    none: ['no yield given', 'The source does not state a yield.'],
  };

  /* ---------------- amounts ---------------- */
  const num = v => { if (typeof v === 'number') return v; const m = String(v || '').match(/[\d.]+/g); return m ? m.map(Number).reduce((a, b) => a + b, 0) / m.length : null; };
  const amountText = a => [a.g != null ? a.g + ' g' : '', a.ml != null ? a.ml + ' mL' : '', a.mol != null ? a.mol + ' mol' : ''].filter(Boolean).join(' · ');
  /* the volume to draw: the stated mL, or a liquid's grams through its density; solids take no volume here */
  function drawnMl(exp, c) {
    if (!c) return 0;
    if (c.ml != null) return num(c.ml) || 0;
    let v = 0;
    for (const it of c.items || []) {
      const s = exp.substances[it.key] || {};
      if (it.ml != null) v += num(it.ml) || 0;
      else if (it.g != null && s.density) v += (num(it.g) || 0) / s.density;
    }
    return v;
  }
  function forView(exp, c) {                              // the data's contents → what Lab3D draws
    if (!c) return null;
    const out = Object.assign({}, c);
    if (c.layers) out.layers = c.layers.map(L => ({ ml: L.ml, color: L.color, cloudy: L.cloudy }));
    else { const v = drawnMl(exp, c); if (v > 0) out.ml = v; else if (!c.solid && !(c.level > 0)) out.level = 0; }
    return out;
  }

  /* ---------------- the page ---------------- */
  function build(host) {
    host.innerHTML = `
      <div class="lit">
        <div class="toolbar lit-top">
          <label class="lit-pick"><span>Experiment</span> <select id="litPick"></select></label>
          <button id="litPlay" class="primary" title="Play the experiment step by step">Play</button>
          <button id="litPrev" title="Previous step">◀</button>
          <button id="litNext" title="Next step">▶</button>
          <button id="litReset" title="Back to the empty bench">Reset</button>
          <label class="lit-speed" title="Playback speed: how long Play stays on each step, and how fast pours and camera moves run"><span>Speed</span> <select id="litSpeed">${SPEEDS.map(v => `<option value="${v}">${v}×</option>`).join('')}</select></label>
          <button id="litWhole" title="See the whole bench">Whole bench</button>
          <span class="lit-yield-badge" id="litBadge"></span>
        </div>
        <div class="lit-main">
          <div class="lit-stage">
            <div id="litScene" class="lit-scene"></div>
            <div id="litTip" class="lit-tip" hidden></div>
            <div class="lit-float" role="group" aria-label="Playback">
              <button id="litFPrev" title="Previous step (←)">◀</button>
              <button id="litFPlay" class="primary" title="Play or pause">Play</button>
              <button id="litFNext" title="Next step (→)">▶</button>
              <span id="litFStep" class="lit-float-step"></span>
              <select id="litFSpeed" title="Playback speed">${SPEEDS.map(v => `<option value="${v}">${v}×</option>`).join('')}</select>
            </div>
            <div class="lit-help">Drag to turn · Shift-drag to move · scroll to zoom · double-click a piece to go to it · click a liquid to see what is in it</div>
          </div>
          <aside class="lit-side">
            <div id="litStep" class="lit-step"></div>
            <div id="litPanel" class="lit-panel"></div>
          </aside>
        </div>
        <div id="litSteps" class="lit-strip" role="list"></div>
        <div id="litSource" class="lit-source"></div>
      </div>`;
    const pick = document.getElementById('litPick');
    const list = window.LAB_EXPERIMENTS || [];
    list.forEach((x, i) => { const o = document.createElement('option'); o.value = String(i); o.textContent = `${i + 1}. ${x.title}`; pick.appendChild(o); });
    pick.addEventListener('change', () => load(+pick.value));
    document.getElementById('litPlay').addEventListener('click', togglePlay);
    document.getElementById('litPrev').addEventListener('click', () => { stop(); goto(S.at - 1, false); });
    document.getElementById('litNext').addEventListener('click', () => { stop(); goto(S.at + 1, true); });
    document.getElementById('litReset').addEventListener('click', () => { stop(); goto(-1, false); });
    /* the same controls float on the bench itself, so they stay at hand while zoomed in on the reaction */
    document.getElementById('litFPlay').addEventListener('click', togglePlay);
    document.getElementById('litFPrev').addEventListener('click', () => { stop(); goto(S.at - 1, false); });
    document.getElementById('litFNext').addEventListener('click', () => { stop(); goto(S.at + 1, true); });
    const speedBoxes = [document.getElementById('litSpeed'), document.getElementById('litFSpeed')];
    try { const v = +localStorage.getItem('lit-speed'); if (SPEEDS.includes(v)) S.speed = v; } catch (e) {}
    speedBoxes.forEach(b => { b.value = String(S.speed); });
    speedBoxes.forEach(speedBox => speedBox.addEventListener('change', () => {
      S.speed = +speedBox.value || 1;
      speedBoxes.forEach(b => { b.value = String(S.speed); });
      try { localStorage.setItem('lit-speed', String(S.speed)); } catch (e) {}
      if (S.view && S.view.setSpeed) S.view.setSpeed(S.speed);
      if (S.playing) { clearTimeout(S.timer); S.timer = setTimeout(S.tick, STEP_MS / S.speed); }   // takes effect at once
    }));
    document.getElementById('litWhole').addEventListener('click', () => { if (S.view) { S.view.focus(null); S.view.highlight(null); } });
    host.addEventListener('keydown', e => {
      if (e.target.tagName === 'SELECT') return;
      if (e.key === 'ArrowRight') { stop(); goto(S.at + 1, true); e.preventDefault(); }
      if (e.key === 'ArrowLeft') { stop(); goto(S.at - 1, false); e.preventDefault(); }
    });
    if (!list.length) { document.getElementById('litStep').innerHTML = '<div class="banner">No experiments are loaded (data/lab_experiments.js).</div>'; return; }
    let first = 0;
    try { const saved = localStorage.getItem('lit-exp'); const k = list.findIndex(x => x.id === saved); if (k >= 0) first = k; } catch (e) {}
    load(first);
  }

  async function ensureView() {
    if (S.view) return S.view;
    if (!S.loading) {
      const box = document.getElementById('litScene');
      box.innerHTML = '<p class="hint lit-loading">Setting up the bench (three.js)…</p>';
      S.loading = (window.THREE ? Promise.resolve() : window.RDKitLoad.script(THREE_SRC)).then(() => {
        if (!window.THREE) throw new Error('three.js did not load');
        box.innerHTML = '';
        S.view = window.Lab3D.create(box, { onPick, onHover });
        S.view.setSpeed(S.speed);
        return S.view;
      }).catch(e => { S.loading = null; box.innerHTML = `<div class="banner">The 3D bench could not start: ${esc(e.message || e)}</div>`; throw e; });
    }
    return S.loading;
  }

  async function load(i) {
    const list = window.LAB_EXPERIMENTS || [];
    stop();
    S.exp = list[i] || list[0];
    S.at = -1; S.pinned = null;
    document.getElementById('litPick').value = String(list.indexOf(S.exp));
    try { localStorage.setItem('lit-exp', S.exp.id); } catch (e) {}
    const yb = YIELD_BADGE[S.exp.yield.kind] || YIELD_BADGE.none;
    const badge = document.getElementById('litBadge');
    badge.className = 'lit-yield-badge ' + S.exp.yield.kind; badge.textContent = 'Yield: ' + yb[0]; badge.title = yb[1] + ' ' + S.exp.yield.text;
    renderSource(); renderStrip(); renderStep(null, -1);
    document.getElementById('litPanel').innerHTML = '';
    try { await ensureView(); } catch (e) { return; }
    S.view.setBench(S.exp.bench);
    goto(-1, false);
  }

  /* the bench after step n (−1: before anything): every step's contents, heat, visibility and readings, in order —
     so stepping back gives exactly what stepping forward did */
  function stateAt(n) {
    const st = { contents: {}, heat: {}, visible: {}, reading: {}, flame: {} };
    for (let k = 0; k <= n && k < S.exp.steps.length; k++) {
      const s = S.exp.steps[k];
      Object.assign(st.contents, s.contents || {}); Object.assign(st.heat, s.heat || {});
      Object.assign(st.visible, s.visible || {}); Object.assign(st.reading, s.reading || {});
      Object.assign(st.flame, s.flame || {});
    }
    return st;
  }
  /* how each added reagent reaches the flask: drops from a pipette, a solid from a weighing boat, a measured liquid
     from a graduated cylinder. Amounts the source leaves open ("as needed") are delivered by their kind too. */
  const TOOL_MS = { pour: 2200, drops: 1700, solid: 1900 };
  function deliveries(step) {
    if (!step || step.pour || !step.vessel || !step.adds || !S.view || !S.view.canDispense(step.vessel)) return [];
    return step.adds.slice(0, 3).map(a => {
      const s = subst(a.key), note = String(a.note || '') + ' ' + String(a.ml || '');
      if (/drop/i.test(note)) return { how: 'drops', key: a.key };
      if (a.ml != null) return { how: 'pour', ml: num(a.ml) || 5, key: a.key };
      if (a.g != null || /solid|scoop|spatula|crystals|teaspoon/i.test(note) || (!s.density && !s.smiles && /leaves|cloves|charcoal|mixture/i.test(s.name))) return { how: 'solid', key: a.key };
      if (s.density) return { how: 'pour', ml: 5, key: a.key };
      return null;
    }).filter(Boolean);
  }
  function goto(n, animate) {
    if (!S.exp || !S.view) return;
    n = Math.max(-1, Math.min(S.exp.steps.length - 1, n));
    const forward = animate && n === S.at + 1;
    S.at = n;
    const st = stateAt(n), step = n >= 0 ? S.exp.steps[n] : null;
    const pour = forward && step && step.pour;
    const plan = forward ? deliveries(step) : [];
    let planMs = 0;
    for (const d of plan) {                                      // one after another, then the vessel fills
      const at = planMs, ms = TOOL_MS[d.how] / S.speed;
      setTimeout(() => { if (S.at === n && S.view) S.view.dispense(step.vessel, d.how, { ml: d.ml, ms, color: d.how === 'solid' ? 0xf4f2ea : 0xdce8ee }); }, at);
      planMs += ms + 120 / S.speed;
    }
    S.planMs = planMs;
    for (const b of S.exp.bench) {
      S.view.setVisible(b.id, b.id in st.visible ? st.visible[b.id] : !b.hidden);
      if (b.kind === 'hose') continue;
      if (HEATERS.has(b.kind)) S.view.setHeat(b.id, st.heat[b.id] || 0);
      if (b.kind === 'burner') S.view.setFlame(b.id, st.flame[b.id] == null ? null : st.flame[b.id]);
      if (GAUGES.has(b.kind)) S.view.setReading(b.id, b.id in st.reading ? st.reading[b.id] : 20);
      const c = st.contents[b.id];
      /* on a pour, the source empties while it is tipped, and the target fills a moment later */
      const k = 1 / S.speed;                              // every animation runs at the chosen speed
      const delay = pour && b.id === pour.to ? 1300 * k : planMs && step && b.id === step.vessel ? Math.max(0, planMs - 900 * k) : 0, ms = forward ? (pour && (b.id === pour.from || b.id === pour.to) ? 1100 : 900) * k : 0;
      if (delay) setTimeout(() => { if (S.at === n) S.view.setContents(b.id, forView(S.exp, c), ms); }, delay);
      else S.view.setContents(b.id, forView(S.exp, c), ms);
    }
    if (pour) S.view.pour(pour.from, pour.to, 2600 / S.speed);
    S.view.highlight(step && step.vessel ? step.vessel : null);
    if (step && (step.frame || step.vessel)) S.view.focus(step.frame || step.vessel, { above: plan.length > 0 && !step.frame }); else if (n === -1) S.view.focus(null);
    renderStep(step, n); renderStrip();
    if (S.pinned) showPiece(S.pinned);
  }
  function togglePlay() {
    if (S.playing) { stop(); return; }
    if (S.at >= S.exp.steps.length - 1) goto(-1, false);
    S.playing = true;
    document.getElementById('litPlay').textContent = 'Pause'; document.getElementById('litFPlay').textContent = 'Pause';
    const tick = S.tick = () => {
      if (!S.playing) return;
      if (S.at >= S.exp.steps.length - 1) { stop(); return; }
      goto(S.at + 1, true);
      S.timer = setTimeout(tick, (STEP_MS + (S.exp.steps[S.at] && S.exp.steps[S.at].pour ? 1200 : 0)) / S.speed + (S.planMs || 0));
    };
    S.timer = setTimeout(tick, 300);
  }
  function stop() {
    S.playing = false; clearTimeout(S.timer);
    for (const id of ['litPlay', 'litFPlay']) { const b = document.getElementById(id); if (b) b.textContent = 'Play'; }
  }

  /* ---------------- words ---------------- */
  const subst = key => (S.exp.substances || {})[key] || { name: key };
  function substTip(s, a) {
    return [s.name, s.cas ? 'CAS ' + s.cas : '', s.mw ? 'M = ' + s.mw + ' g/mol' : '', s.density ? 'density ' + s.density + ' g/mL' : '', s.role || '',
      a && amountText(a) ? 'the source: ' + amountText(a) : '', a && a.note ? a.note : '', s.note || ''].filter(Boolean).join(' · ');
  }
  function amountHTML(a) {
    const s = subst(a.key), q = amountText(a);
    return `<span class="lit-amt" title="${esc(substTip(s, a))}">${esc(s.name)}${q ? ` <b>${esc(q)}</b>` : ''}${a.note ? ` <i>${esc(a.note)}</i>` : ''}</span>`;
  }
  const block = (cls, label, text) => text ? `<div class="lit-b ${cls}"><span class="lit-b-l">${label}</span> ${esc(text)}</div>` : '';
  function renderStep(st, n) {
    const box = document.getElementById('litStep'), x = S.exp;
    if (!st) {
      box.innerHTML = `<h3>${esc(x.title)}</h3><p>${esc(x.summary || '')}</p>` +
        `<p class="lit-tech">${(x.technique || []).map(t => `<span class="lit-chip">${esc(t)}</span>`).join('')}</p>` +
        `<p class="hint">Press <b>Play</b> to run it, or <b>▶</b> (or the → key) to go one step at a time. Every step's amounts are the source's.</p>`;
      return;
    }
    box.innerHTML = `<div class="lit-no">Step ${n + 1} of ${x.steps.length}</div><h3>${esc(st.title)}</h3><p class="lit-text">${esc(st.text)}</p>` +
      (st.adds && st.adds.length ? `<div class="lit-adds"><span class="lit-b-l">Added</span> ${st.adds.map(amountHTML).join(', ')}</div>` : '') +
      (st.time || st.temp ? `<div class="lit-cond">${[st.temp, st.time].filter(Boolean).map(esc).join(' · ')}</div>` : '') +
      block('src', 'The source reports', st.see) +
      block('src yield', 'Reported', st.yieldNote) +
      block('der', 'Worked out', st.derived) +
      block('why', 'Why (the app\'s explanation)', st.why) +
      block('cue', 'Drawn, not in the source', st.cue) +
      block('cue', 'Thermometer', st.readingNote) +
      (st.safety ? `<div class="lit-safety">${esc(st.safety)}</div>` : '');
  }
  function renderStrip() {
    const f = document.getElementById('litFStep');
    if (f) f.textContent = S.at < 0 ? `${S.exp.steps.length} steps` : `Step ${S.at + 1} / ${S.exp.steps.length}`;
    const strip = document.getElementById('litSteps');
    strip.innerHTML = '';
    S.exp.steps.forEach((st, i) => {
      const b = el('button', 'lit-dot' + (i === S.at ? ' on' : '') + (i < S.at ? ' done' : ''), String(i + 1));
      b.title = st.title; b.setAttribute('role', 'listitem');
      b.addEventListener('click', () => { stop(); goto(i, i === S.at + 1); });
      strip.appendChild(b);
    });
  }
  function renderSource() {
    const x = S.exp, s = x.source, yb = YIELD_BADGE[x.yield.kind] || YIELD_BADGE.none;
    document.getElementById('litSource').innerHTML =
      `<div><b>Source:</b> ${esc(s.text)}` + (s.doi ? ` · doi:${esc(s.doi)}` : '') +
      (s.url ? ` · <a href="${esc(s.url)}" target="_blank" rel="noopener">read it as published</a>` : '') +
      (s.kind ? ` <span class="lit-kind">${esc(s.kind)}</span>` : '') + '</div>' +
      `<div><b>Yield</b> <span class="lit-yield-badge ${esc(x.yield.kind)}" title="${esc(yb[1])}">${esc(yb[0])}</span> ${esc(x.yield.text)}</div>` +
      ((x.notes || []).length ? `<details><summary>Notes from the source (${x.notes.length})</summary><ul>${x.notes.map(t => `<li>${esc(t)}</li>`).join('')}</ul></details>` : '') +
      ((x.safety || []).length ? `<details><summary>Safety</summary><ul>${x.safety.map(t => `<li>${esc(t)}</li>`).join('')}</ul></details>` : '') +
      `<p class="hint">Quantities, times, temperatures, observations and yields are the source's; the step wording is this app's own. ` +
      `Liquid heights are drawn from the stated volumes in glassware of true capacity; colours are illustrations, and glassware whose size the source does not give is sized to hold what goes in it.</p>`;
  }

  /* ---------------- clicking and hovering the bench ---------------- */
  function pieceSummary(id, part) {
    const b = (S.exp.bench || []).find(x => x.id === id);
    if (!b) return null;
    if (b.kind === 'hose') return { title: b.type === 'vacuum' ? 'vacuum hose' : b.type === 'gas' ? 'gas hose' : 'water hose (condenser cooling)', body: '' };
    const c = S.view.contentsOf(id);
    return { b, c, title: b.label || id };
  }
  function onHover(id, part, evt) {
    const tip = document.getElementById('litTip');
    if (!id) { tip.hidden = true; return; }
    const p = pieceSummary(id, part); if (!p) { tip.hidden = true; return; }
    const c = p.c, items = (c && c.items) || [];
    tip.innerHTML = `<b>${esc(p.title)}</b>` +
      (c && part === 'liquid' ? `<div>${esc(c.label || '')}</div>` + items.slice(0, 5).map(a => `<div class="lit-tip-i">${esc(subst(a.key).name)}${amountText(a) ? ' — ' + esc(amountText(a)) : ''}</div>`).join('') + (items.length > 5 ? `<div class="lit-tip-i">…</div>` : '') + '<div class="lit-tip-h">click for everything in it</div>'
        : c && c.label ? `<div>${esc(c.label)}</div><div class="lit-tip-h">click for what is in it</div>` : '');
    tip.hidden = false;
    if (evt) {
      const r = document.getElementById('litScene').getBoundingClientRect();
      const x = Math.min(r.width - 260, Math.max(8, evt.clientX - r.left + 14)), y = Math.min(r.height - 40, Math.max(8, evt.clientY - r.top + 14));
      tip.style.left = x + 'px'; tip.style.top = y + 'px';
    }
  }
  function onPick(id) {
    if (!id) { S.pinned = null; document.getElementById('litPanel').innerHTML = ''; return; }
    S.pinned = id; showPiece(id);
  }
  function showPiece(id) {
    const panel = document.getElementById('litPanel'), p = pieceSummary(id);
    if (!p) { panel.innerHTML = ''; return; }
    if (!p.b) { panel.innerHTML = `<h4>${esc(p.title)}</h4>`; return; }
    const c = p.c, items = (c && c.items) || [];
    const has = c && (c.layers || c.ml > 0 || c.level > 0 || c.solid || items.length);
    const bands = c ? [].concat(c.solid || []).filter(x => x && x.label).concat(c.layers || []) : [];   // bottom first: solids under liquids
    panel.innerHTML = `<h4>${esc(p.title)}${p.b.volume ? ` <small>${p.b.volume >= 1000 ? p.b.volume / 1000 + ' L' : p.b.volume + ' mL'}</small>` : ''}</h4>` +
      (has
        ? `<div class="lit-inside"><b>Now:</b> ${esc(c.label || '—')}</div>` +
          (bands.length ? '<div class="lit-layers">' + bands.slice().reverse().map(L => `<div><span class="lit-sw" style="background:#${(L.color || 0xcccccc).toString(16).padStart(6, '0')}"></span>${esc(L.label)}</div>`).join('') + '<div class="hint">top to bottom</div></div>' : '') +
          (items.length ? '<ul class="lit-items">' + items.map(x => `<li>${amountHTML(x)}</li>`).join('') + '</ul>' : '') +
          '<p class="hint">Amounts are the source\'s. Hover a substance for its CAS number, molar mass and role.</p>'
        : '<p class="hint">Nothing in it at this step.</p>');
  }

  /* for tests and the console: the live bench, and stepping by number */
  window.SimLit = { get view() { return S.view; }, get at() { return S.at; }, goto: (n, animate) => goto(n, !!animate), load };

  (window.SimModes = window.SimModes || []).push({
    id: 'literature',
    label: 'Literature',
    blurb: 'Published experiments played out on a 3D bench: the glassware, the amounts weighed out, what the source says you will see, and the yield it reports. Every number comes from the source named under the bench; nothing here is predicted.',
    build,
    activate() { if (S.view) S.view.resize(); },
    deactivate() { stop(); },
  });
})();
