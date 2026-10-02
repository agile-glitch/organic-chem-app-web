/* Precedent — "has this actually been made?" against the recorded-reaction sets.

   Every set has the same shape (built by tools/extract.py or tools/sources.py + build_keys.html):
     molecules [{s: smiles, f: formula, k: canonical key, p: fingerprint hex}]
     reactions [{id, r: [mol], p: mol, a: agents [mol], v: solvents [mol], t: °C, y: %, src, cls, proc, cite?, year?}]
     sources [], classes [], license?, citation?
   The key and fingerprint were produced by the same functions js/chem.js exposes (libKeyHash ∘
   canonicalKey, fpHex ∘ rxnFingerprint), verified identical, so:
     exact(g)    reactions whose product IS the molecule (same canonical key)
     from(g)     reactions that START from it (a reactant with the same key)
     similar(g)  reactions whose product is structurally close (Tanimoto on the fingerprint)
   The built-in set (window.REACTIONS) is always searched; sets loaded through js/sources.js join in.
   Every result carries `_set`, the set it came from, so cards can cite it properly.            */
(() => {
  'use strict';
  const C = window.Chem, NAMES = window.NAMES || { solvents: {}, agents: {} };
  if (!C) return;
  const SETS = [];                                            // {id, title, license, data, byProduct, byReactant, fps}
  function addSet(id, title, license, data) {
    if (!data || !data.molecules || SETS.some(s => s.id === id)) return;
    const s = { id, title, license, data, byProduct: new Map(), byReactant: new Map(), fps: null };
    const M = data.molecules;
    data.reactions.forEach((x, i) => {
      x._set = s;
      const pk = M[x.p].k; if (pk) { if (!s.byProduct.has(pk)) s.byProduct.set(pk, []); s.byProduct.get(pk).push(i); }
      x.r.forEach(j => { const k = M[j].k; if (!k) return; if (!s.byReactant.has(k)) s.byReactant.set(k, []); s.byReactant.get(k).push(i); });
    });
    SETS.push(s);
  }
  if (window.REACTIONS) addSet('db', 'the built-in database (US patents 1976–2001, AstraZeneca ELN, Open Reaction Database)', 'mixed; see README', window.REACTIONS);
  if (window.Sources) { window.Sources.loaded().forEach(id => { const set = window.Sources.SETS.find(x => x.id === id); addSet(id, set.title, set.license, window.Sources.get(id)); }); window.Sources.onLoad((id, data) => { const set = window.Sources.SETS.find(x => x.id === id); addSet(id, set.title, set.license, data); }); }

  const keyOf = g => { try { return C.libKeyHash(C.canonicalKey(g)); } catch (e) { return null; } };
  const heavyOf = (M, i) => (M[i].f.match(/[A-Z][a-z]?\d*/g) || []).reduce((n, t) => n + (t[0] === 'H' && !/^H[a-z]/.test(t) ? 0 : (parseInt(t.replace(/\D/g, ''), 10) || 1)), 0);
  const SOLVENT_LIKE = new Set([...Object.keys(NAMES.solvents), ...Object.keys(NAMES.agents), 'O', 'CO', 'CCO', 'CC(=O)O', 'CCOC(C)=O']);
  /* Data hygiene. Some source records name the wrong thing as the product — the work-up solvent
     (EtOAc, THF) or a small fragment of a much larger reactant. Those are not syntheses of that
     molecule, so they are skipped, and the caller can say how many were. */
  function plausible(x) {
    const M = x._set.data.molecules, p = M[x.p];
    if (SOLVENT_LIKE.has(p.s)) return false;
    const big = Math.max(0, ...x.r.map(i => heavyOf(M, i)));
    if (big >= 10 && heavyOf(M, x.p) < 0.4 * big) return false;
    return true;
  }
  function exact(g, all) {
    const k = keyOf(g); if (!k) return [];
    const rows = SETS.flatMap(s => (s.byProduct.get(k) || []).map(i => s.data.reactions[i]));
    return all ? rows : rows.filter(plausible);
  }
  function skipped(g) { return exact(g, true).length - exact(g).length; }
  function from(g) { const k = keyOf(g); if (!k) return []; return SETS.flatMap(s => (s.byReactant.get(k) || []).map(i => s.data.reactions[i])).filter(plausible); }
  function similar(g, limit = 6, min = 0.55) {
    let q; try { q = C.rxnFingerprint(g); } catch (e) { return []; }
    const k = keyOf(g), best = [];
    for (const s of SETS) {
      const M = s.data.molecules;
      if (!s.fps) s.fps = M.map(m => (m.p ? C.fpFromHex(m.p) : null));
      const seen = new Set();
      s.data.reactions.forEach(x => {
        const m = x.p; if (M[m].k === k || !s.fps[m] || seen.has(m) || !plausible(x)) return;
        seen.add(m);
        const sim = C.fpTanimoto(q, s.fps[m]); if (sim >= min) best.push({ sim, set: s, m });
      });
    }
    return best.sort((a, b) => b.sim - a.sim).slice(0, limit)
      .map(({ sim, set, m }) => ({ sim, rows: set.data.reactions.filter(x => x.p === m && plausible(x)).slice(0, 3) }));
  }

  /* ---- describing a row ---- */
  const molName = (x, i) => { const m = x._set.data.molecules[i]; return NAMES.agents[m.s] || NAMES.solvents[m.s] || (() => { try { return C.displayNameFor(C.parseSmiles(m.s)) || m.f; } catch (e) { return m.f || m.s; } })(); };
  const short = (x, i) => { const m = x._set.data.molecules[i]; return NAMES.agents[m.s] || NAMES.solvents[m.s] || m.f || m.s; };
  const graphOf = (x, i) => { try { return C.parseSmiles(x._set.data.molecules[i].s); } catch (e) { return null; } };
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  /* a citation link: US patents to Google Patents, DOIs to doi.org */
  function citeHtml(cite) {
    if (!cite) return '';
    const pat = cite.match(/^US0*(\d{6,8})/) || cite.match(/^(US\w+)/);
    if (pat && /^US/.test(cite)) return `<a href="https://patents.google.com/patent/US${(cite.match(/^US0*(\d+)/) || [, pat[1]])[1]}" target="_blank" rel="noopener">${esc(cite)}</a>`;
    if (/^10\.\d+\//.test(cite)) return `<a href="https://doi.org/${esc(cite)}" target="_blank" rel="noopener">doi:${esc(cite)}</a>`;
    const ord = cite.match(/^ORD (ord_dataset-([0-9a-f]{2})[0-9a-f]+)$/);
    if (ord) return `<a href="https://github.com/open-reaction-database/ord-data/blob/main/data/${ord[2]}/${ord[1]}.parquet" target="_blank" rel="noopener">Open Reaction Database, ${esc(ord[1].slice(12, 20))}</a>`;
    if (/^https?:/.test(cite)) return `<a href="${esc(cite)}" target="_blank" rel="noopener">source</a>`;
    return esc(cite);
  }

  /* one recorded reaction as a card: reactants → product, then conditions, citation and the procedure */
  function card(x) {
    const D = x._set.data;
    const el = document.createElement('div'); el.className = 'prec';
    const eqn = document.createElement('div'); eqn.className = 'prec-eqn';
    const pic = (i, w) => { const g = graphOf(x, i); const box = document.createElement('div'); box.className = 'prec-mol';
      if (g) box.appendChild(window.MolDraw.svg(g, w || 120, 70)); const nm = document.createElement('div'); nm.className = 'nm'; nm.textContent = molName(x, i); box.appendChild(nm); return box; };
    x.r.forEach((i, k) => { if (k) { const p = document.createElement('span'); p.className = 'plus'; p.textContent = '+'; eqn.appendChild(p); } eqn.appendChild(pic(i)); });
    const ar = document.createElement('span'); ar.className = 'ar'; ar.textContent = '⟶'; eqn.appendChild(ar);
    eqn.appendChild(pic(x.p, 140));
    el.appendChild(eqn);
    const cond = [];
    if (x.a.length) cond.push('<b>Reagents:</b> ' + x.a.map(i => esc(short(x, i))).join(', '));
    if (x.v.length) cond.push('<b>Solvent:</b> ' + x.v.map(i => esc(short(x, i))).join(' / '));
    if (x.t != null) cond.push('<b>T:</b> ' + x.t + ' °C');
    if (x.y != null) cond.push('<b>Yield:</b> ' + x.y + ' %');
    cond.push('<b>Source:</b> ' + esc(D.sources[x.src]) + (D.classes[x.cls] && !/unclassified/.test(D.classes[x.cls]) ? ' · ' + esc(D.classes[x.cls]) : '') +
      (x.cite ? ' · ' + citeHtml(x.cite) : '') + (x.year ? ' (' + x.year + ')' : ''));
    const c = document.createElement('div'); c.className = 'prec-cond'; c.innerHTML = cond.join(' · '); el.appendChild(c);
    if (x.proc) { const d = document.createElement('details'); d.innerHTML = '<summary>Procedure, as recorded</summary><p>' + esc(x.proc) + '</p>'; el.appendChild(d); }
    return el;
  }
  const NOTE = () => 'Searched ' + SETS.map(s => s.data.reactions.length.toLocaleString() + ' reactions in ' + s.title + (s.license ? ' [' + s.license + ']' : '')).join('; ') +
    '. Records are quoted as recorded, not re-verified. No entry here does not mean a reaction is unknown.';

  window.Precedent = { exact, from, similar, skipped, card, citeHtml, get NOTE() { return NOTE(); }, get size() { return SETS.reduce((n, s) => n + s.data.reactions.length, 0); }, molName, plausible, sets: () => SETS };
})();
