/* MolContext — the Molecule tab's "Groups in context": every functional group of the molecule with its surroundings (the Fragment Atlas, js/atlas.js)
   and the reagents that recorded reactions most often used to transform a group like it.

   MolContext.render(host, smiles)

   For each group: where it sits (the reacting atom and what it is bonded to), how crowded it is, its Gasteiger partial charge, its ring and the angle
   strain the ring forces, its typical pKa; and, on a click, the named reagents (Reagents sheet) that most often transformed a group of this kind in
   the 478,000 recorded reactions the Atlas learned from, with how many records stand behind each rate. */
window.MolContext = (() => {
  'use strict';
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const crowd = s => s <= 2 ? 'open' : s === 3 ? 'moderately hindered' : s === 4 ? 'hindered' : 'very hindered';
  const pct = x => (x >= 0.995 ? '>99' : x < 0.005 ? '<1' : String(Math.round(x * 100))) + ' %';
  const short = n => String(n || '').split(' (')[0];
  let seq = 0;

  async function render(host, smiles) {
    const my = ++seq, AT = window.Atlas;
    host.innerHTML = '';
    if (!AT || !smiles) return;
    const R = await window.RDKitLoad();
    const meta = await AT.ready();
    if (!meta || my !== seq) return;
    const RR = window.REACTION_RULES || {};
    const insts = AT.describe(R, smiles);
    // one row per atom: overlapping group definitions (ester / methyl-ethyl ester) are the same place
    const byAtom = new Map();
    insts.forEach(i => { const k = i.atom; if (!byAtom.has(k)) byAtom.set(k, Object.assign({}, i, { names: [short(i.name)], fids: [i.fid] })); else { const c = byAtom.get(k); if (!c.names.includes(short(i.name))) c.names.push(short(i.name)); c.fids.push(i.fid); } });
    const rows = [...byAtom.values()].filter(i => !i.fids.every(f => f === 'F009' || f === 'F108')).sort((a, b) => a.atom - b.atom);
    if (!rows.length) { host.innerHTML = '<p class="mv-hint">No functional group of the Fragments sheet is in this molecule.</p>'; return; }
    let html = '<p class="mv-hint">Each group with its surroundings. Click a group for the reagents that recorded reactions most often used to transform a group like it.</p><table class="mv-tbl"><thead><tr><th>group</th><th>where</th><th>crowding</th><th>charge</th><th>ring</th><th>pKa</th></tr></thead><tbody>';
    rows.forEach((i, k) => {
      const f = (RR.fragments || {})[i.fids[0]] || {};
      const pk = [f.pKa != null ? `${f.pKa}` : '', f.pKaH != null ? `(H⁺: ${f.pKaH})` : ''].filter(Boolean).join(' ');
      html += `<tr class="mv-ctx-row" data-k="${k}" style="cursor:pointer"><td>${esc(i.names.join(' / '))}</td><td><code>${esc(i.where)}</code></td><td>${esc(crowd(i.steric))} (${i.steric})</td>` +
        `<td>${i.charge == null ? '–' : (i.charge >= 0 ? '+' : '') + i.charge.toFixed(2)}</td><td>${i.ring ? i.ring + '-ring' + (i.strain >= 5 ? `, strain ${i.strain.toFixed(0)}°` : '') : 'chain'}</td><td>${esc(pk || '–')}</td></tr>` +
        `<tr class="mv-ctx-more" data-for="${k}" hidden><td colspan="6"></td></tr>`;
    });
    host.innerHTML = html + '</tbody></table>';
    host.querySelectorAll('.mv-ctx-row').forEach(tr => tr.addEventListener('click', async () => {
      const k = +tr.dataset.k, more = host.querySelector(`.mv-ctx-more[data-for="${k}"]`), cell = more.firstChild;
      if (!more.hidden) { more.hidden = true; return; }
      more.hidden = false; cell.innerHTML = '<span class="mv-hint">Looking up recorded reagents…</span>';
      let prof = [];
      try { prof = await AT.reagentProfile(R, rows[k], { top: 8 }); } catch (e) { prof = []; }
      cell.innerHTML = prof.length
        ? `<div class="mv-hint">Reagents most often recorded transforming a group like this (rate over the recorded conditions; records have no temperature or time):</div>` +
          prof.map(p => `<div><b>${esc(p.name)}</b> <span class="mv-hint">${esc(p.role)} · ${pct(p.rateCtx != null ? p.rateCtx : p.rate)} of ${(p.rateCtx != null ? p.nCtx : p.n).toLocaleString()} recorded${p.rateCtx != null ? ' with this surrounding' : ' (the group, any surrounding)'}</span></div>`).join('')
        : '<span class="mv-hint">Not enough recorded reactions of a named reagent with this group.</span>';
    }));
  }
  return { render };
})();
