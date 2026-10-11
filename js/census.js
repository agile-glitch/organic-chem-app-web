/* Validation tab: how predictions are made and what they are built from — the fragment statistics and the functional-group census
   (data/fragment_census.js, written by tools/fg_census.py). The same content is in Reaction_prediction_and_functional_groups.docx (tools/make_report_docx.py). */
(() => {
  'use strict';
  const box = document.getElementById('validCensus');
  const C = window.FRAGMENT_CENSUS;
  if (!box) return;
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  if (!C) { box.innerHTML = '<p class="hint">The census data (data/fragment_census.js) is not installed.</p>'; return; }
  const fmt = n => Number(n).toLocaleString();
  const pct = x => (100 * x).toFixed(x < 0.0005 ? 3 : x < 0.01 ? 2 : 1) + ' %';
  const nf = C.n_fragments, ni = C.n_site + C.n_adv, nt = C.n_templates;

  const tiers = '<ol style="max-width:900px">' +
    '<li><b>Hand-written workbook rules.</b> The 102 textbook reactions of reaction_rules.xlsx. If the reagent and conditions match one, it answers and shows a mechanism where one exists.</li>' +
    '<li><b>Exact look-up in the patent data.</b> If the exact molecules and reagent are in a recorded patent reaction, the recorded product is returned.</li>' +
    '<li><b>The relevance model</b> (the main learned predictor). A neural network trained on about 982,000 recorded patent reactions guesses which of about 39,000 common templates the reaction uses, from the molecules, the reagent and the fragments they contain. ' +
    'The page runs the 100 most probable templates (each also in a looser form), a ranking model scores the products, and a confidence model says how likely the first answer is to be right. ' +
    'The ranking uses the template probabilities, the fragments consumed against how readily the reagent usually takes them, and the interference features (site and neighbouring-group effects).</li>' +
    '<li><b>Older fallback.</b> The earlier reagent-keyed templates answer when the relevance model is not installed or fails.</li></ol>' +
    '<p class="hint" style="max-width:900px">An answer from the relevance model is shown only when its confidence passes a floor (50 % with a reagent that has a reactive role, 70 % with no reagent, 85 % with a reagent the sheet does not classify), ' +
    'and never when the starting materials have fewer than 8 heavy atoms. Otherwise the page says "Unknown". The measured figures are in the card above.</p>';

  const frs = C.fragment_stats.map(r => ({ id: r[0], name: r[1], mean: r[2], share: r[3], pres: r[4] })).sort((a, b) => b.mean - a.mean);
  const fragRows = list => list.map(f => `<tr><td>${esc(f.id)}</td><td>${esc(f.name)}</td><td class="num">${f.mean.toFixed(f.mean < 0.01 ? 4 : 2)}</td><td class="num">${pct(f.share)}</td><td class="num">${f.pres.toFixed(2)}</td></tr>`).join('');
  const fragTable = list => '<table><thead><tr><th>ID</th><th>Fragment</th><th class="num">Average per reaction</th><th class="num">Reactions containing it</th><th class="num">Average when present</th></tr></thead><tbody>' + fragRows(list) + '</tbody></table>';

  const cols = [['group', 'Functional group'],
    ['fragments', `Fragments: entries (of ${nf})`], ['o_fragments', 'Fragments: occurrences'],
    ['interferences', `Interferences: entries (of ${ni})`], ['o_interferences', 'Interferences: occurrences'],
    ['both', `Fragments + interferences: entries (of ${nf + ni})`], ['o_both', 'Fragments + interferences: occurrences'],
    ['tpl_reactant', `Templates, reactant side: entries (of ${fmt(nt)})`], ['o_tpl_reactant', 'Templates, reactant side: occurrences'],
    ['tpl_product', 'Templates, product side: entries'], ['o_tpl_product', 'Templates, product side: occurrences'],
    ['tpl_any', 'Templates, either side: entries'], ['o_tpl_any', 'Templates, either side: occurrences']];
  let sortKey = 'tpl_any', dir = -1, filter = '';
  const censusHTML = () => {
    const rows = C.rows.filter(r => r.group.toLowerCase().includes(filter)).sort((a, b) => sortKey === 'group' ? dir * a.group.localeCompare(b.group) : dir * (a[sortKey] - b[sortKey]));
    return '<div style="overflow-x:auto"><table style="font-size:12px"><thead><tr>' +
      cols.map(([k, t]) => `<th class="${k === 'group' ? '' : 'num'}" data-k="${k}" style="cursor:pointer" title="sort by this column">${esc(t)}${k === sortKey ? (dir < 0 ? ' ▼' : ' ▲') : ''}</th>`).join('') + '</tr></thead><tbody>' +
      rows.map(r => '<tr>' + cols.map(([k]) => k === 'group' ? `<td>${esc(r.group)}</td>` : `<td class="num">${fmt(r[k])}</td>`).join('') + '</tr>').join('') + '</tbody></table></div>' +
      (rows.length ? '' : '<p class="hint">No functional group matches.</p>');
  };

  box.innerHTML =
    '<p style="margin:0 0 10px"><a href="Reaction_prediction_and_functional_groups.docx" download class="primary" style="display:inline-block;padding:6px 14px;border-radius:6px;background:#2f6f5e;color:#fff;text-decoration:none">Download this section as a Word document</a></p>' +
    '<h3 style="margin:14px 0 4px">How a reaction is predicted</h3>' + tiers +
    '<h3 style="margin:18px 0 4px">Fragments in the patent reactions</h3>' +
    `<p class="hint" style="max-width:900px">Counting every fragment in the starting materials of all ${fmt(C.reactions)} recorded reactions (reagents not counted): a reaction has on average ${C.hits_per_reaction.toFixed(1)} fragment hits spread over ${C.distinct_per_reaction.toFixed(1)} different fragments. ` +
    'Fragments that score zero (LiAlH4, DIBAL-H, NaH, Pd, halide, azide and fluoride anions) are reagent-type groups and occur only in reagents.</p>' +
    '<div id="censusFrags">' + fragTable(frs.slice(0, 25)) + '</div>' +
    `<p><button id="censusMore">Show all ${frs.length} fragments</button></p>` +
    '<h3 style="margin:18px 0 4px">How often each functional group appears</h3>' +
    '<p class="hint" style="max-width:900px">Two numbers for each group and list. <b>Entries</b>: in how many items of the list the group appears at all (an item is one fragment, one interference row or one template). ' +
    '<b>Occurrences</b>: how many times it appears in total, counting repeats inside the same item. Example: three fragments contain an alcohol and one of them has two OH groups, so alcohol has 3 entries and 4 occurrences. ' +
    'Entries show how widely a group is covered; occurrences show how heavily it is represented; a group with many template entries but few fragment or interference entries is common in the patent chemistry but weakly described by us. ' +
    `Lists: ${nf} fragments, ${ni} interference rows (${C.n_site} site patterns + ${C.n_adv} advisories), ${fmt(nt)} templates (reactant side = what the reaction consumes, product side = what it makes). ` +
    'These count definitions, not molecules; the fragment averages above say how often a group occurs in the real patent molecules. Click a column heading to sort.</p>' +
    '<p><input id="censusFilter" type="search" placeholder="filter functional groups…" style="padding:4px 8px;width:260px"></p>' +
    '<div id="censusTable">' + censusHTML() + '</div>' +
    '<p class="hint" style="max-width:900px">How the counts were made: each fragment, site pattern and template is a SMARTS pattern; it was turned into a skeleton molecule (elements, charges, aromaticity, bond orders, stated hydrogen counts) and ' +
    `${C.rows.length} functional-group definitions were searched in it (tools/fg_census.py). A template holds only the atoms it changes and their neighbours. Conditions that exclude a group inside a site pattern were ignored; groups outside the list are not counted.</p>`;

  const redraw = () => { document.getElementById('censusTable').innerHTML = censusHTML(); };
  box.addEventListener('click', e => {
    const th = e.target.closest('#censusTable th');
    if (th) { const k = th.dataset.k; if (k === sortKey) dir = -dir; else { sortKey = k; dir = k === 'group' ? 1 : -1; } redraw(); return; }
    if (e.target.id === 'censusMore') {
      const all = e.target.dataset.all !== '1';
      document.getElementById('censusFrags').innerHTML = fragTable(all ? frs : frs.slice(0, 25));
      e.target.dataset.all = all ? '1' : '0'; e.target.textContent = all ? 'Show only the 25 most common' : `Show all ${frs.length} fragments`;
    }
  });
  box.addEventListener('input', e => { if (e.target.id === 'censusFilter') { filter = e.target.value.trim().toLowerCase(); redraw(); } });
})();
