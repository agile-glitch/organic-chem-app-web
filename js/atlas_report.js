/* AtlasReport — the Validation tab's card for the Fragment Atlas: the held-out figures stored with the tables (data/atlas/meta.js), so the number on screen is the number that was measured. */
(() => {
  'use strict';
  const host = document.getElementById('validAtlas');
  if (!host || !window.Atlas) return;
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const L = [['M0', 'the group alone'], ['M1', '+ the kind of reagent'], ['M2', '+ the group\'s surroundings (radius 1, 2, 3)'], ['M3', '+ the named reagents'], ['M4', '+ steric, charge, ring-angle and pKa descriptors']];
  const table = (res, extra) => `<table class="mv-tbl" style="margin:6px 0"><thead><tr><th>model</th><th>log-loss (lower is better)</th><th>AUC</th><th>AUC within one group type</th></tr></thead><tbody>` +
    L.map(([k, t]) => `<tr><td>${esc(t)}</td><td>${res[k].logloss.toFixed(4)}</td><td>${res[k].auc.toFixed(3)}</td><td>${res[k].within_fragment_auc.toFixed(3)}</td></tr>`).join('') + '</tbody></table>';
  window.Atlas.ready().then(M => {
    if (!M || !M.eval) { host.innerHTML = '<p class="hint">The Fragment Atlas tables (data/atlas) are not installed.</p>'; return; }
    const A = M.eval, B = M.evalB, sc = A.site_choice;
    let h = `<p class="hint" style="max-width:900px">Every functional group is read together with its surroundings (radius 1, 2, 3 around its reacting atom) and what is known about it comes from ${M.instances.toLocaleString()} group instances in
      479,035 atom-mapped patent reactions (${M.entries.toLocaleString()} contexts kept). Figures are on reactions held out by molecule: a held-out molecule never occurs in any training reaction
      (${A.train.toLocaleString()} training, ${A.validation.toLocaleString()} validation and ${A.test.toLocaleString()} test instances).</p>`;
    h += '<h3 style="margin:10px 0 2px">Is the group transformed by the reaction?</h3>' + table(A.results);
    if (sc) h += `<p class="hint" style="max-width:900px">Which of several groups of the same type in one molecule reacts (${sc.groups.toLocaleString()} molecules, exactly one transformed): chance ${Math.round(sc.random_top1 * 100)} %, the group and reagent alone ${Math.round(sc.M1 * 100)} %, with the surroundings <b>${Math.round(sc.M3 * 100)} %</b>.</p>`;
    if (B && B.results) h += '<h3 style="margin:10px 0 2px">The reaction happens elsewhere in the molecule: is this group also transformed? (the interference question)</h3>' + table(B.results) +
      `<p class="hint" style="max-width:900px">${B.train.toLocaleString()} training, ${B.validation.toLocaleString()} validation, ${B.test.toLocaleString()} test instances; about 1 % are also transformed. Patent records show only what was tried and worked, so "no recorded case" is shown as <b>no precedent</b>, never as "spared".</p>`;
    host.innerHTML = h;
  });
})();
