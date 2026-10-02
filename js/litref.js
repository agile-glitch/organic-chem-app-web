/* LitRef — "Mentioned together in the literature": for a reaction's reactants and product, the papers in PubChem's
   literature that mention both molecules, each with its title, authors, journal, year, DOI and PubMed link.

   LitRef.button(label, getPairs, panel) → a <button> that fills `panel` when clicked
       getPairs() → [{a: {smiles, name}, b: {smiles, name}}]   usually each reactant with the product
   LitRef.panel() → an empty, hidden panel element to put somewhere on the page

   What it can and cannot say: PubChem links a paper to every compound the paper mentions, so a paper listed here
   mentions both molecules. It may describe this reaction, the reverse one, or something else entirely (both as
   solvents, as metabolites, in an analysis), and the panel says so. The reaction-level evidence is the app's own
   recorded reactions (patents, ORD), shown elsewhere on the page. Looked up only on a click, through PubChem3D's
   request queue (one request at a time, spaced; see js/pubchem3d.js); each answer is kept in this browser and shown
   from there next time, also offline, with the date it was fetched. */
window.LitRef = (() => {
  'use strict';
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const NOTE = 'PubChem links each paper to the compounds it mentions. A paper listed here mentions both molecules: it may ' +
    'describe this reaction, the reverse, or something else (both as solvents, metabolites or analytes). Recorded ' +
    'reactions (patents, the Open Reaction Database) are the reaction-level evidence.';

  function panel() { const p = el('div', 'card litref'); p.hidden = true; return p; }
  function button(label, getPairs, box) {
    const b = el('button', 'litref-btn', esc(label));
    b.title = 'Papers that mention both molecules, from PubChem. Sends the structures (as SMILES) to PubChem (NCBI) over the internet; answers are kept in this browser.';
    b.addEventListener('click', () => show(box, getPairs(), b));
    return b;
  }

  const authors = s => { const a = String(s || '').split(/,\s*/).filter(Boolean); return a.length > 3 ? a.slice(0, 3).join(', ') + ' et al.' : a.join(', '); };
  function article(x) {
    const year = (x.date || '').slice(0, 4);
    const links = [x.doi ? `<a href="https://doi.org/${esc(encodeURI(x.doi))}" target="_blank" rel="noopener">DOI ${esc(x.doi)}</a>` : '',
      x.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(x.pmid)}/" target="_blank" rel="noopener">PubMed ${esc(x.pmid)}</a>` : ''].filter(Boolean).join(' · ');
    return `<li><b>${esc(x.title)}</b>${x.review ? ' <span class="litref-tag">review</span>' : ''}<div class="meta">${esc(authors(x.authors))}` +
      `${x.journal ? ' · <i>' + esc(x.journal) + '</i>' : ''}${year ? ' (' + esc(year) + ')' : ''}</div>${links ? `<div class="meta">${links}</div>` : ''}</li>`;
  }

  async function show(box, pairs, btn) {
    const P = window.PubChem3D;
    box.hidden = false;
    box.innerHTML = '<h3>Mentioned together in the literature <small>PubChem</small></h3><p class="hint">Asking PubChem…</p>';
    box.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    if (btn) btn.disabled = true;
    try {
      if (!P || !P.identify || !P.papersBoth) throw new Error('the PubChem module is not loaded');
      if (!window.RDKitLoad) throw new Error('RDKit is not available');
      const R = await window.RDKitLoad();
      const parts = [];
      for (const { a, b } of pairs) {
        const ia = await P.identify(R, a.smiles), ib = await P.identify(R, b.smiles);
        const head = `<h4>${esc(a.name || a.smiles)} <span class="amp">and</span> ${esc(b.name || b.smiles)}</h4>`;
        if (!ia.ok || !ib.ok) { parts.push(head + `<p class="hint">${esc((!ia.ok ? (a.name || a.smiles) + ': ' + ia.error : '') + (!ia.ok && !ib.ok ? ' ' : '') + (!ib.ok ? (b.name || b.smiles) + ': ' + ib.error : ''))}</p>`); continue; }
        const r = await P.papersBoth(ia.cid, ib.cid);
        if (!r.ok) { parts.push(head + `<p class="hint">${esc(r.error)}</p>`); continue; }
        const all = `https://pubchem.ncbi.nlm.nih.gov/compound/${ia.cid}#section=Chemical-Co-Occurrences-in-Literature`;
        const when = `${r.saved ? 'From your saved record' : 'Fetched'} ${esc(r.fetched)}`;
        parts.push(head + (r.count
          ? `<p><b>${r.count.toLocaleString()}</b> paper${r.count === 1 ? '' : 's'} in PubChem mention both. The ${Math.min(5, r.articles.length)} PubChem ranks most relevant:</p>` +
            `<ol class="litref-list">${r.articles.map(article).join('')}</ol>`
          : '<p>PubChem lists no paper mentioning both. Its lists keep only each compound\'s most frequent co-mentions, so a rare pairing can be missing: this is not proof that no paper exists.</p>') +
          `<p class="hint">${when} · CID ${ia.cid} (${esc(ia.title)}) and CID ${ib.cid} (${esc(ib.title)}) · ` +
          `<a href="${all}" target="_blank" rel="noopener">all co-occurrences on PubChem</a></p>`);
      }
      box.innerHTML = '<h3>Mentioned together in the literature <small>PubChem</small></h3>' +
        `<p class="hint">${esc(NOTE)}</p>` + parts.join('') +
        '<p class="hint">Literature links: NCBI/NLM, from PubMed and publishers\' data in PubChem · ' +
        '<a href="https://www.ncbi.nlm.nih.gov/home/about/policies/" target="_blank" rel="noopener">NCBI policies and disclaimer</a>.</p>';
    } catch (e) {
      box.innerHTML = '<h3>Mentioned together in the literature <small>PubChem</small></h3>' + `<p class="hint">The lookup failed: ${esc(e && e.message || e)}</p>`;
    } finally { if (btn) btn.disabled = false; }
  }

  return { button, panel, show };
})();
