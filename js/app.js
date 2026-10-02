/* Organic Chemistry — app entry point.
   Data available as globals (see README):
     window.COURSE     13 chapters: {chapter, title, concepts[], reactions[], tables[], rules[]}
     window.REACTIONS  14,489 reactions: {id, reactants[], product, agents[], solvents[], tempC, yield, source, class, procedure}
     window.NAMES      {solvents: {smiles: name}, agents: {smiles: name}} */
(() => {
  'use strict';

  /* ---------- tabs ---------- */
  const tabs = document.querySelectorAll('.tab');
  const pages = document.querySelectorAll('.page');

  function showPage(name) {
    tabs.forEach(t => t.classList.toggle('active', t.dataset.page === name));
    pages.forEach(p => p.classList.toggle('active', p.id === 'page-' + name));
    location.hash = name;
  }

  tabs.forEach(t => t.addEventListener('click', () => showPage(t.dataset.page)));

  /* open the tab named in the URL (#reactions), default to the sketcher */
  const names = [...tabs].map(t => t.dataset.page);
  const initial = location.hash.slice(1);
  showPage(names.includes(initial) ? initial : 'sketcher');
})();
