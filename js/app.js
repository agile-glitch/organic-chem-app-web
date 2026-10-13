/* Organic Chemistry — app entry point.
   Data available as globals (see README):
     window.COURSE     13 chapters: {chapter, title, concepts[], reactions[], tables[], rules[]}
     window.REACTIONS  14,489 reactions: {id, reactants[], product, agents[], solvents[], tempC, yield, source, class, procedure}
     window.NAMES      {solvents: {smiles: name}, agents: {smiles: name}}

   app.html runs on its own (a tab bar across the top, one page shown at a time) or inside a tab of the home screen
   (index.html, js/shell.js) as app.html?embed=<page>#<page>. Embedded, it shows just that page; the home screen's tab
   bar covers this page's own (hidden, same height, so nothing below moves). A hand-off to another app (Sketcher's
   "Send to Reactions", Possibilities' "Plan it", Retrosynthesis' "Open in Reactions", …) is passed up to the home
   screen, which opens that app in a tab of its own and makes the same call there. */
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

  const names = [...tabs].map(t => t.dataset.page).filter(Boolean);
  const embed = new URLSearchParams(location.search).get('embed');
  const home = names.includes(embed) && window.parent !== window ? embed : null;

  /* open the tab named in the URL (#reactions), default to the sketcher */
  const initial = home || location.hash.slice(1);
  showPage(names.includes(initial) ? initial : 'sketcher');
  if (!home) return;

  /* ---------- inside a home-screen tab ---------- */
  const post = msg => window.parent.postMessage(Object.assign({ agilechem: true }, msg), '*');
  // the calls one app makes into another, the page that owns each, what the page's own button would also fill in, and
  // (for a call that wants a tab of its own) the new tab's name from the call's arguments
  const HANDOFFS = [
    ['Reactions', 'addStart', 'reactions'], ['Reactions', 'set', 'reactions'], ['Reactions', 'load', 'reactions'],
    ['Retro', 'plan', 'retro', 'retroTarget'],
    ['MolView', 'show', 'mol', 'molInput'],
    ['Draw', 'showMechanism', 'draw', null, args => ({ newTab: true, title: (args[0] && args[0].title) || 'Mechanism' })],
  ];
  const original = {};
  HANDOFFS.forEach(([api, fn, page, input, tabFor]) => {
    const obj = window[api];
    if (!obj || typeof obj[fn] !== 'function') return;
    const own = obj[fn].bind(obj);
    original[api + '.' + fn] = (...args) => {
      const box = input && document.getElementById(input);
      if (box && typeof args[0] === 'string') box.value = args[0];
      return own(...args);
    };
    if (page === home) return;
    obj[fn] = (...args) => {
      showPage(home);              // the caller may already have switched this frame to the other page
      try { post(Object.assign({ type: 'open', app: page, api, fn, args }, tabFor ? tabFor(args) : {})); } catch (e) { return Promise.resolve(false); }
      return Promise.resolve(true);
    };
  });
  window.addEventListener('message', e => {
    const d = e.data;
    if (e.source !== window.parent || !d || !d.agilechem || d.type !== 'call') return;
    const f = original[d.api + '.' + d.fn];
    if (f) { showPage(home); f(...(d.args || [])); }
  });
  if (document.readyState === 'complete') post({ type: 'ready' });
  else window.addEventListener('load', () => post({ type: 'ready' }));
})();
