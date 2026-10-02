/* Simulate tab — two modes, each its own file, loaded before this one:

     literature  real published experiments played out on a 3D bench: the glassware someone actually used, the
                 quantities they actually weighed out, the changes you would actually see, and the yield they
                 actually got. Every experiment carries its citation, and nothing on the bench is invented:
                 quantities, temperatures, times and yields are the source's.          (js/sim_lit.js)
     simulate    things the app works out rather than reads: a reaction played out step by step with its curly
                 arrows (Chem's mechanism engine), the 3D model moving under a force field, bench arithmetic
                 (titration curves from tabulated pKa) and predicted spectra.          (js/sim_model.js)

   Each mode says where its numbers come from and how far to trust them: none of this is a quantum-chemical
   calculation. A mode registers itself as window.SimModes.push({id, label, blurb, build(host), activate(),
   deactivate()}); a mode that is not loaded is not offered, so the app still runs when one is missing. */
(() => {
  'use strict';
  const page = document.getElementById('page-simulate');
  if (!page) return;
  const MODES = window.SimModes || [];
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
  const S = { inited: false, mode: null, built: new Set() };

  function build() {
    page.innerHTML = `
      <div class="sim">
        <div class="toolbar sim-modes" id="simModes" role="tablist"></div>
        <p class="sim-blurb" id="simBlurb"></p>
        <div id="simBody"></div>
      </div>`;
    const bar = document.getElementById('simModes');
    if (!MODES.length) {
      document.getElementById('simBody').innerHTML = '<div class="banner">No simulation mode is loaded.</div>';
      return;
    }
    MODES.forEach(m => {
      const b = document.createElement('button');
      b.className = 'sim-mode'; b.dataset.mode = m.id; b.textContent = m.label;
      b.setAttribute('role', 'tab');
      b.addEventListener('click', () => pick(m.id));
      bar.appendChild(b);
      const host = document.createElement('div');
      host.className = 'sim-pane'; host.id = 'simPane-' + m.id; host.hidden = true;
      document.getElementById('simBody').appendChild(host);
    });
  }
  function pick(id) {
    const m = MODES.find(x => x.id === id) || MODES[0];
    if (!m) return;
    if (S.mode && S.mode !== m.id) {
      const old = MODES.find(x => x.id === S.mode);
      if (old && old.deactivate) { try { old.deactivate(); } catch (e) {} }
    }
    S.mode = m.id;
    document.querySelectorAll('#simModes .sim-mode').forEach(b => {
      const on = b.dataset.mode === m.id;
      b.classList.toggle('active', on); b.setAttribute('aria-selected', String(on));
    });
    document.querySelectorAll('.sim-pane').forEach(p => { p.hidden = p.id !== 'simPane-' + m.id; });
    document.getElementById('simBlurb').innerHTML = m.blurb || '';
    const host = document.getElementById('simPane-' + m.id);
    if (!S.built.has(m.id)) {
      S.built.add(m.id);
      try { m.build(host); } catch (e) { host.innerHTML = `<div class="banner">This mode could not start: ${esc(e.message || e)}</div>`; }
    }
    if (m.activate) { try { m.activate(); } catch (e) {} }
    /* the mode is remembered in this browser, not in the address: js/app.js reopens a tab only when the hash is
       exactly its name, so "#simulate/literature" would have reopened the Sketcher */
    try { localStorage.setItem('sim-mode', m.id); } catch (e) {}
  }
  function init() {
    if (!S.inited) { S.inited = true; build(); }
    let want = S.mode;
    if (!want) { try { want = localStorage.getItem('sim-mode'); } catch (e) { want = null; } }
    pick(MODES.some(m => m.id === want) ? want : (MODES[0] && MODES[0].id));
  }
  const tab = document.querySelector('.tab[data-page="simulate"]');
  if (tab) tab.addEventListener('click', init);
  if (location.hash === '#simulate') init();
  window.Simulate = { init, pick, get mode() { return S.mode; } };
})();
