/* The home screen (index.html). Each app icon opens the app in a new tab: app.html?embed=<page>#<page> in a frame of
   its own, named Draw1, Draw2, Molecule1 … (the lowest free number for that app). Tabs are independent copies of the
   app, so two drawings or two molecules can be open at once. A hand-off from one app to another (js/app.js) arrives
   here as a message and goes to the most recently used tab of the receiving app, or to a new one; the call waits
   until that tab has loaded. index.html#reactions (the old address of a page) opens that app on arrival. */
(() => {
  'use strict';
  const APPS = { draw: 'Draw', sketcher: 'Sketcher', mol: 'Molecule', retro: 'Retrosynthesis', reactions: 'Reactions',
    simulate: 'Simulate', possib: 'Possibilities', valid: 'Validation' };
  const $ = id => document.getElementById(id);
  const strip = $('strip'), frames = $('frames'), home = $('home'), homeTab = $('homeTab');
  const open = [];             // {app, n, name, tab, frame, ready, queue[], used}
  let current = null, clock = 0;

  const CLOSE_ICON = '<svg viewBox="0 0 10 10" aria-hidden="true"><path d="M1.5 1.5l7 7M8.5 1.5l-7 7"/></svg>';

  function openApp(app) {
    if (!APPS[app]) return null;
    let n = 1; while (open.some(t => t.app === app && t.n === n)) n++;
    const t = { app, n, name: APPS[app] + n, ready: false, queue: [], used: 0 };
    t.tab = document.createElement('div');
    t.tab.className = 'tab';
    t.tab.innerHTML = `<button class="label">${t.name}</button><button class="close" title="Close ${t.name}" aria-label="Close ${t.name}">${CLOSE_ICON}</button>`;
    t.tab.querySelector('.label').addEventListener('click', () => show(t));
    t.tab.querySelector('.close').addEventListener('click', () => close(t));
    t.tab.addEventListener('auxclick', e => { if (e.button === 1) { e.preventDefault(); close(t); } });
    // keep the tabs in the order they were opened, the new one at the end
    strip.appendChild(t.tab);
    t.frame = document.createElement('iframe');
    t.frame.title = t.name;
    t.frame.src = `app.html?embed=${app}#${app}`;
    frames.appendChild(t.frame);
    open.push(t);
    show(t);
    t.tab.scrollIntoView({ block: 'nearest', inline: 'nearest' });
    return t;
  }

  function show(t) {
    current = t;
    if (t) t.used = ++clock;
    homeTab.classList.toggle('active', !t);
    home.hidden = !!t;
    open.forEach(o => { o.tab.classList.toggle('active', o === t); o.frame.classList.toggle('active', o === t); });
    document.title = t ? t.name + ' — AgileCHEM' : 'AgileCHEM';
    if (t) { try { t.frame.contentWindow.focus(); } catch (e) { /* not focusable yet */ } }
  }

  function close(t) {
    const i = open.indexOf(t);
    if (i < 0) return;
    open.splice(i, 1);
    t.tab.remove(); t.frame.remove();
    if (current === t) {
      // the most recently used tab left open, else home
      const next = open.slice().sort((a, b) => b.used - a.used)[0] || null;
      show(next);
    }
  }

  function deliver(t) {
    if (!t.ready) return;
    while (t.queue.length) t.frame.contentWindow.postMessage(Object.assign({ agilechem: true, type: 'call' }, t.queue.shift()), '*');
  }

  window.addEventListener('message', e => {
    const d = e.data;
    if (!d || !d.agilechem) return;
    const from = open.find(t => t.frame.contentWindow === e.source);
    if (!from) return;
    if (d.type === 'ready') { from.ready = true; deliver(from); }
    else if (d.type === 'open' && APPS[d.app]) {
      const to = open.filter(t => t.app === d.app).sort((a, b) => b.used - a.used)[0] || openApp(d.app);
      show(to);
      to.queue.push({ api: d.api, fn: d.fn, args: d.args });
      deliver(to);
    }
  });

  homeTab.addEventListener('click', () => show(null));
  document.querySelectorAll('[data-app]').forEach(b => b.addEventListener('click', () => openApp(b.dataset.app)));

  const want = location.hash.slice(1);
  if (APPS[want]) {
    try { history.replaceState(null, '', location.pathname + location.search); } catch (e) { /* the hash stays */ }
    openApp(want);
  }
  window.Shell = { open: openApp, show, close, get tabs() { return open.map(t => t.name); } };
})();
