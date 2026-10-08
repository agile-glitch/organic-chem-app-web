/* The app's version, shown at the right end of the tab bar. Change the two values below when you publish a new version
   (e.g. 1.0.0 -> 1.1.0 for new features, 1.0.1 for fixes) and the date; nothing else needs editing. */
window.APP_VERSION = { version: '1.0.1', date: '2026-10-08' };
(() => {
  const v = window.APP_VERSION, el = document.getElementById('appVersion');
  if (!el) return;
  const d = new Date(v.date + 'T00:00:00');
  el.textContent = 'v' + v.version;
  el.title = 'Version ' + v.version + ', ' + (isNaN(d) ? v.date : d.toLocaleDateString(undefined, { year: 'numeric', month: 'long', day: 'numeric' }));
})();
