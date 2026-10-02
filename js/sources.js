/* Sources — extra public reaction sets, loaded on demand.

   Each set is a script data/sources/<id>.js built by tools/sources.py + tools/build_keys.html, of the
   same shape as data/reactions.js plus `cite` (patent number or DOI) and `year` per reaction and a
   `license` / `citation` for the set. They are large (5–10 MB), so they load only when asked for,
   by script insertion (which works from file:// as well as http). Once loaded they join the
   precedent search on every page and become selectable on the Validation tab.                */
(() => {
  'use strict';
  const SETS = [
    { id: 'uspto', file: 'data/sources/uspto.js', title: 'US patent reactions 1976–2016 (Lowe) — 20,000-reaction sample',
      license: 'CC0 1.0', note: 'Text-mined from patent grants; every reaction carries its patent number. The full set is 1.8 million reactions; this is a fixed random sample.' },
    { id: 'ord', file: 'data/sources/ord.js', title: 'Open Reaction Database — selected datasets',
      license: 'CC-BY-SA 4.0', note: 'Curated, structured reactions with conditions, yields and a DOI or patent per record. Attribution: Open Reaction Database, open-reaction-database.org.' },
  ];
  const loaded = new Map();
  const listeners = [];
  function load(id) {
    const set = SETS.find(s => s.id === id); if (!set) return Promise.reject(new Error('unknown set ' + id));
    if (loaded.has(id)) return Promise.resolve(loaded.get(id));
    return new Promise((resolve, reject) => {
      const s = document.createElement('script'); s.src = set.file;
      s.onload = () => {
        const data = window.SOURCES && window.SOURCES[id];
        if (!data) { reject(new Error('the file loaded but held no data')); return; }
        loaded.set(id, data);
        try { localStorage.setItem('sources.on', JSON.stringify([...loaded.keys()])); } catch (e) {}
        listeners.forEach(fn => fn(id, data)); resolve(data);
      };
      s.onerror = () => reject(new Error(set.file + ' is not built yet — run tools/sources.py ' + id + ', then tools/build_keys.html'));
      document.head.appendChild(s);
    });
  }
  window.Sources = { SETS, load, loaded: () => [...loaded.keys()], get: id => loaded.get(id), onLoad: fn => listeners.push(fn) };
  /* reload the sets that were on last time */
  try { (JSON.parse(localStorage.getItem('sources.on')) || []).forEach(id => load(id).catch(() => {})); } catch (e) {}
})();
