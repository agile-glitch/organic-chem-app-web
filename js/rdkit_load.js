/* RDKitLoad — start RDKit.js once per page, shared by every page that needs it (Reactions, Molecule).

   RDKitLoad()            → Promise<RDKit module>; the first call starts it, later calls get the same promise.
   RDKitLoad.script(src)  → Promise; injects <script src> once per src (a failed load may be retried).
   RDKitLoad.display(smiles, onReady) → the SMILES to SHOW the user: RDKit's canonical form (the standard one, the same however the
                          molecule was typed) once RDKit is running; until then the text given, and RDKit is started and onReady()
                          called when it can answer, so the page can redraw. Only for text people read or copy; lookups keep their own.

   Why a <script> and base64: a page opened from file:// may not fetch() local files, so the WebAssembly comes in
   as a base64 string (vendor/rdkit/RDKit_minimal.wasm.js, written by tools/vendor_rdkit.py) and the page starts it
   itself through Emscripten's `instantiateWasm` hook (this build ignores the `wasmBinary` option).
   RDKit.js: vendor/rdkit/ (BSD-3). */
window.RDKitLoad = (() => {
  'use strict';
  const VENDOR = 'vendor/rdkit/', WASM_JS = VENDOR + 'RDKit_minimal.wasm.js', scripts = new Map();
  let started = null, mod = null;
  const canonCache = new Map(), waiting = new Set();

  function script(src) {
    if (!scripts.has(src)) {
      const p = new Promise((res, rej) => {
        const s = document.createElement('script'); s.src = src;
        s.onload = () => { s.remove(); res(); };
        s.onerror = () => { s.remove(); scripts.delete(src); rej(new Error('could not load ' + src)); };
        document.head.appendChild(s);
      });
      scripts.set(src, p);
    }
    return scripts.get(src);
  }

  function bytesOf(b64) { const bin = atob(b64), u8 = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i); return u8; }

  function RDKitLoad() {
    if (started) return started;
    started = (async () => {
      if (!window.initRDKitModule) await script(VENDOR + 'RDKit_minimal.js');
      if (!window.RDKIT_WASM_B64) await script(WASM_JS);
      if (!window.RDKIT_WASM_B64) throw new Error(WASM_JS + ' held no data');
      /* the 9.8 MB string is freed once decoded; a failed start below forgets the script so a retry reloads it */
      const wasm = bytesOf(window.RDKIT_WASM_B64); window.RDKIT_WASM_B64 = null;
      const R = await new Promise((resolve, reject) => {
        /* Emscripten waits for done() forever, so a failed instantiate must reject here or the page hangs */
        const opts = {
          locateFile: f => VENDOR + f,
          instantiateWasm: (imports, done) => {
            WebAssembly.instantiate(wasm, imports).then(r => done(r.instance, r.module),
              e => reject(new Error('RDKit could not start (WebAssembly): ' + (e && e.message || e))));
            return {};
          },
        };
        try { Promise.resolve(window.initRDKitModule(opts)).then(resolve, reject); } catch (e) { reject(e); }
      });
      try { R.disable_logging && R.disable_logging(); } catch (e) {}
      return R;
    })();
    started.then(R => { mod = R; const fns = [...waiting]; waiting.clear(); fns.forEach(f => { try { f(); } catch (e) { /* the page's own redraw */ } }); },
                 () => { started = null; scripts.delete(WASM_JS); waiting.clear(); });
    return started;
  }

  function display(smiles, onReady) {
    if (!smiles) return smiles;
    if (mod) {
      if (!canonCache.has(smiles)) {
        let c = null;
        try { const m = mod.get_mol(smiles); if (m) { c = m.get_smiles(); m.delete(); } } catch (e) { c = null; }
        if (canonCache.size > 5000) canonCache.clear();
        canonCache.set(smiles, c || smiles);
      }
      return canonCache.get(smiles);
    }
    if (onReady) waiting.add(onReady);
    RDKitLoad().catch(() => { /* RDKit could not start: the text given stays */ });
    return smiles;
  }

  RDKitLoad.script = script;
  RDKitLoad.display = display;
  return RDKitLoad;
})();
