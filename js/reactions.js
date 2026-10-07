(() => {
  'use strict';
  const page = document.getElementById('page-reactions');
  if (!page) return;
  page.style.setProperty('--ax', '50vw');   // arrow centre; relayout() moves these when the content outgrows the window
  page.style.setProperty('--ay', '50vh');
  page.style.setProperty('--rxn-gap', '38px');
  page.style.setProperty('--ah', '44px');   // half the arrow's length: it grows to fit what is written above and below it
  page.innerHTML = `
    <style>
      #rxn-row { position:absolute; width:max-content; left:calc(var(--ax) - var(--ah) - var(--rxn-gap)); top:calc(var(--ay) - 10px); transform:translate(-100%,-50%); display:flex; align-items:stretch; gap:var(--rxn-gap); }
      #rxn-full, #rxn-prod { position:fixed; top:96px; left:16px; width:calc(45.4vw - 10px); box-sizing:border-box; font:14px monospace; color:#aaa; background:#fff; border:1px solid transparent; z-index:2; border-radius:4px; padding:4px 6px; margin-left:-6px; outline:none; }
      #rxn-prod { left:55.15vw; width:auto; right:10px; margin-left:0; color:#000; text-align:right; }
      #rxn-reag { position:fixed; top:96px; left:45.4vw; width:9.75vw; box-sizing:border-box; font:14px monospace; color:#aaa; background:#fff; border:1px solid transparent; border-radius:4px; padding:4px 6px; outline:none; text-align:center; z-index:2; }
      #rxn-codebar { position:fixed; top:44px; left:0; right:0; height:42px; box-sizing:border-box; background:#f6f6f6; border-bottom:1px solid #e0e0e0; z-index:2; }
      #rxn-code { position:absolute; top:8px; left:10px; right:10px; box-sizing:border-box; font:14px monospace; color:#555; background:#fff; border:1px solid #ddd; border-radius:4px; padding:4px 6px; outline:none; }
      #rxn-code:hover { border-color:#ddd; }
      #rxn-gear { position:absolute; top:8px; right:10px; box-sizing:border-box; padding:0; display:flex; align-items:center; justify-content:center;
                  background:#fff; border:1px solid #ddd; border-radius:4px; color:#555; cursor:pointer; }
      #rxn-gear:hover, #rxn-gear.on { border-color:#888; color:#222; }
      #rxn-settings { position:absolute; top:calc(100% + 2px); right:10px; width:300px; background:#fff; border:1px solid #ccc; border-radius:6px;
                      box-shadow:0 4px 14px rgba(0,0,0,.12); padding:6px 0; font:13px Arial,sans-serif; color:#333; z-index:6; }
      .rxn-set-head { padding:4px 12px; font-size:11px; color:#888; text-transform:uppercase; letter-spacing:.04em; }
      .rxn-set-opt { display:block; padding:6px 12px 6px 34px; position:relative; cursor:pointer; }
      .rxn-set-opt:hover { background:#f2f6ff; }
      .rxn-set-opt input { position:absolute; left:12px; top:8px; margin:0; }
      .rxn-set-opt span { display:block; color:#777; font-size:12px; margin-top:2px; }
      #page-reactions.mode-basic .rxn-read, #page-reactions.mode-basic .rxn-mech { display:none; }
      .rxn-set-sep { border-top:1px solid #eee; margin:6px 0; }
      .rxn-set-btn { display:block; width:100%; text-align:left; background:none; border:0; padding:6px 12px; font:13px Arial,sans-serif; color:#333; cursor:pointer; }
      .rxn-set-btn:hover { background:#f2f6ff; }
      .rxn-set-btn span { display:block; color:#777; font-size:12px; margin-top:2px; }
      #rxn-tests { position:fixed; top:96px; right:12px; bottom:12px; width:440px; max-width:calc(100vw - 24px); box-sizing:border-box; background:#fff;
                   border:1px solid #ccc; border-radius:8px; box-shadow:0 6px 24px rgba(0,0,0,.16); z-index:8; display:flex; flex-direction:column; font:13px Arial,sans-serif; color:#333; }
      #rxn-tests[hidden] { display:none; }
      .rxn-t-head { display:flex; align-items:center; gap:10px; padding:10px 12px 6px; }
      .rxn-t-head span { flex:1; color:#555; }
      #rxn-t-close { border:0; background:none; font-size:20px; line-height:1; cursor:pointer; color:#777; }
      .rxn-t-bar { height:6px; margin:0 12px; background:#eee; border-radius:3px; overflow:hidden; }
      #rxn-t-fill { height:100%; width:0; background:#1e7d32; transition:width .2s; }
      #rxn-t-fill.bad { background:#b3261e; }
      .rxn-t-actions { display:flex; gap:6px; padding:8px 12px; }
      .rxn-t-actions button { font:12px Arial,sans-serif; padding:4px 8px; border:1px solid #ccc; background:#fafafa; border-radius:4px; cursor:pointer; }
      .rxn-t-actions button:disabled { opacity:.45; cursor:default; }
      #rxn-t-groups { padding:0 12px 6px; }
      .rxn-t-g { display:flex; justify-content:space-between; padding:2px 0; border-bottom:1px dotted #eee; }
      .rxn-t-g .ok { color:#1e7d32; } .rxn-t-g .bad { color:#b3261e; font-weight:bold; }
      #rxn-t-list { flex:1; overflow:auto; padding:4px 12px 12px; }
      .rxn-t-f { border:1px solid #f0c8c4; background:#fdf5f4; border-radius:6px; padding:6px 8px; margin:6px 0; cursor:pointer; }
      .rxn-t-f:hover { border-color:#b3261e; }
      .rxn-t-f .n { font-weight:bold; }
      .rxn-t-f .d { color:#555; font-size:12px; margin-top:2px; word-break:break-all; }
      .rxn-t-f code { font:11px monospace; }
      .rxn-t-done { color:#1e7d32; padding:6px 0; }
      #rxn-common-btn { position:absolute; top:8px; box-sizing:border-box; padding:0 12px; display:flex; align-items:center; background:#fff; border:1px solid #ddd;
                        border-radius:4px; color:#555; cursor:pointer; font:13px Arial,sans-serif; white-space:nowrap; }
      #rxn-common-btn:hover, #rxn-common-btn.on { border-color:#888; color:#222; }
      #rxn-common { position:fixed; top:96px; right:12px; bottom:12px; width:440px; max-width:calc(100vw - 24px); box-sizing:border-box; background:#fff;
                    border:1px solid #ccc; border-radius:8px; box-shadow:0 6px 24px rgba(0,0,0,.16); z-index:8; display:flex; flex-direction:column; font:13px Arial,sans-serif; color:#333; }
      #rxn-common[hidden] { display:none; }
      .rxn-c-head { display:flex; align-items:center; gap:10px; padding:10px 12px 6px; }
      .rxn-c-head span { flex:1; color:#555; }
      #rxn-c-close { border:0; background:none; font-size:20px; line-height:1; cursor:pointer; color:#777; }
      #rxn-c-q { margin:4px 12px 6px; padding:6px 8px; font:13px Arial,sans-serif; border:1px solid #ccc; border-radius:4px; }
      #rxn-c-q:focus { border-color:#888; outline:none; }
      #rxn-c-list { flex:1; overflow:auto; padding:0 12px 12px; }
      .rxn-c-h { margin:10px 0 4px; font-size:11px; color:#888; text-transform:uppercase; letter-spacing:.04em; display:flex; justify-content:space-between; }
      .rxn-c-item { display:block; width:100%; text-align:left; border:1px solid #e4e4e4; background:#fafafa; border-radius:6px; padding:6px 8px; margin:4px 0; cursor:pointer;
                    font:13px Arial,sans-serif; color:#333; }
      .rxn-c-item:hover, .rxn-c-item:focus { border-color:#888; background:#f2f6ff; outline:none; }
      .rxn-c-item .n { font-weight:bold; }
      .rxn-c-item .i { float:right; color:#999; font:11px monospace; }
      .rxn-c-item .d { display:block; color:#666; font-size:12px; margin-top:2px; }
      .rxn-c-none { color:#777; padding:12px 0; }
      #rxn-code:focus { color:#333; border-color:#888; }
      #rxn-code.bad { color:#b3261e; border-color:#b3261e; }
      #rxn-reag:hover { border-color:#ddd; }
      #rxn-reag:focus { color:#333; border-color:#888; }
      #rxn-reag.bad { color:#b3261e; border-color:#b3261e; }
      #rxn-prod::placeholder { color:#000; opacity:1; }
      #rxn-prod { cursor:default; user-select:text; -webkit-user-select:text; pointer-events:auto; }
      #rxn-prod:hover, #rxn-prod:focus { border-color:transparent; }
      #rxn-full:hover { border-color:#ddd; }
      #rxn-full:focus { color:#333; border-color:#888; }
      #rxn-full.bad { color:#b3261e; border-color:#b3261e; }      #rxn-lower { position:absolute; width:max-content; left:calc(var(--ax) - var(--ah) - var(--rxn-gap)); transform:translateX(-100%); display:flex; flex-direction:column; align-items:flex-end; gap:var(--rxn-gap); }
      #rxn-products { position:absolute; width:max-content; left:calc(var(--ax) + var(--ah) + var(--rxn-gap)); top:calc(var(--ay) - 10px); transform:translateY(-50%); display:none; align-items:stretch; gap:var(--rxn-gap); }
      #rxn-products.on { display:flex; }
      .rxn-menu { max-height:60vh; overflow-y:auto; position:absolute; top:100%; left:50%; transform:translateX(-50%); margin-top:6px; min-width:280px; max-width:420px; background:#fff; border:1px solid #ccc; border-radius:6px; box-shadow:0 4px 14px rgba(0,0,0,.12); font:13px Arial,sans-serif; color:#333; text-align:left; padding:4px 0; z-index:5; cursor:default; }
      #rxn-above .rxn-menu { top:auto; bottom:100%; margin-top:0; margin-bottom:6px; }   /* the reagent list opens upward, above the arrow */
      .rxn-menu-head { padding:7px 12px 3px; font-size:11px; color:#888; text-transform:uppercase; letter-spacing:.03em; }
      .rxn-menu-item { padding:6px 12px; cursor:pointer; }
      .rxn-menu-item:hover { background:#eef3ff; }
      .rxn-menu-item.on::before { content:'✓ '; color:#2e7d32; }
      .rxn-menu-sub { margin-left:8px; color:#888; font-size:12px; }
      .rxn-menu-smi { margin-left:6px; font:12px monospace; color:#666; }
      .rxn-menu-note { padding:4px 12px; color:#999; font-size:12px; }
      .rxn-menu-sep { border-top:1px solid #eee; margin:4px 0; }
      .rxn-mech { margin-top:14px; }
      .rxn-mech-status { font-size:12px; margin:-2px 0 6px; }
      .rxn-mech-status.ok { color:#1e7d32; }
      .rxn-mech-status.bad { color:#a15c00; }
      .rxn-verdict.kind-unknown { color:#a15c00; }
      .rxn-verdict.kind-none, .rxn-verdict.kind-blocked { color:#b3261e; }
      .rxn-read { margin-top:12px; }
      .rxn-read-smi { font:12px monospace; color:#555; margin-right:6px; }
      .rxn-frag { display:inline-block; font-size:11px; color:#1f4e78; background:#eef3fa; border-radius:3px; padding:0 4px; margin:1px 2px 1px 0; }
      .rxn-mech-head { font-weight:bold; color:#333; margin-bottom:6px; }
      .rxn-step { border-top:1px solid #eee; padding:8px 0; display:flex; gap:14px; align-items:flex-start; }
      .rxn-step-text { flex:1 1 300px; min-width:260px; }
      .rxn-step-title { font-weight:bold; color:#333; }
      .rxn-step-line { margin-top:3px; }
      .rxn-step-k { color:#888; }
      .rxn-step-pka.up { color:#a15c00; }
      .rxn-step-pka.down { color:#1e7d32; }
      .rxn-step-note { margin-top:4px; color:#777; font-size:12px; }
      .rxn-step-mol { flex-shrink:0; max-width:40%; text-align:center; align-self:center; }
      .rxn-step-before { flex-shrink:0; align-self:center; max-width:50%; overflow-x:auto; }
      .rxn-step-canvas { display:block; }
      .rxn-step-pics svg { max-width:100%; height:auto; }
      .rxn-step-go { align-self:center; font-size:22px; color:#555; flex-shrink:0; }
      .rxn-step-pics { display:flex; align-items:center; justify-content:center; flex-wrap:wrap; gap:4px; }
      .rxn-step-plus { color:#888; font-size:14px; }
      .rxn-step-mol svg { display:block; margin:0 auto; }
      .rxn-step-smi { font:11px monospace; color:#aaa; word-break:break-all; cursor:pointer; margin-top:2px; }
      .rxn-step-smi:hover { color:#666; }
      #rxn-lower { z-index:1; }
      #rxn-rule { z-index:0; position:absolute; left:var(--rxn-left, 32px); width:calc(100vw - var(--rxn-left, 32px) - 24px); box-sizing:border-box; font:13px Arial,sans-serif; color:#555; }
      .rxn-line { display:flex; align-items:stretch; gap:var(--rxn-gap); }
      .rxn-slot { flex-shrink:0; width:var(--w); height:var(--h); display:flex; flex-direction:column; align-items:center; justify-content:center; box-sizing:border-box; padding:16px; border:2px dashed #888; border-radius:10px; font-size:64px; color:#888; cursor:pointer; user-select:none; }
      .rxn-slot:hover { border-color:#333; color:#333; }
      .rxn-slot.open { cursor:default; font-size:14px; color:#333; border-style:solid; border-color:#333; }
      .rxn-slot.done { border:none; cursor:pointer; font-size:14px; width:auto; height:auto; min-width:0; min-height:var(--h); }
      .rxn-first { --w:280px; --h:235px; }
      .rxn-slot.rxn-locked { cursor:default; }
      .rxn-slot.rxn-locked:not(.done) { visibility:hidden; }      /* no products: nothing is drawn (the rule note says why) */
      .rxn-extra { --w:102px; --h:251px; font-size:40px; border-color:transparent; color:transparent; }
      .rxn-extra.waiting { visibility:hidden; }
      .rxn-extra:hover { border-color:#888; color:#888; }
      .rxn-extra.open { border-color:#333; color:#333; }
      .rxn-extra.done { border:none; --h:235px; }
      .rxn-agent { --w:160px; --h:100px; font-size:36px; border-color:transparent; color:transparent; position:absolute; left:var(--ax); transform:translateX(-50%); }
      .rxn-agent:hover { border-color:#888; color:#888; }
      .rxn-agent.open { border-color:#333; color:#333; }
      .rxn-agent.done { border:none; }
      @media (hover: none) { .rxn-extra:not(.done):not(.open), .rxn-agent:not(.done):not(.open) { border-color:#d0d0d0; color:#c0c0c0; } }
      .rxn-cond-text { font:14px Arial,sans-serif; color:#555; text-align:center; word-break:break-word; max-width:360px; margin:auto; cursor:pointer; }
      #rxn-above { top:calc(var(--ay) - 10px - var(--rxn-gap)); transform:translate(-50%,-100%); }
      #rxn-below { top:calc(var(--ay) + 10px + var(--rxn-gap)); }
      .rxn-slot input { width:100%; min-width:70px; box-sizing:border-box; padding:8px; font:15px monospace; border:1px solid #888; border-radius:4px; }
      .rxn-slot input.bad { border-color:#b3261e; }
      .rxn-msg { margin-top:8px; min-height:18px; font:12px Arial,sans-serif; color:#b3261e; text-align:center; }
      .rxn-slot.done { justify-content:stretch; padding-bottom:8px; }
      .rxn-mol { display:flex; flex-direction:column; align-items:center; flex:1; align-self:stretch; }
      .rxn-mol svg { display:block; margin:auto 0; }
      .rxn-mol-smiles { margin-top:8px; font:13px monospace; color:#aaa; word-break:break-all; text-align:center; cursor:pointer; }
      .rxn-mol-smiles:hover { color:#666; }
    </style>
    <div id="rxn-codebar"><input id="rxn-code" type="text" placeholder="Reaction code: starting materials~reagent~conditions~products" spellcheck="false" autocomplete="off"><button id="rxn-gear" type="button" title="Settings: Basic / Advanced" aria-label="Settings"><svg viewBox="0 0 24 24" width="62%" height="62%" aria-hidden="true"><path fill="currentColor" d="M19.14 12.94c.04-.31.06-.63.06-.94s-.02-.63-.06-.94l2.03-1.58a.5.5 0 0 0 .12-.64l-1.92-3.32a.5.5 0 0 0-.61-.22l-2.39.96a7.03 7.03 0 0 0-1.62-.94l-.36-2.54A.5.5 0 0 0 13.9 2h-3.84a.5.5 0 0 0-.49.42l-.36 2.54c-.59.24-1.13.56-1.62.94l-2.39-.96a.5.5 0 0 0-.61.22L2.67 8.48a.5.5 0 0 0 .12.64l2.03 1.58c-.04.31-.07.63-.07.94s.03.63.07.94l-2.03 1.58a.5.5 0 0 0-.12.64l1.92 3.32c.13.22.39.3.61.22l2.39-.96c.49.38 1.03.7 1.62.94l.36 2.54c.05.24.25.42.49.42h3.84c.24 0 .45-.18.49-.42l.36-2.54c.59-.24 1.13-.56 1.62-.94l2.39.96c.22.08.48 0 .61-.22l1.92-3.32a.5.5 0 0 0-.12-.64l-2.03-1.58zM12 15.6a3.6 3.6 0 1 1 0-7.2 3.6 3.6 0 0 1 0 7.2z"/></svg></button>
      <button id="rxn-common-btn" type="button" aria-expanded="false" aria-controls="rxn-common" title="Browse the common reactions and load one">Common reactions</button>
      <div id="rxn-settings" hidden>
        <div class="rxn-set-head">Show</div>
        <label class="rxn-set-opt"><input type="radio" name="rxn-mode" value="basic"> <b>Basic</b><span>the products and a one-line verdict</span></label>
        <label class="rxn-set-opt"><input type="radio" name="rxn-mode" value="advanced"> <b>Advanced</b><span>also how it was read (fragments), chemoselectivity, stereo, amounts, pKa and the mechanism with arrows</span></label>
        <div class="rxn-set-sep"></div>
        <div class="rxn-set-head">Tests</div>
        <button type="button" class="rxn-set-btn" id="rxn-run-tests"><b>Run all tests</b><span>every reaction in the workbook (its stored reaction code) plus the feature cases in tools/reaction_tests.json, through this page</span></button>
      <div class="rxn-set-sep"></div>
      <div class="rxn-set-head">Reactions</div>
      <button type="button" class="rxn-set-btn" id="rxn-open-common"><b>Common reactions</b><span>browse the known reactions, search them and load one into the code bar</span></button>
      </div></div>
    <div id="rxn-common" hidden role="dialog" aria-label="Common reactions">
      <div class="rxn-c-head"><b>Common reactions</b><span id="rxn-c-sum"></span><button type="button" id="rxn-c-close" title="Close" aria-label="Close">×</button></div>
      <input id="rxn-c-q" type="search" placeholder="Search: name, reagent, group (e.g. Grignard, alkene, NaBH4)" spellcheck="false" autocomplete="off" aria-label="Search the common reactions">
      <div id="rxn-c-list"></div>
    </div>
    <div id="rxn-tests" hidden>
      <div class="rxn-t-head"><b>Tests</b><span id="rxn-t-sum"></span><button type="button" id="rxn-t-close" title="Close">×</button></div>
      <div class="rxn-t-bar"><div id="rxn-t-fill"></div></div>
      <div class="rxn-t-actions"><button type="button" id="rxn-t-again">Run again</button><button type="button" id="rxn-t-failed">Run failures only</button><button type="button" id="rxn-t-stop">Stop</button></div>
      <div id="rxn-t-groups"></div>
      <div id="rxn-t-list"></div>
    </div>
    <input id="rxn-full" type="text" placeholder="SMILES line" spellcheck="false" autocomplete="off">
    <input id="rxn-reag" type="text" placeholder="Reagent SMILES line" spellcheck="false" autocomplete="off">
    <input id="rxn-prod" type="text" placeholder="Product SMILES line" spellcheck="false" autocomplete="off" readonly>
    <div id="rxn-row"></div>
    <div id="rxn-lower"></div>
    <div id="rxn-above" class="rxn-slot rxn-agent" title="Add a molecule">+</div>
    <div id="rxn-below" class="rxn-slot rxn-agent" title="Add a molecule">+</div>
    <div id="rxn-products"></div>
    <div id="rxn-spacer" style="position:absolute;left:0;top:0;visibility:hidden;pointer-events:none;"></div>
    <svg id="rxn-arrow" style="position:absolute;left:var(--ax);top:var(--ay);transform:translate(-50%,-50%);cursor:pointer;" width="88" height="20" viewBox="0 0 88 20" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><line x1="2" y1="10" x2="86" y2="10"/><polyline points="77,3 86,10 77,17"/></svg>`;

  const BOND_PX = 34;        // one bond length on screen, the same for every molecule
  const GAP = 38;            // space between boxes, between rows, and between the last box and the arrow
  const TOP_GAP = 16;        // space kept between the SMILES lines and the highest box
  const LEFT_MARGIN = 32;    // how close to the left edge of the window a row may run before the next molecule drops to a new row
  const ARROW_MIN = 88;      // the arrow's length with nothing written above or below it
  const NOT_SMILES = /[^A-Za-z0-9@+\-\[\]()=#$\/\\%.:*]/g;
  const row1 = document.getElementById('rxn-row'), lower = document.getElementById('rxn-lower');
  row1.style.setProperty('--rxn-gap', GAP + 'px'); lower.style.setProperty('--rxn-gap', GAP + 'px');
  /* Molecules are drawn with the Molecule tab's 2D layout (RDKit, MolView.depict2D), so a SMILES looks the same on both
     pages. The app's graph keeps its atoms; positions, double-bond placement and wedges come from RDKit's drawing (a
     wedge only means what it should with the coordinates it was made for). RDKit numbers atoms in SMILES order, the
     same as the app's parser; when the two do not match one to one (explicit [H] atoms, which RDKit removes), the
     drawing falls back to the app's own layout. */
  let RD = null;
  const loadRD = () => (RD ? Promise.resolve(RD) : window.RDKitLoad ? window.RDKitLoad().then(R => (RD = R), () => null) : Promise.resolve(null));
  function molTabLayout(smiles, graph) {
    if (!RD || !window.MolView || !window.MolView.depict2D) return null;
    const P = window.MolView.depict2D(RD, smiles);
    if (!P || P.atoms.length !== graph.atoms.length || P.bonds.length !== graph.bonds.length) return null;
    const sameEl = (p, a) => p.el === a.element || ((p.el === 'D' || p.el === 'T') && a.element === 'H');
    if (P.atoms.some((p, i) => !sameEl(p, graph.atoms[i]))) return null;
    const idx = new Map(graph.atoms.map((a, i) => [a.id, i]));
    const bonds = [];
    for (const b of graph.bonds) {
      const ia = idx.get(b.a), ib = idx.get(b.b);
      const pb = P.bonds.find(x => (x.a === ia && x.b === ib) || (x.a === ib && x.b === ia));
      if (!pb) return null;
      const nb = Object.assign({}, b, { order: pb.order >= 1 && pb.order <= 3 ? pb.order : b.order });
      delete nb.stereo; delete nb.narrow;
      if (pb.flag === 1 || pb.flag === 6) { nb.stereo = pb.flag === 1 ? 'wedge' : 'dash'; nb.narrow = graph.atoms[pb.a].id; }
      bonds.push(nb);
    }
    // molfile y points up, the drawing's y points down
    return { atoms: graph.atoms.map((a, i) => Object.assign({}, a, { x: P.atoms[i].x, y: -P.atoms[i].y })), bonds, nextId: graph.nextId };
  }

  const esc = t => String(t).replace(/[&<>"]/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[ch]));

  /* a pop-up list under a box: build(add, close) fills it, add(cssClass, html, onClick) adds one row */
  function popMenu(box, build) {
    const menu = document.createElement('div');
    menu.className = 'rxn-menu';
    menu.innerHTML = '<div class="rxn-menu-note">Looking through the reaction rules…</div>';
    box.appendChild(menu);
    const close = () => { menu.remove(); document.removeEventListener('mousedown', outside, true); box._menu = null; };
    const outside = e => { if (!menu.contains(e.target) && !box.contains(e.target)) close(); };
    document.addEventListener('mousedown', outside, true);
    box._menu = close;
    const add = (cls, html, act) => {
      const d = document.createElement('div'); d.className = cls; d.innerHTML = html;
      if (act) d.addEventListener('click', e => { e.stopPropagation(); act(); });
      menu.appendChild(d);
    };
    Promise.resolve(build(add, close, () => { menu.innerHTML = ''; })).catch(() => {});
  }

  function fallbackCopy(text, done) {
    const ta = document.createElement('textarea');
    ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
    document.body.appendChild(ta); ta.select();
    try { if (document.execCommand('copy')) done(); } catch (e) { /* nothing more to try */ }
    ta.remove();
  }

  // Boxes are filled in order and packed into rows that all end at the same right edge, just left of the arrow.
  // When a row is full the next box starts a new row below. The last box is always empty and sits at the end of the last row.
  let queued = false;
  function relayout() {
    if (queued) return;
    queued = true;
    const run = () => {
      if (!queued) return;
      queued = false;
      const vw = document.documentElement.clientWidth, vh = document.documentElement.clientHeight;
      try { sizeGear(); } catch (e) { /* the gear is set up a little later during start-up */ }
      // the arrow spans the widest reagent or condition written over/under it, so neither can reach the molecules
      const above = document.getElementById('rxn-above'), below = document.getElementById('rxn-below');
      const aw = Math.max(ARROW_MIN, ...[above, below].map(b => b.classList.contains('done') ? b.offsetWidth : 0));
      const ah = aw / 2;
      page.style.setProperty('--ah', ah + 'px');
      const arrow = document.getElementById('rxn-arrow');
      arrow.setAttribute('width', aw); arrow.setAttribute('viewBox', `0 0 ${aw} 20`);
      arrow.querySelector('line').setAttribute('x2', aw - 2);
      arrow.querySelector('polyline').setAttribute('points', `${aw - 11},3 ${aw - 2},10 ${aw - 11},17`);
      // nothing may sit under the tab bar or the two SMILES lines: everything starts at least TOP_GAP below them
      const topLimit = fullInput.getBoundingClientRect().bottom + TOP_GAP;
      const maxW = vw / 2 - ah - GAP - LEFT_MARGIN;
      const rows = [[]]; let used = 0;
      slots.forEach(s => {
        const w = s.box.offsetWidth, cur = rows[rows.length - 1];
        const need = used + (cur.length ? GAP : 0) + w;
        if (cur.length && need > maxW) { rows.push([s]); used = w; }
        else { cur.push(s); used = need; }
      });
      rows.forEach((r, i) => {
        let host;
        if (i === 0) host = row1;
        else {
          host = lower.children[i - 1];
          if (!host) { host = document.createElement('div'); host.className = 'rxn-line'; lower.appendChild(host); }
        }
        r.forEach((s, j) => { if (host.children[j] !== s.box) host.insertBefore(s.box, host.children[j] || null); });
      });
      while (lower.children.length > rows.length - 1) lower.lastChild.remove();
      // Put the arrow where everything fits: at the window centre, or further right/down when the content is too big.
      const rowsW = Math.max(...[row1, ...lower.children].map(e => e.offsetWidth));
      const ha = above.offsetHeight;
      const hp = products.classList.contains('on') ? products.offsetHeight : 0;
      page.style.setProperty('--ax', Math.max(vw / 2, rowsW + LEFT_MARGIN + GAP + ah) + 'px');
      page.style.setProperty('--ay', Math.max(vh / 2, topLimit + 10 + Math.max(row1.offsetHeight, hp) / 2, topLimit + 10 + GAP + ha) + 'px');
      lower.style.top = (row1.getBoundingClientRect().bottom + window.scrollY + GAP) + 'px';
      // The rule note, how it was read and the mechanism use the full width under the whole reaction.
      const note = document.getElementById('rxn-rule');
      if (note) {
        let low = 0;
        // a lower row holding only the invisible empty "add" box doesn't count (it stays clickable above the text)
        [row1, ...[...lower.children].filter(h => h.querySelector('.done')), products, above, below].forEach(e => {
          const r = e.getBoundingClientRect();
          if (r.width || r.height) low = Math.max(low, r.bottom + window.scrollY);
        });
        note.style.setProperty('--rxn-left', LEFT_MARGIN + 'px');
        note.style.width = Math.max(300, vw - LEFT_MARGIN - 24) + 'px';   // the visible width (100vw would run under the scrollbar)
        note.style.top = (low + GAP) + 'px';
      }
      // Give the page room to scroll to the far edge of everything on it.
      let right = 0, bottom = 0;
      [row1, lower, products, document.getElementById('rxn-above'), document.getElementById('rxn-below'), document.getElementById('rxn-rule')].forEach(e => {
        if (!e) return;
        const r = e.getBoundingClientRect();
        if (r.width || r.height) { right = Math.max(right, r.right + window.scrollX); bottom = Math.max(bottom, r.bottom + window.scrollY); }
      });
      const spacer = document.getElementById('rxn-spacer');
      spacer.style.width = (right + 40) + 'px'; spacer.style.height = (bottom + 40) + 'px';
    };
    requestAnimationFrame(run);
    setTimeout(run, 60);   // also runs when the window is hidden and frames are paused
  }
  const watcher = new ResizeObserver(relayout);
  window.addEventListener('resize', relayout);

  /* Typing a SMILES into a box, drawing the molecule, click-to-copy. Used by every box on the page.
     Click a placed molecule to edit its SMILES; clear the text and press Enter (or click away) to remove it; Esc
     cancels the edit. baseCls is the box's class list, or a function giving it (a box's place in a row can change).
     A locked box can't be typed into or edited: only code fills it (the product boxes, filled by the reaction rules). */
  function wireBox(box, baseCls, onFilled, onCleared, menuFor, locked) {
    const cls = extra => (typeof baseCls === 'function' ? baseCls() : baseCls) + (extra ? ' ' + extra : '');
    let cur = null, curGraph = null;                // the molecule placed in this box

    function reset() { box.className = cls(); box.title = locked ? 'Products come from the reaction rules: click the arrow' : 'Add a molecule'; box.textContent = locked ? '' : '+'; relayout(); }

    async function showMolecule(smiles, graph) {
      const R = await loadRD();                      // the Molecule tab's layout needs RDKit (loading starts when a box opens)
      const typed = smiles;
      if (R) {                                       // one spelling everywhere: the box, the SMILES lines and the code bar
        let m = null; try { m = R.get_mol(smiles); } catch (e) { m = null; }
        const c = m ? m.get_smiles() : null; if (m) m.delete();
        if (c && c !== smiles) { let g2 = null; try { g2 = window.Chem.parseSmiles(c); } catch (e) { g2 = null; } if (g2) { smiles = c; graph = g2; } }
      }
      cur = smiles; curGraph = graph;
      box.className = cls('done'); box.title = locked ? 'Drawn by the reaction rules' : 'Click to edit · clear the SMILES to remove';
      box.innerHTML = '<div class="rxn-mol"></div>';
      const wrap = box.firstChild;
      const laid = molTabLayout(smiles, graph);
      const drawing = laid ? window.MolDraw.svg(laid, 3000, 3000, { layout: false }) : window.MolDraw.svg(graph, 3000, 3000);
      wrap.appendChild(drawing);
      const bb = drawing.getBBox(), pad = 6, k = BOND_PX / 30;
      drawing.setAttribute('viewBox', `${bb.x - pad} ${bb.y - pad} ${bb.width + 2 * pad} ${bb.height + 2 * pad}`);
      drawing.setAttribute('width', (bb.width + 2 * pad) * k);
      drawing.setAttribute('height', (bb.height + 2 * pad) * k);
      const label = document.createElement('div');
      label.className = 'rxn-mol-smiles'; label.textContent = smiles;
      label.style.width = drawing.getAttribute('width') + 'px';
      label.title = (typed !== smiles ? `Typed as ${typed} (written here in the standard form) · ` : '') + 'Click to copy';
      label.addEventListener('click', e => {
        e.stopPropagation();                        // copying is not editing
        const done = () => { label.textContent = 'Copied'; setTimeout(() => { label.textContent = smiles; }, 900); };
        if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(smiles).then(done, () => fallbackCopy(smiles, done));
        else fallbackCopy(smiles, done);
      });
      wrap.appendChild(label);
      relayout();
      if (onFilled) onFilled(smiles, graph);
    }

    function openInput() {
      loadRD();
      box.className = cls('open'); box.title = '';
      box.innerHTML = '<input type="text" placeholder="SMILES, e.g. C=CC" spellcheck="false" autocomplete="off"><div class="rxn-msg"></div>';
      const input = box.querySelector('input'), msg = box.querySelector('.rxn-msg');
      let finished = false;
      input.value = cur || '';
      input.focus(); input.select();
      function commit() {
        if (finished) return;
        const s = input.value.trim();
        if (!s) {                                   // empty: an empty box stays empty, a placed molecule is removed
          finished = true;
          const had = cur !== null;
          cur = null; curGraph = null;
          if (had && onCleared) onCleared(); else reset();
          return;
        }
        let g;
        try { g = window.Chem.parseSmiles(s); } catch (err) { g = null; }
        if (!g) { input.className = 'bad'; msg.textContent = 'Not a valid SMILES'; return; }
        finished = true; showMolecule(s, g);
      }
      function cancel() {
        finished = true;
        if (cur !== null) showMolecule(cur, curGraph); else reset();
      }
      input.addEventListener('input', () => { input.value = input.value.replace(NOT_SMILES, ''); input.className = ''; msg.textContent = ''; });
      input.addEventListener('keydown', e => { if (e.key === 'Enter') commit(); else if (e.key === 'Escape') cancel(); });
      input.addEventListener('blur', commit);
    }

    box.addEventListener('click', () => {
      if (locked || box.classList.contains('open')) return;
      if (box._menu) { box._menu(); return; }
      if (cur === null && menuFor) popMenu(box, (add, close, empty) => menuFor(add, close, empty, openInput));
      else openInput();
    });
    return { show: showMolecule, get: () => cur, clear: () => { cur = null; curGraph = null; reset(); } };
  }

  // The conditions box: plain text instead of a SMILES. Click a filled box to edit it; clearing the text empties the box.
  function wireText(box, baseCls, optionsFor) {
    const cls = extra => baseCls + (extra ? ' ' + extra : '');
    let text = '';
    function reset() { text = ''; box.className = cls(); box.title = 'Choose conditions'; box.textContent = '+'; relayout(); updateCode(); }
    function show() {
      box.className = cls('done'); box.title = 'Click to change the conditions';
      box.innerHTML = '<div class="rxn-cond-text"></div>';
      box.firstChild.textContent = text;
      relayout(); updateCode();
    }
    function openInput() {
      box.className = cls('open'); box.title = '';
      box.innerHTML = '<input type="text" placeholder="e.g. heat, hv, 25 C, CH2Cl2" spellcheck="false" autocomplete="off">';
      const input = box.querySelector('input');
      input.value = text; input.focus(); input.select();
      let finished = false;
      function commit() {
        if (finished) return;
        finished = true;
        text = input.value.trim();
        if (text) show(); else reset();
      }
      input.addEventListener('keydown', e => { if (e.key === 'Enter') commit(); });
      input.addEventListener('blur', commit);
    }
    /* Clicking opens a list: the conditions the matching rules need (optionsFor), the other conditions the rules know,
       then "Type your own…". Picking a condition adds it to the box (a comma list); picking it again removes it. */
    let menu = null;
    const parts = () => splitConds(text).map(x => x.trim()).filter(Boolean);
    function closeMenu() { if (menu) { menu.remove(); menu = null; document.removeEventListener('mousedown', outside, true); } }
    function outside(e) { if (menu && !menu.contains(e.target) && !box.contains(e.target)) closeMenu(); }
    function setParts(list) { text = list.join(', '); if (text) show(); else reset(); }
    function toggle(c) {
      const p = parts(), k = p.findIndex(x => x.toLowerCase() === c.toLowerCase());
      if (k >= 0) p.splice(k, 1); else p.push(c);
      closeMenu(); setParts(p);
    }
    async function openMenu() {
      menu = document.createElement('div');
      menu.className = 'rxn-menu';
      menu.innerHTML = '<div class="rxn-menu-note">Finding the conditions…</div>';
      box.appendChild(menu);
      document.addEventListener('mousedown', outside, true);
      const opts = optionsFor ? await optionsFor() : { matched: [], free: [], other: [] };
      if (!menu) return;
      menu.innerHTML = '';
      const have = parts().map(x => x.toLowerCase());
      const add = (cls, html, act) => {
        const d = document.createElement('div'); d.className = cls; d.innerHTML = html;
        if (act) d.addEventListener('click', e => { e.stopPropagation(); act(); });
        menu.appendChild(d);
      };
      if ((opts.suggest || []).length) {
        add('rxn-menu-head', 'Suggested by the fragments in these molecules');
        opts.suggest.forEach(sg => {
          const items = splitConds(sg.text).map(x => x.trim()).filter(Boolean);
          const on = items.every(x => have.includes(x.toLowerCase()));
          add('rxn-menu-item' + (on ? ' on' : ''), `${esc(sg.text)}<span class="rxn-menu-sub">${esc(sg.why)}</span>`, () => {
            const p = parts();
            items.forEach(x => { if (!p.some(y => y.toLowerCase() === x.toLowerCase())) p.push(x); });
            closeMenu(); setParts(p);
          });
        });
      }
      add('rxn-menu-head', 'Needed by reactions that match');
      if (!opts.matched.length && !opts.free.length) add('rxn-menu-note', esc(opts.why || 'No rule matches these molecules and this reagent yet.'));
      opts.matched.forEach(o => add('rxn-menu-item' + (have.includes(o.cond.toLowerCase()) ? ' on' : ''), `${esc(o.label || o.cond)}<span class="rxn-menu-sub">${esc(o.rules.join(', '))}</span>`, () => toggle(o.cond)));
      opts.free.forEach(r => add('rxn-menu-note', `${esc(r)} needs no special conditions`));
      if (opts.other.length) {
        add('rxn-menu-head', 'Other conditions the rules know');
        opts.other.forEach(c => add('rxn-menu-item' + (have.includes(c) ? ' on' : ''), esc(c), () => toggle(c)));
      }
      if ((opts.solvents || []).length) {
        add('rxn-menu-head', 'Solvents (reflux runs at the boiling point)');
        const solv = (window.REACTION_RULES || {}).solvents || {};
        opts.solvents.forEach(c => add('rxn-menu-item' + (have.includes(c.toLowerCase()) ? ' on' : ''), `${esc(c)}<span class="rxn-menu-sub">bp ${esc(fmtT(solv[c].bp))}</span>`, () => toggle(c)));
      }
      add('rxn-menu-sep', '');
      add('rxn-menu-item', 'Type your own…', () => { closeMenu(); openInput(); });
      if (text) add('rxn-menu-item', 'Clear', () => { closeMenu(); setParts([]); });
    }
    box.addEventListener('click', () => { if (box.classList.contains('open')) return; if (menu) closeMenu(); else openMenu(); });
    return { get: () => text, set: t => { closeMenu(); setParts(splitConds(t).map(x => x.trim()).filter(Boolean)); } };
  }

  /* A row of molecule boxes: the first is the visible dashed box, the rest are invisible until hovered. There is
     always exactly one empty box at the end (and at least two boxes); removing a molecule closes the gap. */
  function makeList(host, title, onChange, locked) {
    const list = [];
    const base = slot => 'rxn-slot ' + (list.indexOf(slot) === 0 ? 'rxn-first' : 'rxn-extra') + (locked ? ' rxn-locked' : '');
    function addBox() {
      const box = document.createElement('div');
      box.title = title; box.textContent = locked ? '' : '+';
      host.appendChild(box);
      const slot = { box, smiles: null, graph: null };
      list.push(slot);
      watcher.observe(box);
      slot.api = wireBox(box, () => base(slot),
        (smiles, graph) => { slot.smiles = smiles; slot.graph = graph; tidy(); onChange(); },
        () => { remove(slot); }, null, locked);
    }
    function remove(slot) {
      watcher.unobserve(slot.box); slot.box.remove();
      list.splice(list.indexOf(slot), 1);
      tidy(); onChange();
    }
    function tidy() {
      // filled boxes in the same order as the SMILES line; a box being edited and the empty boxes keep their places
      const filled = list.filter(x => x.smiles && x.graph && !x.box.classList.contains('open'));
      if (filled.length > 1) {
        let sorted = filled; try { sorted = sortByRules(filled); } catch (e) { sorted = filled; }
        let k = 0;
        for (let i = 0; i < list.length; i++) if (filled.includes(list[i])) list[i] = sorted[k++];
        list.forEach(x => host.appendChild(x.box));    // DOM order = list order (the left rows are re-packed by relayout)
      }
      if (locked) { if (!list.length) addBox(); }
      else {
        if (!list.length || list[list.length - 1].smiles) addBox();
        if (list.length === 1) addBox();
      }
      list.forEach((s, i) => {
        if (s.box.classList.contains('open')) return;
        s.box.className = base(s) + (s.smiles ? ' done' : i === 1 && !list[0].smiles ? ' waiting' : '');
      });
      relayout();
    }
    tidy();
    list.setAll = async smilesList => {                  // replace the row with these molecules
      for (const s of list.slice()) { watcher.unobserve(s.box); s.box.remove(); }
      list.length = 0;
      tidy();
      for (const smi of smilesList) {
        let g; try { g = window.Chem.parseSmiles(smi); } catch (e) { continue; }
        if (!list.some(x => !x.smiles)) addBox();         // (an unlocked row adds the next empty box by itself)
        await list.find(x => !x.smiles).api.show(smi, g);
      }
      onChange();
    };
    return list;
  }

  const fullInput = document.getElementById('rxn-full');
  // Order: most carbon atoms first; then most aromatic atoms; then most C-C, then C=C, then C#C bonds; then the order added.
  // Aromatic atoms are read from the SMILES text (the parser turns rings into alternating single/double bonds and forgets
  // which atoms were aromatic). Bonds are counted on that parsed form, so each benzene ring counts as 3 C-C and 3 C=C.
  function aromaticAtoms(smiles) {
    let n = 0;
    for (let i = 0; i < smiles.length; i++) {
      const ch = smiles[i];
      if (ch === '[') {
        const end = smiles.indexOf(']', i), m = smiles.slice(i + 1, end).match(/^\d*([A-Za-z])/);
        if (m && m[1] >= 'a' && m[1] <= 'z') n++;
        i = end;
      } else if ('bcnops'.includes(ch)) n++;
      else if ((ch === 'C' && smiles[i + 1] === 'l') || (ch === 'B' && smiles[i + 1] === 'r')) i++;
    }
    return n;
  }
  function carbonStats(g) {
    const isC = id => g.atoms.find(a => a.id === id).element === 'C';
    const bonds = order => g.bonds.filter(b => b.order === order && isC(b.a) && isC(b.b)).length;
    return { carbons: g.atoms.filter(a => a.element === 'C').length, cc: bonds(1), dbl: bonds(2), trp: bonds(3) };
  }
  // Last tie-break before "order added": heaviest first, by amu with implicit hydrogens included. Molecules with no
  // carbon already sit after all the carbon ones, so this orders them among themselves too.
  function amu(g, smiles) {
    let total = 0;
    for (const [, el, n] of window.Chem.formula(g).matchAll(/([A-Z][a-z]?)(\d*)/g)) total += window.MolData.element(el).aw * (n ? +n : 1);
    // the app's SMILES parser drops isotope labels, so they are read from the text: [2H], [13CH3] … each weigh their
    // mass number (within 0.02 u of the true isotope mass) instead of the standard atomic weight
    for (const [, iso, sym] of String(smiles || '').matchAll(/\[(\d+)([A-Z][a-z]?|[a-z]{1,2})/g)) {
      const el = sym.charAt(0).toUpperCase() + sym.slice(1), e = window.MolData.element(el);
      if (e && e.aw) total += +iso - e.aw;
    }
    return total;
  }
  function sortByRules(list) {                      // the SMILES-line order, used by the lines AND the boxes
    return list.filter(s => s.smiles)
      .map((s, i) => { const st = carbonStats(s.graph); return { item: s, s: s.smiles, i, arom: aromaticAtoms(s.smiles), ...st, mass: amu(s.graph, s.smiles) }; })
      .sort((a, b) => b.carbons - a.carbons || b.arom - a.arom || b.cc - a.cc || b.dbl - a.dbl || b.trp - a.trp || b.mass - a.mass || a.i - b.i)
      .map(x => x.item);
  }
  function orderedSmiles(list) { return sortByRules(list).map(x => x.smiles).join('.'); }
  function updateFull() { fullInput.value = orderedSmiles(slots); updateCode(); }

  // Clicking into the substrate line selects all of its text. The first mouseup is cancelled so it doesn't undo that.
  let justFocused = false, viaMouse = false;
  fullInput.addEventListener('mousedown', () => { viaMouse = document.activeElement !== fullInput; });
  fullInput.addEventListener('focus', () => { if (fullInput.value) { fullInput.select(); justFocused = viaMouse; } viaMouse = false; });
  fullInput.addEventListener('mouseup', e => { if (justFocused) { e.preventDefault(); justFocused = false; } });
  fullInput.addEventListener('blur', () => { justFocused = false; });
  // Enter: the starting materials become exactly what the line says, one box per dot-separated SMILES
  fullInput.addEventListener('keydown', async e => {
    if (e.key !== 'Enter') return;
    e.preventDefault();
    const parts = fullInput.value.split('.').map(x => x.trim()).filter(Boolean);
    const bad = parts.find(p => { try { return !window.Chem.parseSmiles(p); } catch (err) { return true; } });
    if (bad) { fullInput.classList.add('bad'); fullInput.title = 'Not a valid SMILES: ' + bad; return; }
    fullInput.classList.remove('bad'); fullInput.title = '';
    fullInput.blur();
    await slots.setAll(parts);
  });
  fullInput.addEventListener('input', () => { fullInput.classList.remove('bad'); fullInput.title = ''; });

  // The product line shrinks its text (down to 6px) so the whole SMILES stays visible in its box.
  const prodInput = document.getElementById('rxn-prod');
  const measure = document.createElement('canvas').getContext('2d');
  function fitLine(el) {                                // shrink the text (down to 6px) so all of it shows
    const avail = el.clientWidth - 14;
    if (avail <= 0) return;
    measure.font = '14px monospace';
    const w = measure.measureText(el.value || el.placeholder).width;
    el.style.fontSize = Math.max(6, Math.min(14, 14 * avail / (w || 1))) + 'px';
  }
  const fitProduct = () => fitLine(prodInput);
  prodInput.addEventListener('input', fitProduct);
  prodInput.addEventListener('click', () => { if (prodInput.value) prodInput.select(); });   // one click selects the whole line, ready to copy
  new ResizeObserver(fitProduct).observe(prodInput);

  const slots = makeList(row1, 'Add a molecule', updateFull);

  /* Clicking the arrow shows the product area and runs the reaction rules (data/reaction_rules.js, generated from
     reaction_rules.xlsx and molecule_data.xlsx by tools/rules_to_js.py). A rule is tried when its reagents are the
     molecule above the arrow (matched by InChIKey); its reaction SMARTS is run in RDKit on the molecules on the left,
     in every order. The first rule that gives a product fills the product boxes, followed by its by-products. */
  const ruleNote = document.createElement('div');
  ruleNote.id = 'rxn-rule'; ruleNote.style.display = 'none';
  page.appendChild(ruleNote);
  const keyOf = (R, smi) => { let m = null; try { m = R.get_mol(smi); return m ? R.get_inchikey_for_inchi(m.get_inchi()) : null; } catch (e) { return null; } finally { if (m) m.delete(); } };
  function runRule(R, smarts, subs) {
    let rx = null;
    try { rx = R.get_rxn(smarts); } catch (e) { rx = null; }
    if (!rx) return null;
    const n = smarts.split('>>')[0].split('.').length;
    const orders = [];
    const pick = (cur) => { if (cur.length === n) { orders.push(cur.slice()); return; } subs.forEach((_, i) => { if (!cur.includes(i)) { cur.push(i); pick(cur); cur.pop(); } }); };
    pick([]);
    try {
      for (const order of orders) {
        const ml = new R.MolList(), mols = order.map(i => R.get_mol(subs[i]));
        mols.forEach(m => ml.append(m));
        let out = null;
        try {
          const res = rx.run_reactants(ml, 10);
          if (res.size()) {
            const pl = res.get(0); out = [];
            for (let b = 0; b < pl.size(); b++) { const pm = pl.at(b); try { out.push(pm.get_smiles()); } catch (e) {} pm.delete(); }
            pl.delete();
          }
          res.delete();
        } catch (e) { out = null; }
        ml.delete(); mols.forEach(m => m.delete());
        if (out && out.length) return out.map(smi => { const m = R.get_mol(smi); const c = m ? m.get_smiles() : smi; if (m) m.delete(); return c; });
      }
      return null;
    } finally { rx.delete(); }
  }
  /* Reading the conditions box (sheets "Conditions" and "Solvents" in reaction_rules.xlsx).
     Each comma item can hold temperatures (80 °C, 353 K, 176 °F), named conditions (reflux, ice bath, anhydrous …)
     and solvents (ethanol, THF …), matched as whole words, longest phrase first ("dry ice" before "dry").
     Then: reflux runs at the solvent's boiling point; a temperature counts for every condition whose range holds it
     (80 °C -> heat, 0 °C -> ice bath); a condition counts as the ones it implies (reflux -> heat). */
  /* the conditions box is a comma list, but a comma between digits is part of a name (1,4-dioxane, 1,2-dichloroethane) */
  const splitConds = t => String(t || '').split(/,(?!\d)/);
  const toCelsius = (v, unit) => unit === 'k' ? v - 273.15 : unit === 'f' ? (v - 32) * 5 / 9 : v;   // unit: c, k or f
  const TEMP_RE = /(?<![\d.])(-?\d+(?:\.\d+)?)\s*°?\s*(c|k|f)(?![a-z])/g;
  const inRange = (t, lo, hi) => (lo == null || t >= lo) && (hi == null || t <= hi);
  const FILLER = new Set(['in', 'at', 'under', 'with', 'and', 'of', 'the', 'a', 'to', 'then']);
  function readItem(item) {
    const RR = window.REACTION_RULES || {}, conds = RR.conditions || {}, solv = RR.solvents || {};
    let low = ' ' + String(item).toLowerCase().replace(/[;:()/+]+/g, ' ').replace(/\s+/g, ' ').trim() + ' ';
    const r = { temps: [], keys: [], solvents: [], equiv: null };
    const EQ_RE = /(?<![\d.])(\d+(?:\.\d+)?)\s*(?:equivalents?|equiv|eq)\.?(?![a-z])/g;
    for (const m of low.matchAll(EQ_RE)) r.equiv = +m[1];
    low = low.replace(EQ_RE, ' ');
    for (const m of low.matchAll(TEMP_RE)) r.temps.push(toCelsius(+m[1], m[2]));
    low = low.replace(TEMP_RE, ' ');
    const phrases = [];
    for (const [k, c] of Object.entries(conds)) for (const w of c.words) phrases.push({ w, list: r.keys, k });
    for (const [k, v] of Object.entries(solv)) for (const w of v.words) phrases.push({ w, list: r.solvents, k });
    phrases.sort((a, b) => b.w.length - a.w.length);
    for (const ph of phrases) {
      const pat = ' ' + ph.w + ' ';
      if (!low.includes(pat)) continue;
      low = low.split(pat).join('  ');
      if (!ph.list.includes(ph.k)) ph.list.push(ph.k);
    }
    r.rest = low.split(' ').filter(w => w && !FILLER.has(w)).join(' ');
    return r;
  }
  function readConditions(text) {
    const RR = window.REACTION_RULES || {}, conds = RR.conditions || {}, solv = RR.solvents || {};
    const P = { temps: [], keys: new Set(), solvents: [], refluxT: null, equiv: null };
    for (const item of splitConds(text)) {
      const r = readItem(item);
      P.temps.push(...r.temps); r.keys.forEach(k => P.keys.add(k)); if (r.equiv != null) P.equiv = r.equiv; r.solvents.forEach(k => { if (!P.solvents.includes(k)) P.solvents.push(k); });
    }
    if (P.keys.has('reflux') && P.solvents.length) { P.refluxT = Math.min(...P.solvents.map(k => solv[k].bp)); P.temps.push(P.refluxT); }
    for (const [k, c] of Object.entries(conds))
      if ((c.minC != null || c.maxC != null) && P.temps.some(t => inRange(t, c.minC, c.maxC))) P.keys.add(k);
    for (let grew = true; grew;) {
      grew = false;
      for (const k of [...P.keys]) for (const i of (conds[k] && conds[k].implies) || []) if (!P.keys.has(i)) { P.keys.add(i); grew = true; }
    }
    return P;
  }
  /* one required item: a named condition, or a temperature range that a known temperature (typed, or reflux in a
     solvent) must fall in; a condition whose whole range fits also counts (ice bath for "<= 5 °C"), plain "heat" doesn't */
  function requirementMet(q, P) {
    const conds = (window.REACTION_RULES || {}).conditions || {};
    if (q.cond) return P.keys.has(q.cond);
    if (q.solvent) return P.solvents.includes(q.solvent);
    if (P.temps.some(t => inRange(t, q.tmin, q.tmax))) return true;
    return [...P.keys].some(k => { const c = conds[k]; return c && c.minC != null && c.maxC != null && inRange(c.minC, q.tmin, q.tmax) && inRange(c.maxC, q.tmin, q.tmax); });
  }
  const reqLabel = q => q.cond || q.solvent || q.label;
  const fmtT = t => minusSign(String(Math.round(t))) + ' °C';
  const minusSign = x => String(x).replace('-', '−');
  /* what the conditions list offers: the conditions needed by the rules that match the molecules on the left and the
     reagent above the arrow, the rules among them that need none, and every other condition the rules know */
  async function reagentMenu(add, close, empty, typeOwn) {
    const R = await loadRD(), RR = window.REACTION_RULES;
    empty();
    add('rxn-menu-head', 'Reagents that react with these starting materials');
    const subs = slots.filter(x => x.smiles).map(x => x.smiles);
    if (!R || !RR || !RR.rules) add('rxn-menu-note', 'The reaction rules are not loaded.');
    else if (!subs.length) add('rxn-menu-note', 'Add the starting molecules on the left first.');
    else {
      const byReagent = new Map(), free = [], byAB = new Map();
      for (const rule of RR.rules) {
        if (!runRule(R, rule.smarts, subs)) continue;
        if (rule.acidBase) {
          for (const e of (RR.acidsBases || []).filter(x => x.role === rule.acidBase.role && acidBaseCheck(rule.acidBase, x).ok)) {
            if (!byAB.has(e.smiles)) byAB.set(e.smiles, { e, rules: [] });
            byAB.get(e.smiles).rules.push(`${rule.id} ${rule.name}`);
          }
          continue;
        }
        if (!rule.reagents.length) { free.push(`${rule.id} ${rule.name}`); continue; }
        for (const k of rule.reagents) { if (!byReagent.has(k)) byReagent.set(k, []); byReagent.get(k).push(`${rule.id} ${rule.name}`); }
      }
      for (const { e, rules } of byAB.values()) {
        add('rxn-menu-item', `${esc(e.name)} <span class="rxn-menu-smi">${esc(e.smiles)}</span><span class="rxn-menu-sub">${esc(e.role)}, pKa ${esc(minus(e.pKa))} · ${esc(rules.join(', '))}</span>`, async () => {
          close();
          let g = null; try { g = window.Chem.parseSmiles(e.smiles); } catch (err) { g = null; }
          if (g) await aboveApi.show(e.smiles, g);
        });
      }
      if (!byReagent.size && !free.length && !byAB.size) add('rxn-menu-note', 'No rule in reaction_rules.xlsx reacts with these starting materials.');
      for (const [k, rules] of byReagent) {
        const m = RR.molecules[k];
        add('rxn-menu-item', `${esc(m.name)} <span class="rxn-menu-smi">${esc(m.smiles)}</span><span class="rxn-menu-sub">${esc(rules.join(', '))}</span>`, async () => {
          close();
          let g = null; try { g = window.Chem.parseSmiles(m.smiles); } catch (e) { g = null; }
          if (g) await aboveApi.show(m.smiles, g);
        });
      }
      free.forEach(r => add('rxn-menu-note', `${esc(r)} needs no reagent`));
    }
    add('rxn-menu-sep', '');
    add('rxn-menu-item', 'Type a SMILES…', () => { close(); typeOwn(); });
  }
  async function conditionOptions() {
    const R = await loadRD(), RR = window.REACTION_RULES, res = { matched: [], free: [], other: [] };
    if (!R || !RR || !RR.rules) { res.why = 'The reaction rules are not loaded.'; return res; }
    const subs = slots.filter(x => x.smiles).map(x => x.smiles);
    const reagent = aboveApi.get(), key = reagent ? keyOf(R, reagent) : null;
    const need = new Map();
    const offer = (cond, label, who) => { if (!need.has(cond)) need.set(cond, { label, rules: [] }); need.get(cond).rules.push(who); };
    if (!subs.length) res.why = 'Add the starting molecules on the left first.';
    else for (const rule of RR.rules) {
      if (rule.reagents.length && !rule.reagents.includes(key)) continue;
      if (!runRule(R, rule.smarts, subs)) continue;
      if (!(rule.required || []).length) { res.free.push(`${rule.id} ${rule.name}`); continue; }
      const who = `${rule.id} ${rule.name}` + (rule.required.length > 1 ? ' (one of its options)' : '');
      for (const alt of rule.required) {             // each alternative set is one choice ("reflux + anhydrous")
        const named = alt.filter(q => !q.label).map(q => q.cond || q.solvent);   // conditions and solvents
        const ranges = alt.filter(q => q.label);                                 // temperatures
        if (!ranges.length) { offer(named.join(', '), named.join(' + '), who); continue; }
        // a temperature range: offer a temperature inside it and reflux in each solvent that boils inside it
        for (const q of ranges) {
          const lo = q.tmin, hi = q.tmax;
          const mid = lo != null && hi != null ? Math.round((lo + hi) / 2 / 5) * 5 : lo != null ? lo + 30 : hi;
          offer([fmtT(mid), ...named].join(', '), [fmtT(mid), ...named].join(' + '), `${who}: needs ${q.label}`);
          Object.entries(RR.solvents || {}).filter(([, v]) => inRange(v.bp, lo, hi)).sort((a, b) => a[1].bp - b[1].bp)
            .forEach(([k, v]) => offer([`reflux in ${k}`, ...named].join(', '), [`reflux in ${k} (${fmtT(v.bp)})`, ...named].join(' + '), `${who}: needs ${q.label}`));
        }
      }
    }
    if (subs.length && !need.size && !res.free.length) res.why = 'No rule matches these molecules' + (reagent ? ' with this reagent.' : ' yet: add the reagent above the arrow first.');
    res.matched = [...need].map(([cond, o]) => ({ cond, label: o.label, rules: o.rules }));
    res.other = Object.keys(RR.conditions || {}).filter(c => !need.has(c));
    res.solvents = Object.keys(RR.solvents || {});
    // suggestions from fragments (sheet Condition_Suggestions, plus what Fragment_Conditions requires)
    if (subs.length) {
      const frags = new Set(subs.concat(reagent ? [reagent] : []).flatMap(x => fragmentsOf(R, x)));
      const seen = new Set();
      res.suggest = [];
      const push = (text, why) => { const k = text.toLowerCase(); if (!seen.has(k)) { seen.add(k); res.suggest.push({ text, why }); } };
      for (const fc of RR.fragmentConditions || [])
        if (frags.has(fc.frag) && (fc.requires || []).length)
          fc.requires.forEach(alt => push(alt.map(q => q.cond || q.solvent || (q.tmax != null && q.tmin == null ? `${Math.round(q.tmax)} °C` : q.label)).join(', '),
                                          `${fragName(fc.frag)} requires it (Fragment_Conditions)`));
      for (const cs of RR.conditionSuggestions || []) if (frags.has(cs.frag)) push(cs.suggest, `${fragName(cs.frag)}: ${cs.why}`);
    }
    return res;
  }
  /* The reaction's identifying code: "starting materials~reagent~conditions~products". The molecule parts are RDKit
     canonical SMILES in the SMILES-line order, so the same reaction always gives the same text however its molecules
     were typed; the conditions are the conditions box's comma list. ~ never occurs in a SMILES, so it separates the parts. */
  /* Conditions are written in one standard form so the same conditions always give the same code: each item becomes
     its condition's name from the Conditions sheet (refluxing -> reflux, Δ or hot -> heat, rt -> room temperature), a
     solvent its name from the Solvents sheet (EtOH -> ethanol, CH2Cl2 -> dichloromethane); a temperature becomes "N C" (no degree sign)
     (whole degrees; K and °F converted); anything else is kept as typed.
     Temperatures come first, then the words A-Z. */
  function condCode(t) {
    const temps = new Set(), words = new Map();
    const addWord = w => { if (!words.has(w.toLowerCase())) words.set(w.toLowerCase(), w); };
    for (const part of splitConds(String(t || '').replace(/~/g, '')).map(x => x.trim().replace(/\s+/g, ' ')).filter(Boolean)) {
      const r = readItem(part);
      if (r.equiv != null) addWord(`${+r.equiv} equiv`);
      if (!r.temps.length && !r.keys.length && !r.solvents.length) { if (r.equiv == null) addWord(part); continue; }   // unknown: kept as typed
      r.temps.forEach(v => temps.add(Math.round(v)));
      r.keys.forEach(addWord); r.solvents.forEach(addWord);
      if (r.rest) addWord(r.rest);                       // anything left over (e.g. "2 h")
    }
    return [...[...temps].sort((a, b) => a - b).map(v => v + ' C'),   // no degree sign in codes (plain text, easy to type)
            ...[...words.keys()].sort().map(k => words.get(k))].join(', ');
  }
  function reactionCode(R, subs, reagent, products, conditions) {
    const canonList = list => list.map(smi => { const m = R.get_mol(smi); const c = m ? m.get_smiles() : smi; if (m) m.delete(); return c; })
      .map(smi => { let g = null; try { g = window.Chem.parseSmiles(smi); } catch (e) { g = null; } return { smiles: smi, graph: g }; })
      .filter(x => x.graph);
    return [orderedSmiles(canonList(subs)), reagent ? orderedSmiles(canonList([reagent])) : '', condCode(conditions), orderedSmiles(canonList(products))].join('~');
  }
  /* Fragments (sheet "Fragments" in molecule_data.xlsx): each starting material is read as the fragments its
     SMARTS patterns find in it. A rule is only tried when every fragment it consumes ('Fragments consumed' in
     reaction_rules.xlsx) is present; its product must then contain the fragment(s) it forms. */
  const fragQ = {};
  function fragmentsOf(R, smiles) {
    const F = (window.REACTION_RULES || {}).fragments || {}, found = [];
    let m = null;
    try { m = R.get_mol(smiles); } catch (e) { m = null; }
    if (!m) return found;
    try {
      for (const [id, f] of Object.entries(F)) {
        if (!(id in fragQ)) { try { fragQ[id] = R.get_qmol(f.smarts); } catch (e) { fragQ[id] = null; } }
        const q = fragQ[id];
        if (q && Object.keys(JSON.parse(m.get_substruct_match(q) || '{}')).length) found.push(id);
      }
    } finally { m.delete(); }
    return found;
  }
  const fragName = id => { const f = ((window.REACTION_RULES || {}).fragments || {})[id]; return f ? `${id} ${f.name}` : id; };
  const itemsText = items => items.map(alt => alt.map(fragName).join(' or ')).join(' + ');
  /* Conditions that come from fragments (sheet "Fragment_Conditions"): for every fragment found in the starting
     materials or the reagent, its 'requires' must be met, nothing it 'forbids' may be in the conditions box, and no
     'forbidden fragment' may be in any molecule present. A rule listed under 'except' uses the fragment up and is exempt. */
  function fragmentConditionCheck(rule, present, P) {
    const problems = [], passed = [];
    const fcs = (window.REACTION_RULES || {}).fragmentConditions || [];
    for (const fc of fcs) {
      if (!present.has(fc.frag) || (fc.except || []).includes(rule.id)) continue;
      const bad = [];
      if ((fc.requires || []).length && !fc.requires.some(alt => alt.every(q => requirementMet(q, P))))
        bad.push(`needs ${fc.requires.map(alt => alt.map(reqLabel).join(' and ')).join(', or ')}`);
      const hit = (fc.forbids || []).filter(q => q.label ? requirementMet(q, P) : q.solvent ? P.solvents.includes(q.solvent) : P.keys.has(q.cond));
      if (hit.length) bad.push(`can't be used with ${hit.map(reqLabel).join(', ')}`);
      const ff = (fc.forbidFrags || []).filter(id => present.has(id));
      if (ff.length) bad.push(`can't be used with ${ff.map(fragName).join(', ')} present`);
      const what = fragName(fc.frag);
      if (bad.length) problems.push(`${what} ${bad.join('; ')} (${fc.reason})`);
      else passed.push(what + ': ' + [fc.requires && fc.requires.length ? 'has ' + fc.requires.map(alt => alt.map(reqLabel).join(' + ')).join(' or ') : '',
        fc.forbids && fc.forbids.length ? 'no ' + fc.forbids.map(reqLabel).join(', ') : '',
        fc.forbidFrags && fc.forbidFrags.length ? 'no ' + fc.forbidFrags.join(', ') : ''].filter(Boolean).join('; ') + ' ✓');
    }
    return { problems, passed };
  }
  /* Acid/base by pKa ('Acid/base needed' in reaction_rules.xlsx, sheet Acids_Bases): a rule that says
     'base: F_CFJ46WPI, ΔpKa ≥ -5' takes any base from the sheet whose conjugate acid pKa minus F_CFJ46WPI's pKa is at least -5;
     'acid: F_WH4T6PP5, ΔpKa ≥ -5' takes any acid where the pKa of the protonated fragment minus the acid's pKa is at least -5.
     K ≈ 10^ΔpKa for that proton transfer. */
  const abKeys = new Map();
  function acidBaseEntry(R, role, key) {
    if (!key) return null;
    for (const e of ((window.REACTION_RULES || {}).acidsBases || [])) {
      if (e.role !== role) continue;
      if (!abKeys.has(e.smiles)) abKeys.set(e.smiles, keyOf(R, e.smiles));
      if (abKeys.get(e.smiles) === key) return e;
    }
    return null;
  }
  /* ---- pKa of each fragment IN ITS MOLECULE (tools/pka.py computes the same) ----
     Every acidic and basic SITE of a molecule (one atom: the H-bearing O/N/S, else C, of an acidic fragment with a
     typical pKa <= 20; the N of a basic fragment with a typical pKaH >= 0) gets, in this order:
       1. the molecule's measured value (data/pka.js: IUPAC Digitized pKa Dataset, water, 20-30 °C), matched to the
          sites by the assignment closest to the expected values (so an amino acid's zwitterion values land right);
       2. else a Hammett estimate (phenol, benzoic acid, anilinium; meta/para substituents on a plain benzene ring):
          pKa = pKa0 - rho * sum(sigma); 4-nitrophenol 9.99 - 2.23 * 1.27 = 7.16 (measured 7.15);
       3. else the fragment's typical value (Fragments sheet), saying why there is no estimate. */
  const SIGMA = [   // [name, SMARTS (first atom = the one bonded to the ring), sigma_m, sigma_p, sigma_p-]  Hansch, Leo, Taft, Chem. Rev. 1991
    ['nitro', '[N+:1](=O)[O-]', 0.71, 0.78, 1.27], ['cyano', '[C:1]#N', 0.56, 0.66, 1.00], ['trifluoromethyl', '[C:1](F)(F)F', 0.43, 0.54, 0.65],
    ['methylsulfonyl', '[S:1](=O)(=O)[#6]', 0.60, 0.72, 1.13], ['sulfonamide', '[S:1](=O)(=O)[NX3]', 0.53, 0.60, 0.94],
    ['formyl', '[CH1:1]=O', 0.35, 0.42, 1.03], ['acetyl (ketone)', '[C:1](=O)[#6]', 0.38, 0.50, 0.84], ['carboxylic acid', '[C:1](=O)[OH]', 0.37, 0.45, 0.77],
    ['ester', '[C:1](=O)O[#6]', 0.37, 0.45, 0.75], ['amide', '[C:1](=O)[NX3]', 0.28, 0.36, 0.61], ['trifluoromethoxy', '[OX2;+0:1]C(F)(F)F', 0.38, 0.35, 0.35],
    ['fluoro', '[F:1]', 0.34, 0.06, -0.03], ['chloro', '[Cl:1]', 0.37, 0.23, 0.19], ['bromo', '[Br:1]', 0.39, 0.23, 0.25], ['iodo', '[I:1]', 0.35, 0.18, 0.27],
    ['trifluoromethylthio', '[SX2;+0:1]C(F)(F)F', 0.40, 0.50, 0.57], ['methylthio', '[SX2;+0:1][CX4;!$(C~[!#6;!#1])]', 0.15, 0.00, 0.06],
    ['acetamido', '[NX3;H1;+0:1]C(=O)', 0.21, 0.00, -0.46], ['dimethylamino', '[NX3;+0:1]([CX4])[CX4]', -0.16, -0.83, -0.12],
    ['amino', '[NX3;H2;+0:1]', -0.16, -0.66, -0.15], ['hydroxy', '[OX2;H1;+0:1]', 0.12, -0.37, -0.37], ['phenoxy', '[OX2;+0:1]c', 0.25, -0.03, -0.10],
    ['alkoxy', '[OX2;+0:1][CX4]', 0.12, -0.27, -0.26], ['tert-butyl', '[CX4;!$(C~[!#6;!#1]):1](C)(C)C', -0.10, -0.20, -0.13],
    ['methyl', '[CH3;!$(C~[!#6;!#1]):1]', -0.07, -0.17, -0.17], ['alkyl', '[CX4;!$(C~[!#6;!#1]):1]', -0.07, -0.15, -0.15],
    ['phenyl', '[c:1]', 0.06, -0.01, 0.02], ['vinyl', '[CX3:1]=[CX3]', 0.06, -0.04, -0.08]];
  const HAMMETT = {
    F_EBFEFTNC: { kind: 'acid', smarts: '[c:1][OX2H1]', pKa0: 9.99, rho: 2.23, para: 'sigma-', name: 'phenol' },
    F_WH4T6PP5: { kind: 'acid', smarts: '[c:1][CX3](=O)[OX2H1]', pKa0: 4.20, rho: 1.00, para: 'sigma', name: 'benzoic acid' },
    F_UY2ANYI2: { kind: 'base', smarts: '[c:1][NX3;H2;+0]', pKa0: 4.60, rho: 2.89, para: 'sigma-', name: 'anilinium' } };
  const pkaQ = {}, pkaCache = new Map();
  const qmolOf = (R, sm) => { if (!(sm in pkaQ)) { try { pkaQ[sm] = R.get_qmol(sm); } catch (e) { pkaQ[sm] = null; } } return pkaQ[sm]; };
  const EL = { 1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I' };
  function moleculePkas(R, smiles) {
    if (pkaCache.has(smiles)) return pkaCache.get(smiles);
    let m = null; try { m = R.get_mol(smiles); } catch (e) { m = null; }
    let out = { acid: [], base: [] };
    if (m) { try { out = pkaOfMol(R, m, smiles); } finally { m.delete(); } }
    pkaCache.set(smiles, out);
    return out;
  }
  /* the sites of an RDKit molecule, atom numbers in ITS order (the Molecule page labels its atoms with them); key =
     its canonical SMILES, for the measured values */
  function pkaOfMol(R, m, key) {
    const out = { acid: [], base: [] };
    {
      const J = JSON.parse(m.get_json()).molecules[0], ext = (J.extensions || []).find(x => x.name === 'rdkitRepresentation') || {};
      const atoms = J.atoms.map(a => ({ el: EL[a.z == null ? 6 : a.z] || '?', h: a.impHs || 0, chg: a.chg || 0, nb: [] }));
      (J.bonds || []).forEach(b => { atoms[b.atoms[0]].nb.push(b.atoms[1]); atoms[b.atoms[1]].nb.push(b.atoms[0]);
        if ((b.bo == null ? 1 : b.bo) !== 1) { atoms[b.atoms[0]].multi = true; atoms[b.atoms[1]].multi = true; } });
      (J.atoms || []).forEach((a, i) => { atoms[i].nb.forEach(j => { if (atoms[j].el === 'H') atoms[i].h++; }); });
      const aromatic = new Set(ext.aromaticAtoms || []), rings = ext.atomRings || [];
      const ringCount = new Map(); rings.forEach(r => r.forEach(i => ringCount.set(i, (ringCount.get(i) || 0) + 1)));
      const F = (window.REACTION_RULES || {}).fragments || {};
      const sites = { acid: new Map(), base: new Map() };
      for (const [fid, f] of Object.entries(F)) {
        if (!(f.pKa != null && f.pKa <= 20) && !(f.pKaH != null && f.pKaH >= 0)) continue;
        const q = qmolOf(R, f.smarts); if (!q) continue;
        for (const mt of matchesOf(m, q)) {
          if (f.pKa != null && f.pKa <= 20) {
            const het = mt.atoms.filter(i => ['O', 'N', 'S'].includes(atoms[i].el) && atoms[i].h > 0), car = mt.atoms.filter(i => atoms[i].el === 'C' && atoms[i].h > 0 && !atoms[i].multi && !aromatic.has(i));   // sp3 C: the alpha C, not a formyl C-H
            const site = het.length ? het[0] : car.length ? car[0] : null;
            if (site != null) { const s_ = sites.acid.get(site) || { atom: site, fragments: [], typical: f.pKa }; s_.fragments.push(fid); s_.typical = Math.min(s_.typical, f.pKa); sites.acid.set(site, s_); }
          }
          if (f.pKaH != null && f.pKaH >= 0) {
            const ns = mt.atoms.filter(i => atoms[i].el === 'N' && atoms[i].chg === 0);
            if (ns.length) { const s_ = sites.base.get(ns[0]) || { atom: ns[0], fragments: [], typical: f.pKaH }; s_.fragments.push(fid); s_.typical = Math.max(s_.typical, f.pKaH); sites.base.set(ns[0], s_); }
          }
        }
      }
      const hammett = fid => {                           // [{value|null, text, siteAtom}] per occurrence of the fragment's site
        const h = HAMMETT[fid], res = [];
        for (const mt of matchesOf(m, qmolOf(R, h.smarts))) {
          const site = mt.atoms[0], own = new Set(mt.atoms);
          const siteAtom = mt.atoms.slice(1).find(i => (h.kind === 'acid' && atoms[i].el === 'O' && atoms[i].h > 0) || (h.kind === 'base' && atoms[i].el === 'N'));
          const rs = rings.filter(r => r.includes(site) && r.length === 6);
          if (rs.length !== 1 || rs[0].some(i => ringCount.get(i) !== 1)) { res.push({ value: null, text: 'fused ring: no estimate', siteAtom }); continue; }
          const ring = rs[0];
          if (ring.some(i => !aromatic.has(i) || atoms[i].el !== 'C')) { res.push({ value: null, text: 'not a benzene ring: no estimate', siteAtom }); continue; }
          const pos = i => { const d = Math.abs(ring.indexOf(i) - ring.indexOf(site)); return Math.min(d, 6 - d); };
          let total = 0, ok = true, why = ''; const parts = [];
          for (const i of ring) {
            if (i === site) continue;
            for (const j of atoms[i].nb) {
              if (ring.includes(j) || own.has(j) || atoms[j].el === 'H') continue;
              if (pos(i) === 1) { ok = false; why = 'ortho substituent'; break; }
              let kind = null;
              for (const [name, sm, sMeta, sPara, sMinus] of SIGMA) {
                const q = qmolOf(R, sm);
                if (q && matchesOf(m, q).some(x => x.atoms[0] === j)) {
                  const paraMinus = h.para === 'sigma-' && sMinus > sPara;     // only resonance-withdrawing para groups
                  kind = [name, pos(i) === 2 ? sMeta : (paraMinus ? sMinus : sPara), pos(i) === 3 && paraMinus]; break;
                }
              }
              if (!kind) { ok = false; why = 'a substituent with no sigma value'; break; }
              total += kind[1];
              parts.push(`${pos(i) === 2 ? 'meta' : 'para'}-${kind[0]} (${pos(i) === 2 ? 'σm' : kind[2] ? 'σp⁻' : 'σp'} ${kind[1] >= 0 ? '+' : ''}${kind[1].toFixed(2)})`);
            }
            if (!ok) break;
          }
          if (!ok) { res.push({ value: null, text: `${why}: no estimate`, siteAtom }); continue; }
          const v = Math.round((h.pKa0 - h.rho * total) * 100) / 100;
          res.push({ value: v, text: `Hammett: ${h.name} ${h.pKa0} − ${h.rho} × (${parts.length ? parts.join(' + ') : 'no substituent'}) = ${v}`, siteAtom });
        }
        return res;
      };
      /* symmetry classes (tools/pka.py symmetry_classes): element, heavy degree, H count, charge, aromatic, refined by
         the neighbours' classes; equivalent sites (acetone's two CH3) become one site with all its atoms */
      const cls = (() => {
        const deg = atoms.map(a => a.nb.filter(j => atoms[j].el !== 'H').length);
        const renum = labs => { const u = [...new Set(labs)].sort(); const o = new Map(u.map((l, k) => [l, k])); return labs.map(l => o.get(l)); };
        let c = renum(atoms.map((a, i) => `${a.el},${deg[i]},${a.h},${a.chg},${aromatic.has(i) ? 1 : 0}`));
        for (let it = 0; it < atoms.length; it++) {
          const nx = renum(atoms.map((a, i) => `${c[i]}|` + a.nb.filter(j => atoms[j].el !== 'H').map(j => c[j]).sort((x, y) => x - y).join(',')));
          if (new Set(nx).size === new Set(c).size) break;
          c = nx;
        }
        return c;
      })();
      for (const kind of ['acid', 'base']) {
        const merged = new Map();
        for (const a of [...sites[kind].keys()].sort((x, y) => x - y)) {
          const s_ = sites[kind].get(a), g = merged.get(cls[a]);
          if (g) { g.atoms.push(a); g.fragments.push(...s_.fragments); g.typical = kind === 'acid' ? Math.min(g.typical, s_.typical) : Math.max(g.typical, s_.typical); }
          else merged.set(cls[a], Object.assign({}, s_, { atoms: [a], fragments: s_.fragments.slice() }));
        }
        sites[kind] = new Map([...merged.values()].map(g => [g.atom, g]));
      }
      for (const kind of ['acid', 'base']) {
        for (const s_ of sites[kind].values()) {
          s_.pKa = s_.typical; s_.source = 'typical value for ' + [...new Set(s_.fragments)].sort().join(', ');
          for (const fid of s_.fragments) {
            if (!HAMMETT[fid] || HAMMETT[fid].kind !== kind) continue;
            for (const e of hammett(fid)) {
              if (!s_.atoms.includes(e.siteAtom)) continue;
              if (e.value != null) { s_.pKa = e.value; s_.source = 'estimated: ' + e.text; } else s_.source += ` (${e.text})`;
            }
          }
        }
        out[kind] = [...sites[kind].values()].sort((a, b) => kind === 'base' ? b.pKa - a.pKa : a.pKa - b.pKa);
      }
      // measured values -> sites, the assignment closest to the expected values (acid pKa and pKaH on one scale)
      const M = ((window.PKA_MEASURED || {}).molecules || {})[key];
      // an aldehyde/ketone measured in water is partly hydrated: with no O-H/N-H/S-H of its own, a measured acid value
      // below 15 is the hydrate's O-H (acetaldehyde 13.57); it goes on the carbonyl O (tools/pka.py does the same)
      const measAcid = M ? (M.acid || []).slice() : [];
      let hydrate = null;
      const co = matchesOf(m, qmolOf(R, '[CX3;$([CH1](=O)[#6]),$([CH2]=O),$(C(=O)([#6])[#6])]=[OX1]'));
      if (co.length && measAcid.length && Math.min(...measAcid) < 15 && !out.acid.some(x => ['O', 'N', 'S'].includes(atoms[x.atom].el))) {
        const v = Math.min(...measAcid); measAcid.splice(measAcid.indexOf(v), 1);
        hydrate = { atom: co[0].atoms[1], atoms: [co[0].atoms[1]], fragments: [atoms[co[0].atoms[0]].h > 0 ? 'F_ADOAMFFM' : 'F_CFJ46WPI'], typical: v, pKa: v,
          source: 'measured (IUPAC Digitized pKa Dataset) in water, where the C=O is partly hydrated: the O–H of the hydrate R2C(OH)2' };
      }
      // a measured pKaH below 0 is the protonation of a very weak base (C=O, ether: acetone -7.2), not a basic site here
      const groups = out.acid.concat(out.base), meas = measAcid.concat(M ? (M.base || []).filter(v => v >= 0) : []);
      // one slot per atom: a group of equivalent sites (succinic acid's two COOH) takes as many measured values as atoms
      const all = groups.flatMap(g => g.atoms.map(() => g));
      if (meas.length && meas.length <= all.length) {
        let pairs = null;
        if (all.length <= 8) {
          let best = null;
          const perm = (cur, used) => {
            if (cur.length === meas.length) { const c = cur.reduce((t, j, k) => t + Math.abs(all[j].pKa - meas[k]), 0); if (best === null || c < best.c) best = { c, p: cur.slice() }; return; }
            for (let j = 0; j < all.length; j++) if (!used[j]) { used[j] = true; cur.push(j); perm(cur, used); cur.pop(); used[j] = false; }
          };
          perm([], []);
          pairs = best.p.map((j, k) => [j, meas[k]]);
        } else {
          const order = all.map((_, j) => j).sort((a, b) => all[a].pKa - all[b].pKa), ms = meas.slice().sort((a, b) => a - b);
          pairs = ms.map((v, k) => [order[k], v]);
        }
        const got = new Map();
        for (const [j, v] of pairs) { const g = all[j]; if (!got.has(g)) got.set(g, []); got.get(g).push(v); }
        for (const [g, vals] of got) {
          vals.sort((a, b) => out.base.includes(g) ? b - a : a - b);
          const before = (g.source.startsWith('estimated') ? 'the estimate was ' : 'the typical value was ') + g.pKa;
          g.source = 'measured (IUPAC Digitized pKa Dataset)' + (Math.abs(vals[0] - g.pKa) > 0.005 ? `; ${before}` : '') + vals.slice(1).map(v => `; then ${v} for the next equivalent site`).join('');
          g.pKa = vals[0];
        }
      } else if (meas.length) groups.forEach(s_ => { s_.source += `; the molecule has ${meas.length} measured values for ${all.length} site(s), not assigned`; });
      if (hydrate) out.acid.push(hydrate);
    }
    return out;
  }
  window.MoleculePka = { forMol: (R, M) => { try { return pkaOfMol(R, M, M.get_smiles()); } catch (e) { return null; } } };   // the Molecule page
  /* ---- Steric and electronic effects of each fragment IN ITS MOLECULE (tools/effects.py computes the same) ----
     Each fragment occurrence is read at its reacting atom: the C of a C=O / C=N / C≡N; else a C carrying a halogen;
     else an N, O or S (one with an H first); else a C of a C=C / C≡C; else the first atom of the match.
       steric score = the reacting atom's heavy neighbours (not counting a doubly bonded O/S partner) + the most of its
         neighbours' own further heavy neighbours (for an aromatic neighbour only the ortho substituents count, the ring
         itself is flat): <= 2 open, 3 moderately hindered, 4 hindered, >= 5 very hindered
         (n-butylamine N 2 open; tert-butylamine N 4 hindered; diisopropylethylamine N 5 very hindered; pinacolone C=O 5).
       electronic: the Gasteiger-Marsili partial charge on the reacting atom (MolInfo); for an atom on, or bonded to, a
         benzene ring, the meta/para substituents' Hammett sigma (same table as the pKa estimates) and their sum
         (> 0 electron-poor ring, < 0 electron-rich); whether it is conjugated; a carbon's oxidation state. */
  /* Taft parameters of a group R on the reacting atom (CH3 = 0): sigma* (polar: > 0 pulls electrons, < 0 pushes) and Es (steric:
     more negative = bulkier). Taft, in Steric Effects in Organic Chemistry (1956); Hansch, Leo, Taft, Chem. Rev. 91, 165 (1991);
     Hansch & Leo, Substituent Constants for Correlation Analysis (1979). null = no reliable value. Same table as tools/effects.py.
     [name, SMARTS (first atom = the one bonded to the reacting atom), sigma*, Es]; the first match wins. */
  const TAFT = [
    ['methyl', '[CH3X4]', 0.00, 0.00], ['ethyl', '[CH2X4;D2][CH3]', -0.10, -0.07], ['n-propyl', '[CH2X4;D2][CH2X4;D2][CH3]', -0.115, -0.36],
    ['n-butyl', '[CH2X4;D2][CH2X4;D2][CH2X4;D2][CH3]', -0.13, -0.39], ['isopropyl', '[CH1X4;D3]([CH3])[CH3]', -0.19, -0.47],
    ['isobutyl', '[CH2X4;D2][CH1X4;D3]([CH3])[CH3]', -0.125, -0.93], ['sec-butyl', '[CH1X4;D3]([CH3])[CH2X4;D2][CH3]', -0.21, -1.13],
    ['tert-butyl', '[CX4;D4]([CH3])([CH3])[CH3]', -0.30, -1.54], ['neopentyl', '[CH2X4;D2][CX4;D4]([CH3])([CH3])[CH3]', -0.165, -1.74],
    ['cyclohexyl', '[CH1X4;D3;R1]1[CH2;D2][CH2;D2][CH2;D2][CH2;D2][CH2;D2]1', -0.15, -0.79],
    ['benzyl', '[CH2X4;D2][c;D3]1[cH][cH][cH][cH][cH]1', 0.215, -0.38],
    ['trifluoromethyl', '[CX4;D4](F)(F)F', 2.61, -1.16], ['trichloromethyl', '[CX4;D4](Cl)(Cl)Cl', 2.65, -2.06],
    ['dichloromethyl', '[CH1X4;D3](Cl)Cl', 1.94, -1.54], ['fluoromethyl', '[CH2X4;D2]F', 1.10, -0.24], ['chloromethyl', '[CH2X4;D2]Cl', 1.05, -0.24],
    ['bromomethyl', '[CH2X4;D2]Br', 1.00, -0.27], ['iodomethyl', '[CH2X4;D2]I', 0.85, -0.37],
    ['hydroxymethyl', '[CH2X4;D2][OX2H1]', 0.555, -0.03], ['methoxymethyl', '[CH2X4;D2][OX2;D2][CH3]', 0.52, -0.19],
    ['phenyl', '[c;D3]1[cH][cH][cH][cH][cH]1', 0.60, -1.01],
    ['aryl (taken as phenyl; its ring substituents are in the Hammett reading)', '[c;D3]', 0.60, -1.01],
    ['fluoro', '[F]', null, 0.78], ['chloro', '[Cl]', null, 0.27], ['bromo', '[Br]', null, 0.08], ['iodo', '[I]', null, -0.16],
    ['hydroxy', '[OX2H1]', null, 0.69], ['methoxy', '[OX2;D2][CH3]', null, 0.69],
  ];
  // a carbon group not in the table: Es estimated from its branching (sigma* left out: it depends on what the group carries)
  const TAFT_BRANCH = { 2: ['a CH2R group (Es taken as n-propyl)', -0.36], 3: ['a CHR2 group (Es taken as isopropyl)', -0.47], 4: ['a CR3 group (Es taken as tert-butyl)', -1.54] };
  function effectsOfMol(R, m, record) {
    const RR = window.REACTION_RULES || {}, F = RR.fragments || {};
    const A = record.atoms, n = record.n;
    const bondOrder = (i, j) => { const e = A[i].neighbours.find(x => x.atom === j); return e ? (e.aromatic ? 1.5 : e.order) : 0; };
    const heavy = i => A[i].neighbours.filter(x => x.atom < n && !A[x.atom].isH).map(x => x.atom);
    const ringOf = i => (record.rings || []).find(r => r.length === 6 && r.includes(i) && r.every(k => A[k].aromatic));
    const pick = atoms => {
      const isC = i => A[i].el === 'C';
      const multi = (i, els) => heavy(i).some(j => els.includes(A[j].el) && bondOrder(i, j) >= 2);
      return atoms.find(i => isC(i) && multi(i, ['O', 'N', 'S']))
        ?? atoms.find(i => isC(i) && heavy(i).some(j => ['F', 'Cl', 'Br', 'I'].includes(A[j].el)))
        ?? atoms.find(i => ['N', 'O', 'S'].includes(A[i].el) && A[i].hCount > 0)
        ?? atoms.find(i => ['N', 'O', 'S'].includes(A[i].el))
        ?? atoms.find(i => isC(i) && heavy(i).some(j => A[j].el === 'C' && bondOrder(i, j) >= 2 && bondOrder(i, j) !== 1.5))
        ?? atoms[0];
    };
    const steric = x => {
      const partner = heavy(x).filter(j => !(['O', 'S'].includes(A[j].el) && bondOrder(x, j) === 2));
      let k = 0, why = [];
      for (const j of partner) {
        let kk;
        if (A[j].aromatic) { const ring = ringOf(j); kk = ring ? ring.filter(r => r !== j && Math.min(Math.abs(ring.indexOf(r) - ring.indexOf(j)), 6 - Math.abs(ring.indexOf(r) - ring.indexOf(j))) === 1)
          .reduce((t, o) => t + heavy(o).filter(q => !ring.includes(q)).length, 0) : 0; }
        else kk = heavy(j).filter(q => q !== x).length;
        k = Math.max(k, kk);
      }
      const score = partner.length + k;
      const cls = score <= 2 ? 'open' : score === 3 ? 'moderately hindered' : score === 4 ? 'hindered' : 'very hindered';
      return { score, cls, text: `${cls} (score ${score}: ${partner.length} heavy neighbour${partner.length === 1 ? '' : 's'}; the most crowded one carries ${k} more)` };
    };
    const electronic = x => {
      const a = A[x], parts = [];
      if (a.gasteiger != null) parts.push(`partial charge ${a.gasteiger >= 0 ? '+' : ''}${a.gasteiger.toFixed(2)}`);
      const site = a.aromatic ? x : heavy(x).find(j => A[j].aromatic);
      let sigma = null;
      const ring = site != null ? ringOf(site) : null;
      if (ring && ring.every(k => A[k].el === 'C')) {
        let t = 0; const named = [];
        for (const r of ring) {
          if (r === site) continue;
          const d = Math.min(Math.abs(ring.indexOf(r) - ring.indexOf(site)), 6 - Math.abs(ring.indexOf(r) - ring.indexOf(site)));
          for (const j of heavy(r)) {
            if (ring.includes(j) || j === x) continue;
            if (d === 1) { named.push('ortho substituent'); continue; }
            const hit = SIGMA.find(([name, sm]) => { const q = qmolOf(R, sm); return q && matchesOf(m, q).some(mt => mt.atoms[0] === j); });
            if (!hit) { named.push(`${d === 2 ? 'meta' : 'para'} group with no σ value`); continue; }
            const v = d === 2 ? hit[2] : hit[3]; t += v;
            named.push(`${d === 2 ? 'meta' : 'para'}-${hit[0]} (σ ${v >= 0 ? '+' : ''}${v.toFixed(2)})`);
          }
        }
        sigma = Math.round(t * 100) / 100;
        parts.push(named.length ? `ring: ${named.join(', ')}; Σσ ${sigma >= 0 ? '+' : ''}${sigma.toFixed(2)} → ${sigma > 0.1 ? 'electron-poor ring' : sigma < -0.1 ? 'electron-rich ring' : 'about neutral ring'}` : 'ring: no other substituent');
      }
      const conj = (record.bonds || []).some(b => (b.a === x || b.b === x) && b.conjugated);
      if (conj) parts.push(a.hybridization && a.hybridization.promoted === 'amide' ? 'lone pair in an amide (delocalised onto the C=O)'
        : ['N', 'O', 'S'].includes(a.el) ? 'lone pair conjugated (delocalised: less available)' : 'conjugated');
      if (a.el === 'C' && a.oxidationState && a.oxidationState.value != null) parts.push(`oxidation state ${a.oxidationState.text}`);
      return { charge: a.gasteiger, sigma, conjugated: conj, text: parts.join('; ') };
    };
    const tmatch = new Map();                         // TAFT index -> set of first atoms, computed once per molecule
    const firstAtoms = k => { if (!tmatch.has(k)) { const q = qmolOf(R, TAFT[k][1]); tmatch.set(k, new Set(q ? matchesOf(m, q).map(mt => mt.atoms[0]) : [])); } return tmatch.get(k); };
    const fmt = v => (v >= 0 ? '+' : '') + v.toFixed(2);
    const f0 = v => fmt(v).replace('+0.00', '0.00');
    const taft = x => {                               // the groups on the reacting atom, not counting its own ring or a =O/=S partner
      const xr = (record.rings || []).filter(r => r.includes(x));
      const groups = [];
      for (const j of heavy(x)) {
        if (['O', 'S'].includes(A[j].el) && bondOrder(x, j) === 2) continue;
        if (xr.some(r => r.includes(j))) continue;
        let hit = null;
        for (let k = 0; k < TAFT.length; k++) if (firstAtoms(k).has(j)) { hit = [TAFT[k][0], TAFT[k][2], TAFT[k][3], false]; break; }
        if (!hit && A[j].el === 'C' && !A[j].aromatic && A[j].neighbours.every(e => !e.aromatic && e.order === 1) && TAFT_BRANCH[heavy(j).length]) {
          const [nm, es] = TAFT_BRANCH[heavy(j).length]; hit = [nm, null, es, true];
        }
        if (hit) groups.push({ atom: j, name: hit[0], sigmaStar: hit[1], Es: hit[2], estimated: hit[3] });
      }
      const ss = groups.filter(g => g.sigmaStar != null).map(g => g.sigmaStar), es = groups.filter(g => g.Es != null).map(g => g.Es);
      const sSum = ss.length ? Math.round(ss.reduce((a, b) => a + b, 0) * 1000) / 1000 : null, eMin = es.length ? Math.min(...es) : null;
      const st = es.length ? 'Taft Es: ' + groups.filter(g => g.Es != null).map(g => `${g.name} ${f0(g.Es)}`).join(', ') + `; bulkiest ${f0(eMin)}` : '';
      const el = ss.length ? 'Taft σ*: ' + groups.filter(g => g.sigmaStar != null).map(g => `${g.name} ${f0(g.sigmaStar)}`).join(', ') + `; Σσ* ${f0(sSum)} → `
        + (sSum >= 0.5 ? 'its groups pull electrons away' : sSum <= -0.1 ? 'its groups push electrons in' : 'about like methyl groups') : '';
      return { groups, sigmaStar: sSum, EsMin: eMin, stericText: st, electronicText: el };
    };
    const out = new Map();
    for (const [fid, f] of Object.entries(F)) {
      const q = qmolOf(R, f.smarts); if (!q) continue;
      for (const mt of matchesOf(m, q)) {
        const x = pick(mt.atoms.filter(i => i < n));
        if (x == null) continue;
        const key = fid + ':' + x;
        if (out.has(key)) continue;
        const t = taft(x), st = steric(x), el = electronic(x);
        if (t.stericText) st.text += '; ' + t.stericText;
        if (t.electronicText) el.text = (el.text ? el.text + '; ' : '') + t.electronicText;
        out.set(key, { fragment: fid, name: f.name, atom: x, steric: st, electronic: el, taft: t });
      }
    }
    return [...out.values()].sort((p, q) => p.atom - q.atom || p.fragment.localeCompare(q.fragment));
  }
  const effectsCache = new Map();
  function moleculeEffects(R, smiles) {                 // the Reactions page: per starting material, cached
    if (effectsCache.has(smiles)) return effectsCache.get(smiles);
    let m = null, out = [];
    try { m = R.get_mol(smiles); if (m && window.MolInfo) out = effectsOfMol(R, m, window.MolInfo.analyse(R, m)); } catch (e) { out = []; }
    finally { if (m) m.delete(); }
    effectsCache.set(smiles, out);
    return out;
  }
  window.MoleculeEffects = { forMol: (R, M, record) => { try { return effectsOfMol(R, M, record); } catch (e) { return null; } } };   // the Molecule page

  /* the pKa of a fragment in a given molecule ('acid': its acidic H; 'base': its protonated form), or null */
  function fragmentPkaIn(R, smiles, fid, kind) {
    const s_ = moleculePkas(R, smiles)[kind].find(x => x.fragments.includes(fid));
    return s_ ? { value: s_.pKa, source: s_.source } : null;
  }
  function acidBaseCheck(ab, e, own) {
    const F = ((window.REACTION_RULES || {}).fragments || {})[ab.frag] || {};
    // the fragment's pKa in THIS molecule when known (measured or estimated), else its typical value
    const fp = own && own.value != null ? own.value : (ab.role === 'base' ? F.pKa : F.pKaH);
    const d = Math.round((ab.role === 'base' ? e.pKa - fp : fp - e.pKa) * 10) / 10;
    const ok = d >= ab.min;
    const text = ab.role === 'base'
      ? `${e.name} (its conjugate acid, ${e.conj}, pKa ${minus(e.pKa)}) takes the H of ${fragName(ab.frag)} (pKa ${minus(fp)})`
      : `${e.name} (pKa ${minus(e.pKa)}) protonates ${fragName(ab.frag)} (its conjugate acid has pKa ${minus(fp)})`;
    const srcNote = own && own.value != null ? ` [${own.source}]` : ' [typical value for the fragment]';
    return { ok, d, text: `${text}${srcNote}: ΔpKa ${d > 0 ? '+' : ''}${minus(d)}, K ≈ 10^${minus(d)} — needs ΔpKa ≥ ${minus(ab.min)}` };
  }
  /* SN1 / SN2 / E1 / E2 for alkyl halides, decided from fragments (sheets Halide_Pathways, Nucleophiles, Solvents in
     reaction_rules.xlsx): substrate class from F_KNBZHWVH methyl / F_OSWUXBJP 1° / F_4KF4VNAH 2° / F_YUNA7OOT 3°, reagent class from the
     Nucleophiles sheet (no reagent + a protic solvent = solvolysis, the solvent is the nucleophile), solvent type and heat.
     The first matching table row gives the pathway; the app's mechanism engine (Chem.buildMechanism) builds that
     pathway's steps, arrows and product. */
  const HALIDE_CLASS = [['F_KNBZHWVH', 'methyl'], ['F_OSWUXBJP', '1°'], ['F_4KF4VNAH', '2°'], ['F_YUNA7OOT', '3°']];
  function halidePathway(R, subs, read, reagent, reagentKey, P, condText, reagentFrags) {
    const RR = window.REACTION_RULES || {}, table = RR.halidePathways || [], nucs = RR.nucleophiles || [];
    if (!table.length || !window.Chem.buildMechanism) return null;
    let sub = null, cls = null, clsFrag = null;
    for (const x of read) { const hit = HALIDE_CLASS.find(([id]) => x.frags.includes(id)); if (hit) { sub = x; clsFrag = hit[0]; cls = hit[1]; break; } }
    if (!sub) return null;
    // the nucleophile / base: the reagent if the Nucleophiles sheet knows it, else a protic solvent (solvolysis)
    let nuc = null, via = 'reagent';
    if (reagentKey) nuc = nucs.find(n => { if (!abKeys.has(n.smiles)) abKeys.set(n.smiles, keyOf(R, n.smiles)); return abKeys.get(n.smiles) === reagentKey; }) || null;
    if (!nuc && !reagent && read.length === 1) {           // solvolysis: the halide alone in a protic solvent (with a second
      const prot = P.solvents.find(k => (RR.solvents[k] || {}).type === 'protic' && nucs.some(n => n.name === k));   // starting material,
      if (prot) { nuc = nucs.find(n => n.name === prot); via = 'solvent'; }   // e.g. an amine, the table does not apply)
    }
    if (!nuc) return null;
    const types = P.solvents.map(k => (RR.solvents[k] || {}).type);
    const solvType = types.includes('polar aprotic') ? 'polar aprotic' : types.includes('protic') ? 'protic' : (types[0] || 'none');
    const heat = P.keys.has('heat');
    const row = table.find(h => (h.substrate === cls) && (h.reagent === 'any' || h.reagent === nuc.class) &&
      (h.solvent === 'any' || h.solvent === solvType) && (h.heat === 'any' || (h.heat === 'yes') === heat));
    if (!row) return null;
    const head = `Halide_Pathways (reaction_rules.xlsx): ${fragName(clsFrag)} + ${nuc.name} (${nuc.class}${via === 'solvent' ? ', as the solvent' : ''})` +
      `, ${solvType} solvent${heat ? ', heat' : ''} → `;
    const reading = { subs: read, reagent: reagent ? { smiles: reagent, frags: reagentFrags } : null,
      pathway: `${head}${row.major === 'none' ? 'no reaction' : row.major + (row.minor ? ' (major), ' + row.minor + ' (minor)' : '')}. ${row.why}` };
    if (row.major === 'none') return { reading, kind: 'none', note: `No reaction (known, Halide_Pathways sheet): ${row.why}` };
    // build the pathway with the app's engine
    let g0 = null; try { g0 = window.Chem.parseSmiles(sub.smiles); } catch (e) { g0 = null; }
    if (!g0) return null;
    const lg0 = g0.atoms.find(a => ['Cl', 'Br', 'I'].includes(a.element) && g0.bonds.some(b => (b.a === a.id || b.b === a.id)));
    const bd = lg0 && g0.bonds.find(b => b.a === lg0.id || b.b === lg0.id);
    if (!bd) return null;
    const lg = { id: lg0.id };
    let cId = bd.a === lg.id ? bd.b : bd.a;
    const rk = window.Chem.resolveReagent(nuc.smiles);
    if (!rk || !rk.key) return null;
    let mech = null, rearr = [], shiftSteps = [], firstStep = null;
    if (row.major === 'SN1' || row.major === 'E1') {
      const cr = carbocationRoute(R, sub.smiles, { smarts: '[CX4][Cl,Br,I]', cation: 1, leaving: 2, then: 'rerun' });
      if (cr) {
        // step 1 (ionisation) from the real starting material, then the shift(s), then the rest on the rearranged skeleton
        try { const m1 = window.Chem.buildMechanism(g0, cId, lg.id, rk.key, { heat, force: row.major }); firstStep = m1 && m1.steps[0]; } catch (e) { firstStep = null; }
        let g1 = null; try { g1 = window.Chem.parseSmiles(cr.smiles); } catch (e) { g1 = null; }
        const lg1 = g1 && g1.atoms.find(a => ['Cl', 'Br', 'I'].includes(a.element));
        const bd1 = lg1 && g1.bonds.find(b => b.a === lg1.id || b.b === lg1.id);
        if (bd1) { g0 = g1; rearr = cr.shifts; shiftSteps = cr.steps; cId = bd1.a === lg1.id ? bd1.b : bd1.a; lg.id = lg1.id; }
      }
    }
    try { mech = window.Chem.buildMechanism(g0, cId, lg.id, rk.key, { heat, force: row.major, forceMinor: row.minor || null }); } catch (e) { mech = null; }
    if (!mech || !mech.steps || !mech.steps.length) return { reading, note: `No reaction drawn: the mechanism engine could not build ${row.major} here.` };
    const smi = g => { try { return window.Chem.toSmiles(g); } catch (e) { return ''; } };
    const canon = x => { const m = R.get_mol(x); const c = m ? m.get_smiles() : x; if (m) m.delete(); return c; };
    const last = mech.steps[mech.steps.length - 1];
    const nucCanon = canon(nuc.smiles);
    const products = smi(last.graph).split('.').filter(Boolean).map(canon).filter(x => x !== nucCanon);
    const engineSteps = mech.steps.slice(0, -1).map((st, i) => ({
      step: i + 1, type: st.title.replace(/^(Step \d+|One step) — /, ''), note: st.note,
      graph: st.graph, gArrows: st.arrows || [], arrowLabels: (st.arrows || []).map(a => a.label).filter(Boolean),
      smiles: smi(mech.steps[i + 1].graph).split('.').filter(Boolean).map(canon).filter(x => x !== nucCanon).join('.'),
    }));
    // with a rearrangement: the original ionisation, the shift(s), then the rearranged route after its own ionisation
    const steps = shiftSteps.length
      ? [].concat(firstStep ? [{ step: 1, type: firstStep.title.replace(/^(Step \d+|One step) — /, ''), note: firstStep.note, graph: firstStep.graph, gArrows: firstStep.arrows || [],
                                arrowLabels: (firstStep.arrows || []).map(a => a.label).filter(Boolean), smiles: shiftSteps[0].smiles }] : [],
                 shiftSteps, engineSteps.slice(1).map((st, i) => Object.assign({}, st, { step: i + 2 })))
      : engineSteps;
    reading.rearr = rearr;
    const stereo = [];
    if (row.major === 'SN2' && sub.smiles.includes('@')) {
      const anion = String(nuc.smiles).split('.').find(x => x.includes('-'));
      if (anion && nuc.atom) {
        const outs = runRuleAll(R, `[C@:1]-[Cl,Br,I].[${nuc.atom}-:2]>>[C@@:1]-[${nuc.atom}+0:2]`, [sub.smiles, anion]);
        if (outs.length) outs[0].products.forEach(sp => { const k = stripStereo(R, sp); const i = products.findIndex(x => stripStereo(R, x) === k); if (i >= 0) products[i] = sp; });
        stereo.push('SN2: inversion at the stereocentre (backside attack)');
      }
    } else if ((row.major === 'SN1' || row.major === 'E1') && sub.smiles.includes('@')) {
      // the carbocation is flat: the configuration is lost (the engine copies the starting material's, which is wrong here)
      for (let i = 0; i < products.length; i++) if (/@/.test(products[i])) products[i] = stripStereo(R, products[i]);
      stereo.push(`${row.major}: the flat carbocation is attacked from both faces, so the stereocentre is racemised`);
    }
    products.forEach(x => { const t = stereoText(R, x, sub.smiles.includes('@')); if (t) stereo.push(`${x}: ${t}`); });
    reading.stereo = stereo;
    const prodRead = products.map(x => ({ smiles: x, frags: fragmentsOf(R, x) }));
    reading.prods = prodRead;
    reading.rule = { id: row.major, consumes: [[clsFrag]], forms: [] };
    reading.formsOk = true;
    return { products, reading, mechanism: steps, code: reactionCode(R, subs, reagent, products, condText),
      note: `${row.major} · alkyl halide pathway from the Halide_Pathways sheet (reaction_rules.xlsx)` + (row.minor ? ` · minor: ${row.minor}` : '') +
        (heat ? ' · with heat' : '') };
  }
  /* every distinct outcome of a rule's SMARTS on these molecules (all sites, all orders), not just the first */
  function runRuleAll(R, smarts, subs) {
    let rx = null;
    try { rx = R.get_rxn(smarts); } catch (e) { rx = null; }
    if (!rx) return [];
    const n = smarts.split('>>')[0].split('.').length, orders = [], seen = new Set(), outs = [];
    const pick = cur => { if (cur.length === n) { orders.push(cur.slice()); return; } subs.forEach((_, i) => { if (!cur.includes(i)) { cur.push(i); pick(cur); cur.pop(); } }); };
    pick([]);
    const canon = smi => { const m = R.get_mol(smi); const c = m ? m.get_smiles() : null; if (m) m.delete(); return c; };
    try {
      for (const order of orders) {
        const ml = new R.MolList(), mols = order.map(i => R.get_mol(subs[i]));
        mols.forEach(m => ml.append(m));
        try {
          const res = rx.run_reactants(ml, 50);
          for (let k = 0; k < res.size(); k++) {
            const pl = res.get(k), out = [];
            for (let b = 0; b < pl.size(); b++) { const pm = pl.at(b); try { out.push(pm.get_smiles()); } catch (e) {} pm.delete(); }
            pl.delete();
            const cs = out.map(canon);
            if (!cs.length || cs.some(x => !x)) continue;
            const key = cs.slice().sort().join('.');
            if (!seen.has(key)) { seen.add(key); outs.push({ products: cs, used: order.map(i => subs[i]) }); }
          }
          res.delete();
        } catch (e) { /* this order does not fit */ }
        ml.delete(); mols.forEach(m => m.delete());
      }
    } finally { rx.delete(); }
    return outs;
  }
  /* how many times each fragment occurs (chemoselectivity: the fragment whose count drops is the one that reacted) */
  function fragmentCounts(R, list) {
    const F = (window.REACTION_RULES || {}).fragments || {}, counts = {};
    for (const smiles of list) {
      let m = null; try { m = R.get_mol(smiles); } catch (e) { m = null; }
      if (!m) continue;
      try {
        for (const [id, f] of Object.entries(F)) {
          if (!(id in fragQ)) { try { fragQ[id] = R.get_qmol(f.smarts); } catch (e) { fragQ[id] = null; } }
          if (!fragQ[id]) continue;
          const n = matchesOf(m, fragQ[id]).length;
          if (n) counts[id] = (counts[id] || 0) + n;
        }
      } finally { m.delete(); }
    }
    return counts;
  }
  function consumedFragment(R, used, products, group) {
    const a = fragmentCounts(R, used), b = fragmentCounts(R, products);
    const dropped = Object.keys(a).filter(id => (b[id] || 0) < a[id]);
    const ranks = ((window.REACTION_RULES || {}).reactivity || {})[group] || {};
    const ranked = dropped.filter(id => ranks[id]).sort((x, y) => ranks[x].rank - ranks[y].rank);
    return ranked[0] || dropped[0] || null;
  }
  /* Directing groups (sheet Directing_Groups): for an aromatic substitution, which ring position reacts. Each outcome's
     new attachment point is found by matching the starting ring in the product; it is then ortho / meta / para to each
     group already on the ring, and each group votes with its weight (strongest activator decides; ortho 0.6 for
     crowding; a position between two groups 0.1). Friedel–Crafts fails if any group blocks it. */
  const dirQ = {};
  /* all substructure matches as a list ([{atoms, bonds}, ...]); RDKit gives '{}' when there are none */
  function matchesOf(mol, q) {
    let r = null; try { r = JSON.parse(mol.get_substruct_matches(q) || '[]'); } catch (e) { r = null; }
    return Array.isArray(r) ? r : (r && r.atoms ? [r] : []);
  }
  function aromaticSubstitution(R, rule, outcomes) {
    const RR = window.REACTION_RULES, DG = RR.directing || [];
    const sub = (outcomes[0].used || []).find(x => /c/.test(x.replace(/Cl/g, '')));
    if (!sub) return { outcomes };
    let sm = null; try { sm = R.get_mol(sub); } catch (e) { sm = null; }
    if (!sm) return { outcomes };
    const groups = [];
    let aromatic = new Set();
    try {
      for (const d of DG) {
        if (!(d.smarts in dirQ)) { try { dirQ[d.smarts] = R.get_qmol(d.smarts); } catch (e) { dirQ[d.smarts] = null; } }
        if (!dirQ[d.smarts]) continue;
        const ms = matchesOf(sm, dirQ[d.smarts]);
        ms.forEach(x => { const ring = x.atoms[0]; if (!groups.some(g => g.ring === ring && g.weight >= d.weight)) groups.push(Object.assign({ ring }, d)); });
      }
      if (!dirQ.__a) dirQ.__a = R.get_qmol('a');
      aromatic = new Set(matchesOf(sm, dirQ.__a).map(x => x.atoms[0]));
    } finally { sm.delete(); }
    // a group already claimed by a stronger one on the same ring carbon is dropped
    const gs = groups.filter(g => !groups.some(h => h !== g && h.ring === g.ring && h.weight > g.weight));
    if (rule.eas.fc) {
      const block = gs.filter(g => g.blocksFC);
      if (block.length) return { blocked: `Friedel–Crafts fails on this ring: ${block.map(g => `${g.name} (${g.activation}${g.name.startsWith('amino') ? '; the N binds AlCl3' : ''})`).join(', ')}` };
    }
    if (!gs.length) return { outcomes, info: '' };
    let sg = null; try { sg = window.Chem.parseSmiles(sub); } catch (e) { sg = null; }
    if (!sg) return { outcomes };
    const nb = i => sg.bonds.filter(b => b.a === i + 1 || b.b === i + 1).map(b => (b.a === i + 1 ? b.b : b.a) - 1);
    const ringDist = (a, b) => {                         // steps around the aromatic ring: 1 ortho, 2 meta, 3 para
      const seen = new Map([[a, 0]]), q = [a];
      while (q.length) { const x = q.shift(); if (x === b) return seen.get(x); if (seen.get(x) >= 3) continue;
        for (const y of nb(x)) if (aromatic.has(y) && !seen.has(y)) { seen.set(y, seen.get(x) + 1); q.push(y); } }
      return null;
    };
    const REL = { 1: 'ortho', 2: 'meta', 3: 'para' };
    let qm = null; try { qm = R.get_qmol(sub); } catch (e) { qm = null; }
    const scored = [];
    for (const o of outcomes) {
      const prod = o.products.find(x => { const pm = R.get_mol(x); if (!pm) return false; let ok = false; try { ok = qm && Object.keys(JSON.parse(pm.get_substruct_match(qm) || '{}')).length > 0; } finally { pm.delete(); } return ok; });
      if (!prod) continue;
      const pm = R.get_mol(prod); let all = [];
      try { all = matchesOf(pm, qm).map(x => x.atoms); } finally { pm.delete(); }
      let pg = null; try { pg = window.Chem.parseSmiles(prod); } catch (e) { pg = null; }
      if (!pg || !all.length) continue;
      const pn = i => pg.bonds.filter(b => b.a === i + 1 || b.b === i + 1).map(b => (b.a === i + 1 ? b.b : b.a) - 1);
      /* the right match keeps every starting-material atom's connections, except ONE ring atom that gained the new
         group (a looser match could put toluene's CH3 on the new acetyl carbon) */
      let at = -1;
      for (const match of all) {
        const changed = match.map((pi, i) => pn(pi).length !== nb(i).length ? i : -1).filter(i => i >= 0);
        if (changed.length === 1 && aromatic.has(changed[0]) && pn(match[changed[0]]).length === nb(changed[0]).length + 1) { at = changed[0]; break; }
      }
      if (at < 0) continue;
      let score = 0; const rels = [];
      let orthoCount = 0;
      for (const g of gs) {
        const d = ringDist(at, g.ring);
        if (!d) continue;
        if (d === 1) orthoCount++;
        rels.push(`${REL[d]} to ${g.name}`);
        if (g.directs === 'o/p' && d === 1) score += g.weight * 0.6;
        if (g.directs === 'o/p' && d === 3) score += g.weight;
        if (g.directs === 'm' && d === 2) score += g.weight;
      }
      if (orthoCount >= 2) score *= 0.1;
      scored.push({ o, score, rels });
    }
    if (qm) qm.delete();
    if (!scored.length) return { outcomes };
    scored.sort((a, b) => b.score - a.score);
    const best = scored[0], minor = scored.slice(1).filter(x => x.score > 0 && x.score >= best.score * 0.4);
    const lead = gs.slice().sort((a, b) => b.weight - a.weight)[0];
    const info = `${gs.map(g => `${g.name}: ${g.activation}, directs ${g.directs}`).join('; ')} → the new group goes ${best.rels.join(', ') || 'to the ring'} (major)` +
      (minor.length ? `; minor: ${[...new Set(minor.map(x => x.rels.join(', ')))].join(' / ')}` : '') + `. ${lead.why}`;
    return { outcomes: [best.o], info };
  }
  /* Can this rule run here? 'skip' (wrong reagent, fragments or no SMARTS match), 'short' (a condition, fragment
     condition or acid/base strength is missing: the message says which) or 'ok' with all its outcomes. */
  function evalRule(R, rule, subs, reagentKey, P, presentAll) {
    const RR = window.REACTION_RULES;
    let abEntry = null;
    if (rule.acidBase) { abEntry = acidBaseEntry(R, rule.acidBase.role, reagentKey); if (!abEntry) return { status: 'skip' }; }
    else if (rule.reagents.length && !rule.reagents.includes(reagentKey)) return { status: 'skip' };
    const present = new Set(subs.flatMap(x => fragmentsOf(R, x)));
    if ((rule.consumes || []).length && !rule.consumes.every(alt => alt.some(id => present.has(id)))) return { status: 'skip' };
    const outcomes = runRuleAll(R, rule.smarts, subs);
    if (!outcomes.length) return { status: 'skip' };
    const alts = rule.required || [];
    if (alts.length && !alts.some(alt => alt.every(q => requirementMet(q, P)))) {
      const noTemp = alts.some(alt => alt.some(q => q.label)) && !P.temps.length;
      return { status: 'short', msg: `${rule.id} · ${rule.name} needs ${alts.map(alt => alt.map(reqLabel).join(' and ')).join(', or ')}` +
        (noTemp ? ' (no temperature is known: give one, e.g. 80 °C, or reflux with a solvent, e.g. reflux in ethanol)'
          : P.temps.length ? ` (these conditions are at ${P.refluxT != null ? `${fmtT(P.refluxT)}: reflux in ${P.solvents.join(' / ')}` : P.temps.map(fmtT).join(', ')})` : '') +
        ' (choose it in the conditions box under the arrow)' };
    }
    const fcheck = fragmentConditionCheck(rule, presentAll, P);
    if (fcheck.problems.length) return { status: 'short', msg: `${rule.id} · ${rule.name}: ${fcheck.problems.join('; ')}` };
    if (P.keys.has('catalytic') && typeof rule.equiv === 'number')
      return { status: 'short', msg: `${rule.id} · ${rule.name}: the reagent is used up (${rule.equiv} equiv per reaction), so a catalytic amount runs out after a few turnovers; give an amount (e.g. ${rule.equiv} equiv) or "excess"` };
    let own = null;                                       // the fragment's pKa in the molecule that reacts
    if (abEntry) for (const x of (outcomes[0] && outcomes[0].used) || subs) { own = fragmentPkaIn(R, x, rule.acidBase.frag, rule.acidBase.role === 'base' ? 'acid' : 'base'); if (own) break; }
    const abCheck = abEntry ? acidBaseCheck(rule.acidBase, abEntry, own) : null;
    if (abCheck && !abCheck.ok) return { status: 'short', msg: `${rule.id} · ${rule.name}: ${abEntry.name} is too weak ${rule.acidBase.role === 'acid' ? 'an acid' : 'a base'} (${abCheck.text})` };
    let easInfo = '';
    if (rule.eas) {
      const ea = aromaticSubstitution(R, rule, outcomes);
      if (ea.blocked) return { status: 'short', msg: `${rule.id} · ${rule.name}: ${ea.blocked}` };
      if (ea.outcomes && ea.outcomes !== outcomes) { const keep = ea.outcomes.slice(); outcomes.length = 0; keep.forEach(x => outcomes.push(x)); }
      easInfo = ea.info || '';
    }
    return { status: 'ok', rule, outcomes, fcheck, abCheck, easInfo };
  }
  /* Chemoselectivity (Reactivity_Order sheet): among the candidates, the first rule (sheet order) sets the group; every
     outcome of every rule in that group is ranked by the fragment it consumes, and the most reactive one wins. */
  function chooseCandidate(R, cands) {
    const first = cands[0], group = first.rule.group;
    if (!group) return { cand: first, outcome: first.outcomes[0], group: '', frag: null, others: [] };
    const ranks = (window.REACTION_RULES.reactivity || {})[group] || {};
    const options = [];
    cands.filter(c => c.rule.group === group).forEach((c, ci) => c.outcomes.forEach((o, oi) => {
      const frag = consumedFragment(R, o.used, o.products, group);
      options.push({ cand: c, outcome: o, frag, rank: frag && ranks[frag] ? ranks[frag].rank : 50, ci, oi });
    }));
    options.sort((a, b) => a.rank - b.rank || a.ci - b.ci || a.oi - b.oi);
    const best = options[0];
    const others = [...new Set(options.filter(o => o.rank > best.rank && o.frag).map(o => o.frag))];
    return { cand: best.cand, outcome: best.outcome, group, frag: best.frag, others };
  }
  /* ---- Stereochemistry ----
     A rule's 'Stereo outcome' says what happens; its stereo templates build it. Anti/syn additions across a
     1,2-disubstituted C=C use the cis template for a ring or Z alkene and the trans template for an E alkene (an acyclic
     C=C with no / or \ has no defined geometry: no stereo is drawn). The chiral-starting-material template is used when
     the starting material has a stereocentre (SN2 inversion). Every product is then labelled from RDKit's CIP tags. */
  const stQ = {};
  const qOf = (R, sm) => { if (!(sm in stQ)) { try { stQ[sm] = R.get_qmol(sm); } catch (e) { stQ[sm] = null; } } return stQ[sm]; };
  function stripStereo(R, smi) {
    const m = R.get_mol(String(smi).replace(/@+/g, '').replace(/[\/\\]/g, ''));
    const c = m ? m.get_smiles() : smi; if (m) m.delete(); return c;
  }
  function stereoTags(R, smi) {
    const m = R.get_mol(smi); if (!m) return { CIP_atoms: [], CIP_bonds: [] };
    try { return JSON.parse(m.get_stereo_tags() || '{}'); } catch (e) { return { CIP_atoms: [], CIP_bonds: [] }; } finally { m.delete(); }
  }
  // geometry of the one 1,2-disubstituted C=C (or ring-fused epoxide) in a molecule: 'cis', 'trans' or null
  function siteGeometry(R, smi, epoxide) {
    const m = R.get_mol(smi); if (!m) return null;
    try {
      if (epoxide) {
        const ring = matchesOf(m, qOf(R, '[#6][CH1;R2]1O[CH1;R2]1[#6]'));
        return ring.length ? 'cis' : null;
      }
      const all = matchesOf(m, qOf(R, '[#6][CH1]=[CH1][#6]'));
      const sites = [...new Set(all.map(x => [x.atoms[1], x.atoms[2]].sort((a, b) => a - b).join('-')))];
      if (sites.length !== 1) return null;
      if (matchesOf(m, qOf(R, '[CH1;R]=[CH1;R]')).length) return 'cis';           // a C=C in a ring is cis
      const [a, b] = sites[0].split('-').map(Number);
      let tags = {}; try { tags = JSON.parse(m.get_stereo_tags() || '{}'); } catch (e) { tags = {}; }
      const bt = (tags.CIP_bonds || []).find(x => (x[0] === a && x[1] === b) || (x[0] === b && x[1] === a));
      return bt ? (bt[2] === '(Z)' ? 'cis' : 'trans') : null;
    } finally { m.delete(); }
  }
  function applyStereo(R, rule, outcome) {
    const st = rule.stereo || {};
    let tmpl = '', why = '', sub = null;
    if (st.chiral && outcome.used.some(x => x.includes('@'))) {
      tmpl = st.chiral; sub = outcome.used.find(x => x.includes('@')); why = `${st.outcome} at the stereocentre`;
    } else if (st.cis || st.trans) {
      const epo = /\]1\[O/.test(st.cis);
      sub = outcome.used.find(x => siteGeometry(R, x, epo));
      if (!sub) return st.cis ? { info: `${st.outcome}; the C=C geometry is not defined here (write it with / or \\, or use a ring), so no stereo is drawn` } : null;
      const geo = siteGeometry(R, sub, epo);
      tmpl = geo === 'cis' ? st.cis : st.trans;
      if (!tmpl) return null;
      why = `${st.outcome} to a ${geo === 'cis' ? (epo ? 'cis (ring-fused) epoxide' : 'cis (Z or ring) C=C') : 'trans (E) C=C'}`;
    } else return null;
    const outs = runRuleAll(R, tmpl, [sub]);
    if (!outs.length) return null;
    const prods = outcome.products.slice();
    outs[0].products.forEach(sp => { const k = stripStereo(R, sp); const i = prods.findIndex(x => stripStereo(R, x) === k); if (i >= 0) prods[i] = sp; });
    return { products: prods, info: why };
  }
  function stereoText(R, smi, chiralStart) {
    const t = stereoTags(R, smi), atoms = t.CIP_atoms || [], bonds = t.CIP_bonds || [];
    const set = atoms.filter(a => a[1] !== '(?)'), unset = atoms.filter(a => a[1] === '(?)');
    const out = [];
    if (set.length) {
      const lab = '(' + set.map(a => a[1].replace(/[()]/g, '')).join(',') + ')';
      const mir = String(smi).replace(/@@/g, '!').replace(/@/g, '@@').replace(/!/g, '@');
      const meso = set.length > 1 && stripStereo(R, smi) && (() => { const a = R.get_mol(mir), b = R.get_mol(smi); const r = a && b && a.get_smiles() === b.get_smiles(); if (a) a.delete(); if (b) b.delete(); return r; })();
      out.push(meso ? `meso ${lab}: achiral, it is its own mirror image`
        : chiralStart ? `single enantiomer ${lab}`
        : `racemic: ${lab} is shown and its mirror image forms equally (the starting materials are achiral)`);
    }
    if (unset.length) out.push(`${unset.length} new stereocentre${unset.length > 1 ? 's' : ''} formed with no control: racemic (both configurations, 50:50)`);
    bonds.forEach(b => out.push(`${b[2]} C=C`));
    return out.join('; ');
  }
  /* ---- Carbocation rearrangements (sheet Rearrangements) ----
     Form the cation (protonate the C=C, or let the leaving group go), then look at each sp3 carbon next door: a 1,2-hydride
     shift, a 1,2-methyl/alkyl shift or a ring expansion (strained 4-membered ring) is made only if the new cation scores
     higher (stability + bonus). Repeat until nothing helps, then finish: add the nucleophile, eliminate (Zaitsev), or put
     the leaving group back on the new carbon so the reaction can be re-run on the rearranged skeleton. Returns null when
     no shift happens. */
  function carbocationRoute(R, smiles, spec) {
    const RA = (window.REACTION_RULES || {}).rearrangements || {}, ST = RA.stability || {}, BO = RA.bonus || {}, PR = RA.priority || {};
    const q = qOf(R, spec.smarts); if (!q) return null;
    const m = R.get_mol(smiles); if (!m) return null;
    let hit = null; try { hit = matchesOf(m, q)[0]; } finally { m.delete(); }
    if (!hit) return null;
    let g = null; try { g = window.Chem.parseSmiles(smiles); } catch (e) { g = null; }
    if (!g) return null;
    const A = i => g.atoms[i];                            // RDKit and the app both number atoms in SMILES order
    const cat0 = A(hit.atoms[spec.cation - 1]);
    if (!cat0 || cat0.element !== 'C') return null;
    const bondsOf = a => g.bonds.filter(b => b.a === a.id || b.b === a.id);
    const atomById = id => g.atoms.find(x => x.id === id);
    const other = (b, a) => atomById(b.a === a.id ? b.b : b.a);
    const bondBetweenAtoms = (a, b) => g.bonds.find(x => (x.a === a.id && x.b === b.id) || (x.a === b.id && x.b === a.id));
    const used = a => bondsOf(a).reduce((t, b) => t + b.order, 0);
    const hOf = a => a.element === 'C' ? Math.max(0, (a.charge > 0 ? 3 : 4) - used(a)) : 0;
    const carbonNb = a => bondsOf(a).map(b => other(b, a)).filter(x => x && x.element === 'C');
    const piNext = a => carbonNb(a).some(n => bondsOf(n).some(b => b.order === 2 && other(b, n) !== a && other(b, n).element === 'C'));
    const deg = a => carbonNb(a).length;
    const score = a => [ST['methyl cation'], ST['1° cation'], ST['2° cation'], ST['3° cation']][Math.min(3, deg(a))] + (piNext(a) ? (BO['benzylic or allylic'] || 0) : 0);
    const cls = a => ['methyl', '1°', '2°', '3°'][Math.min(3, deg(a))] + (piNext(a) ? ' benzylic/allylic' : '');
    const ringSize = (a, b) => {                          // smallest ring through the a–b bond (null if none)
      const seen = new Map([[a.id, 0]]), qq = [a];
      while (qq.length) { const x = qq.shift(); for (const bd of bondsOf(x)) { const y = other(bd, x); if ((x === a && y === b) || seen.has(y.id)) continue;
        seen.set(y.id, seen.get(x.id) + 1); if (y === b) return seen.get(y.id) + 1; if (seen.get(y.id) < 6) qq.push(y); } }
      return null;
    };
    // form the cation
    let lgEl = null;
    if (spec.leaving) {
      const lg = A(hit.atoms[spec.leaving - 1]); lgEl = lg.element;
      g.bonds = g.bonds.filter(b => b.a !== lg.id && b.b !== lg.id); g.atoms = g.atoms.filter(x => x !== lg);
    } else {
      const o = A(hit.atoms[1]); const bd = bondBetweenAtoms(cat0, o); if (bd) bd.order = 1;    // H+ adds to the other C=C carbon
    }
    cat0.charge = 1;
    const smi = () => { try { const x = window.Chem.toSmiles(g); const mm = R.get_mol(x); const c = mm ? mm.get_smiles() : x; if (mm) mm.delete(); return c; } catch (e) { return ''; } };
    const snapshot = () => JSON.parse(JSON.stringify({ atoms: g.atoms, bonds: g.bonds, nextId: g.nextId }));
    let c = cat0; const shifts = [], steps = [];
    for (let iter = 0; iter < 4; iter++) {
      const cur = score(c), cands = [];
      for (const n of carbonNb(c)) {
        if (bondsOf(n).some(b => b.order > 1) || n.charge) continue;          // only an sp3 neighbour can shift a group
        // 1,2-hydride shift: the charge simply moves (the H counts follow)
        if (hOf(n) > 0) { c.charge = 0; n.charge = 1; cands.push({ kind: '1,2-hydride shift', n, score: score(n), pri: PR['1,2-hydride shift'] || 1 }); n.charge = 0; c.charge = 1; }
        // 1,2-alkyl shift, or ring expansion when the moving bond is part of a small ring that does not contain c
        for (const mv of carbonNb(n)) {
          if (mv === c || bondsOf(mv).some(b => b.order > 1)) continue;
          const rs = ringSize(n, mv), inRingWithC = ringSize(c, n);
          const ring = rs && rs <= 4 && !inRingWithC;
          if (rs && !ring) continue;                                         // other ring bonds don't migrate here
          const bd = bondBetweenAtoms(n, mv);
          const swap = () => { if (bd.a === n.id) bd.a = c.id; else bd.b = c.id; };
          const back = () => { if (bd.a === c.id) bd.a = n.id; else bd.b = n.id; };
          swap(); c.charge = 0; n.charge = 1;
          const sc = score(n) + (ring ? (BO['ring expansion of a 4-membered ring'] || 0) : 0);
          n.charge = 0; c.charge = 1; back();
          const kind = ring ? 'ring expansion' : (deg(mv) === 1 && hOf(mv) === 3 ? '1,2-methyl shift' : '1,2-alkyl shift');
          cands.push({ kind, n, mv, score: sc, pri: PR[kind] || 3, ringFrom: ring ? rs : null });
        }
      }
      const better = cands.filter(x => x.score > cur + 1e-9).sort((a, b) => b.score - a.score || a.pri - b.pri);
      if (!better.length) break;
      const best = better[0], before = cls(c);
      // drawing of this step: the cation with the moving H (made explicit) or group, and its arrow
      const pic = snapshot();
      const arrows = [];
      if (best.kind === '1,2-hydride shift') {
        const nn = pic.atoms.find(x => x.id === best.n.id), cc = pic.atoms.find(x => x.id === c.id);
        const h = { id: pic.nextId++, element: 'H', charge: 0, x: nn.x + (nn.x - cc.x) * 0.3 + (cc.y - nn.y) * 0.6, y: nn.y + (nn.y - cc.y) * 0.3 + (nn.x - cc.x) * 0.6 };
        pic.atoms.push(h); pic.bonds.push({ a: nn.id, b: h.id, order: 1 });
        arrows.push({ from: { kind: 'bond', a: nn.id, b: h.id }, to: { kind: 'bond', a: cc.id, b: h.id }, label: 'the C–H bond (H with its electron pair) moves to the C⁺' });
      } else {
        arrows.push({ from: { kind: 'bond', a: best.n.id, b: best.mv.id }, to: { kind: 'bond', a: c.id, b: best.mv.id },
                      label: best.kind === 'ring expansion' ? 'a ring C–C bond moves to the C⁺ (the ring grows by one)' : 'the C–C bond (the group with its electron pair) moves to the C⁺' });
      }
      // make the shift for real
      if (best.kind === '1,2-hydride shift') { c.charge = 0; best.n.charge = 1; }
      else { const bd = bondBetweenAtoms(best.n, best.mv); if (bd.a === best.n.id) bd.a = c.id; else bd.b = c.id; c.charge = 0; best.n.charge = 1; }
      const after = cls(best.n);
      const text = `${best.kind}: ${before} cation → ${after} cation` + (best.ringFrom ? ` (a ${best.ringFrom}-membered ring becomes ${best.ringFrom + 1}-membered, releasing ring strain)` : '');
      shifts.push(text);
      steps.push({ step: 'shift', type: `carbocation rearrangement: ${best.kind}`, note: `Happens as soon as the cation forms, before the nucleophile or base can react: ${text}. A shift only happens when it gives a more stable cation (Rearrangements sheet).`,
                   graph: pic, gArrows: arrows, arrowLabels: arrows.map(a => a.label), smiles: smi() });
      c = best.n;
    }
    if (!shifts.length) return null;
    c.charge = 0;
    if (spec.then === 'add' || spec.then === 'rerun') {
      const el = spec.then === 'add' ? spec.add : lgEl;
      const x = { id: g.nextId++, element: el, charge: 0, x: c.x + 30, y: c.y + 30 };
      g.atoms.push(x); g.bonds.push({ a: c.id, b: x.id, order: 1 });
    } else if (spec.then === 'eliminate') {
      const betas = carbonNb(c).filter(n => hOf(n) > 0 && bondsOf(n).every(b => b.order === 1)).sort((a, b) => deg(b) - deg(a));
      if (betas.length) bondBetweenAtoms(c, betas[0]).order = 2;                 // Zaitsev: the more substituted alkene
    }
    const out = smi();
    return out ? { smiles: out, shifts, steps } : null;
  }
  /* ---- Unknown vs no reaction (sheet Known_No_Reaction) ---- */
  function knownNoReaction(R, reagentKey, presentSubs) {
    const RR = window.REACTION_RULES, hits = [];
    for (const k of RR.knownNone || []) {
      if (!presentSubs.has(k.frag)) continue;
      let ok = false, who = k.name;
      if (k.reagent === 'any nucleophile') {
        const n = reagentKey && (RR.nucleophiles || []).find(x => { if (!abKeys.has(x.smiles)) abKeys.set(x.smiles, keyOf(R, x.smiles)); return abKeys.get(x.smiles) === reagentKey; });
        ok = !!n; if (n) who = n.name;
      }
      else if (k.reagent === 'any acid' || k.reagent === 'any base') ok = !!acidBaseEntry(R, k.reagent.slice(4), reagentKey);
      else { if (!abKeys.has(k.reagent)) abKeys.set(k.reagent, keyOf(R, k.reagent)); ok = !!reagentKey && abKeys.get(k.reagent) === reagentKey; }
      if (ok) hits.push(Object.assign({}, k, { name: who }));
    }
    return hits;
  }
  /* ---- Learned reactions (learned_reactions.xlsx -> data/learned/, built by tools/learned_to_js.py) ----
     Reactions recorded in patents (data/rxnkb), each checked when it was learned: its template (a reaction SMARTS
     that software extracted from the recorded reaction) gives the recorded product. The rules in reaction_rules.xlsx
     come first (they explain: fragments, conditions, pKa, stereo, mechanism); a rule's answer is then compared with a
     recorded reaction of the same molecules (predict). When no rule answers, the learned reactions are used:
       looked up: the same starting materials + reagent as a learned reaction -> its recorded product;
       predicted: new molecules -> the templates recorded with this reagent (or, for a new reagent, with its parts)
                  are run on them, and the product most of them agree on wins (each template weighted by how many
                  learned reactions used it with this reagent). */
  /* the learned data is cut into small files (tools/learned_to_js.py); each is loaded once, when a reaction needs it */
  const LEARNED = { meta: undefined, parts: new Map() };
  const loadScriptOnce = src => new Promise((res, rej) => {
    const el = document.createElement('script'); el.src = src;
    el.onload = () => { el.remove(); res(); }; el.onerror = () => { el.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(el);
  });
  async function unpackGz(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return JSON.parse(await new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).text());
  }
  async function learnedMeta() {
    if (LEARNED.meta === undefined) {
      try { await loadScriptOnce('data/learned/learned_meta.js'); LEARNED.meta = window.LEARNED_META || null; }
      catch (e) { LEARNED.meta = null; }               // no data/learned: the page works with the workbook rules only
    }
    return LEARNED.meta;
  }
  function learnedPart(name) {                         // 'rx/00012', 't/003', 'x/AB', 'g/AB', 'p/AB' (a missing file = null)
    if (!LEARNED.parts.has(name)) LEARNED.parts.set(name, (async () => {
      try { await loadScriptOnce('data/learned/' + name + '.js'); } catch (e) { return null; }
      const P = (window.LEARNED_PARTS || {})[name];
      if (!P) return null;
      delete window.LEARNED_PARTS[name];
      try { return await unpackGz(P.gz); }
      catch (e) { LEARNED.parts.delete(name); return null; }       // e.g. no DecompressionStream: the rules still answer
    })());
    return LEARNED.parts.get(name);
  }
  const hasCarbon = smi => { let g = null; try { g = window.Chem.parseSmiles(smi); } catch (e) { g = null; } return !!g && g.atoms.some(a => a.element === 'C'); };
  /* runRuleAll for a learned template: the molecules are the starting materials plus the reagent parts poolExtras
     allows (OH-, Br2, HNO3, H2O, MeO-, MeI ...: a recorded reaction counts them as reactants); every outcome must use
     at least one starting material (tools/learn_reactions.py run_template checks the learned reactions the same way) */
  function runLearned(R, smarts, subs, extras) {
    let rx = null;
    try { rx = R.get_rxn(smarts); } catch (e) { rx = null; }
    if (!rx) return [];
    const pool = subs.concat(extras.filter(x => !subs.includes(x)));
    const n = smarts.split('>>')[0].split('.').length, orders = [], seen = new Set(), outs = [];
    const pick = cur => { if (cur.length === n) { if (Math.min(...cur) < subs.length) orders.push(cur.slice()); return; } pool.forEach((_, i) => { if (!cur.includes(i)) { cur.push(i); pick(cur); cur.pop(); } }); };
    if (n <= pool.length) pick([]);
    const canon = smi => { const m = R.get_mol(smi); const c = m ? m.get_smiles() : null; if (m) m.delete(); return c; };
    try {
      for (const order of orders) {
        const ml = new R.MolList(), mols = order.map(i => R.get_mol(pool[i]));
        if (mols.some(m => !m)) { mols.forEach(m => m && m.delete()); ml.delete(); continue; }
        mols.forEach(m => ml.append(m));
        try {
          const res = rx.run_reactants(ml, 50);
          for (let k = 0; k < res.size(); k++) {
            const pl = res.get(k), out = [];
            for (let b = 0; b < pl.size(); b++) { const pm = pl.at(b); try { out.push(pm.get_smiles()); } catch (e) {} pm.delete(); }
            pl.delete();
            const cs = out.map(canon);
            if (!cs.length || cs.some(x => !x)) continue;
            const key = cs.slice().sort().join('.');
            if (!seen.has(key)) { seen.add(key); outs.push({ products: cs, used: order.map(i => pool[i]).filter(x => subs.includes(x)) }); }
          }
          res.delete();
        } catch (e) { /* this order does not fit */ }
        ml.delete(); mols.forEach(m => m.delete());
      }
    } finally { rx.delete(); }
    return outs;
  }
  /* By-products of a learned reaction: its leaving pieces (the starting-material atoms the product does not contain,
     H added where a bond broke) looked up in the Byproducts sheet of reaction_rules.xlsx, which says what they really
     leave as (a methyl ester hydrolysed by hydroxide loses CH3 + O as methanol; a Boc group leaves as CO2 + isobutene).
     A row may need a reagent part ("when the reagent contains"). No row: the pieces themselves, said so. */
  function learnedByproducts(leaving, reagent) {
    if (!leaving) return { list: [], why: '' };
    const rows = (window.REACTION_RULES || {}).byproducts || [];
    const parts = new Set(String(reagent || '').split('.').filter(Boolean));
    const fits = r => r.pieces === leaving && (!r.when || r.when.split('.').every(x => parts.has(x)));
    const hit = rows.filter(fits).sort((a, b) => (b.when ? 1 : 0) - (a.when ? 1 : 0))[0];
    if (hit) return { list: hit.byproducts.filter(Boolean), why: `${hit.text || hit.byproducts.join(' + ')} (Byproducts sheet: ${hit.why})` };
    const pieces = leaving.split('.').filter(x => x && x !== '?');      // raw pieces are described, not drawn as products
    if (!pieces.length) return { list: [], why: '' };
    return { list: [], why: `not in the Byproducts sheet yet; the starting-material atoms the product does not contain are ${pieces.join(' + ')} (H added where a bond broke)` };
  }
  /* Selectivity (tools/learn_selectivity.py): each template's record when it could have reacted (wins, times applicable),
     learned from the recorded reactions; the products are ranked by the method that predicted held-out reactions best */
  async function learnedWeights() {
    const W = await learnedPart('w/winrates');
    return W || { method: 'count', win: {} };
  }
  function weightOf(method, n, win) {
    const wr = win ? (win[0] + 1) / (win[1] + 2) : 0.5;
    return method === 'count*winrate' ? n * wr : method === 'count*winrate^2' ? n * wr * wr : method === 'winrate-sum' ? wr : method === 'winrate-max' ? wr : n;
  }
  /* Keys (tools/learn_reactions.py computes the same): a reagent is keyed by its distinct parts (a record lists each
     species once, K2CO3 as CO3(2-) + one K+; a student writes both K+); a recorded reaction is keyed by every distinct
     molecule, starting materials and reagent parts together, so a reagent typed in a starting-material box (or the
     other way round) still finds it. '' = no reagent; null = a reagent RDKit cannot key (never treated as 'no reagent'). */
  const reagentParts = reagent => [...new Set(String(reagent || '').split('.').filter(Boolean))].sort();
  // a molecule counts as many times as written (O.O, two waters, is a different reaction from O); a single-atom ion
  // once (records list Na+, K+, Cl- once whatever the salt) -- learn_reactions.reagent_counted
  const singleIon = x => /^\[[A-Z][a-z]?[+-]\d*\]$/.test(x);
  const reagentCounted = reagent => { const p = String(reagent || '').split('.').filter(Boolean); return p.filter(x => !singleIon(x)).concat([...new Set(p.filter(singleIon))].sort()).sort(); };
  const learnedReagentKey = (R, reagent) => reagent ? (keyOf(R, reagentCounted(reagent).join('.')) || null) : '';
  const learnedLookupKey = (R, subs, reagent) => { const s = [...new Set(subs)].sort(); return keyOf(R, s.concat(reagentCounted(reagent).filter(x => !s.includes(x))).sort().join('.')) || null; };
  /* reagent parts that may take part in a learned template: carbon-free pieces (OH-, Br2, HNO3, H2O) and the parts the
     Reagents sheet says transfer atoms (nucleophiles such as MeO-, alkylating agents such as MeI ...); bases like Et3N
     or carbonate are left out so they never turn up as false nucleophiles (learn_reactions.pool_extras: same rule) */
  const TRANSFER_ROLES = new Set(['nucleophile', 'electrophile / alkylating agent', 'organometallic reagent', 'halogenating agent',
    'protecting-group reagent', 'reactant (inorganic)', 'nucleophile (sulfur ylide precursor)']);
  function poolExtras(reagent) {
    const RG = (window.REACTION_RULES || {}).reagents || {};
    // as many copies of each part as the reagent lists: two waters are written O.O and may be used twice, one O once
    return String(reagent || '').split('.').filter(Boolean).sort()
      .filter(x => !hasCarbon(x) || (RG[x] && (TRANSFER_ROLES.has(RG[x].role) || TRANSFER_ROLES.has(RG[x].also))));
  }
  const MIN_OBSERVATIONS = 2;                            // a template must be seen at least twice to predict (tools/learn_selectivity.py)
  let currentRun = 0;                                    // the arrow's run number: an older run never writes over a newer one
  /* hand control back to the page for a moment (lets it draw and take clicks). One channel kept for the whole page: a
     MessageChannel made per call can be garbage-collected before its message arrives, and the wait would never end;
     a timer would be slowed to once a second (or a minute) when the tab is hidden. */
  const breathe = (() => {
    const ch = new MessageChannel(), waiting = [];
    ch.port1.onmessage = () => { const r = waiting.shift(); if (r) r(); };
    return () => new Promise(r => { waiting.push(r); ch.port2.postMessage(0); });
  })();
  /* Recorded reactions with exactly these molecules (learned_reactions.xlsx), most recorded first */
  async function learnedLookup(R, subs, reagent) {
    const meta = await learnedMeta();
    if (!meta || !subs.length) return null;
    const key = learnedLookupKey(R, subs, reagent);
    if (!key) return null;
    const X = await learnedPart('x/' + key.slice(0, 2));
    const hits = X && X[key];
    if (!hits || !hits.length) return null;
    const rows = [];
    for (const k of hits) {
      const c = Math.floor(k / meta.rx_chunk), part = await learnedPart('rx/' + String(c).padStart(5, '0'));
      if (part && part[k - c * meta.rx_chunk]) rows.push(part[k - c * meta.rx_chunk]);
    }
    rows.sort((a, b) => b[6] - a[6]);
    return rows.length ? rows : null;
  }
  const learnedSource = r => [r[7], r[8], r[9] !== '' && r[9] != null ? `${r[9]}% yield` : ''].filter(Boolean).join(', ') + (r[6] > 1 ? ` (recorded ${r[6]}×)` : '');
  async function learnedReaction(R, subs, reagent, reagentKey, condText, read, reagentFrags) {
    const meta = await learnedMeta();
    if (!meta) return null;
    const myRun = currentRun;
    const pad = (k, n) => String(k).padStart(n, '0');
    const shard = key => key ? key.slice(0, 2) : '__';
    const rxRow = async k => { const c = Math.floor(k / meta.rx_chunk), part = await learnedPart('rx/' + pad(c, 5)); return part ? part[k - c * meta.rx_chunk] || null : null; };
    const tmpl = async t => {
      const c = Math.floor(t / meta.tmpl_chunk), part = await learnedPart('t/' + pad(c, 3)), j = t - c * meta.tmpl_chunk;
      return part && part.smarts[j] ? { id: part.id[j], smarts: part.smarts[j], recorded: part.recorded[j], leaving: (part.leaving || [])[j] || '' } : null;
    };
    const fragDiff = (used, prods) => {
      const a = fragmentCounts(R, used), b = fragmentCounts(R, prods);
      return { consumes: Object.keys(a).filter(id => (b[id] || 0) < a[id]).sort().map(id => [id]),
               forms: Object.keys(b).filter(id => b[id] > (a[id] || 0)).sort().map(id => [id]) };
    };
    const finish = (products, used, id, learnedText, note, leaving, leavingIsUsual) => {
      const prodRead = products.map(smiles => ({ smiles, frags: fragmentsOf(R, smiles) }));
      const fd = fragDiff(used, products);
      const bp = learnedByproducts(leaving, reagent);
      if (bp.why) learnedText += `. By-products${leavingIsUsual ? ' (from the leaving pieces this template usually has)' : ''}: ${bp.why}`;
      return { products: products.concat(bp.list.filter(x => !products.includes(x))), kind: 'learned', note,
        mechNote: 'Learned reactions have no mechanism: a recorded reaction gives only the starting materials and the product, never the arrows. If this reaction belongs in reaction_rules.xlsx, its mechanism can be added there.',
        reading: { subs: read, reagent: reagent ? { smiles: reagent, frags: reagentFrags } : null, learned: learnedText,
                   rule: { id, consumes: fd.consumes, forms: fd.forms }, prods: prodRead, formsOk: true },
        code: reactionCode(R, subs, reagent, products, condText) };
    };
    // looked up: recorded reactions with exactly these molecules (the most recorded one first)
    const rows = await learnedLookup(R, subs, reagent);
    if (rows) {
      const r = rows[0], T = await tmpl(r[2]), products = r[1].split('~')[3].split('.').filter(Boolean);
      const also = [...new Set(rows.slice(1).map(x => x[1].split('~')[3]))].filter(x => x !== r[1].split('~')[3]);
      return finish(products, subs, r[0],
        `looked up: these molecules are learned reaction ${r[0]} (recorded: ${learnedSource(r)}; its code: ${r[1]}); template ${T ? T.id : '?'}` +
          (T ? `, recorded ${T.recorded.toLocaleString()}× in all` : '') +
          (also.length ? `; the same molecules are also recorded giving: ${also.join('; ')}` : ''),
        `${r[0]} · learned reaction (learned_reactions.xlsx), looked up: recorded in ${learnedSource(r)} · template ${T ? T.id : '?'} · conditions recorded: ${r[1].split('~')[2] || 'none'} · no mechanism`, r[5], false);
    }
    // predicted: never with no reagent at all. A recorded reaction without one is usually a record that left the
    // reagent out (a hydrogenation written without its H2/Pd), so it would teach the wrong thing; those are only looked up.
    if (!reagent) return { miss: 'with no reagent the page only looks up recorded reactions (recorded reactions without a reagent are often ones whose reagent was not written down), and these molecules are not one of them' };
    /* Which templates to try, in turn (each step only when no template of the step before fits these molecules;
       measured on held-out reactions by tools/selectivity_experiments.py, this order predicts best):
         1. the templates learned with this exact reagent (skipped if RDKit cannot key the reagent);
         2. those learned with its parts (lone ions such as Na+ left out unless nothing else is there);
         3. those learned with reagents that do the same jobs (Reagents sheet roles: K2CO3 ~ Cs2CO3, "base"). */
    const rk = learnedReagentKey(R, reagent);
    const stages = [];
    {
      /* this reagent; and this reagent + a starting material that is itself a known reagent (MeI typed as a starting
         material: records often list it WITH the reagents, "MeI + K2CO3"); counts summed (learn_selectivity.exact_list) */
      const RG = (window.REACTION_RULES || {}).reagents || {}, keys = rk ? [rk] : [], also = [];
      for (const x of subs) if (RG[x]) { const k = learnedReagentKey(R, reagent + '.' + x); if (k) { keys.push(k); also.push(RG[x].abbreviation || RG[x].abbr || RG[x].name || x); } }
      const sums = new Map();
      for (const k of [...new Set(keys)]) {
        const G = await learnedPart('g/' + shard(k));
        for (const [t, n, ex] of (G && G[k]) || []) { const o = sums.get(t); if (o) o[0] += n; else sums.set(t, [n, ex]); }
      }
      if (sums.size) stages.push([[...sums.entries()].map(([t, [n, ex]]) => [t, n, ex]).sort((a, b) => b[1] - a[1]),
        'this reagent' + (also.length ? ` (also counting records that list ${also.join(', ')} with the reagent)` : ''), false]);
    }
    {
      /* only records whose reagent is PART of this one: every combination of this reagent's essential parts (lone
         ions such as Na+ left out unless nothing else is there) is looked up, so Br2 + HBr is never borrowed for HBr
         alone, while a record with just the base still counts for base + catalyst (learn_reactions.subset_keys) */
      const parts = reagentParts(reagent);
      const single = x => /^\[[A-Za-z]+[+-]\d*\]$/.test(x);
      const use = parts.filter(x => !single(x)).length ? parts.filter(x => !single(x)) : parts;
      const keys = [...new Set(use.map(x => keyOf(R, x) || '').filter(Boolean))].sort();
      const sum = new Map(), low = keys.length <= 8 ? 1 : keys.length - 2;
      const combos = [];
      const pick = (start, cur) => { if (cur.length >= low) combos.push(cur.slice()); for (let i = start; i < keys.length; i++) { cur.push(keys[i]); pick(i + 1, cur); cur.pop(); } };
      pick(0, []);
      combos.sort((a, b) => b.length - a.length);         // most parts first, each size in itertools order (ties match Python)
      for (const c of combos) {
        const sk = c.join(' '), P = await learnedPart('s/' + shard(sk));
        for (const [t, n, ex] of (P && P[sk]) || []) { const o = sum.get(t) || { t, n: 0, ex, size: 0 }; o.n += n; o.size = Math.max(o.size, c.length); sum.set(t, o); }
      }
      if (sum.size) stages.push([[...sum.values()].sort((a, b) => b.size - a.size || b.n - a.n).map(o => [o.t, o.n, o.ex]), `reagents made of parts of this one (${use.join(', ')})`, true]);
    }
    const sig = roleSignature(reagent);
    if (sig) { const RL = await learnedPart('r/' + roleShard(sig)); if (RL && RL[sig]) stages.push([RL[sig], `reagents with the same role${sig.includes('|') ? 's' : ''} (${sig.split('|').join(' + ')}; Reagents sheet)`, false]); }
    if (!stages.length) return { miss: `no learned reaction uses this reagent${rk ? '' : ' (RDKit cannot key it)'}, any of its parts or a reagent with the same roles` };
    const extras = poolExtras(reagent);
    const W = await learnedWeights();
    const stats = { stages: stages.map(x => x[1]), load: 0, run: 0, slowest: null, tried: 0 };
    window.ReactionsLearnedStats = stats;
    const report = [];
    let score = new Map(), how = '', list = [], summed = false, tried = 0;
    for ([list, how, summed] of stages) {
      score = new Map(); tried = 0;
      // at least two observations: a template seen only once with this reagent (one record, possibly a recording
      // error or a one-off) is not used to predict for other molecules; it can still be looked up exactly
      const once = list.filter(x => x[1] < MIN_OBSERVATIONS).length;
      list = list.filter(x => x[1] >= MIN_OBSERVATIONS);
      const want = list.length;                          // every template learned with it: no limit on how many are tried
      for (const [t, n, ex] of list) {
        if (tried % 10 === 0) { if (myRun === currentRun) ruleNote.textContent = `Running the learned templates for ${how}: ${tried} of ${want}…`; await breathe(); }
        const ta = Date.now();
        const T = await tmpl(t);
        stats.load += Date.now() - ta;
        if (!T) continue;
        tried++; stats.tried++;
        const w = weightOf(W.method, n, W.win[t]);
        const tb = Date.now(), outs = runLearned(R, T.smarts, subs, extras), dt = Date.now() - tb;
        stats.run += dt;
        if (!stats.slowest || dt > stats.slowest.ms) stats.slowest = { id: T.id, ms: dt };
        for (const o of outs) {
          const key = o.products.slice().sort().join('.');
          const e = score.get(key) || { products: o.products, used: o.used, w: 0, best: null };
          e.w = W.method === 'winrate-max' ? Math.max(e.w, w + 1e-6 * n) : e.w + w;
          if (!e.best || w > e.best.w) e.best = { t, n, ex, id: T.id, w, leaving: T.leaving, win: W.win[t] };
          score.set(key, e);
        }
      }
      report.push(`${how}: tried all ${tried.toLocaleString()} seen at least ${MIN_OBSERVATIONS}×` + (once ? ` (${once.toLocaleString()} seen only once left out)` : ''));
      if (score.size) break;                             // this step found products: the later steps are not needed
    }
    if (!score.size) return { miss: `no learned template that was tried fits these molecules (${report.join('; ')})` };
    const ranked = [...score.values()].sort((a, b) => b.w - a.w), top = ranked[0];
    const ex = await rxRow(top.best.ex);
    const wrOf = e => e.best.win ? (e.best.win[0] + 1) / (e.best.win[1] + 2) : 0.5;
    const share = W.method === 'winrate-max'
      ? `its best template wins ${Math.round(100 * wrOf(top))}% of the time when it could react` + (ranked.length > 1 ? `; next: ${ranked[1].products.join('.')} (${Math.round(100 * wrOf(ranked[1]))}%)` : '')
      : `${Math.round(100 * top.w / ranked.reduce((t, e) => t + e.w, 0))}% of the weight` + (ranked.length > 1 ? `; next: ${ranked[1].products.join('.')}` : '');
    const wr = top.best.win ? `; when it could react it won ${top.best.win[0]} of ${top.best.win[1]} times in the recorded reactions (selectivity, tools/learn_selectivity.py)` : '';
    const learnedN = `learned ${top.best.n}×${summed ? ' in total' : ''} with ${how}`;
    return finish(top.products, top.used, top.best.id,
      `predicted: ${report.join('; ')}; ${ranked.length} different product(s); this one: ${share} (ranking: ${W.method}). ` +
        `Leading template ${top.best.id}, ${learnedN}${wr}` + (ex ? `; example ${ex[0]}: ${ex[1]} (${learnedSource(ex)})` : ''),
      `Predicted from learned reactions (learned_reactions.xlsx): template ${top.best.id}, ${learnedN}` + (ex ? ` · example ${ex[0]}` : '') + ' · no mechanism',
      top.best.leaving, true);
  }
  /* the jobs a reagent's parts do (Reagents sheet), spectator ions and solvents left out; none if a part is not in the
     sheet (tools/learned_to_js.py role_signature builds the index the same way) */
  const SKIP_ROLES = new Set(['counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent']);
  const roleShard = sig => { let t = 0; for (const ch of sig) t += ch.codePointAt(0); return (t % 256).toString(16).padStart(2, '0'); };   // learned_to_js.role_shard
  function roleSignature(reagent) {
    const RG = (window.REACTION_RULES || {}).reagents || {}, sig = new Set();
    for (const p of reagentParts(reagent)) {
      const r = RG[p];
      if (!r) return null;
      if (!SKIP_ROLES.has(r.role)) sig.add(r.role);
    }
    return sig.size ? [...sig].sort().join('|') : null;
  }
  /* The workbook rules explain; a recorded reaction of these exact molecules shows what one patent got. They are
     compared (a recorded reaction has no temperature or amounts, so its conditions are never compared):
       a record gives the rule's product -> the rule's answer, citing the record;
       same compound, and the record has MORE stereo (stereocentres / C=C geometry) -> the recorded stereo;
       a different product -> the recorded product, saying what the rule would have given (often a condition the rule
       does not model: SOCl2 IN METHANOL gives the methyl ester);
       a rule result that depends on what the student set (no reaction without heat / too weak a base, excess, 2 equiv)
       is never replaced, only annotated with the record. */
  async function predict() {
    const R0 = await loadRD();
    const subs0 = slots.filter(x => x.smiles).map(x => x.smiles), reagent0 = aboveApi.get(), cond0 = belowApi.get();
    const res = await predictRules();
    if (!res || res.kind === 'learned' || res.kind === 'unknown' || !R0) return res;
    let rows = null;
    try { rows = await learnedLookup(R0, subs0, reagent0); } catch (e) { rows = null; }
    if (!rows) return res;
    const R = R0, subs = subs0, reagent = reagent0, P = readConditions(cond0);
    const fullKey = x => keyOf(R, x) || '';
    const stereoMarks = x => (x.match(/@|\/|\\/g) || []).length;
    // InChIKey: first block = the compound (connectivity), last letter = protonation; stereo sits in the middle block
    const sameCompound = (a, b) => { const ka = fullKey(a), kb = fullKey(b); return ka && kb && ka.slice(0, 14) === kb.slice(0, 14); };
    const sameButStereo = (a, b) => sameCompound(a, b) && fullKey(a).slice(-1) === fullKey(b).slice(-1);
    const got = res.products || [];
    const recProducts = rows.map(r => r[1].split('~')[3].split('.').filter(Boolean));
    const hit = rows.find((r, i) => got.length && recProducts[i].every(x => got.includes(x)));
    if (hit) { res.note += ` · also recorded: ${hit[0]} (${learnedSource(hit)})`; return res; }
    const r = rows[0], want = recProducts[0];
    const studentSet = res.kind === 'blocked' || P.equiv != null || P.keys.has('excess') || P.keys.has('catalytic') || /applied \d+×|limiting reagent/.test(res.note || '');
    if (got.length && want.every(x => got.some(g => sameCompound(g, x)))) {          // the same compound(s), maybe as a salt
      if (!studentSet && want.every(x => { const g = got.find(y => sameButStereo(y, x)); return g && stereoMarks(x) > stereoMarks(g); })) {
        res.products = got.map(g => want.find(x => sameButStereo(x, g)) || g);
        res.note += ` · stereochemistry from the recorded reaction ${r[0]} (${learnedSource(r)})`;
        res.code = reactionCode(R, subs, reagent, res.products, cond0);
      } else res.note += ` · also recorded: ${r[0]} (${learnedSource(r)}; written there as ${want.join('.')}${want.some(x => got.some(g => sameCompound(g, x) && !sameButStereo(g, x))) ? ', a protonated (salt) form' : ''})`;
      return res;
    }
    if (studentSet) { res.note += ` · recorded under unstated conditions: ${r[0]} gives ${want.join('.')} (${learnedSource(r)})`; return res; }
    const read = subs.map(smiles => ({ smiles, frags: fragmentsOf(R, smiles) }));
    let lr = null;
    try { lr = await learnedReaction(R, subs, reagent, null, cond0, read, reagent ? fragmentsOf(R, reagent) : []); } catch (e) { lr = null; }
    if (!lr || !lr.products) return res;
    const ruleSaid = got.length ? `the workbook rule (${res.note.split(' · ')[0]}) would give ${got.join('.')}` : `the workbook said: ${res.note}`;
    lr.note += ` · ${ruleSaid}; these exact molecules are recorded giving ${want.join('.')} (the record's conditions were not compared)`;
    if (lr.reading) lr.reading.learned += `. Compared with reaction_rules.xlsx: ${ruleSaid}. The recorded product is shown because a real reaction of exactly these molecules gave it; the record has no temperature or amounts, so the difference may come from conditions the rule does not model (e.g. the solvent taking part)`;
    return lr;
  }
  async function predictRules() {
    const R = await loadRD(), RR = window.REACTION_RULES;
    if (!R) return { note: 'RDKit could not start, so no rule can be run.' };
    if (!RR || !RR.rules) return { note: 'The reaction rules are not loaded (data/reaction_rules.js: run tools/rules_to_js.py).' };
    const subs = slots.filter(x => x.smiles).map(x => x.smiles);
    if (!subs.length) return { note: 'Add the starting molecules on the left first.' };
    const reagent = aboveApi.get(), reagentKey = reagent ? keyOf(R, reagent) : null;
    const condText = belowApi.get(), P = readConditions(condText);
    const short = [];                                   // rules that would have fired but lack a condition
    const read = subs.map(smiles => ({ smiles, frags: fragmentsOf(R, smiles) }));
    const reagentFrags = reagent ? fragmentsOf(R, reagent) : [];
    const presentAll = new Set([...read.flatMap(x => x.frags), ...reagentFrags]);   // the reagent counts for fragment conditions
    const cands = [];
    for (const rule of RR.rules) {
      const ev = evalRule(R, rule, subs, reagentKey, P, presentAll);
      if (ev.status === 'short') short.push(ev.msg);
      else if (ev.status === 'ok') cands.push(ev);
    }
    if (cands.length) {
      let pick = chooseCandidate(R, cands);
      const rule = pick.cand.rule, fcheck = pick.cand.fcheck, abCheck = pick.cand.abCheck;
      const ranks = (RR.reactivity || {})[pick.group] || {};
      const selection = [];
      if (pick.cand.easInfo) selection.push('directing groups (Directing_Groups): ' + pick.cand.easInfo);
      if (pick.group && pick.frag) {
        selection.push(`group "${pick.group}" (Reactivity_Order): ${fragName(pick.frag)} reacts first` +
          (ranks[pick.frag] && ranks[pick.frag].why ? ` (${ranks[pick.frag].why})` : ''));
      }
      // the reagent's main product, and any by-products
      let main = pick.outcome.products.slice(), extra = rule.byproducts.map(k => RR.molecules[k].smiles);
      const stereo = [], chiralStart = subs.some(x => x.includes('@'));
      if (!P.keys.has('excess') && !(P.equiv != null && typeof rule.equiv === 'number' && P.equiv >= 2 * rule.equiv)) {
        const st = applyStereo(R, rule, pick.outcome);
        if (st && st.products) main = st.products;
        if (st && st.info) stereo.push(st.info);
      }
      if (!stereo.length && rule.stereo && rule.stereo.outcome) stereo.push(rule.stereo.outcome);
      let shiftSteps = [], rearr = [];
      if (rule.cation && !P.keys.has('excess')) {
        const csub = pick.outcome.used.find(x => { const mm = R.get_mol(x); if (!mm) return false; const ok = matchesOf(mm, qOf(R, rule.cation.smarts)).length > 0; mm.delete(); return ok; });
        const cr = csub ? carbocationRoute(R, csub, rule.cation) : null;
        if (cr) {
          rearr = cr.shifts; shiftSteps = cr.steps;
          if (rule.cation.then === 'rerun') {                // run the reaction again on the rearranged starting material
            const ev2 = evalRule(R, rule, subs.map(x => x === csub ? cr.smiles : x), reagentKey, P, presentAll);
            if (ev2.status === 'ok') { const p2 = chooseCandidate(R, [ev2]); main = p2.outcome.products.slice(); if (p2.cand.easInfo) selection.push('directing groups (Directing_Groups): ' + p2.cand.easInfo); }
          } else main = [cr.smiles];
          stereo.length = 0;
        }
      }
      const steps = [rule.id];
      /* Amounts (Reagent equivalents column): with no amount, one reaction (as before); 'excess' = keep going until
         nothing in reach is left; a number = as many reactions as that amount pays for, most reactive site first (the
         same group's rules, or the same rule again on its product); less than one reaction's worth = only that fraction
         of the starting material reacts and the rest is left over. */
      const eqOf = r => (typeof r.equiv === 'number' ? r.equiv : 0);          // catalytic / not stated: not used up
      const excess = P.keys.has('excess'), have = excess ? Infinity : P.equiv;
      const amounts = [], rname = (RR.molecules[reagentKey] || {}).name || (reagent ? 'the reagent' : '');
      let leftover = null, budget = have == null ? 0 : have - eqOf(rule);
      if (have != null && have !== Infinity && eqOf(rule) > 0 && have < eqOf(rule) - 1e-9) {
        leftover = pick.outcome.used.slice();                // the starting materials that don't get to react
        const frac = Math.round(have / eqOf(rule) * 100);
        amounts.push(`${have} equiv of ${rname}: one reaction needs ${eqOf(rule)} equiv, so only about ${frac}% of the starting material reacts; the rest is left over (the reagent is the limiting reagent)`);
      } else if (have != null && have !== Infinity && eqOf(rule) > 0) {
        amounts.push(`${have} equiv of ${rname}; ${rule.id} uses ${eqOf(rule)} equiv per reaction`);
      }
      const canRepeat = !leftover && (excess || (have != null && budget > 1e-9));
      if (canRepeat) {
        for (let i = 0; i < 6; i++) {
          const big = main.slice().sort((a, b) => b.length - a.length)[0];
          const pool = pick.group ? RR.rules.filter(r => r.group === pick.group) : [rule];
          const more = pool.map(r => evalRule(R, r, [big], reagentKey, P, presentAll)).filter(e => e.status === 'ok');
          if (!more.length) break;
          const nx = chooseCandidate(R, more), need = eqOf(nx.cand.rule);
          if (!excess && budget < need - 1e-9) { if (budget > 1e-9) amounts.push(`not enough ${rname} left for another ${nx.cand.rule.id} (needs ${need} equiv, ${+budget.toFixed(2)} left)`); break; }
          budget -= need;
          selection.push(`${excess ? 'excess' : 'more reagent'}: then ${nx.frag ? fragName(nx.frag) : 'the same reaction again'} (${nx.cand.rule.id})`);
          main = main.filter(x => x !== big).concat(nx.outcome.products);
          extra = extra.concat(nx.cand.rule.byproducts.map(k => RR.molecules[k].smiles));
          steps.push(nx.cand.rule.id);
        }
        if (!excess && have != null && steps.length > 1) amounts.push(`${steps.length} reactions used ${+(have - budget).toFixed(2)} of the ${have} equiv`);
      } else if (pick.others.length) {
        selection.push(`left untouched (less reactive): ${pick.others.map(fragName).join(', ')}. Add more reagent (e.g. "2 equiv") or "excess" to react them too`);
      }
      if (leftover) extra = extra.concat(leftover);
      // 'X reacts first' only means something when there was a choice
      if (steps.length === 1 && !pick.others.length) { const k = selection.findIndex(t => t.startsWith('group ')); if (k >= 0) selection.splice(k, 1); }
      const products = main.concat(extra);
      const prodRead = main.map(smiles => ({ smiles, frags: fragmentsOf(R, smiles) }));
      const made = new Set(prodRead.flatMap(x => x.frags));
      main.forEach(x => { const t = stereoText(R, x, chiralStart); if (t) stereo.push(`${x}: ${t}`); });
      const ranAt = P.refluxT != null ? ` · ran at ${fmtT(P.refluxT)} (reflux in ${P.solvents.join(' / ')})`
        : P.temps.length ? ` · ran at ${P.temps.map(fmtT).join(', ')}` : '';
      const reading = { subs: read, reagent: reagent ? { smiles: reagent, frags: reagentFrags } : null, fragChecks: fcheck.passed,
        acidBase: abCheck ? abCheck.text + ' ✓' : '', selection, stereo, rearr, amounts, rule, prods: prodRead,
        formsOk: (rule.forms || []).every(alt => alt.some(id => made.has(id))) };
      return { products, reading, mechanism: shiftSteps.concat(rule.mechanism || []), code: reactionCode(R, subs, reagent, products, condText),
        note: `${rule.id} · ${rule.name} (reaction_rules.xlsx)` + (steps.length > 1 ? ` · applied ${steps.length}×: ${steps.join(' → ')}` : '') + (leftover ? ' · partial (limiting reagent)' : '') +
          (rule.conditions ? ` · conditions: ${rule.conditions}` : '') + ranAt };
    }
    const hp = halidePathway(R, subs, read, reagent, reagentKey, P, condText, reagentFrags);
    if (hp) return hp;
    const reading = { subs: read, reagent: reagent ? { smiles: reagent, frags: reagentFrags } : null };
    if (short.length) return { reading, kind: 'blocked', note: 'No reaction under these conditions: ' + short.join('; ') + '.' };
    const presentSubs = new Set(read.flatMap(x => x.frags));
    const kn = knownNoReaction(R, reagentKey, presentSubs);
    if (kn.length) return { reading, kind: 'none', note: 'No reaction (known, Known_No_Reaction sheet): ' + kn.map(k => `${k.name || 'this reagent'} + ${fragName(k.frag)}: ${k.why}`).join('; ') };
    const learned = await learnedReaction(R, subs, reagent, reagentKey, condText, read, reagentFrags);
    if (learned && learned.products) return learned;
    // unknown: say so, and what the app does know about these fragments and this reagent
    const fr = [...presentSubs].filter(id => !['F_7GV2Y2WW'].includes(id) || presentSubs.size === 1);
    const forFrags = RR.rules.filter(r => (r.consumes || []).some(alt => alt.some(id => presentSubs.has(id)))).map(r => `${r.id} ${r.name}`);
    const forReagent = reagentKey ? RR.rules.filter(r => r.reagents.includes(reagentKey) || (r.acidBase && acidBaseEntry(R, r.acidBase.role, reagentKey))).map(r => `${r.id} ${r.name}`) : [];
    return { reading, kind: 'unknown', note: `Unknown: nothing in reaction_rules.xlsx covers ${fr.length ? fr.map(fragName).join(', ') : 'these molecules'} ${reagent ? 'with this reagent' : 'with no reagent'}, so the app cannot say whether they react (this is not the same as "no reaction").` +
      (learned && learned.miss ? ` Learned reactions (learned_reactions.xlsx): ${learned.miss}.` : '') +
      (forFrags.length ? ` Reactions it knows for these fragments: ${forFrags.slice(0, 8).join('; ')}${forFrags.length > 8 ? ' …' : ''}.` : '') +
      (forReagent.length ? ` This reagent is used in: ${forReagent.slice(0, 6).join('; ')}.` : '') };
  }
  /* The rule's mechanism (sheet "Mechanism" in reaction_rules.xlsx): one row per curved-arrow step, saying where the
     electrons come from and go to. Proton transfers show the donor's pKa and the pKa of the acid formed:
     ΔpKa = formed − donor, K ≈ 10^ΔpKa (above 1 the proton moves; below 1 only a small amount does). Each step draws
     the intermediate (or product) it makes. */
  const minus = x => String(x).replace('-', '−');
  function smallDrawing(smiles) {
    let g = null; try { g = window.Chem.parseSmiles(smiles); } catch (e) { g = null; }
    if (!g) return null;
    const laid = molTabLayout(smiles, g);
    const svg = laid ? window.MolDraw.svg(laid, 3000, 3000, { layout: false }) : window.MolDraw.svg(g, 3000, 3000);
    return svg;
  }
  /* The structure at the start of a step with its curved arrows drawn in red (Chem.render draws them: full head =
     an electron pair, fishhook = one electron). 'before' is atom-mapped SMILES; 'arrows' uses the Mechanism sheet's
     notation: 5>1 (lone pair on 5 to atom 1), 1-2>3 (bond 1–2 electrons onto atom 3), 1-2>2-4 (to a new bond 2–4),
     ~ in front for one electron. Each dot-separated species is laid out on its own and set side by side. */
  function arrowDrawing(before, arrows) {
    const G = { atoms: [], bonds: [], nextId: 1 }, byMap = new Map();
    let cursor = 0;
    for (const piece of before.split('.')) {
      let g = null; try { g = window.Chem.parseSmiles(piece); } catch (e) { g = null; }
      if (!g) return null;
      const laid = molTabLayout(piece, g) || g;
      const at = laid.atoms;
      const lens = laid.bonds.map(b => { const p = at.find(a => a.id === b.a), q = at.find(a => a.id === b.b); return Math.hypot(p.x - q.x, p.y - q.y); }).filter(x => x > 0);
      const k = lens.length ? 40 / (lens.reduce((s, x) => s + x, 0) / lens.length) : 1;
      const xs = at.map(a => a.x * k), ys = at.map(a => a.y * k);
      const minX = Math.min(...xs), maxX = Math.max(...xs), midY = (Math.min(...ys) + Math.max(...ys)) / 2;
      const ids = new Map();
      at.forEach((a, i) => {
        const n = Object.assign({}, a, { id: G.nextId++, x: xs[i] - minX + cursor, y: ys[i] - midY });
        ids.set(a.id, n.id); G.atoms.push(n);
        if (a.mapNum) byMap.set(a.mapNum, n);
      });
      laid.bonds.forEach(b => G.bonds.push(Object.assign({}, b, { a: ids.get(b.a), b: ids.get(b.b) })));
      cursor += (maxX - minX) + 95;                       // room between species for the arrows to cross
    }
    G.arrows = [];
    for (const t of String(arrows || '').split(';').map(x => x.trim()).filter(Boolean)) {
      const m = t.match(/^(~)?(\d+)(?:-(\d+))?>(\d+)(?:-(\d+))?$/);
      if (!m) continue;
      const A = byMap.get(+m[2]), B = m[3] ? byMap.get(+m[3]) : null, C = byMap.get(+m[4]), Dd = m[5] ? byMap.get(+m[5]) : null;
      if (!A || !C || (m[3] && !B) || (m[5] && !Dd)) continue;
      pushArrow(G, A, B, C, Dd, !!m[1]);
    }
    return arrowCanvas(G);
  }
  /* one curved arrow: from the lone pair on A (or the A–B bond) to atom C (or to a new C–D bond) */
  function pushArrow(G, A, B, C, Dd, fish) {
    const mid = (p, q) => ({ x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 });
    const s = B ? mid(A, B) : { x: A.x, y: A.y };
    let e = Dd ? mid(C, Dd) : { x: C.x, y: C.y };
    if (B && !Dd && (C === A || C === B)) {
      /* a bond's electrons moving onto one of its own atoms: curl out to the side of the bond (away from the
         rest of the molecule) and land beside the atom, so the short arrow is easy to see */
      const dx = B.x - A.x, dy = B.y - A.y, L = Math.hypot(dx, dy) || 1;
      let nx = -dy / L, ny = dx / L;
      const near = G.atoms.filter(a => a !== A && a !== B && Math.hypot(a.x - s.x, a.y - s.y) < 70);
      if (near.length) {
        const cx = near.reduce((t, a) => t + a.x, 0) / near.length, cy = near.reduce((t, a) => t + a.y, 0) / near.length;
        if ((cx - s.x) * nx + (cy - s.y) * ny > 0) { nx = -nx; ny = -ny; }
      }
      const st = { x: s.x + nx * 5, y: s.y + ny * 5 };
      const en = { x: C.x + (s.x - C.x) * 0.2 + nx * 15, y: C.y + (s.y - C.y) * 0.2 + ny * 15 };
      const side = (en.x - st.x) * ny - (en.y - st.y) * nx;          // bulge outward, away from the bond
      G.arrows.push({ x1: st.x, y1: st.y, x2: en.x, y2: en.y, fromBond: true, autoSide: false,
                      bend: side > 0 ? 14 : -14, type: fish ? 'fishhook' : 'full' });
      return;
    }
    if (!Dd) {                                            // stop short of the atom label
      const dx = s.x - e.x, dy = s.y - e.y, L = Math.hypot(dx, dy) || 1;
      e = { x: e.x + dx / L * 11, y: e.y + dy / L * 11 };
    }
    if (!B && fish) {                                     // a single electron leaves from beside its atom
      const dx = e.x - s.x, dy = e.y - s.y, L = Math.hypot(dx, dy) || 1;
      s.x += dx / L * 10; s.y += dy / L * 10;
    }
    G.arrows.push({ x1: s.x, y1: s.y, x2: e.x, y2: e.y, fromBond: !!B, type: fish ? 'fishhook' : 'full' });
  }
  function arrowCanvas(G) {
    const xs = G.atoms.map(a => a.x), ys = G.atoms.map(a => a.y);
    const S = 0.9, pad = 30;
    const canvas = document.createElement('canvas');
    canvas.className = 'rxn-step-canvas';
    canvas.style.width = Math.round((Math.max(...xs) - Math.min(...xs)) * S + pad * 2) + 'px';
    canvas.style.height = Math.round((Math.max(...ys) - Math.min(...ys)) * S + pad * 2) + 'px';
    return { canvas, draw: () => window.Chem.renderFit(canvas, G, { pad, maxScale: S }) };
  }
  /* a step from the app's SN1/SN2/E1/E2 engine (Chem.buildMechanism): its own graph, arrows given as
     {from: {kind:'lp', atom} | {kind:'bond', a, b}, to: {kind:'atom', atom} | {kind:'bond', a, b}, half} */
  function graphArrowDrawing(graph, arrows) {
    const lens = graph.bonds.map(b => { const p = graph.atoms.find(a => a.id === b.a), q = graph.atoms.find(a => a.id === b.b); return p && q ? Math.hypot(p.x - q.x, p.y - q.y) : 0; }).filter(x => x > 0);
    const k = lens.length ? 40 / (lens.reduce((t, x) => t + x, 0) / lens.length) : 1;
    const G = { atoms: graph.atoms.map(a => Object.assign({}, a, { x: a.x * k, y: a.y * k })), bonds: graph.bonds.map(b => Object.assign({}, b)), nextId: graph.nextId, arrows: [] };
    const at = id => G.atoms.find(a => a.id === id);
    for (const ar of arrows || []) {
      const f = ar.from || {}, t = ar.to || {};
      const A = at(f.kind === 'bond' ? f.a : f.atom), B = f.kind === 'bond' ? at(f.b) : null;
      const C = at(t.kind === 'bond' ? t.a : t.atom), Dd = t.kind === 'bond' ? at(t.b) : null;
      if (A && C) pushArrow(G, A, B, C, Dd, !!ar.half);
    }
    return arrowCanvas(G);
  }
  /* the acidic and basic sites of a molecule with their pKa in this molecule and where each number comes from */
  function pkaLine(smiles) {
    if (!RD) return '';
    let P = null; try { P = moleculePkas(RD, smiles); } catch (e) { P = null; }
    if (!P || (!P.acid.length && !P.base.length)) return '';
    const item = (s_, kind) => `${esc([...new Set(s_.fragments)].map(fragName).join(' / '))} ${kind === 'acid' ? 'pKa' : 'pKaH'} ${esc(minus(s_.pKa))} <span class="rxn-step-k">(${esc(s_.source)})</span>`;
    return `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">pKa in this molecule:</span> ${P.acid.map(s_ => item(s_, 'acid')).concat(P.base.map(s_ => item(s_, 'base'))).join('; ')}</div>`;
  }
  /* steric and electronic effects of each fragment of a starting material, at its reacting atom */
  function effectsLine(smiles) {
    if (!RD) return '';
    let E = []; try { E = moleculeEffects(RD, smiles); } catch (e) { E = []; }
    const seen = new Set(), items = [];
    for (const e of E) {
      if (e.fragment === 'F_7GV2Y2WW') continue;                         // every arene C-H: too many to list
      const k = e.fragment; if (seen.has(k)) continue; seen.add(k);
      items.push(`${esc(fragName(e.fragment))}: ${esc(e.steric.cls)}; ${esc(e.electronic.text)}`);
    }
    return items.length ? `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Steric and electronic effects in this molecule:</span> ${items.join(' · ')}</div>` : '';
  }
  function showReading(rd) {
    if (!rd) return;
    const wrap = document.createElement('div');
    wrap.className = 'rxn-read';
    const fr = list => list.length ? list.map(id => `<span class="rxn-frag">${esc(fragName(id))}</span>`).join(' ') : '<span class="rxn-step-k">no fragment from the sheet</span>';
    let html = '<div class="rxn-mech-head">How it was read (fragments, molecule_data.xlsx)</div>';
    html += rd.subs.map(x => `<div class="rxn-step-line"><span class="rxn-read-smi">${esc(x.smiles)}</span> ${fr(x.frags)}</div>` + pkaLine(x.smiles) + effectsLine(x.smiles)).join('');
    if (rd.reagent && rd.reagent.frags.length) html += `<div class="rxn-step-line"><span class="rxn-step-k">Reagent:</span> <span class="rxn-read-smi">${esc(rd.reagent.smiles)}</span> ${fr(rd.reagent.frags)}</div>`;
    const RG = (window.REACTION_RULES || {}).reagents || {};
    if (rd.reagent) {                                      // what each reagent part is (Reagents sheet, reaction_rules.xlsx)
      const named = rd.reagent.smiles.split('.').map(x => RG[x] ? `${RG[x].abbr || RG[x].name} (${RG[x].role}${RG[x].also ? ', ' + RG[x].also : ''})` : null).filter(Boolean);
      if (named.length) html += `<div class="rxn-step-line"><span class="rxn-step-k">Reagent parts (Reagents sheet):</span> ${esc([...new Set(named)].join('; '))}</div>`;
    }
    if ((rd.amounts || []).length) html += rd.amounts.map(t => `<div class="rxn-step-line rxn-step-pka ${/only about|not enough/.test(t) ? 'up' : 'down'}"><span class="rxn-step-k">Amounts:</span> ${esc(t)}</div>`).join('');
    if ((rd.rearr || []).length) html += rd.rearr.map(t => `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Rearrangement:</span> ${esc(t)}</div>`).join('');
    if ((rd.stereo || []).length) html += rd.stereo.map(t => `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Stereochemistry:</span> ${esc(t)}</div>`).join('');
    if ((rd.selection || []).length) html += rd.selection.map(t => `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Chemoselectivity:</span> ${esc(t)}</div>`).join('');
    if (rd.pathway) html += `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Pathway:</span> ${esc(rd.pathway)}</div>`;
    if (rd.learned) html += `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Learned:</span> ${esc(rd.learned)}</div>`;
    if (rd.acidBase) html += `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">${esc(rd.rule.acidBase.role === 'base' ? 'Base' : 'Acid')} strength (pKa):</span> ${esc(rd.acidBase)}</div>`;
    if ((rd.fragChecks || []).length) html += rd.fragChecks.map(c => `<div class="rxn-step-line rxn-step-pka down"><span class="rxn-step-k">Fragment conditions:</span> ${esc(c)}</div>`).join('');
    if (rd.rule) {
      if (!rd.pathway) html += `<div class="rxn-step-line"><span class="rxn-step-k">Rule:</span> ${esc(rd.rule.id)} consumes ${esc(itemsText(rd.rule.consumes || []) || 'no listed fragment')}` +
        ` → forms ${esc(itemsText(rd.rule.forms || []) || 'no new functional group (C–C/C–H only)')}</div>`;
      html += rd.prods.map(x => `<div class="rxn-step-line"><span class="rxn-read-smi">${esc(x.smiles)}</span> ${fr(x.frags)}</div>`).join('');
      if (!rd.formsOk) html += '<div class="rxn-step-line rxn-step-pka up">The product does not contain the fragment this rule should form: check the rule in reaction_rules.xlsx.</div>';
    }
    wrap.innerHTML = html;
    ruleNote.appendChild(wrap);
  }
  function showMechanism(steps, hasProducts, status, note) {
    if (!steps || !steps.length) {
      if (hasProducts) {
        const w = document.createElement('div'); w.className = 'rxn-mech';
        w.innerHTML = '<div class="rxn-mech-head">Mechanism: electron flow</div><div class="rxn-step-note">' + esc(note || 'The mechanism for this reaction is not in the workbook yet (Mechanism sheet). The product is still predicted from the Reactions row.') + '</div>';
        ruleNote.appendChild(w);
      }
      return;
    }
    const wrap = document.createElement('div');
    wrap.className = 'rxn-mech';
    wrap.innerHTML = '<div class="rxn-mech-head">Mechanism: electron flow</div>' + (status ? `<div class="rxn-mech-status ${/^verified/.test(status) ? 'ok' : 'bad'}">${/^verified/.test(status) ? '✓ ' : ''}${esc(status)}${/^verified/.test(status) ? ': every step\'s arrows were applied (bonds, charges, lone pairs) and give exactly that step\'s result' : ''}</div>` : '');
    ruleNote.appendChild(wrap);                          // in the page first, so the drawings can be measured
    for (const st of steps) {
      const row = document.createElement('div'); row.className = 'rxn-step';
      const lines = [`<div class="rxn-step-title">${st.step}. ${esc(st.type)}</div>`];
      if (st.from || st.to) lines.push(`<div class="rxn-step-line"><span class="rxn-step-k">Electrons:</span> ${esc(st.from)} → ${esc(st.to)}</div>`);
      if (st.arrowLabels && st.arrowLabels.length) lines.push(`<div class="rxn-step-line"><span class="rxn-step-k">Arrows:</span> ${esc(st.arrowLabels.join('; '))}</div>`);
      if (st.made || st.broken) lines.push(`<div class="rxn-step-line"><span class="rxn-step-k">Bond made:</span> ${esc(st.made || '—')} · <span class="rxn-step-k">broken:</span> ${esc(st.broken || '—')}</div>`);
      if (st.dpka != null) {
        const up = st.dpka < 0;
        lines.push(`<div class="rxn-step-line rxn-step-pka ${up ? 'up' : 'down'}"><span class="rxn-step-k">pKa:</span> ${esc(st.donor)} (${minus(st.donorPka)}) gives H⁺ → ${esc(st.formed)} (${minus(st.formedPka)}) · ΔpKa ${st.dpka > 0 ? '+' : ''}${minus(+st.dpka.toFixed(1))}, K ≈ 10<sup>${minus(+st.dpka.toFixed(1))}</sup> · ${up ? 'uphill (small amount)' : 'downhill (favoured)'}</div>`);
      }
      if (st.lg) lines.push(`<div class="rxn-step-line"><span class="rxn-step-k">Leaving group:</span> ${esc(st.lg)}` + (st.lgPka != null ? ` (pKa of its conjugate acid ${minus(st.lgPka)}: the lower, the better it leaves)` : '') + '</div>');
      if (st.note) lines.push(`<div class="rxn-step-note">${esc(st.note)}</div>`);
      row.innerHTML = `<div class="rxn-step-text">${lines.join('')}</div>`;
      // the start of the step with its red curved arrows, then (→) what the step makes
      let drawArrows = null;
      if (st.graph) {
        const ad = graphArrowDrawing(st.graph, st.gArrows);
        const box = document.createElement('div'); box.className = 'rxn-step-before';
        box.appendChild(ad.canvas);
        row.appendChild(box);
        row.insertAdjacentHTML('beforeend', '<div class="rxn-step-go">→</div>');
        drawArrows = ad.draw;
      } else if (st.before) {
        const ad = arrowDrawing(st.before, st.arrows);
        if (ad) {
          const box = document.createElement('div'); box.className = 'rxn-step-before';
          box.appendChild(ad.canvas);
          row.appendChild(box);
          row.insertAdjacentHTML('beforeend', '<div class="rxn-step-go">→</div>');
          drawArrows = ad.draw;
        }
      }
      if (st.smiles) {
        const mol = document.createElement('div'); mol.className = 'rxn-step-mol';
        row.appendChild(mol);
        const pics = document.createElement('div'); pics.className = 'rxn-step-pics';
        mol.appendChild(pics);
        const svgs = [];
        st.smiles.split('.').forEach((piece, i) => {
          if (i) pics.insertAdjacentHTML('beforeend', '<span class="rxn-step-plus">+</span>');
          const svg = smallDrawing(piece);
          if (svg) { pics.appendChild(svg); svgs.push(svg); }
        });
        const smi = document.createElement('div'); smi.className = 'rxn-step-smi'; smi.textContent = st.smiles; smi.title = 'Click to copy';
        smi.addEventListener('click', () => {
          const done = () => { smi.textContent = 'Copied'; setTimeout(() => { smi.textContent = st.smiles; }, 900); };
          if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(st.smiles).then(done, () => fallbackCopy(st.smiles, done));
          else fallbackCopy(st.smiles, done);
        });
        mol.appendChild(smi);
        wrap.appendChild(row);
        for (const svg of svgs) {                        // fit each drawing once it is in the page (getBBox needs that)
          const bb = svg.getBBox(), pad = 6, k = 0.6 * BOND_PX / 30;
          svg.setAttribute('viewBox', `${bb.x - pad} ${bb.y - pad} ${bb.width + 2 * pad} ${bb.height + 2 * pad}`);
          svg.setAttribute('width', (bb.width + 2 * pad) * k);
          svg.setAttribute('height', (bb.height + 2 * pad) * k);
        }
      } else wrap.appendChild(row);
      if (drawArrows) drawArrows();                       // the canvas is in the page now, so it has a size
    }
    relayout();
  }
  const products = document.getElementById('rxn-products');
  products.style.setProperty('--rxn-gap', GAP + 'px');
  const updateProdLine = () => { prodInput.value = products.classList.contains('on') ? orderedSmiles(pslots) : ''; fitProduct(); updateCode(); };
  const pslots = makeList(products, 'Products come from the reaction rules: click the arrow', updateProdLine, true);
  document.getElementById('rxn-arrow').addEventListener('click', async () => {
    const on = !products.classList.contains('on');
    if (!on) currentRun++;                               // closing the arrow cancels a run still going
    products.classList.toggle('on', on);
    ruleNote.style.display = on ? '' : 'none';
    if (on) {
      const myRun = ++currentRun;
      ruleNote.textContent = 'Running the reaction rules…';
      let res;
      try { res = await predict(); } catch (e) { res = { note: 'Something went wrong while running the rules: ' + e.message }; }
      if (myRun !== currentRun || !products.classList.contains('on')) return;     // a newer run (or a closed arrow) took over
      ruleNote.innerHTML = '';
      const verdict = document.createElement('div');
      verdict.className = 'rxn-verdict' + (res.kind ? ' kind-' + res.kind : '');
      verdict.textContent = res.note;
      ruleNote.appendChild(verdict);
      showReading(res.reading);
      showMechanism(res.mechanism, !!(res.products && res.products.length), res.reading && res.reading.rule && res.reading.rule.mechStatus, res.mechNote);
      await pslots.setAll(res.products || []);           // no reaction: no products left over from an earlier run
    }
    updateProdLine(); relayout();
  });
  const above = document.getElementById('rxn-above'), below = document.getElementById('rxn-below');
  above.style.setProperty('--rxn-gap', GAP + 'px'); below.style.setProperty('--rxn-gap', GAP + 'px');
  above.title = 'Add reagents';
  below.title = 'Add conditions';
  const reagLine = document.getElementById('rxn-reag');
  const setReagLine = smi => { reagLine.value = smi || ''; reagLine.classList.remove('bad'); reagLine.title = ''; fitLine(reagLine); updateCode(); };
  const aboveApi = wireBox(above, 'rxn-slot rxn-agent',      // reagents: molecules, typed as SMILES
    smi => setReagLine(smi),
    () => { aboveApi.clear(); setReagLine(''); },
    reagentMenu);
  reagLine.addEventListener('keydown', async e => {
    if (e.key !== 'Enter') return;
    e.preventDefault();
    const t = reagLine.value.trim();
    if (!t) { aboveApi.clear(); setReagLine(''); reagLine.blur(); return; }
    let g = null; try { g = window.Chem.parseSmiles(t); } catch (err) { g = null; }
    if (!g) { reagLine.classList.add('bad'); reagLine.title = 'Not a valid SMILES: ' + t; return; }
    reagLine.blur();
    await aboveApi.show(t, g);
  });
  reagLine.addEventListener('input', () => { reagLine.classList.remove('bad'); reagLine.title = ''; fitLine(reagLine); });
  reagLine.addEventListener('focus', () => { if (reagLine.value) reagLine.select(); });
  new ResizeObserver(() => fitLine(reagLine)).observe(reagLine);
  const belowApi = wireText(below, 'rxn-slot rxn-agent', conditionOptions);     // conditions: free text (heat, hv, 25 C, solvent ...), not molecules
  watcher.observe(products); watcher.observe(above); watcher.observe(below);
  relayout();

  /* The reaction code bar (its own strip under the tabs, full width). It always shows the page's reaction as
     "starting materials~reagent~conditions~products", built with reactionCode. Typing or pasting a code and pressing
     Enter sets the boxes on the left, the reagent above the arrow and the conditions under it (an empty part clears
     that part). The products part is ignored: products are only drawn by clicking the arrow, which runs the rules. */
  const codeBar = document.getElementById('rxn-code');
  /* Settings: Basic shows the products and a one-line verdict; Advanced also shows how it was read, the checks and
     the mechanism. Both use the same workbook data. The gear is a square as tall as the code bar. */
  const gear = document.getElementById('rxn-gear'), settings = document.getElementById('rxn-settings');
  function sizeGear() {
    const h = codeBar.getBoundingClientRect().height;       // exact (offsetHeight rounds to whole pixels)
    if (!h) return;
    gear.style.width = gear.style.height = h + 'px';
    const cb = document.getElementById('rxn-common-btn');                 // the Common reactions button sits left of the gear
    cb.style.height = h + 'px'; cb.style.right = (10 + h + 6) + 'px';
    codeBar.style.right = (10 + h + 6 + cb.getBoundingClientRect().width + 6) + 'px';
  }
  new ResizeObserver(sizeGear).observe(codeBar);
  sizeGear();
  // the page may still be hidden at start-up (no size yet), so try again once it can be measured
  [0, 150, 600].forEach(ms => setTimeout(sizeGear, ms));
  window.addEventListener('load', sizeGear);
  window.addEventListener('hashchange', () => setTimeout(sizeGear, 50));
  let mode = 'advanced';
  try { mode = localStorage.getItem('rxn-mode') || 'advanced'; } catch (e) { mode = 'advanced'; }
  function setMode(v) {
    mode = v === 'basic' ? 'basic' : 'advanced';
    page.classList.toggle('mode-basic', mode === 'basic');
    settings.querySelectorAll('input[name="rxn-mode"]').forEach(r => { r.checked = r.value === mode; });
    gear.title = 'Settings: ' + (mode === 'basic' ? 'Basic' : 'Advanced') + ' mode';
    try { localStorage.setItem('rxn-mode', mode); } catch (e) { /* not stored */ }
    relayout();
  }
  setMode(mode);
  gear.addEventListener('click', e => { e.stopPropagation(); settings.hidden = !settings.hidden; gear.classList.toggle('on', !settings.hidden); });
  settings.addEventListener('change', e => { if (e.target.name === 'rxn-mode') setMode(e.target.value); });
  document.addEventListener('mousedown', e => { if (!settings.hidden && !settings.contains(e.target) && !gear.contains(e.target)) { settings.hidden = true; gear.classList.remove('on'); } }, true);
  let codeSeq = 0;
  async function updateCode() {
    const seq = ++codeSeq, R = await loadRD();
    if (seq !== codeSeq || document.activeElement === codeBar) return;    // a newer update came in, or the user is typing
    const subs = slots.filter(x => x.smiles).map(x => x.smiles);
    const prods = products.classList.contains('on') ? pslots.filter(x => x.smiles).map(x => x.smiles) : [];
    const reagent = aboveApi.get(), cond = belowApi.get();
    codeBar.value = !subs.length && !reagent && !cond && !prods.length ? ''
      : R ? reactionCode(R, subs, reagent, prods, cond) : [fullInput.value, reagLine.value, condCode(cond), prodInput.value].join('~');
    codeBar.classList.remove('bad'); codeBar.title = '';
  }
  codeBar.addEventListener('focus', () => { if (codeBar.value) codeBar.select(); });
  codeBar.addEventListener('input', () => { codeBar.classList.remove('bad'); codeBar.title = ''; });
  codeBar.addEventListener('blur', () => updateCode());
  let codeBarBusy = Promise.resolve();                 // filling the boxes from a code (the test runner waits for it)
  codeBar.addEventListener('keydown', e => {
    if (e.key === 'Escape') { codeBar.blur(); return; }
    if (e.key !== 'Enter') return;
    e.preventDefault();
    codeBarBusy = fillFromCode();
  });
  async function fillFromCode() {
    const t = codeBar.value.trim();
    const parts = t.split('~').map(x => x.trim());
    if (!t || parts.length > 4) { codeBar.classList.add('bad'); codeBar.title = 'A reaction code is starting materials~reagent~conditions~products (at most three ~)'; return; }
    while (parts.length < 4) parts.push('');
    const cond = parts.splice(2, 1)[0];
    const [subs, reag] = parts.slice(0, 2).map(x => x.split('.').map(y => y.trim()).filter(Boolean));   // the products part is ignored
    const bad = [...subs, ...reag].find(q => { try { return !window.Chem.parseSmiles(q); } catch (err) { return true; } });
    if (bad) { codeBar.classList.add('bad'); codeBar.title = 'Not a valid SMILES: ' + bad; return; }
    codeBar.blur();
    await slots.setAll(subs);
    if (reag.length) { const r = reag.join('.'); await aboveApi.show(r, window.Chem.parseSmiles(r)); }
    else { aboveApi.clear(); setReagLine(''); }
    belowApi.set(cond);
    products.classList.remove('on');                    // products only ever come from running the rules (the arrow)
    ruleNote.style.display = 'none';
    await pslots.setAll([]);
    updateProdLine(); relayout();
  }

  /* ---- Run all tests (gear menu) ----
     Every reaction in the workbook is run from its stored reaction code (starting materials~reagent~conditions~) and must
     give the products recorded in that code, fire its own rule and (if it has a mechanism) show it as verified. The
     feature cases in tools/reaction_tests.json check each decision layer. Everything goes through this page exactly as a
     click would; whatever was on the page is put back afterwards. */
  const testsPanel = document.getElementById('rxn-tests');
  let testRun = null, lastResults = [];
  async function loadTestCases() {
    const RR = window.REACTION_RULES || {}, cases = [];
    for (const r of RR.rules || []) {
      const parts = String(r.code || '').split('~');
      if (parts.length < 4 || !parts[0]) continue;
      cases.push({ group: 'reactions (workbook)', name: `${r.id} ${r.name}`, code: parts.slice(0, 3).join('~') + '~',
                   expect: r.id, products: parts[3].split('.').filter(Boolean), absent: [], mech: (r.mechanism || []).length > 0 });
    }
    try {
      const T = await (await fetch('tools/reaction_tests.json', { cache: 'no-store' })).json();
      for (const f of T.features || []) cases.push({ group: f.feature, name: f.name, code: f.code, expect: f.expect, products: f.products || [], absent: f.absent || [], mech: false });
    } catch (e) { /* no feature file: only the workbook reactions are tested */ }
    try {                                                // a sample of the learned reactions (tools/learned_to_js.py writes it)
      const T = await (await fetch('tools/learned_tests.json', { cache: 'no-store' })).json();
      for (const f of T.learned || []) cases.push({ group: f.feature, name: f.name, code: f.code, expect: f.expect, products: f.products || [], absent: f.absent || [], mech: false });
    } catch (e) { /* no learned reactions built */ }
    return cases;
  }
  const waitFor = async (test, ms) => { const t0 = Date.now(); while (Date.now() - t0 < ms) { if (test()) return true; await new Promise(r => setTimeout(r, 50)); } return false; };
  async function runOneCase(R, t) {
    const canon = x => { const m = R.get_mol(x); if (!m) return x; const c = m.get_smiles(); m.delete(); return c; };
    if (products.classList.contains('on')) { document.getElementById('rxn-arrow').dispatchEvent(new MouseEvent('click', { bubbles: true })); await new Promise(r => setTimeout(r, 80)); }
    codeBar.value = t.code;
    codeBar.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    await codeBarBusy;                                                // the boxes are filled
    await new Promise(r => setTimeout(r, 150));
    document.getElementById('rxn-arrow').dispatchEvent(new MouseEvent('click', { bubbles: true }));
    await waitFor(() => { const v = ruleNote.firstChild; return v && v.textContent && !v.textContent.startsWith('Running'); }, 60000);
    await new Promise(r => setTimeout(r, 250));                       // drawings settle
    const note = (ruleNote.firstChild && ruleNote.firstChild.textContent) || '';
    const got = prodInput.value.split('.').filter(Boolean);
    const exp = t.expect + (/^(R\d{3}|L\d{6}|SN1|SN2|E1|E2)$/.test(t.expect) ? ' ' : '');
    const why = [];
    if (!note.startsWith(exp) && !(/^L\d+$/.test(t.expect) && note.includes(t.expect))) why.push(`expected "${t.expect}…", got "${note.slice(0, 90)}"`);
    const missing = t.products.filter(p => !got.includes(canon(p))), extra = t.absent.filter(p => got.includes(canon(p)));
    if (missing.length) why.push(`missing product ${missing.join(', ')}; page gave ${got.join('.') || 'nothing'}`);
    if (extra.length) why.push(`should not contain ${extra.join(', ')}`);
    if (t.mech) { const st = (ruleNote.querySelector('.rxn-mech-status') || {}).textContent || ''; if (!st.startsWith('✓')) why.push('mechanism not shown as verified'); }
    return { t, ok: !why.length, why };
  }
  function renderTests(results, total, running) {
    const done = results.length, fails = results.filter(x => !x.ok);
    document.getElementById('rxn-t-sum').textContent = running ? `running ${done} of ${total}…` : `${done - fails.length} of ${done} passed` + (done < total ? ` (stopped at ${done} of ${total})` : '');
    const fill = document.getElementById('rxn-t-fill');
    fill.style.width = (total ? 100 * done / total : 0) + '%'; fill.classList.toggle('bad', fails.length > 0);
    const groups = {};
    results.forEach(x => { const g = groups[x.t.group] = groups[x.t.group] || [0, 0]; g[1]++; if (x.ok) g[0]++; });
    document.getElementById('rxn-t-groups').innerHTML = Object.entries(groups).map(([g, [ok, n]]) =>
      `<div class="rxn-t-g"><span>${esc(g)}</span><span class="${ok === n ? 'ok' : 'bad'}">${ok} / ${n}</span></div>`).join('');
    const list = document.getElementById('rxn-t-list');
    list.innerHTML = fails.length ? fails.map((x, i) => `<div class="rxn-t-f" data-i="${i}"><div class="n">✗ ${esc(x.t.name)}</div><div class="d"><code>${esc(x.t.code)}</code></div>${x.why.map(w => `<div class="d">${esc(w)}</div>`).join('')}</div>`).join('')
      : (running ? '' : '<div class="rxn-t-done">All passed ✓</div>');
    list.querySelectorAll('.rxn-t-f').forEach(el => el.addEventListener('click', async () => {      // load that case on the page
      const x = fails[+el.dataset.i]; await runOneCase(await loadRD(), x.t);
    }));
    document.getElementById('rxn-t-again').disabled = document.getElementById('rxn-t-failed').disabled = !!running;
    document.getElementById('rxn-t-failed').disabled = !!running || !fails.length;
    document.getElementById('rxn-t-stop').disabled = !running;
  }
  async function runTests(onlyFailed) {
    if (testRun) return;
    const R = await loadRD();
    if (!R) return;
    settings.hidden = true; gear.classList.remove('on'); testsPanel.hidden = false;
    commonPanel.hidden = true; commonBtn.classList.remove('on'); commonBtn.setAttribute('aria-expanded', 'false');
    let cases = await loadTestCases();
    if (onlyFailed) { const names = new Set(lastResults.filter(x => !x.ok).map(x => x.t.name)); cases = cases.filter(c => names.has(c.name)); }
    const saved = codeBar.value, wasOn = products.classList.contains('on');
    testRun = { stop: false };
    const results = [];
    renderTests(results, cases.length, true);
    for (const t of cases) {
      if (testRun.stop) break;
      let r; try { r = await runOneCase(R, t); } catch (e) { r = { t, ok: false, why: ['error: ' + e.message] }; }
      results.push(r); renderTests(results, cases.length, true);
    }
    testRun = null;
    lastResults = onlyFailed ? lastResults.map(x => results.find(y => y.t.name === x.t.name) || x) : results;
    renderTests(onlyFailed ? lastResults : results, onlyFailed ? lastResults.length : cases.length, false);
    // put the page back as it was
    if (products.classList.contains('on')) document.getElementById('rxn-arrow').dispatchEvent(new MouseEvent('click', { bubbles: true }));
    if (saved) { codeBar.value = saved; codeBar.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true })); await codeBarBusy; if (wasOn) document.getElementById('rxn-arrow').dispatchEvent(new MouseEvent('click', { bubbles: true })); }
  }
  document.getElementById('rxn-run-tests').addEventListener('click', () => runTests(false));
  document.getElementById('rxn-t-again').addEventListener('click', () => runTests(false));
  document.getElementById('rxn-t-failed').addEventListener('click', () => runTests(true));
  document.getElementById('rxn-t-stop').addEventListener('click', () => { if (testRun) testRun.stop = true; });
  document.getElementById('rxn-t-close').addEventListener('click', () => { if (testRun) testRun.stop = true; testsPanel.hidden = true; });
  window.ReactionsTests = { run: runTests, results: () => lastResults, running: () => !!testRun };

  /* Common reactions: the workbook's reactions (window.REACTION_RULES.rules) as a searchable list. Clicking one puts its
     reaction code into the code bar exactly as typing it and pressing Enter would, then clicks the arrow, so the page
     predicts the products itself from the rules (the stored products are left out, as the test runner leaves them out).
     The families come from the first fragment each rule consumes (FAMILIES maps fragment names to a family); a rule
     whose fragment is in no family is listed under "Other". */
  const commonBtn = document.getElementById('rxn-common-btn'), commonPanel = document.getElementById('rxn-common');
  const commonQ = document.getElementById('rxn-c-q'), commonList = document.getElementById('rxn-c-list'), commonSum = document.getElementById('rxn-c-sum');
  const FAMILIES = [
    ['Alkenes and alkynes', ['alkene', 'alkyne', 'conjugated diene', 'allylic CH2']],
    ['Alkyl halides', ['alkyl bromide', 'alkyl chloride', 'alkyl iodide', 'alkyl halide', 'methyl halide', 'vicinal dibromide', 'benzylic bromide', 'halohydrin']],
    ['Alcohols, ethers and epoxides', ['alcohol', '1,2-diol', 'epoxide', 'ether', 'acetal']],
    ['Aldehydes and ketones', ['aldehyde', 'ketone', 'enone', 'oxime', 'imine']],
    ['Carboxylic acids and their derivatives', ['carboxylic acid', 'ester', 'amide', 'acyl chloride', 'carboxylic anhydride', 'nitrile']],
    ['Aromatic compounds', ['arene', 'aryl', 'phenol', 'nitroarene', 'alkyl arene']],
    ['Organometallic reagents', ['Grignard', 'aluminium hydride']],
  ];
  const RRC = window.REACTION_RULES || {};
  const SPECTATORS = new Set(['counter-ion / spectator ion', 'salt (work-up or additive)']);
  const familyOf = rule => {
    const first = ((rule.consumes || [])[0] || [])[0];
    const nm = first && RRC.fragments && RRC.fragments[first] ? RRC.fragments[first].name : '';
    const f = FAMILIES.find(([, pre]) => pre.some(p => nm.startsWith(p)));
    return f ? f[0] : 'Other';
  };
  const reagentText = reagentPart => {                          // the reagent part of the code as names (abbreviations), spectator ions left out
    const out = [];
    for (const c of String(reagentPart || '').split('.').filter(Boolean)) {
      const e = (RRC.reagents || {})[c];
      if (e && SPECTATORS.has(e.role)) continue;
      const t = e ? (e.abbr || e.name) : (c === 'O' ? 'H2O' : c);
      if (!out.includes(t)) out.push(t);
    }
    return out.join(' + ');
  };
  const COMMON = (RRC.rules || []).filter(r => { const p = String(r.code || '').split('~'); return p.length >= 4 && p[0]; }).map(r => {
    const parts = r.code.split('~'), fam = familyOf(r), reag = reagentText(parts[1]), cond = parts[2] || '';
    const frs = [].concat(...(r.consumes || []), ...(r.forms || [])).map(id => ((RRC.fragments || {})[id] || {}).name || '');
    return { r, parts, fam, sub: [reag, cond].filter(Boolean).join(' · '), text: [r.id, r.name, fam, reag, cond, ...frs].join(' ').toLowerCase() };
  });
  const COMMON_ORDER = FAMILIES.map(f => f[0]).concat('Other');
  function renderCommon() {
    const q = commonQ.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    const hits = COMMON.filter(c => q.every(t => c.text.includes(t)));
    commonSum.textContent = hits.length === COMMON.length ? COMMON.length + ' reactions' : hits.length + ' of ' + COMMON.length;
    commonList.textContent = '';
    if (!hits.length) { const d = document.createElement('div'); d.className = 'rxn-c-none'; d.textContent = 'No reaction matches. Try a shorter word.'; commonList.appendChild(d); return; }
    for (const fam of COMMON_ORDER) {
      const g = hits.filter(c => c.fam === fam);
      if (!g.length) continue;
      const h = document.createElement('div'); h.className = 'rxn-c-h';
      const a = document.createElement('span'); a.textContent = fam; const b = document.createElement('span'); b.textContent = g.length;
      h.append(a, b); commonList.appendChild(h);
      for (const c of g) {
        const it = document.createElement('button'); it.type = 'button'; it.className = 'rxn-c-item'; it.dataset.id = c.r.id;
        const i = document.createElement('span'); i.className = 'i'; i.textContent = c.r.id;
        const n = document.createElement('span'); n.className = 'n'; n.textContent = c.r.name;
        const d = document.createElement('span'); d.className = 'd'; d.textContent = c.sub;
        it.append(i, n, d); commonList.appendChild(it);
      }
    }
  }
  function toggleCommon(open) {
    if (testRun) return;                                        // the tests own the page while they run
    const show = open === undefined ? commonPanel.hidden : open;
    commonPanel.hidden = !show; commonBtn.classList.toggle('on', show); commonBtn.setAttribute('aria-expanded', String(show));
    if (show) { settings.hidden = true; gear.classList.remove('on'); testsPanel.hidden = true; renderCommon(); commonQ.focus(); commonQ.select(); }
  }
  let commonBusy = false;
  async function loadCommon(id) {
    const c = COMMON.find(x => x.r.id === id);
    if (!c || commonBusy) return;
    commonBusy = true;
    try {
      toggleCommon(false);
      const arrow = document.getElementById('rxn-arrow');
      if (products.classList.contains('on')) { arrow.dispatchEvent(new MouseEvent('click', { bubbles: true })); await new Promise(r => setTimeout(r, 80)); }
      codeBar.value = c.parts.slice(0, 3).join('~') + '~';
      codeBar.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      await codeBarBusy;                                        // the boxes are filled
      await new Promise(r => setTimeout(r, 150));
      arrow.dispatchEvent(new MouseEvent('click', { bubbles: true }));   // run the rules
    } finally { commonBusy = false; }
  }
  commonBtn.addEventListener('click', e => { e.stopPropagation(); toggleCommon(); });
  document.getElementById('rxn-open-common').addEventListener('click', () => toggleCommon(true));   // the same list, from the settings menu
  document.getElementById('rxn-c-close').addEventListener('click', () => toggleCommon(false));
  commonQ.addEventListener('input', renderCommon);
  commonQ.addEventListener('keydown', e => {
    if (e.key === 'Escape') { toggleCommon(false); commonBtn.focus(); }
    else if (e.key === 'Enter') { const first = commonList.querySelector('.rxn-c-item'); if (first) loadCommon(first.dataset.id); }
    else if (e.key === 'ArrowDown') { const first = commonList.querySelector('.rxn-c-item'); if (first) { e.preventDefault(); first.focus(); } }
  });
  commonList.addEventListener('click', e => { const b = e.target.closest('.rxn-c-item'); if (b) loadCommon(b.dataset.id); });
  commonList.addEventListener('keydown', e => {
    if (e.key === 'Escape') { toggleCommon(false); commonBtn.focus(); return; }
    if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
    const items = [...commonList.querySelectorAll('.rxn-c-item')], k = items.indexOf(document.activeElement);
    if (k < 0) return;
    e.preventDefault();
    if (e.key === 'ArrowUp' && k === 0) commonQ.focus(); else (items[k + (e.key === 'ArrowDown' ? 1 : -1)] || items[k]).focus();
  });
  window.ReactionsCommon = { list: () => COMMON.map(c => ({ id: c.r.id, name: c.r.name, family: c.fam, sub: c.sub })), load: loadCommon, open: () => toggleCommon(true) };

  window.Reactions = { set: () => {} };
})();
