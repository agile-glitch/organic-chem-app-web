/* Lab3D — a 3D laboratory bench: glassware, heat sources, hoses, and what is inside the glassware.

   The Simulate tab's Literature mode plays recorded experiments on this bench (js/sim_lit.js supplies them; nothing
   here knows any chemistry). Everything is built from three.js primitives, so no model files are downloaded: a
   flask is a lathe of its own profile, and a liquid is the same profile filled to a height.

   Scale is real, in centimetres. Every vessel is sized so that its capacity (up to the neck) is its nominal volume —
   a "1 L" round-bottom flask holds 1000 mL below its neck — and a liquid given in mL fills it to the height those
   mL actually reach, found by integrating the vessel's own profile. Colours are illustrative; amounts are not.

     Lab3D.create(container, opts) → view
       view.setBench(items)          the apparatus: [{id, kind, label, x, z, y?, volume?, on?, tilt?, ...}]
                                     on: 'otherId' or 'otherId:port' stands the piece on that port (a condenser on
                                     a flask's neck, a funnel in a side neck); hoses: {kind:'hose', from, to, type}
                                     with from/to 'id:port' or [x, y, z], type 'water' | 'vacuum' | 'gas'
       view.setContents(id, c, ms?)  what is in a vessel, animated over ms:
                                     {ml | level, color, label, items[], layers:[{ml, color, label}], solid:{depth, color}
                                      (or a list of them, bottom first: a packed column), cloudy, bubbling, stirring,
                                      dripping, fuming, coat:{color} (crystals on a culture dish's lid), plate:{front,
                                      spots:[{lane, rf, color}]} (a TLC plate in its tank)}
       view.setHeat(id, w)           0 (off) … 1 (full): a burner's flame, a hotplate's, mantle's or oil bath's glow;
                                     a microwave oven shuts its door and lights up
       view.setFlame(id, color)      a burner flame's colour (null: the ordinary blue) — a Beilstein test's green
       view.setReading(id, °C)       a thermometer's (or melting-point apparatus's) red column
       view.pour(fromId, toId, ms)   lift, tip and pour one vessel into another (levels are set separately)
       view.focus(id | [ids] | null) glide the camera to one piece (or frame several), or back to the whole bench
       view.highlight(id | null)     ring the piece the current step is about
       view.resize() / view.dispose()
       opts.onPick(id, part)         a click: part 'liquid' (what is inside) or 'glass' (the piece itself)
       opts.onHover(id, part)        the same for hovering (null when nothing is under the pointer)

   Kinds — glassware: rbflask, threeneck, receiver, erlenmeyer, filterflask (side-arm), beaker, cylinder (graduated),
   sepfunnel, droppingfunnel, testtube, falcon, vial, buchner, hirsch, watchglass, conicalvial, craig, funnel, culturedish
   (with its lid), column (chromatography), tlc (developing tank; plate: true stands a plate in it); fittings:
   condenser, distillhead, vigreux, thermometer, dryingtube, adapter, aircondenser, rod (a glass rod, or a copper
   wire); equipment: burner (Bunsen), hotplate, stirplate, mantle, oilbath, waterbath, steambath, icebath (beaker: a
   glass beaker of ice), alblock, microwave, meltemp (melting-point apparatus), rotovap, labjack, stand, clamp, ring,
   corkring, rack (test tubes), overheadstirrer, tap (water tap with aspirator), vacuumtrap, desiccator; and hose
   (type water | vacuum | gas | nylon; taut: true runs it straight). A port 'id:ml30' is the point on a vessel's axis
   at the height 30 mL fills it to. */
window.Lab3D = (() => {
  'use strict';
  /* A dark bench with mid-grey walls: most solids that come out of these reactions are white and the glass is clear,
     so both read best against something dark, while the walls stay light enough to frame the glass. */
  const COL = { glass: 0xeaf4f8, steel: 0xb6bdc2, bench: 0x141618, wall: 0x6f767b, floor: 0x121416, dark: 0x3a4045,
    cork: 0xc79b62, oil: 0xc9a227, clampBlue: 0x3a6fb0 };
  const HOSE = { water: 0x3b82c4, vacuum: 0x1a1a1a, gas: 0xd08a2a, nylon: 0xf4f1e8 };
  const WALL_T = 0.14;                                         // glass wall thickness, cm

  /* ---------------- profiles: [[radius, height], …] bottom to top, in cm ---------------- */
  /* at a flat floor (a horizontal first segment) the radius is the floor's full width, not the axis point it starts
     from: a liquid band that starts on the floor is a cylinder, not a cone */
  function radiusAt(profile, y) {
    if (y < profile[0][1]) return profile[0][0];
    for (let i = 1; i < profile.length; i++) {
      const [r0, y0] = profile[i - 1], [r1, y1] = profile[i];
      if (y <= y1) { const t = y1 === y0 ? 1 : (y - y0) / (y1 - y0); return r0 + (r1 - r0) * t; }
    }
    return profile[profile.length - 1][0];
  }
  /* the volume (mL = cm³) inside the wall below height y */
  function volumeBelow(profile, y) {
    const y0 = profile[0][1], dy = 0.02;
    let v = 0;
    for (let h = y0; h < y; h += dy) { const r = Math.max(0, radiusAt(profile, Math.min(y, h + dy / 2)) - WALL_T); v += Math.PI * r * r * Math.min(dy, y - h); }
    return v;
  }
  /* the height a given volume reaches (bisection on volumeBelow) */
  function heightFor(profile, ml, capTop) {
    let lo = profile[0][1], hi = capTop;
    if (!(ml > 0)) return lo;
    if (volumeBelow(profile, hi) <= ml) return hi;
    for (let k = 0; k < 40; k++) { const mid = (lo + hi) / 2; if (volumeBelow(profile, mid) < ml) lo = mid; else hi = mid; }
    return (lo + hi) / 2;
  }
  /* scale a unit profile so that the volume inside its wall, below its neck, is `ml`. The wall does not scale with
     the vessel, so the scale is found by a few rounds of correction rather than in one step. */
  function scaled(unit, ml) {
    const at = s => unit.profile.map(([r, y]) => [r * s, y * s]);
    let s = 1;
    for (let k = 0; k < 6; k++) {
      const cap = volumeBelow(at(s), unit.neck * s);
      if (!(cap > 0)) { s *= 2; continue; }
      s *= Math.cbrt(ml / cap);
    }
    const out = Object.assign({}, unit);
    out.profile = at(s);
    for (const k of ['height', 'neck', 'body', 'stem', 'plate']) if (typeof unit[k] === 'number') out[k] = unit[k] * s;
    out.scale = s;
    return out;
  }

  const SHAPES = {
    rbflask(v = 250) {                                          // sphere with a neck; the neck joint does not scale up
      const p = [], neckR = 0.27, ang = Math.acos(neckR);
      for (let a = -Math.PI / 2; a < ang; a += Math.PI / 30) p.push([Math.cos(a), 1 + Math.sin(a)]);
      const yTop = 1 + Math.sin(ang);
      p.push([neckR, yTop], [neckR, yTop + 0.55], [neckR * 1.2, yTop + 0.6]);
      return scaled({ profile: p, neck: yTop, height: yTop + 0.6, body: 1 }, v);
    },
    threeneck(v = 500) { const s = SHAPES.rbflask(v); s.threeNeck = true; return s; },
    receiver(v = 100) { return SHAPES.rbflask(v); },
    erlenmeyer(v = 250) {                                       // a cone: wide base, narrow neck (a 250 mL one is ~8.5 × 13.5 cm)
      const p = [[0, 0], [1, 0], [1, 0.15], [0.26, 2.4], [0.26, 3.0], [0.31, 3.06]];
      return scaled({ profile: p, neck: 2.4, height: 3.06, body: 1 }, v);
    },
    filterflask(v = 500) { const s = SHAPES.erlenmeyer(v); s.sidearm = true; return s; },
    beaker(v = 250) {
      return scaled({ profile: [[0, 0], [1, 0], [1, 2.6], [1.05, 2.68]], neck: 2.6, height: 2.68, body: 1 }, v);
    },
    cylinder(v = 100) {                                         // a graduated cylinder, tall and thin, on a foot
      const s = scaled({ profile: [[0, 0.4], [1, 0.4], [1, 12], [1.12, 12.3]], neck: 12, height: 12.3, body: 1 }, v);
      s.foot = true; s.profile.unshift([0, 0]); return s;
    },
    sepfunnel(v = 250) {                                        // pear-shaped, with a stopcock below
      const p = [[0, 0], [0.1, 0.05], [0.1, 1.0], [0.4, 1.45], [0.85, 2.1], [1, 2.8], [1, 3.3], [0.7, 3.9], [0.3, 4.3], [0.3, 4.75], [0.36, 4.85]];
      const s = scaled({ profile: p, neck: 4.3, height: 4.85, body: 1, stem: 1.0 }, v); s.tap = true; return s;
    },
    droppingfunnel(v = 100) {                                   // a pressure-equalising addition funnel
      const p = [[0, 0], [0.1, 0.05], [0.1, 1.1], [0.55, 1.5], [1, 1.9], [1, 4.4], [0.35, 4.9], [0.35, 5.4], [0.42, 5.5]];
      const s = scaled({ profile: p, neck: 4.9, height: 5.5, body: 1, stem: 1.1 }, v); s.tap = true; return s;
    },
    testtube() {
      const p = [], R = 0.8;
      for (let a = -Math.PI / 2; a <= 0; a += Math.PI / 16) p.push([R * Math.cos(a), R + R * Math.sin(a)]);
      p.push([R, 15], [R * 1.1, 15.2]);
      return { profile: p, neck: 15, height: 15.2, body: R };
    },
    falcon() {                                                   // a 50 mL conical centrifuge ("Falcon") tube: cone tip, straight body
      return { profile: [[0, 0], [0.22, 0.05], [1.32, 2.0], [1.38, 10.6], [1.45, 10.8]], neck: 10.6, height: 10.8, body: 1.38 };
    },
    vial() { return { profile: [[0, 0], [1.3, 0], [1.3, 5], [0.95, 5.5], [0.95, 6.2]], neck: 5.5, height: 6.2, body: 1.3 }; },
    buchner() {                                                 // porcelain filter funnel: stem, cone, perforated plate
      return { profile: [[0, 0], [0.55, 0], [0.55, 3.2], [1.1, 3.6], [4.2, 4.2], [4.3, 7.4], [4.45, 7.6]], neck: 7.4, height: 7.6,
        body: 4.3, plate: 4.2, porcelain: true };
    },
    watchglass(v, spec = {}) {                                   // r: its radius (default 4.4 cm)
      const k = (spec.r || 4.4) / 4.4;
      return { profile: [[0, 0], [2 * k, 0.14 * k], [4.2 * k, 0.55 * k], [4.4 * k, 0.6 * k]], neck: 0.55 * k, height: 0.6 * k, body: 4.4 * k };
    },
    conicalvial(v = 5) {                                         // a microscale conical vial (Reacti-vial style)
      const p = [[0, 0], [0.25, 0.3], [1, 1.25], [1, 2.6], [0.82, 2.85], [0.82, 3.3], [0.9, 3.4]];
      return scaled({ profile: p, neck: 2.6, height: 3.4, body: 1 }, v);
    },
    hirsch() {                                                   // a small porcelain filter funnel for microscale work
      return { profile: [[0, 0], [0.3, 0], [0.3, 1.8], [0.55, 2.0], [1.7, 2.3], [1.8, 3.6], [1.9, 3.7]], neck: 3.6, height: 3.7, body: 1.8, plate: 2.3, porcelain: true };
    },
    funnel(v = 50) {                                             // a glass filter / powder funnel: stem and cone
      const p = [[0, 0], [0.18, 0], [0.18, 1.6], [0.35, 1.8], [1, 3.2], [1.04, 3.3]];
      return scaled({ profile: p, neck: 3.2, height: 3.3, body: 1 }, v);
    },
    craig() {                                                    // a Craig tube for microscale recrystallisation
      return { profile: [[0, 0], [0.45, 0.08], [0.5, 0.4], [0.5, 4.4], [0.62, 4.6], [0.62, 5.2]], neck: 4.4, height: 5.2, body: 0.5 };
    },
    tlc(v, spec = {}) {                                          // a developing tank (a jar with a lid); r and h in cm
      const R = spec.r || 3.4, H = spec.h || 11;
      return { profile: [[0, 0], [R, 0], [R, H], [R * 1.04, H + 0.2]], neck: H, height: H + 0.2, body: R, lid: true, tlcPlate: !!spec.plate };
    },
    culturedish() {                                              // a 100 × 15 mm culture (Petri) dish, its larger half on top
      return { profile: [[0, 0], [4.7, 0], [4.7, 1.45], [4.75, 1.5]], neck: 1.45, height: 1.72, body: 4.75, dishLid: 5.05 };
    },
    column() {                                                   // a chromatography column: stopcock, frit, tube, reservoir
      return { profile: [[0, 0], [0.12, 0.05], [0.12, 2.2], [0.6, 2.6], [0.6, 16], [1.25, 17.6], [1.25, 21], [1.33, 21.2]],
        neck: 21, height: 21.2, body: 0.6, stem: 2.2, plate: 2.7, tap: true, open: true, scale: 1 };
    },
  };

  function create(container, opts = {}) {
    const THREE = window.THREE;
    if (!THREE || !THREE.WebGLRenderer) throw new Error('Lab3D: three.js (window.THREE) is not loaded');
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x3c4247);
    const camera = new THREE.PerspectiveCamera(38, 1, 0.5, 600);
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    container.appendChild(renderer.domElement);
    renderer.domElement.style.cssText = 'display:block;width:100%;height:100%;touch-action:none;cursor:grab';
    renderer.domElement.setAttribute('aria-label', '3D laboratory bench');

    /* filmic tone mapping and soft shadows */
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.NeutralToneMapping; renderer.toneMappingExposure = 1.1;
    renderer.localClippingEnabled = true;                       // a tipped glass's liquid is cut by a level plane
    renderer.shadowMap.enabled = true; renderer.shadowMap.type = THREE.PCFShadowMap;           // (r185 folds PCFSoft into PCF; key.shadow.radius softens it)

    /* an environment for reflections: a lab room made on the fly (grey walls, strip lights in the ceiling, a bright
       window), rendered once into a prefiltered map. Glass, liquids and steel reflect it, which is most of what makes
       glass read as glass. No image files: the page must run from file://. */
    {
      const env = new THREE.Scene(), pm = new THREE.PMREMGenerator(renderer);
      const basic = c => new THREE.MeshBasicMaterial({ color: c, side: THREE.BackSide });
      const room = new THREE.Mesh(new THREE.BoxGeometry(20, 8, 20), basic(0x9aa1a6)); room.position.y = 3; env.add(room);
      const glow = (w, d, x, y, z, k, rx = Math.PI / 2) => { const m = new THREE.Mesh(new THREE.PlaneGeometry(w, d), new THREE.MeshBasicMaterial({ color: new THREE.Color(k, k, k), side: THREE.DoubleSide }));
        m.position.set(x, y, z); m.rotation.x = rx; env.add(m); };
      for (const x of [-6, -2, 2, 6]) glow(1.6, 14, x, 6.9, 0, 7);             // ceiling strip lights
      glow(9, 4.5, 0, 3.5, -9.9, 3.2, 0);                                       // a window behind the bench
      const bench = new THREE.Mesh(new THREE.PlaneGeometry(20, 20), new THREE.MeshBasicMaterial({ color: 0x24262a }));
      bench.rotation.x = -Math.PI / 2; bench.position.y = -0.9; env.add(bench);
      scene.environment = pm.fromScene(env, 0.035).texture;
      pm.dispose();
    }
    /* an evenly lit room: most of the light is the sky and the environment, one soft overhead lamp for shape, and a
       gentle fill from the front so nothing falls into shadow. No hard, dramatic key light. */
    scene.add(new THREE.HemisphereLight(0xffffff, 0x6c757b, 1.25));
    const key = new THREE.DirectionalLight(0xfff8ef, 1.15); key.position.set(20, 44, 26); scene.add(key);
    key.castShadow = true; key.shadow.mapSize.set(2048, 2048); key.shadow.bias = -0.0004; key.shadow.normalBias = 0.03; key.shadow.radius = 6;
    key.shadow.intensity = 0.62;                                // soft shadows that place things, not deep dramatic ones
    Object.assign(key.shadow.camera, { left: -60, right: 60, top: 45, bottom: -45, near: 5, far: 150 }); key.shadow.camera.updateProjectionMatrix();
    const rim = new THREE.DirectionalLight(0xe6edff, 0.5); rim.position.set(-24, 18, -20); scene.add(rim);
    const fill = new THREE.DirectionalLight(0xffffff, 0.55); fill.position.set(-10, 16, 34); scene.add(fill);

    /* the room: a black epoxy bench top with a soft sheen, tiled grey walls, a dark floor */
    const benchTop = new THREE.Mesh(new THREE.BoxGeometry(96, 3, 52), new THREE.MeshPhysicalMaterial({ color: COL.bench, roughness: 0.38, metalness: 0,
      clearcoat: 0.25, clearcoatRoughness: 0.5, envMapIntensity: 0.3 }));
    benchTop.position.y = -1.5; benchTop.receiveShadow = true; scene.add(benchTop);
    /* one small noise image, made on the fly, used as a bump and roughness map: powders read as grains rather than
       smooth resin, steel as brushed rather than mirror-polished */
    const grain = (() => {
      const c = document.createElement('canvas'); c.width = c.height = 128; const g = c.getContext('2d');
      const img = g.createImageData(128, 128);
      for (let i = 0; i < 128 * 128; i++) { const v = 120 + Math.floor(Math.random() * 135); img.data[i * 4] = img.data[i * 4 + 1] = img.data[i * 4 + 2] = v; img.data[i * 4 + 3] = 255; }
      g.putImageData(img, 0, 0);
      const t = new THREE.CanvasTexture(c); t.wrapS = t.wrapT = THREE.RepeatWrapping; t.repeat.set(6, 6); return t;
    })();
    const brushed = (() => {                                     // fine lengthwise streaks, for steel
      const c = document.createElement('canvas'); c.width = c.height = 128; const g = c.getContext('2d');
      g.fillStyle = '#b4b4b4'; g.fillRect(0, 0, 128, 128);
      for (let i = 0; i < 260; i++) { const v = 140 + Math.floor(Math.random() * 80); g.strokeStyle = `rgb(${v},${v},${v})`;
        g.beginPath(); g.moveTo(Math.random() * 128, 0); g.lineTo(Math.random() * 128, 128); g.stroke(); }
      const t = new THREE.CanvasTexture(c); t.wrapS = t.wrapT = THREE.RepeatWrapping; t.repeat.set(2, 2); return t;
    })();
    const tiles = (() => {                                        // wall tiles, drawn on a canvas (no image files)
      const c = document.createElement('canvas'); c.width = c.height = 256; const g = c.getContext('2d');
      g.fillStyle = '#9ba1a5'; g.fillRect(0, 0, 256, 256);
      g.fillStyle = '#8a9094'; for (let i = 0; i <= 256; i += 128) { g.fillRect(i - 1.5, 0, 3, 256); g.fillRect(0, i - 1.5, 256, 3); }
      const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace; t.wrapS = t.wrapT = THREE.RepeatWrapping; t.anisotropy = 4;
      return t;
    })();
    const wallTex = (w, h) => { const t = tiles.clone(); t.needsUpdate = true; t.repeat.set(w / 30, h / 30); return t; };
    const wallMat = new THREE.MeshStandardMaterial({ color: 0xd2d7d9, map: wallTex(96, 60), roughness: 0.7, side: THREE.DoubleSide });
    const sideMat = new THREE.MeshStandardMaterial({ color: 0xd2d7d9, map: wallTex(52, 60), roughness: 0.7, side: THREE.DoubleSide });
    const back = new THREE.Mesh(new THREE.PlaneGeometry(96, 60), wallMat); back.position.set(0, 27, -26); back.receiveShadow = true; scene.add(back);
    for (const sx of [-1, 1]) {
      const side = new THREE.Mesh(new THREE.PlaneGeometry(52, 60), sideMat); side.receiveShadow = true;
      side.rotation.y = sx * Math.PI / 2; side.position.set(sx * 48, 27, 0); scene.add(side);
    }
    const floor = new THREE.Mesh(new THREE.PlaneGeometry(200, 200), new THREE.MeshStandardMaterial({ color: COL.floor, roughness: 1 }));
    floor.rotation.x = -Math.PI / 2; floor.position.y = -30; scene.add(floor);

    const group = new THREE.Group(); scene.add(group);
    const mats = {
      /* clear glass: the body barely tints what is behind it, and the reflections of the room do the work. The inner
         (BackSide) shell brightens the rims, where a real vessel's glass is thickest. */
      glass: new THREE.MeshPhysicalMaterial({ color: COL.glass, roughness: 0.02, metalness: 0, clearcoat: 1, clearcoatRoughness: 0.03,
        ior: 1.52, specularIntensity: 1, envMapIntensity: 2.4, transparent: true, opacity: 0.15, side: THREE.DoubleSide, depthWrite: false }),
      glassEdge: new THREE.MeshPhysicalMaterial({ color: 0xdff0f6, roughness: 0.06, metalness: 0, envMapIntensity: 2.2,
        transparent: true, opacity: 0.2, side: THREE.BackSide, depthWrite: false }),
      steel: new THREE.MeshStandardMaterial({ color: 0xb6bdc2, roughness: 0.34, metalness: 0.88, roughnessMap: brushed, envMapIntensity: 1.2 }),
      dark: new THREE.MeshStandardMaterial({ color: COL.dark, roughness: 0.7 }),
      porcelain: new THREE.MeshStandardMaterial({ color: 0xf6f4f0, roughness: 0.42, bumpMap: grain, bumpScale: 0.006, envMapIntensity: 0.8 }),
      blue: new THREE.MeshStandardMaterial({ color: COL.clampBlue, roughness: 0.6 }),
      cork: new THREE.MeshStandardMaterial({ color: COL.cork, roughness: 0.9 }),
      white: new THREE.MeshStandardMaterial({ color: 0xf7f7f6, roughness: 0.45, bumpMap: grain, bumpScale: 0.004 }),
      red: new THREE.MeshStandardMaterial({ color: 0xc0392b, roughness: 0.5 }),
      print: new THREE.MeshBasicMaterial({ color: 0xf6f6f2, transparent: true, opacity: 0.8, depthWrite: false, side: THREE.DoubleSide }),
      frost: new THREE.MeshStandardMaterial({ color: 0xf2f4f5, roughness: 1, transparent: true, opacity: 0.38, depthWrite: false, side: THREE.DoubleSide }),
    };
    const lathe = (profile, seg = 48) => new THREE.LatheGeometry(profile.map(([r, y]) => new THREE.Vector2(Math.max(r, 1e-3), y)), seg);
    const cyl = (r, h, m, open = false, seg = 20) => new THREE.Mesh(new THREE.CylinderGeometry(r, r, h, seg, 1, open), m);

    const items = new Map();          // id → {spec, root, shape, ports, liquid[], ...}
    const hoses = [];
    let anims = [];                   // running animations: {t0, ms, step(t), done?}

    /* ---------------- building the pieces ---------------- */
    function glassMesh(geo, id) {
      const m = new THREE.Mesh(geo, mats.glass); m.renderOrder = 3; m.userData = { id, part: 'glass' };
      const edge = new THREE.Mesh(geo, mats.glassEdge); edge.renderOrder = 2; edge.scale.setScalar(1.004); edge.userData = { id, part: 'glass' };
      const g = new THREE.Group(); g.add(edge, m); return g;
    }
    function buildVessel(spec) {
      const shape = (SHAPES[spec.kind] || SHAPES.beaker)(spec.volume, spec);
      const root = new THREE.Group();
      const ports = { neck: new THREE.Vector3(0, shape.height, 0) };
      const body = shape.porcelain
        ? (() => { const m = new THREE.Mesh(lathe(shape.profile), new THREE.MeshStandardMaterial({ color: 0xf4f2ee, roughness: 0.35, side: THREE.DoubleSide }));
          m.userData = { id: spec.id, part: 'body' }; return m; })()
        : glassMesh(lathe(shape.profile), spec.id);
      root.add(body);
      if (shape.threeNeck) {                                    // two side necks, leaning out at 24°
        const R = shape.body, tilt = 0.42;
        for (const [name, sx] of [['neckL', -1], ['neckR', 1]]) {
          const g = new THREE.Group();
          const tube = glassMesh(new THREE.CylinderGeometry(0.27 * R, 0.3 * R, 0.6 * R, 20, 1, true), spec.id);
          tube.position.y = 0.3 * R; g.add(tube);
          g.position.set(sx * 0.62 * R, 1.6 * R, 0); g.rotation.z = -sx * tilt;
          root.add(g);
          ports[name] = new THREE.Vector3(sx * (0.62 * R + 0.6 * R * Math.sin(tilt)), 1.6 * R + 0.6 * R * Math.cos(tilt), 0);
          ports[name].tilt = -sx * tilt;
        }
      }
      if (shape.sidearm) {                                      // the filter flask's hose barb, just below the neck
        const y = shape.neck + 0.2 * (shape.height - shape.neck), x0 = radiusAt(shape.profile, y);
        const arm = glassMesh(new THREE.CylinderGeometry(0.35, 0.35, 3.2, 14), spec.id);
        arm.rotation.z = Math.PI / 2; arm.position.set(x0 + 1.5, y, 0); root.add(arm);
        ports.vacuum = new THREE.Vector3(x0 + 3.1, y, 0);
      }
      if (shape.foot) {
        const foot = new THREE.Mesh(new THREE.CylinderGeometry(shape.body * 2.2, shape.body * 2.4, 0.5, 6), mats.blue);
        foot.position.y = 0.25; root.add(foot);
        for (let k = 1; k < 10; k++) {                          // graduation marks
          const mark = new THREE.Mesh(new THREE.BoxGeometry(0.03, 0.03, shape.body * 0.8), mats.white);
          mark.position.set(shape.body + 0.02, shape.profile[1][1] + (shape.neck - shape.profile[1][1]) * k / 10, 0); root.add(mark);
        }
      }
      if (shape.lid) { const lid = glassMesh(new THREE.CylinderGeometry(shape.body * 1.1, shape.body * 1.1, 0.4, 32), spec.id); lid.position.y = shape.height + 0.2; root.add(lid); }
      if (shape.dishLid) {                                      // the culture dish's larger half, upside down over it
        const R = shape.dishLid, lid = glassMesh(lathe([[R, 0.2], [R, shape.height - 0.08], [R - 0.12, shape.height], [0, shape.height]]), spec.id);
        root.add(lid);
      }
      let tlcPlate = null;
      if (shape.tlcPlate) {                                     // a TLC plate leaning in the tank: 2.5 cm wide, origin 1 cm up
        const H = Math.min(7.5, shape.neck - 1.2), W = 2.5;
        tlcPlate = new THREE.Group();
        const face = new THREE.Mesh(new THREE.BoxGeometry(W, H, 0.06), new THREE.MeshStandardMaterial({ color: 0xf7f7f2, roughness: 0.9 }));
        face.position.y = H / 2; tlcPlate.add(face);
        const origin = new THREE.Mesh(new THREE.BoxGeometry(W * 0.9, 0.025, 0.07), new THREE.MeshBasicMaterial({ color: 0x9aa0a4 }));
        origin.position.y = 1; tlcPlate.add(origin);             // the pencil line
        tlcPlate.position.set(0, 0.1, -0.2); tlcPlate.rotation.x = -0.1;
        root.add(tlcPlate);
        tlcPlate.userData.size = { W, H, origin: 1 };
      }
      if (shape.tap) {                                          // the stopcock: a barrel across the stem and a handle
        const y = shape.stem * 0.6;
        const barrel = glassMesh(new THREE.CylinderGeometry(0.32, 0.32, 2.2, 14), spec.id);
        barrel.rotation.z = Math.PI / 2; barrel.position.y = y; root.add(barrel);
        const handle = new THREE.Mesh(new THREE.BoxGeometry(0.35, 0.35, 1.8), mats.blue); handle.position.set(1.15, y, 0); root.add(handle);
        if (!shape.open && !spec.open) {                        // a stoppered funnel (a column, or a funnel stirred open, has none)
          const stopper = new THREE.Mesh(new THREE.CylinderGeometry(0.4 * shape.scale, 0.3 * shape.scale, 0.35 * shape.scale, 16), mats.blue);
          stopper.position.y = shape.height + 0.4; root.add(stopper);
        }
        ports.tip = new THREE.Vector3(0, 0, 0);
      }
      if (spec.cap) {                                           // a sealed vessel: a crimped aluminium cap, or a coloured screw cap
        const rTop = shape.profile[shape.profile.length - 1][0];
        const capMat = spec.cap === true ? mats.steel : new THREE.MeshStandardMaterial({ color: new THREE.Color(spec.cap), roughness: 0.55 });
        const cap = new THREE.Mesh(new THREE.CylinderGeometry(rTop * 1.18, rTop * 1.18, 0.7, 24), capMat);
        cap.position.y = shape.height + 0.2; cap.castShadow = true; root.add(cap);
      }
      if (!shape.porcelain && !shape.dishLid && !shape.lid) {   // the fire-polished bead round the rim
        const rTop = shape.profile[shape.profile.length - 1][0];
        const bead = glassMesh(new THREE.TorusGeometry(rTop, Math.min(0.1, Math.max(0.04, rTop * 0.07)), 8, 40), spec.id);
        bead.rotation.x = Math.PI / 2; bead.position.y = shape.height; root.add(bead);
      }
      const decal = m => { m.userData = { id: spec.id, part: 'glass' }; m.renderOrder = 4; root.add(m); return m; };
      if (/^(beaker|erlenmeyer|filterflask)$/.test(spec.kind) && spec.volume >= 20) {   // white printed graduations
        const cap = spec.volume, step = cap <= 50 ? 10 : cap <= 150 ? 25 : cap <= 300 ? 50 : 100;
        for (let v = step, k = 1; v <= cap * 0.81; v += step, k++) {
          const y = heightFor(shape.profile, v, shape.height), r = radiusAt(shape.profile, y);
          const m = decal(new THREE.Mesh(new THREE.BoxGeometry(k % 2 ? r * 0.16 : r * 0.3, 0.04, 0.02), mats.print));
          m.position.set(-r * 0.25, y, r + 0.012); m.rotation.y = -0.25;
        }
        if (spec.kind === 'beaker') {                          // the white marking patch
          const r = shape.body, m = decal(new THREE.Mesh(new THREE.PlaneGeometry(r * 0.5, shape.height * 0.12), mats.print));
          m.position.set(r * 0.3, shape.height * 0.78, r * 0.96); m.rotation.y = 0.3;
        }
      }
      if (/^(rbflask|threeneck|receiver)$/.test(spec.kind)) {  // the frosted ground-glass joint on the neck
        const rn = radiusAt(shape.profile, (shape.neck + shape.height) / 2), h = (shape.height - shape.neck) * 0.85;
        const j = decal(new THREE.Mesh(new THREE.CylinderGeometry(rn + 0.012, rn + 0.012, h, 24, 1, true), mats.frost));
        j.position.y = shape.height - h / 2 - 0.05;
      }
      return { root, shape, ports, tlcPlate };
    }
    function buildCondenser(spec) {                             // Liebig: inner tube in a water jacket, two spouts
      const L = spec.length || 20, root = new THREE.Group(), ports = {};
      const inner = glassMesh(new THREE.CylinderGeometry(0.62, 0.62, L, 20, 1, true), spec.id); inner.position.y = L / 2; root.add(inner);
      const jacket = glassMesh(new THREE.CylinderGeometry(1.55, 1.55, L - 5, 24, 1, true), spec.id); jacket.position.y = L / 2; root.add(jacket);
      const water = new THREE.Mesh(new THREE.CylinderGeometry(1.4, 1.4, L - 5.2, 24, 1, true),
        new THREE.MeshPhysicalMaterial({ color: 0x7fb8e0, transparent: true, opacity: 0.18, depthWrite: false, side: THREE.DoubleSide }));
      water.position.y = L / 2; water.renderOrder = 1; root.add(water);
      /* water spouts near each end. end0 is the root end, end1 the far end. Water always goes in at the LOWER end and
         out at the upper one: for an upright reflux condenser that is end0 → end1 ('in' and 'out'); for a distillation
         condenser, which slopes down from the still head, it is end1 → end0. */
      for (const [name, y] of [['end0', 3.6], ['end1', L - 3.6]]) {
        const spout = glassMesh(new THREE.CylinderGeometry(0.3, 0.3, 2.6, 12), spec.id);
        spout.rotation.z = Math.PI / 2; spout.position.set(2.6, y, 0); root.add(spout);
        ports[name] = new THREE.Vector3(3.9, y, 0);
      }
      ports.in = ports.end0; ports.out = ports.end1;
      ports.top = new THREE.Vector3(0, L, 0);
      return { root, ports, water };
    }
    function buildDistillHead(spec) {                           // a still head: vertical tube, side arm down at 75°
      const root = new THREE.Group(), ports = {};
      const up = glassMesh(new THREE.CylinderGeometry(0.62, 0.62, 7, 20, 1, true), spec.id); up.position.y = 3.5; root.add(up);
      /* the side arm leaves at 5.2 cm and slopes 15° below the horizontal, away from the flask */
      const tilt = -(Math.PI / 2 + 0.26), dx = Math.sin(-tilt), dy = Math.cos(-tilt);   // (0.966, −0.259)
      const arm = glassMesh(new THREE.CylinderGeometry(0.55, 0.55, 7, 16, 1, true), spec.id);
      arm.rotation.z = tilt; arm.position.set(3.5 * dx, 5.2 + 3.5 * dy, 0); root.add(arm);
      ports.side = new THREE.Vector3(7 * dx, 5.2 + 7 * dy, 0); ports.side.tilt = tilt;
      ports.junction = new THREE.Vector3(0, 5.2, 0);
      ports.top = new THREE.Vector3(0, 7, 0);
      return { root, ports };
    }
    function buildVigreux(spec) {
      const L = spec.length || 18, root = new THREE.Group();
      const tube = glassMesh(new THREE.CylinderGeometry(0.8, 0.8, L, 20, 1, true), spec.id); tube.position.y = L / 2; root.add(tube);
      for (let k = 1; k < 8; k++) for (let a = 0; a < 3; a++) {        // the indentations
        const dent = new THREE.Mesh(new THREE.ConeGeometry(0.25, 0.6, 8), mats.glass);
        const ang = a * 2 * Math.PI / 3 + k * 0.6;
        dent.position.set(Math.cos(ang) * 0.55, L * k / 8, Math.sin(ang) * 0.55); dent.rotation.z = Math.PI / 2; root.add(dent);
      }
      return { root, ports: { top: new THREE.Vector3(0, L, 0) } };
    }
    function buildThermometer(spec) {
      const L = spec.length || 26, root = new THREE.Group();
      const stem = glassMesh(new THREE.CylinderGeometry(0.28, 0.28, L, 12), spec.id); stem.position.y = L / 2; root.add(stem);
      const bulb = new THREE.Mesh(new THREE.SphereGeometry(0.34, 12, 8), mats.red); bulb.position.y = 0.3; root.add(bulb);
      const col = new THREE.Mesh(new THREE.CylinderGeometry(0.08, 0.08, 1, 8), mats.red); col.position.y = 0.3; root.add(col);
      return { root, column: col, length: L };
    }
    function buildDryingTube(spec) {
      const root = new THREE.Group();
      const tube = glassMesh(new THREE.CylinderGeometry(0.9, 0.9, 7, 16), spec.id); tube.position.y = 4; root.add(tube);
      const fill = new THREE.Mesh(new THREE.CylinderGeometry(0.78, 0.78, 5, 16), new THREE.MeshStandardMaterial({ color: 0xf2f2f2, roughness: 1 }));
      fill.position.y = 4; root.add(fill);
      const neck = glassMesh(new THREE.CylinderGeometry(0.3, 0.3, 1.2, 12), spec.id); neck.position.y = 0.3; root.add(neck);
      return { root };
    }
    function buildAdapter(spec) {
      const root = new THREE.Group();
      const t = glassMesh(new THREE.CylinderGeometry(0.62, 0.45, 5, 16, 1, true), spec.id); t.position.y = 2.5; root.add(t);
      return { root, ports: { top: new THREE.Vector3(0, 5, 0) } };
    }
    function buildBurner(spec) {                                // Bunsen: base, barrel, air collar, flame
      const root = new THREE.Group();
      const base = new THREE.Mesh(new THREE.CylinderGeometry(2.6, 3.0, 0.8, 28), mats.dark); base.position.y = 0.4; root.add(base);
      const barrel = cyl(0.62, 9, mats.steel); barrel.position.y = 5.2; root.add(barrel);
      const collar = cyl(0.8, 1.2, mats.steel); collar.position.y = 1.8; root.add(collar);
      const inlet = cyl(0.3, 2.6, mats.steel); inlet.rotation.z = Math.PI / 2; inlet.position.set(1.8, 1.6, 0); root.add(inlet);
      const outer = new THREE.Mesh(new THREE.ConeGeometry(0.6, 4.4, 20, 1, true), new THREE.MeshBasicMaterial({ color: 0x4f82d8, transparent: true, opacity: 0.55, depthWrite: false }));
      outer.position.y = 12; outer.visible = false; root.add(outer);
      const inner = new THREE.Mesh(new THREE.ConeGeometry(0.3, 1.7, 16), new THREE.MeshBasicMaterial({ color: 0x9fd8ff, transparent: true, opacity: 0.9, depthWrite: false }));
      inner.position.y = 10.6; inner.visible = false; root.add(inner);
      return { root, flame: [outer, inner], ports: { gas: new THREE.Vector3(3.1, 1.6, 0), top: new THREE.Vector3(0, 9.7, 0) } };
    }
    function buildHotplate(spec, stir) {
      const root = new THREE.Group();
      const body = new THREE.Mesh(new THREE.BoxGeometry(14, 4, 16), mats.dark); body.position.y = 2; root.add(body);
      const plate = new THREE.Mesh(new THREE.CylinderGeometry(6, 6, 0.4, 36), new THREE.MeshStandardMaterial({ color: 0x5b5f62, roughness: 0.5 }));
      plate.position.y = 4.2; root.add(plate);
      for (const [k, x] of [[0, -3], [1, 3]]) {                 // the two knobs: heat and, on a stirrer, speed
        if (k === 1 && !stir) continue;
        const knob = cyl(0.9, 0.8, mats.steel); knob.rotation.x = Math.PI / 2; knob.position.set(x, 2, 8.3); root.add(knob);
      }
      const lamp = new THREE.Mesh(new THREE.SphereGeometry(0.35, 10, 8), new THREE.MeshStandardMaterial({ color: 0x3a1210, emissive: 0x000000 }));
      lamp.position.set(stir ? 0 : 3, 3, 8.05); root.add(lamp);
      return { root, plate, lamp, ports: { top: new THREE.Vector3(0, 4.4, 0) } };
    }
    function buildMantle(spec) {
      const root = new THREE.Group(), R = spec.fits || 5;
      const bowl = new THREE.Mesh(lathe([[R + 0.4, 0], [R + 1.4, 0], [R + 1.8, R * 0.9], [R + 0.5, R * 0.95], [R * 0.2, 0.9]]),
        new THREE.MeshStandardMaterial({ color: 0x4a4f53, roughness: 0.9, side: THREE.DoubleSide }));
      root.add(bowl);
      const box = new THREE.Mesh(new THREE.BoxGeometry(3.2, 2.2, 3.2), mats.dark); box.position.set(R + 4, 1.1, 2); root.add(box);
      const lamp = new THREE.Mesh(new THREE.SphereGeometry(0.3, 10, 8), new THREE.MeshStandardMaterial({ color: 0x3a1210, emissive: 0x000000 }));
      lamp.position.set(R + 4, 1.6, 3.65); root.add(lamp);
      return { root, plate: bowl, lamp, ports: { top: new THREE.Vector3(0, 0.9, 0) } };
    }
    /* a crystallising dish on a stirrer hotplate, filled with oil (an oil bath), water (a water bath) or sand */
    function buildOilBath(spec, fill = 'oil') {
      const hp = buildHotplate(spec, true), R = spec.fits || 6.5, D = spec.depth || 5.5;
      const dish = glassMesh(lathe([[0, 0], [R, 0], [R, D], [R * 1.03, D + 0.2]]), spec.id); dish.position.y = 4.4; hp.root.add(dish);
      const look = { oil: [COL.oil, 0.55, 0.2], water: [0x8fc0e0, 0.38, 0.2], sand: [0xcdb68a, 0.97, 0.95] }[fill];
      const oil = new THREE.Mesh(lathe([[0, 0.1], [R - WALL_T, 0.1], [R - WALL_T, D * 0.76], [0, D * 0.76]]),
        new THREE.MeshPhysicalMaterial({ color: look[0], transparent: true, opacity: look[1], roughness: look[2], depthWrite: fill === 'sand' }));
      oil.position.y = 4.4; oil.renderOrder = 1; hp.root.add(oil);
      hp.ports.top = new THREE.Vector3(0, 4.4 + D * 0.33, 0);   // a flask stands with its bulb in the oil (water, sand)
      if (fill === 'oil') hp.oil = oil;
      return hp;
    }
    function buildIceBath(spec) {
      const root = new THREE.Group();
      /* a metal bowl of ice (or, noIce, of water); glass: true, a glass dish instead, so a small flask in it stays in
         sight; or, beaker: 250, a glass beaker of that size packed with ice */
      const bowlProfile = spec.beaker ? SHAPES.beaker(spec.beaker).profile : (() => { const R = spec.fits || 8; return [[0, 0], [R, 0], [R, 6], [R * 1.04, 6.2]]; })();
      const R = bowlProfile[1][0], top = bowlProfile[bowlProfile.length - 2][1];
      const bowl = spec.beaker || spec.glass ? glassMesh(lathe(bowlProfile), spec.id)
        : new THREE.Mesh(lathe(bowlProfile), new THREE.MeshStandardMaterial({ color: 0xbfc6cb, roughness: 0.6, side: THREE.DoubleSide }));
      root.add(bowl);
      const fill = spec.beaker ? top * 0.7 : 4.6;
      const ice = new THREE.Mesh(lathe([[0, 0.2], [R - 0.3, 0.2], [R - 0.3, fill], [0, fill]]),
        new THREE.MeshPhysicalMaterial({ color: spec.noIce ? 0x8fbfe0 : 0xd7eef8, roughness: 0.3, transparent: true, opacity: spec.noIce ? 0.45 : 0.7, depthWrite: false }));
      ice.renderOrder = 1; root.add(ice);
      const cube = Math.min(1.4, R / 3);
      for (let k = 0; k < (spec.noIce || spec.beaker ? 0 : 14); k++) {         // (a beaker of ice has no room for loose cubes)
        const c = new THREE.Mesh(new THREE.BoxGeometry(cube, cube, cube), new THREE.MeshPhysicalMaterial({ color: 0xeaf6fb, roughness: 0.15, transparent: true, opacity: 0.85 }));
        const a = k * 2.4, r = (k % 3 + 1) * (R - cube) / 3.4;
        c.position.set(Math.cos(a) * r, fill - 0.1, Math.sin(a) * r); c.rotation.set(a, a * 1.3, 0); root.add(c);
      }
      return { root, ports: { top: new THREE.Vector3(0, spec.beaker ? 0.3 : 1.2, 0) } };
    }
    function buildAlBlock(spec) {                               // an aluminium heating block for microscale vials
      const root = new THREE.Group();
      const block = new THREE.Mesh(new THREE.BoxGeometry(7, 3.2, 7), new THREE.MeshStandardMaterial({ color: 0xb9c0c4, roughness: 0.35, metalness: 0.8 }));
      block.position.y = 1.6; root.add(block);
      const hole = new THREE.Mesh(new THREE.CylinderGeometry(1.15, 1.15, 0.05, 20), mats.dark); hole.position.y = 3.22; root.add(hole);
      /* a vial stands in a hole 1.3 cm deep: its conical tip is in the metal, and what is in it stays in sight above */
      return { root, plate: block, noGlow: true, ports: { top: new THREE.Vector3(0, 1.9, 0) } };
    }
    function buildAirCondenser(spec) {                          // a plain glass tube: air cools the vapour
      const L = spec.length || 12, root = new THREE.Group();
      const t = glassMesh(new THREE.CylinderGeometry(0.45, 0.45, L, 16, 1, true), spec.id); t.position.y = L / 2; root.add(t);
      return { root, ports: { top: new THREE.Vector3(0, L, 0) } };
    }
    function buildSteamBath(spec) {                             // a copper pot with concentric rings; wisps of steam when on
      const root = new THREE.Group();
      const pot = new THREE.Mesh(lathe([[0, 0], [7.5, 0], [7.5, 8], [7.9, 8.3]]), new THREE.MeshStandardMaterial({ color: 0xb87333, roughness: 0.45, metalness: 0.6, side: THREE.DoubleSide }));
      root.add(pot);
      for (const r of [6.2, 4.6, 3.0]) { const ring = new THREE.Mesh(new THREE.TorusGeometry(r, 0.35, 8, 36), mats.steel); ring.rotation.x = Math.PI / 2; ring.position.y = 8.3; root.add(ring); }
      return { root, plate: pot, steam: true, ports: { top: new THREE.Vector3(0, 6.8, 0) } };
    }
    function buildRotovap(spec) {                               // a rotary evaporator: bath, tilted flask, condenser, motor
      const root = new THREE.Group();
      const bath = new THREE.Mesh(lathe([[0, 0], [8, 0], [8, 6], [8.3, 6.3]]), new THREE.MeshStandardMaterial({ color: 0x9aa4aa, roughness: 0.4, metalness: 0.5, side: THREE.DoubleSide }));
      root.add(bath);
      const water = new THREE.Mesh(lathe([[0, 0.2], [7.8, 0.2], [7.8, 4.8], [0, 4.8]]), new THREE.MeshPhysicalMaterial({ color: 0x8fc0e0, transparent: true, opacity: 0.4, depthWrite: false }));
      water.renderOrder = 1; root.add(water);
      const post = cyl(0.8, 30, mats.steel); post.position.set(-12, 15, 0); root.add(post);
      const motor = new THREE.Mesh(new THREE.BoxGeometry(6, 5, 6), mats.dark); motor.position.set(-7, 18, 0); motor.rotation.z = -0.6; root.add(motor);
      const coil = glassMesh(new THREE.CylinderGeometry(2.2, 2.2, 14, 20, 1, true), spec.id); coil.position.set(-11, 29, 0); root.add(coil);
      return { root, ports: { flask: new THREE.Vector3(-1.5, 5, 0), top: new THREE.Vector3(-1.5, 5, 0) } };
    }
    function buildLabJack(spec) {                               // a scissor jack: raises a bath up around a clamped flask
      const root = new THREE.Group(), H = spec.height || 8, W = spec.width || 14;
      const top = new THREE.Mesh(new THREE.BoxGeometry(W, 0.6, W), mats.steel); top.position.y = H - 0.3; root.add(top);
      const base = new THREE.Mesh(new THREE.BoxGeometry(W, 0.6, W), mats.steel); base.position.y = 0.3; root.add(base);
      for (const sz of [-1, 1]) for (const sx of [-1, 1]) {
        const bar = cyl(0.2, Math.hypot(W * 0.8, H - 1.2), mats.steel, false, 8);
        bar.position.set(0, H / 2, sz * W * 0.35); bar.rotation.z = sx * Math.atan2(W * 0.8, H - 1.2); root.add(bar);
      }
      return { root, ports: { top: new THREE.Vector3(0, H, 0) } };
    }
    function buildStand(spec) {
      const root = new THREE.Group(), H = spec.height || 60;
      const base = new THREE.Mesh(new THREE.BoxGeometry(12, 1, 20), mats.dark); base.position.set(0, 0.5, 5); root.add(base);
      const rod = cyl(0.5, H, mats.steel, false, 16); rod.position.set(0, H / 2 + 1, 0); root.add(rod);
      return { root };
    }
    /* a clamp: a boss on the stand's rod and an arm reaching to the piece it holds */
    function buildClamp(spec) {
      const root = new THREE.Group();
      const boss = new THREE.Mesh(new THREE.BoxGeometry(1.6, 1.6, 1.6), mats.steel); root.add(boss);
      return { root, boss };
    }
    function buildRing(spec) {
      const root = new THREE.Group();
      const ring = new THREE.Mesh(new THREE.TorusGeometry(spec.r || 4.5, 0.25, 8, 36), mats.steel); ring.rotation.x = Math.PI / 2; root.add(ring);
      return { root, ports: { top: new THREE.Vector3(0, 0, 0) } };
    }
    function buildTap(spec) {                                   // a lab water tap, with an aspirator on its outlet
      const root = new THREE.Group();
      const post = cyl(0.9, 14, mats.steel); post.position.y = 7; root.add(post);
      const arm = cyl(0.7, 7, mats.steel); arm.rotation.z = Math.PI / 2; arm.position.set(3.3, 13.5, 0); root.add(arm);
      const spout = cyl(0.6, 3, mats.steel); spout.position.set(6.6, 12, 0); root.add(spout);
      const handle = new THREE.Mesh(new THREE.BoxGeometry(3, 0.6, 0.6), mats.red); handle.position.set(0, 15, 0); root.add(handle);
      const asp = cyl(0.45, 4, mats.steel); asp.position.set(6.6, 8.6, 0); root.add(asp);
      const side = cyl(0.25, 2, mats.steel); side.rotation.z = Math.PI / 2; side.position.set(7.6, 9.2, 0); root.add(side);
      const sink = new THREE.Mesh(new THREE.BoxGeometry(14, 0.6, 12), mats.steel); sink.position.set(6, 0.3, 0); root.add(sink);
      return { root, ports: { water: new THREE.Vector3(6.6, 10.6, 0), vacuum: new THREE.Vector3(8.6, 9.2, 0), drain: new THREE.Vector3(6.6, 1, 2) } };
    }
    function buildVacuumTrap(spec) {                            // a heavy-walled trap flask between aspirator and filter flask
      const v = buildVessel(Object.assign({}, spec, { kind: 'filterflask', volume: spec.volume || 250 }));
      const stopper = new THREE.Mesh(new THREE.CylinderGeometry(1.0, 0.8, 1.4, 16), mats.cork);
      stopper.position.y = v.shape.height + 0.3; v.root.add(stopper);
      v.ports.in = v.ports.neck.clone().add(new THREE.Vector3(0, 1.2, 0));
      return v;
    }
    function buildDesiccator(spec) {
      const root = new THREE.Group(), R = 9;
      const body = glassMesh(lathe([[0, 0], [R * 0.8, 0], [R, 5], [R, 12], [R * 1.08, 12.4]]), spec.id); root.add(body);
      const lid = glassMesh(lathe([[R * 1.08, 0], [R, 1], [R * 0.7, 4.6], [0.9, 6], [0, 6]]), spec.id); lid.position.y = 12.4; root.add(lid);
      const knob = glassMesh(new THREE.SphereGeometry(1.2, 16, 10), spec.id); knob.position.y = 19; root.add(knob);
      const agent = new THREE.Mesh(lathe([[0, 0.3], [R * 0.75, 0.3], [R * 0.8, 3.2], [0, 3.2]]), new THREE.MeshStandardMaterial({ color: 0xf1f1ee, roughness: 1 }));
      root.add(agent);
      const plate = new THREE.Mesh(new THREE.CylinderGeometry(R * 0.95, R * 0.95, 0.3, 32), mats.porcelain); plate.position.y = 5; root.add(plate);
      return { root, ports: { top: new THREE.Vector3(0, 5.2, 0) } };
    }

    function buildCorkRing(spec) {                              // a round-bottom flask stands on one of these
      const root = new THREE.Group(), R = spec.r || 3.4;
      const ring = new THREE.Mesh(new THREE.TorusGeometry(R, 0.9, 10, 32), mats.cork); ring.rotation.x = Math.PI / 2; ring.position.y = 0.9; root.add(ring);
      return { root, ports: { top: new THREE.Vector3(0, 1.3, 0) } };          // a flask's bulb rests in the ring
    }
    /* an overhead (mechanical) stirrer: a motor on the rod, a shaft down the flask's centre neck, a paddle at the
       bottom. setHeat(id, w) sets its speed. */
    function buildOverheadStirrer(spec) {
      const root = new THREE.Group(), L = spec.shaft || 14;
      const motor = new THREE.Mesh(new THREE.BoxGeometry(4.5, 7, 4.5), mats.dark); motor.position.y = 9; root.add(motor);
      const chuck = cyl(0.6, 2, mats.steel); chuck.position.y = 4.6; root.add(chuck);
      const spinner = new THREE.Group(); root.add(spinner);
      const shaft = cyl(0.22, L + 3.6, mats.steel, false, 10); shaft.position.y = 3.6 - (L + 3.6) / 2; spinner.add(shaft);
      const paddle = new THREE.Mesh(new THREE.BoxGeometry(spec.paddle || 5, 1.2, 0.18), new THREE.MeshStandardMaterial({ color: 0xe8e8e8, roughness: 0.5 }));
      paddle.position.y = -L + 0.6; spinner.add(paddle);
      return { root, spinner };
    }

    /* a microwave oven (46 × 27 × 36 cm): a cavity with a glass turntable, a control panel, and a shut door with a
       dark see-through window. setHeat(id, w) lights the cavity while it runs. */
    function buildMicrowave(spec) {
      const root = new THREE.Group(), W = 46, D = 36;
      const shell = new THREE.MeshStandardMaterial({ color: 0xd4d8db, roughness: 0.5, metalness: 0.25 });
      const inside = new THREE.MeshStandardMaterial({ color: 0xb9bdc0, roughness: 0.7 });
      const box = (w, h, d, x, y, z, m, parent = root) => { const b = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), m); b.position.set(x, y, z); parent.add(b); return b; };
      box(W, 3, D, 0, 1.5, 0, shell);                                  // floor
      box(W, 4, D, 0, 25, 0, shell);                                   // roof
      box(W, 20, 2, 0, 13, -D / 2 + 1, inside);                        // back
      box(2, 20, D - 2, -W / 2 + 1, 13, 1, shell);                     // left side
      box(14, 20, D - 2, W / 2 - 7, 13, 1, shell);                     // the control-panel block
      box(8, 12, 0.3, W / 2 - 7, 15, D / 2 + 0.05, mats.dark);
      box(5, 1.6, 0.1, W / 2 - 7, 19, D / 2 + 0.3, new THREE.MeshBasicMaterial({ color: 0x1d2b22 }));
      for (let k = 0; k < 9; k++) box(1.4, 0.9, 0.2, W / 2 - 9 + (k % 3) * 2, 16 - Math.floor(k / 3) * 1.5, D / 2 + 0.25, mats.steel);
      const turntable = new THREE.Mesh(new THREE.CylinderGeometry(12, 12, 0.3, 40), mats.glass);
      turntable.position.set(-6, 3.35, 1); turntable.userData = { id: spec.id, part: 'glass' }; root.add(turntable);
      const lamp = box(4, 0.2, 3, -6, 22.9, 1, new THREE.MeshStandardMaterial({ color: 0xe8e2c8, emissive: 0x000000 }));
      /* the door: a frame round the window, across the cavity's open front */
      const frame = new THREE.MeshStandardMaterial({ color: 0x3a3f43, roughness: 0.6 });
      const x0 = -W / 2, dw = 32, dh = 24, z = D / 2 + 0.7;
      for (const [w, h, x, y] of [[dw, 3, dw / 2, 2.5], [dw, 3, dw / 2, dh - 0.5], [4, dh, 2, dh / 2 + 1], [6, dh, dw - 3, dh / 2 + 1]]) box(w, h, 1.2, x0 + x, y, z, frame);
      const win = new THREE.Mesh(new THREE.BoxGeometry(dw - 10, dh - 6, 0.3),
        new THREE.MeshPhysicalMaterial({ color: 0x20262a, transparent: true, opacity: 0.3, roughness: 0.1, depthWrite: false }));
      win.position.set(x0 + 4 + (dw - 10) / 2, dh / 2 + 1, z); win.renderOrder = 3; win.userData = { id: spec.id, part: 'glass' }; root.add(win);
      return { root, cavityLamp: lamp, window: win, ports: { top: new THREE.Vector3(-6, 3.5, 1) } };
    }
    /* a test-tube rack: a base and a top plate with a row of holes; ports slot0, slot1, … stand a tube in each hole */
    function buildRack(spec) {
      const root = new THREE.Group(), n = spec.n || 4, gap = 2.6, W = n * gap + 1.6, ports = {};
      const mat = new THREE.MeshStandardMaterial({ color: 0xe9e4d6, roughness: 0.6 });
      const part = (w, h, d, x, y) => { const b = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), mat); b.position.set(x, y, 0); root.add(b); };
      part(W, 0.6, 5, 0, 0.3); part(W, 0.5, 5, 0, 6.2); part(0.5, 6, 5, -W / 2 + 0.25, 3.2); part(0.5, 6, 5, W / 2 - 0.25, 3.2);
      for (let k = 0; k < n; k++) {
        const x = -W / 2 + 0.8 + gap / 2 + k * gap;
        const hole = new THREE.Mesh(new THREE.CylinderGeometry(1.0, 1.0, 0.06, 20), mats.dark); hole.position.set(x, 6.47, 0); root.add(hole);
        ports['slot' + k] = new THREE.Vector3(x, 0.6, 0);
      }
      return { root, ports };
    }
    /* a glass stirring rod, or (metal: 'copper') a length of copper wire, from its tip (the root) along +y */
    function buildRod(spec) {
      const L = spec.length || 20, root = new THREE.Group();
      if (spec.metal) {
        const m = new THREE.Mesh(new THREE.CylinderGeometry(0.07, 0.07, L, 8), new THREE.MeshStandardMaterial({ color: 0xb87333, roughness: 0.35, metalness: 0.85 }));
        m.position.y = L / 2; root.add(m);
      } else { const g = glassMesh(new THREE.CylinderGeometry(0.25, 0.25, L, 12), spec.id); g.position.y = L / 2; root.add(g); }
      return { root, ports: { tip: new THREE.Vector3(0, 0, 0), mid: new THREE.Vector3(0, L / 2, 0), top: new THREE.Vector3(0, L, 0) } };
    }
    /* a melting-point apparatus (Mel-Temp style): a heated block in a case, a lens to watch the capillary through, a
       thermometer in the top (its scale runs 0–400 °C), and a lamp that lights while it heats */
    function buildMeltTemp(spec) {
      const root = new THREE.Group();
      const body = new THREE.Mesh(new THREE.BoxGeometry(12, 13, 15), new THREE.MeshStandardMaterial({ color: 0x5d6d7e, roughness: 0.55 }));
      body.position.y = 6.5; root.add(body);
      const lens = cyl(1.3, 3, mats.dark); lens.rotation.x = Math.PI / 2 - 0.5; lens.position.set(-2, 10.5, 8.2); root.add(lens);
      const knob = cyl(1, 0.8, mats.steel); knob.rotation.x = Math.PI / 2; knob.position.set(3, 4, 7.9); root.add(knob);
      const lamp = new THREE.Mesh(new THREE.SphereGeometry(0.35, 10, 8), new THREE.MeshStandardMaterial({ color: 0x3a1210, emissive: 0x000000 }));
      lamp.position.set(3, 7, 7.6); root.add(lamp);
      const L = 18, tg = new THREE.Group(); tg.position.set(2, 12.5, -2); root.add(tg);
      const stem = glassMesh(new THREE.CylinderGeometry(0.28, 0.28, L, 12), spec.id); stem.position.y = L / 2; tg.add(stem);
      const col = new THREE.Mesh(new THREE.CylinderGeometry(0.08, 0.08, 1, 8), mats.red); col.position.y = 0.3; tg.add(col);
      return { root, lamp, column: col, length: L, scaleMax: 400, ports: { top: new THREE.Vector3(0, 13, 0) } };
    }

    /* a Claisen head adapter: a straight tube with a side arm branching up at 35°, both capped with red septa */
    function buildClaisen(spec) {
      const root = new THREE.Group(), H = spec.length || 6, R = 0.42;
      const main = glassMesh(new THREE.CylinderGeometry(R, R, H, 16, 1, true), spec.id); main.position.y = H / 2; root.add(main);
      const tilt = -0.61, L = 4, y0 = 2.2;                           // the side arm leans away to +x
      const arm = glassMesh(new THREE.CylinderGeometry(R, R, L, 16, 1, true), spec.id);
      arm.rotation.z = tilt; arm.position.set(Math.sin(-tilt) * L / 2, y0 + Math.cos(tilt) * L / 2, 0); root.add(arm);
      const side = new THREE.Vector3(Math.sin(-tilt) * L, y0 + Math.cos(tilt) * L, 0); side.tilt = tilt;
      const septum = (at, rz) => { const c = cyl(R + 0.12, 0.5, mats.red); c.position.copy(at); c.rotation.z = rz; root.add(c); };
      septum(new THREE.Vector3(0, H + 0.25, 0), 0);
      septum(side.clone().add(new THREE.Vector3(Math.sin(-tilt) * 0.25, Math.cos(tilt) * 0.25, 0)), tilt);
      return { root, ports: { top: new THREE.Vector3(0, H + 0.5, 0), side: side } };
    }
    /* a 1 mL syringe, needle down: the root is the needle's tip */
    function buildSyringe(spec) {
      const root = new THREE.Group(), N = 3.5, B = 6.5;
      const needle = cyl(0.03, N, mats.steel, false, 6); needle.position.y = N / 2; root.add(needle);
      const hub = cyl(0.14, 0.5, mats.steel, false, 10); hub.position.y = N + 0.25; root.add(hub);
      const barrel = glassMesh(new THREE.CylinderGeometry(0.3, 0.3, B, 14, 1, true), spec.id); barrel.position.y = N + 0.5 + B / 2; root.add(barrel);
      const flange = new THREE.Mesh(new THREE.BoxGeometry(1.6, 0.12, 0.5), mats.white); flange.position.y = N + 0.5 + B; root.add(flange);
      const plunger = cyl(0.12, B * 0.8, mats.white, false, 8); plunger.position.y = N + 0.5 + B * 0.6 + B * 0.4; root.add(plunger);
      const thumb = cyl(0.45, 0.12, mats.white, false, 16); thumb.position.y = N + 0.5 + B * 1.4; root.add(thumb);
      return { root, ports: { tip: new THREE.Vector3(0, 0, 0), top: new THREE.Vector3(0, N + 0.5 + B * 1.4, 0) } };
    }

    /* a 254 nm viewing cabinet: a dark hood over a tray, with the lamp behind a baffle. setHeat(id, w) switches it on,
       and the inside glows the violet-white that makes TLC spots show up. */
    function buildUVLamp(spec) {
      const root = new THREE.Group(), W = 20, H = 16, D = 14;
      const shell = new THREE.MeshStandardMaterial({ color: 0x4a5055, roughness: 0.6 });
      const box = (w, h, d, x, y, z, m) => { const b = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), m); b.position.set(x, y, z); root.add(b); return b; };
      box(W, 0.8, D, 0, 0.4, 0, shell);                          // the tray a plate lies on
      box(W, H, 0.8, 0, H / 2, -D / 2, shell);                   // back
      box(0.8, H, D, -W / 2, H / 2, 0, shell); box(0.8, H, D, W / 2, H / 2, 0, shell);
      box(W, 0.8, D, 0, H, 0, shell);                            // hood
      const tube = box(W * 0.8, 1, 1, 0, H - 1.6, -D / 4, new THREE.MeshStandardMaterial({ color: 0xe8e6f2, emissive: 0x000000 }));
      const knob = cyl(0.8, 0.6, mats.steel); knob.rotation.x = Math.PI / 2; knob.position.set(W / 2 - 2.5, 2.4, D / 2 + 0.2); root.add(knob);
      return { root, lamp: tube, uv: true, ports: { top: new THREE.Vector3(0, 0.8, 0) } };
    }
    const BUILD = { alblock: buildAlBlock, uvlamp: buildUVLamp, aircondenser: buildAirCondenser, steambath: buildSteamBath, rotovap: buildRotovap,
      microwave: buildMicrowave, rack: buildRack, rod: buildRod, meltemp: buildMeltTemp, waterbath: s => buildOilBath(s, 'water'),
      sandbath: s => buildOilBath(s, 'sand'), claisen: buildClaisen, syringe: buildSyringe,
      corkring: buildCorkRing, overheadstirrer: buildOverheadStirrer, labjack: buildLabJack, condenser: buildCondenser, distillhead: buildDistillHead, vigreux: buildVigreux, thermometer: buildThermometer,
      dryingtube: buildDryingTube, adapter: buildAdapter, burner: buildBurner, hotplate: s => buildHotplate(s, false),
      stirplate: s => buildHotplate(s, true), mantle: buildMantle, oilbath: buildOilBath, icebath: buildIceBath, stand: buildStand,
      clamp: buildClamp, ring: buildRing, tap: buildTap, vacuumtrap: buildVacuumTrap, desiccator: buildDesiccator };

    /* a port's position in bench coordinates, 'id:port' or [x, y, z] */
    function portWorld(ref) {
      if (Array.isArray(ref)) return new THREE.Vector3(ref[0], ref[1], ref[2]);
      const [id, name] = String(ref).split(':');
      const it = items.get(id); if (!it) return null;
      const ml = /^ml([\d.]+)$/.exec(name || '');
      const p = ml && it.shape ? new THREE.Vector3(0, heightFor(it.shape.profile, +ml[1], it.shape.height), 0)
        : (it.ports && it.ports[name || 'neck']) || (it.ports && it.ports.top) || new THREE.Vector3();
      it.root.updateMatrixWorld(true);
      const w = it.root.localToWorld(p.clone()); w.tilt = p.tilt || 0; return w;
    }
    function disposeTree(o) {
      o.traverse(n => { if (n.geometry) n.geometry.dispose(); if (n.material && !Object.values(mats).includes(n.material)) { (Array.isArray(n.material) ? n.material : [n.material]).forEach(m => m.dispose()); } });
    }
    function buildHose(spec) {
      const a = portWorld(spec.from), b = portWorld(spec.to);
      if (!a || !b) return null;
      const mid = a.clone().lerp(b, 0.5), d = a.distanceTo(b);
      mid.y = spec.taut ? mid.y : Math.max(0.6, Math.min(a.y, b.y) - 0.35 * d);   // a hose hangs down between its ends
      const pts = spec.taut ? [a, mid, b]
        : [a, a.clone().lerp(mid, 0.5).setY((a.y + mid.y) / 2 - 0.05 * d), mid, b.clone().lerp(mid, 0.5).setY((b.y + mid.y) / 2 - 0.05 * d), b];
      if (spec.via) spec.via.forEach((v, k) => pts.splice(1 + k, 0, new THREE.Vector3(v[0], v[1], v[2])));
      const curve = new THREE.CatmullRomCurve3(pts);
      const mesh = new THREE.Mesh(new THREE.TubeGeometry(curve, 64, spec.radius || 0.42, 10, false),
        new THREE.MeshStandardMaterial({ color: HOSE[spec.type] || HOSE.water, roughness: 0.55, transparent: spec.type === 'water', opacity: spec.type === 'water' ? 0.85 : 1 }));
      mesh.userData = { id: spec.id, part: 'hose' }; mesh.castShadow = true;
      return mesh;
    }

    function setBench(list) {
      for (const it of items.values()) { group.remove(it.root); disposeTree(it.root); }
      for (const h of hoses) { group.remove(h); h.geometry.dispose(); h.material.dispose(); }
      for (const d of drips) { d.parent.remove(d.mesh); if (d.ring) d.parent.remove(d.ring); } drips = [];
      items.clear(); hoses.length = 0; anims = [];
      if (ring) { group.remove(ring); ring = null; }
      const hoseSpecs = [];
      for (let spec of list || []) {                         // let: a self-sizing lab jack gets its height filled in
        if (spec.kind === 'hose') { hoseSpecs.push(spec); continue; }
        /* 'around': under a piece placed earlier. A bath goes round its bottom; a lab jack is made just tall enough to
           fill the gap from the bench up to it (a receiver below a sloping condenser, or the bath around it). */
        const aroundRef = spec.around && items.get(spec.around);
        let aroundBox = null;
        if (aroundRef) { aroundRef.root.updateMatrixWorld(true); aroundBox = boundsOf([aroundRef.root]); }
        const sink = spec.sink != null ? spec.sink : spec.kind === 'labjack' ? 0 : 1.2;
        if (aroundBox && spec.kind === 'labjack') spec = Object.assign({}, spec, { height: Math.max(0.6, aroundBox.min.y - sink), width: spec.width || Math.max(8, Math.min(16, aroundBox.max.x - aroundBox.min.x + 1)) });
        const make = BUILD[spec.kind];
        const built = make ? make(spec) : buildVessel(spec);
        built.root.userData = { id: spec.id };
        /* every mesh of the piece is pickable: clear glass keeps its 'glass' part (it can be seen, and clicked,
           through), everything else is an opaque 'body' that stops a click */
        built.root.traverse(o => { if (o.isMesh && !(o.userData && o.userData.id)) o.userData = { id: spec.id, part: 'body' }; });
        built.root.traverse(o => { if (o.isMesh && !(o.material && o.material.transparent)) { o.castShadow = true; o.receiveShadow = true; } });
        const under = spec.under ? portWorld(spec.under) : null;
        const base = spec.on ? portWorld(spec.on) : null;
        if (aroundBox) {
          built.root.position.set(aroundRef.root.position.x, spec.kind === 'labjack' ? 0 : aroundBox.min.y - sink, aroundRef.root.position.z);
        } else if (under) {
          const neck = (built.ports && built.ports.neck) || new THREE.Vector3(0, built.shape ? built.shape.height : 0, 0);
          built.root.position.copy(under).sub(neck).add(new THREE.Vector3(spec.dx || 0, spec.dy || 0, spec.dz || 0));
        } else if (base) {
          built.root.position.copy(base).add(new THREE.Vector3(spec.dx || 0, spec.dy || 0, spec.dz || 0));
          built.root.rotation.z = spec.tilt != null ? spec.tilt : base.tilt;
        } else {
          built.root.position.set(spec.x || 0, spec.y || 0, spec.z || 0);
          if (spec.tilt) built.root.rotation.z = spec.tilt;
        }
        if (spec.rotY) built.root.rotation.y = spec.rotY;
        group.add(built.root);
        const it = Object.assign({ spec, contents: null, liquid: [], home: built.root.position.clone(), homeRot: built.root.rotation.clone() }, built);
        items.set(spec.id, it);
      }
      /* support rings reach from their stand's rod */
      for (const it of items.values()) {
        if (it.spec.kind !== 'ring' || !it.spec.stand) continue;
        const stand = items.get(it.spec.stand); if (!stand) continue;
        const dx = stand.root.position.x - it.root.position.x, dz = stand.root.position.z - it.root.position.z, len = Math.hypot(dx, dz);
        if (len < 0.5) continue;
        const R = it.spec.r || 4.5, arm = cyl(0.25, len - R, mats.steel, false, 8);
        arm.rotation.z = Math.PI / 2; arm.position.x = R + (len - R) / 2;
        const holder = new THREE.Group(); holder.add(arm); holder.rotation.y = -Math.atan2(dz, dx); it.root.add(holder);
      }
      /* clamps reach from their stand to what they hold */
      for (const it of items.values()) {
        if (it.spec.kind !== 'clamp') continue;
        const stand = items.get(it.spec.stand), held = items.get(it.spec.holds);
        if (!stand || !held) continue;
        const at = portWorld(it.spec.holds + ':' + (it.spec.at || 'neck')) || held.root.position.clone();
        at.y -= it.spec.below || 1.2;
        const rodX = stand.root.position.x, rodZ = stand.root.position.z;
        it.root.position.set(rodX, at.y, rodZ);
        const dx = at.x - rodX, dz = at.z - rodZ, len = Math.hypot(dx, dz);
        const arm = cyl(0.3, len, mats.steel); arm.rotation.z = Math.PI / 2; arm.position.x = len / 2;
        const jaw = new THREE.Mesh(new THREE.TorusGeometry(1.1, 0.28, 8, 20, Math.PI * 1.4), mats.blue); jaw.position.x = len; jaw.rotation.x = Math.PI / 2;
        const holder = new THREE.Group(); holder.add(arm, jaw); holder.rotation.y = -Math.atan2(dz, dx); it.root.add(holder);
      }
      for (const h of hoseSpecs) { const m = buildHose(h); if (m) { m.visible = !h.hidden; group.add(m); hoses.push(m); } }
      for (const it of items.values()) if (it.spec.hidden) it.root.visible = false;
      /* a lagged (foil-wrapped) piece */
      for (const it of items.values()) if (it.spec.foil) {
        const box = boundsOf([it.root]), s = box.getSize(new THREE.Vector3());
        const wrap = new THREE.Mesh(new THREE.CylinderGeometry(Math.max(s.x, s.z) / 2 + 0.15, Math.max(s.x, s.z) / 2 + 0.15, s.y * 0.8, 20, 1, true),
          new THREE.MeshStandardMaterial({ color: 0xdfe3e6, roughness: 0.28, metalness: 0.35, side: THREE.DoubleSide }));
        wrap.position.y = s.y * 0.5; it.root.add(wrap);
      }
      frame(true);
      render();
    }

    /* ---------------- what is inside ---------------- */
    function bandGeometry(profile, y0, y1, meniscus) {
      const pts = [new THREE.Vector2(1e-3, y0), new THREE.Vector2(Math.max(1e-3, radiusAt(profile, y0) - WALL_T), y0)];
      for (const [r, y] of profile) if (y > y0 && y < y1) pts.push(new THREE.Vector2(Math.max(1e-3, r - WALL_T), y));
      const rTop = Math.max(1e-3, radiusAt(profile, y1) - WALL_T);
      pts.push(new THREE.Vector2(rTop, y1));
      if (meniscus) {                                            // a liquid climbs the glass: the surface dips a little in the middle
        const dip = Math.min(0.14, 0.05 * rTop, 0.3 * (y1 - y0));
        pts.push(new THREE.Vector2(rTop * 0.82, y1 - dip * 0.55), new THREE.Vector2(rTop * 0.45, y1 - dip * 0.92), new THREE.Vector2(1e-3, y1 - dip));
      } else pts.push(new THREE.Vector2(1e-3, y1));
      return new THREE.LatheGeometry(pts, 48);
    }
    /* front faces only: a liquid's far wall, seen from inside through its own top surface, would otherwise blend over
       the near faces in index order and show as a bright wedge */
    /* a solid in a vessel: matt, finely grained, and it does not let the light through the way a liquid does */
    function solidMat(color) {
      return new THREE.MeshStandardMaterial({ color: new THREE.Color(color == null ? 0xf2f0ea : color), roughness: 0.92, metalness: 0,
        bumpMap: grain, bumpScale: 0.012, envMapIntensity: 0.65, side: THREE.FrontSide });
    }
    function liquidMat(color, cloudy) {
      return new THREE.MeshPhysicalMaterial({ color: new THREE.Color(color == null ? 0xbcd6e2 : color), roughness: cloudy ? 0.8 : 0.12,
        ior: 1.33, clearcoat: cloudy ? 0 : 0.8, clearcoatRoughness: 0.08, envMapIntensity: cloudy ? 0.45 : 1.15,
        specularIntensity: cloudy ? 0.3 : 1,
        transparent: true, opacity: cloudy ? 0.96 : 0.74, depthWrite: false, side: THREE.FrontSide });
    }
    /* the heights of each band, from the contents: mL when given (through the vessel's own profile), else a level */
    function bandsFor(it, c) {
      const sh = it.shape, floor = sh.plate != null ? sh.plate : sh.profile[0][1], cap = sh.neck;
      const toY = (ml, level) => (ml != null ? heightFor(sh.profile, ml, sh.height) : floor + (cap - floor) * Math.max(0, Math.min(1, level || 0)));
      const out = [];
      let y = floor, below = 0;
      for (const s of [].concat(c.solid || [])) {                // a solid, or a stack of them (alumina, sand, a band)
        if (!(s.depth > 0)) continue;
        out.push({ y0: y, y1: y + s.depth, color: s.color || 0xf6f5f0, cloudy: true, solid: true }); y += s.depth;
      }
      const baseVol = volumeBelow(sh.profile, y);                // liquid stands on whatever is below it (the plate, a solid)
      const layers = c.layers && c.layers.length ? c.layers : (c.ml != null || c.level > 0 ? [{ ml: c.ml, level: c.level, color: c.color }] : []);
      for (const L of layers) {                                  // layers are listed bottom first; their mL add up
        below += L.ml != null ? L.ml : 0;
        const top = L.ml != null ? Math.max(y + 0.05, toY(baseVol + below)) : Math.max(y + 0.05, toY(null, L.level));
        out.push({ y0: y, y1: Math.min(top, sh.height), color: L.color != null ? L.color : c.color, cloudy: L.cloudy != null ? L.cloudy : c.cloudy });
        y = Math.min(top, sh.height);
      }
      return out;
    }
    function drawBands(it, bands) {
      for (const m of it.liquid) { it.root.remove(m); m.geometry.dispose(); m.material.dispose(); }
      it.liquid = [];
      bands.forEach((b, k) => {
        if (!(b.y1 > b.y0 + 1e-3)) return;
        const m = new THREE.Mesh(bandGeometry(it.shape.profile, b.y0, b.y1, k === bands.length - 1 && !b.solid),
          b.solid ? solidMat(b.color) : liquidMat(b.color, b.cloudy));
        m.renderOrder = 1; m.userData = { id: it.spec.id, part: 'liquid', layer: k };
        if (it.pouring) m.visible = false;                       // while it pours, the level-cut liquid stands in for the bands
        it.root.add(m); it.liquid.push(m);
      });
      it.bands = bands;
      it.surface = bands.length ? bands[bands.length - 1].y1 : null;
    }
    let drips = [];                    // falling condensate: {owner, parent, mesh, ring?, origin, fall, n}
    const onTopOf = id => [...items.values()].filter(o => o.spec.on && String(o.spec.on).split(':')[0] === id);
    function makeDrips(owner, parent, origin, fall, ring) {
      const n = 3, mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(0.1, 8, 6),
        new THREE.MeshPhysicalMaterial({ color: 0xe8f2f6, roughness: 0.05, transparent: true, opacity: 0.85, depthWrite: false }), n);
      mesh.renderOrder = 3; parent.add(mesh);
      let ringMesh = null;
      if (ring) {                                                // the reflux ring: where the vapour condenses
        ringMesh = new THREE.Mesh(new THREE.TorusGeometry(0.5, 0.14, 8, 28), new THREE.MeshBasicMaterial({ color: 0xf4f8fa, transparent: true, opacity: 0.35, depthWrite: false }));
        ringMesh.rotation.x = Math.PI / 2; ringMesh.position.copy(origin); ringMesh.renderOrder = 3; parent.add(ringMesh);
      }
      drips.push({ owner, parent, mesh, ring: ringMesh, origin, fall, n });
    }
    function attachCondensate(it) {
      for (const up of onTopOf(it.spec.id)) {
        const k = up.spec.kind;
        if (k === 'condenser' || k === 'aircondenser') {         // refluxing: drops run back down inside the tube
          const L = up.spec.length || (k === 'condenser' ? 20 : 12);
          makeDrips(it.spec.id, up.root, new THREE.Vector3(0, L * 0.45, 0), L * 0.45, true);
        }
        let head = k === 'distillhead' ? up : null;               // distilling: drops fall from the condenser's far end
        if (k === 'vigreux') head = onTopOf(up.spec.id).find(o => o.spec.kind === 'distillhead') || null;
        const cond = head && onTopOf(head.spec.id).find(o => o.spec.kind === 'condenser' && /:side$/.test(o.spec.on));
        if (cond) { cond.root.updateMatrixWorld(true); makeDrips(it.spec.id, group, cond.root.localToWorld(cond.ports.top.clone()), 2.2, false); }
      }
    }
    function clearEffects(it) {
      drips = drips.filter(d => { if (d.owner !== it.spec.id) return true; d.parent.remove(d.mesh); d.mesh.geometry.dispose(); d.mesh.material.dispose();
        if (d.ring) { d.parent.remove(d.ring); d.ring.geometry.dispose(); d.ring.material.dispose(); } return false; });
      for (const k of ['bubbles', 'stirbar', 'drops', 'fumes', 'coat']) if (it[k]) { it.root.remove(it[k].mesh); it[k].mesh.geometry.dispose(); it[k].mesh.material.dispose(); it[k] = null; }
      if (it.marks) { for (const m of it.marks) { m.parent.remove(m); m.geometry.dispose(); m.material.dispose(); } it.marks = null; }
    }
    function setContents(id, c, ms = 0) {
      const it = items.get(id);
      if (!it || !it.shape) return;
      if (it.contents) { it.lastContents = it.contents; it.lastSurface = it.surface; }   // a pour still draws what was just in it
      it.contents = c || null;
      clearEffects(it);
      const target = c ? bandsFor(it, c) : [];
      const from = it.bands || [];
      if (ms > 0 && (from.length || target.length)) {
        /* animate the heights: every band grows or shrinks from where it was (a new band from the old surface) */
        const start = from.map(b => Object.assign({}, b)), top0 = it.surface != null ? it.surface : (target[0] ? target[0].y0 : 0);
        anims.push({ t0: performance.now(), ms, step: t => {
          const e = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
          const bands = target.map((b, k) => {
            const s = start[k] || { y0: Math.min(b.y0, top0), y1: Math.min(b.y0, top0) };
            return Object.assign({}, b, { y0: s.y0 + (b.y0 - s.y0) * e, y1: s.y1 + (b.y1 - s.y1) * e });
          });
          drawBands(it, bands);
        }, done: () => { drawBands(it, target); addEffects(it, c); } });
      } else { drawBands(it, target); addEffects(it, c); }
      render();
    }
    function addEffects(it, c) {
      if (c && c.coat && it.shape.dishLid) {                     // crystals grown on the underside of the lid
        const mesh = new THREE.Mesh(new THREE.CylinderGeometry(it.shape.body * 0.72, it.shape.body * 0.72, 0.1, 36),
          new THREE.MeshStandardMaterial({ color: new THREE.Color(c.coat.color || 0xf0a040), roughness: 0.8 }));
        mesh.position.y = it.shape.height - 0.13; mesh.userData = { id: it.spec.id, part: 'liquid' }; it.root.add(mesh);
        it.coat = { mesh };
      }
      if (c && c.plate && it.tlcPlate) {                         // a developed (or developing) TLC plate
        const { W, H, origin } = it.tlcPlate.userData.size, run = H - origin - 0.6, front = Math.max(0, Math.min(1, c.plate.front || 0));
        const lanes = c.plate.lanes || 3, marks = [];
        const add = m => { it.tlcPlate.add(m); marks.push(m); };
        if (front > 0) {
          const f = new THREE.Mesh(new THREE.BoxGeometry(W * 0.96, 0.03, 0.07), new THREE.MeshBasicMaterial({ color: 0x6f8fa8 }));
          f.position.y = origin + front * run; add(f);
          const wet = new THREE.Mesh(new THREE.BoxGeometry(W * 0.98, origin + front * run - 0.1, 0.065),
            new THREE.MeshBasicMaterial({ color: 0xdfe7ee, transparent: true, opacity: 0.45, depthWrite: false }));
          wet.position.y = (origin + front * run - 0.1) / 2 + 0.05; add(wet);
        }
        for (const sp of c.plate.spots || []) {
          const x = -W / 2 + (W / lanes) * ((sp.lane || 0) + 0.5), y = origin + (sp.rf || 0) * front * run;
          const d = new THREE.Mesh(new THREE.CylinderGeometry(0.16, 0.16, 0.08, 14), new THREE.MeshBasicMaterial({ color: new THREE.Color(sp.color || 0x5b4a6b) }));
          d.rotation.x = Math.PI / 2; d.position.set(x, y, 0); add(d);
        }
        it.marks = marks;
      }
      if (!c || it.surface == null) return;
      const floor = it.bands.length ? it.bands[0].y0 : 0, rMax = Math.max(0.3, radiusAt(it.shape.profile, (floor + it.surface) / 2) - WALL_T - 0.2);
      if (c.bubbling) attachCondensate(it);
      if (c.bubbling) {
        const n = 26, mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(0.12, 10, 8), new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.05, transparent: true, opacity: 0.4, depthWrite: false }), n);
        mesh.renderOrder = 2; it.root.add(mesh);
        /* a handful of nucleation sites on the floor, each releasing bubbles of its own size and rate — not an even
           curtain of identical spheres */
        const sites = Array.from({ length: 5 }, (_, j) => ({ a: j * 2.39 + 0.4, r: rMax * Math.sqrt((j + 0.4) / 5) * 0.75 }));
        it.bubbles = { mesh, floor, top: it.surface, seeds: Array.from({ length: n }, (_, i) => {
          const st = sites[i % sites.length];
          return { x: Math.cos(st.a) * st.r, z: Math.sin(st.a) * st.r, t: (i * 0.618) % 1, s: 0.55 + (i % 7) * 0.09, sc: 0.55 + (i % 5) * 0.22, w: (i % 3) * 0.05 + 0.04, ph: i * 1.7 };
        }) };
      }
      if (c.stirring) {
        const mesh = new THREE.Mesh(new THREE.CapsuleGeometry(0.28, Math.min(2.4, rMax * 1.2), 4, 10), new THREE.MeshStandardMaterial({ color: 0xf8f8f8, roughness: 0.3 }));
        mesh.rotation.z = Math.PI / 2; mesh.position.y = floor + 0.35; it.root.add(mesh);
        it.stirbar = { mesh };
      }
      if (c.fuming) {                                            // a column of fumes above the neck, in the given colour
        const n = 18, mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(0.9, 10, 8),
          new THREE.MeshBasicMaterial({ color: new THREE.Color(c.fuming === true ? 0x9a4a1c : c.fuming), transparent: true, opacity: 0.28, depthWrite: false }), n);
        mesh.renderOrder = 4; it.root.add(mesh);
        it.fumes = { mesh, y0: it.shape.height, seeds: Array.from({ length: n }, (_, i) => ({ t: i / n, a: i * 2.1 })) };
      }
      if (c.dripping && it.shape.tap) {
        const mesh = new THREE.Mesh(new THREE.SphereGeometry(0.16, 10, 8), new THREE.MeshBasicMaterial({ color: new THREE.Color(c.color || 0xcccccc) }));
        it.root.add(mesh); it.drops = { mesh };
      }
    }

    function setVisible(id, on) {
      const it = items.get(id);
      if (it) it.root.visible = !!on;
      else for (const h of hoses) if (h.userData.id === id) h.visible = !!on;
      render();
    }
    /* a thermometer's red column, for a reading in °C (its scale runs 0–250 °C along the stem; 0–400 °C on a
       melting-point apparatus) */
    function setReading(id, degC) {
      const it = items.get(id); if (!it || !it.column) return;
      const L = it.length || 26, h = Math.max(0.3, Math.min(1, (+degC || 0) / (it.scaleMax || 250))) * (L - 1);
      it.column.scale.y = h; it.column.position.y = 0.3 + h / 2; it.reading = degC; render();
    }
    function setHeat(id, w) {
      const it = items.get(id); if (!it) return;
      it.heat = Math.max(0, Math.min(1, +w || 0));
      if (it.flame) it.flame.forEach(f => { f.visible = it.heat > 0; f.scale.set(1, 0.5 + 0.7 * it.heat, 1); });
      if (it.steam) {
        if (it.heat > 0 && !it.fumes) {
          const n = 14, mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(0.9, 10, 8), new THREE.MeshBasicMaterial({ color: 0xf2f4f5, transparent: true, opacity: 0.15, depthWrite: false }), n);
          mesh.renderOrder = 4; it.root.add(mesh);
          it.fumes = { mesh, y0: 8.3, r0: 7.2, rise: 6, grow: 1.2, seeds: Array.from({ length: n }, (_, i) => ({ t: i / n, a: i * 2.3 })) };   // light wisps: a flask on the bath stays visible
        } else if (!(it.heat > 0) && it.fumes) { it.root.remove(it.fumes.mesh); it.fumes.mesh.geometry.dispose(); it.fumes.mesh.material.dispose(); it.fumes = null; }
      } else if (it.plate && !it.noGlow) {
        /* nothing on a lab bench glows at these temperatures (metal shows red only above about 500 °C): a faint warm
           tint, and the heater's own indicator lamp, say that it is on */
        it.plate.material.emissive = new THREE.Color(0xff5a1f); it.plate.material.emissiveIntensity = 0.12 * it.heat;
      }
      if (it.cavityLamp) {                                       // a microwave oven running: the cavity lit behind the window
        it.cavityLamp.material.emissive = new THREE.Color(it.heat > 0 ? 0xfff1c0 : 0x000000); it.cavityLamp.material.emissiveIntensity = it.heat > 0 ? 1.4 : 0;
        it.window.material.color = new THREE.Color(it.heat > 0 ? 0x8a7a40 : 0x20262a); it.window.material.opacity = it.heat > 0 ? 0.22 : 0.3;
      }
      if (it.uv && it.lamp) {                                    // a 254 nm cabinet: the inside lights violet-white
        it.lamp.material.emissive = new THREE.Color(it.heat > 0 ? 0xb9a8ff : 0x000000); it.lamp.material.emissiveIntensity = it.heat > 0 ? 2.2 : 0;
      } else if (it.lamp) { it.lamp.material.color = new THREE.Color(it.heat > 0 ? 0xff3b2f : 0x3a1210); it.lamp.material.emissive = new THREE.Color(it.heat > 0 ? 0xff2a1a : 0x000000); it.lamp.material.emissiveIntensity = it.heat > 0 ? 1.2 : 0; }
      if (it.oil) it.oil.material.color = new THREE.Color(COL.oil).lerp(new THREE.Color(0xe0a93a), it.heat);
      render();
    }

    /* a burner flame's colour: a Beilstein test turns it green; null is the ordinary blue cone */
    function setFlame(id, color) {
      const it = items.get(id); if (!it || !it.flame) return;
      it.flame[0].material.color = new THREE.Color(color == null ? 0x4f82d8 : color);
      it.flame[1].material.color = new THREE.Color(color == null ? 0x9fd8ff : new THREE.Color(color).lerp(new THREE.Color(0xffffff), 0.5));
      render();
    }

    /* pour: lift the vessel over the other one, tip it, hold, and put it back */
    function pour(fromId, toId, ms = 2600) {
      const a = items.get(fromId), b = items.get(toId);
      if (!a || !b) return;
      const home = a.home.clone(), homeRot = a.homeRot.clone();
      /* the target's mouth: a vessel's rim, or the top of anything else (a condenser, a funnel) */
      const topY = b.shape ? b.root.position.y + b.shape.height : boundsOf([b.root]).max.y;
      const mouth = new THREE.Vector3(b.root.position.x, topY, b.root.position.z);
      /* the pouring lip is the rim on the side facing the target; the glass tips about that lip, which stays put just
         above the middle of the target's mouth, so the stream falls straight in */
      const side = mouth.x >= home.x ? 1 : -1, tipMax = 1.8;
      const H = a.shape ? a.shape.height : 6, rim = a.shape ? radiusAt(a.shape.profile, H) : 1;
      const lipLocal = new THREE.Vector3(side * rim, H, 0);
      const rotOf = th => new THREE.Matrix4().makeRotationZ(-side * th);
      const offset = th => lipLocal.clone().applyMatrix4(rotOf(th));          // lip relative to the glass's base
      /* lift the lip high enough that no part of the tipped glass dips below the target's rim */
      let low = 0;
      if (a.shape) for (const [r, y] of a.shape.profile) for (const sx of [-1, 1]) {
        const d = new THREE.Vector3(sx * r, y, 0).sub(lipLocal).applyMatrix4(rotOf(tipMax));
        low = Math.min(low, d.y);
      }
      const anchor = mouth.clone(); anchor.y = topY + Math.max(1.2, 0.8 - low);
      const baseAt = th => anchor.clone().sub(offset(th));
      const start = baseAt(0);
      let stream = null, splash = null;
      /* the liquid stays level as the glass tips: the glass's whole inside is drawn with the liquid's colour and cut by a
         horizontal plane. The plane stays where the surface was until the lip comes down to it; from then on it follows
         the lip, so what is left is the liquid still below the lip — and it runs out as the glass goes over */
      const held = a.contents || a.lastContents, surf = a.contents ? a.surface : a.lastSurface;
      let level = null, cut = null, plane = null;
      if (held && surf != null && a.shape) {
        a.pouring = true; a.liquid.forEach(m => { m.visible = false; });
        const floor = a.shape.plate != null ? a.shape.plate : a.shape.profile[0][1];
        const col = held.color != null ? held.color : (held.layers && held.layers[0] && held.layers[0].color);
        const mat = liquidMat(col, held.cloudy); mat.side = THREE.DoubleSide;
        plane = new THREE.Plane(new THREE.Vector3(0, -1, 0), 0); mat.clippingPlanes = [plane];
        cut = new THREE.Mesh(bandGeometry(a.shape.profile, floor, a.shape.height - 0.03, false), mat);
        cut.renderOrder = 1; a.root.add(cut);
      }
      const lipW = new THREE.Vector3();
      const setLevel = (t, th) => {
        if (!cut) return;
        a.root.updateMatrixWorld(true); lipW.copy(lipLocal).applyMatrix4(a.root.matrixWorld);
        if (th === 0 && t < 0.5) level = a.root.position.y + surf;          // upright: the surface where it was
        else if (t < 0.8) level = Math.min(level, lipW.y - 0.05);         // tipped: never above the lip
        else level = -1e4;                                                // poured out
        plane.constant = level;
      };
      const ease = x => x < 0.5 ? 2 * x * x : 1 - Math.pow(-2 * x + 2, 2) / 2;
      const arc = (p, q, u) => {                                 // up first, across, then down: a quadratic curve over the top
        const c = p.clone().lerp(q, 0.5); c.y = Math.max(p.y, q.y) + 5;
        return p.clone().multiplyScalar((1 - u) * (1 - u)).add(c.multiplyScalar(2 * u * (1 - u))).add(q.clone().multiplyScalar(u * u));
      };
      anims.push({ t0: performance.now(), ms, step: t => {
        let th = 0, pos;
        if (t < 0.3) pos = arc(home, start, ease(t / 0.3));                              // lift, carry over in an arc
        else if (t < 0.8) { th = tipMax * ease(Math.min(1, (t - 0.3) / 0.15)); pos = baseAt(th); }   // tip about the lip
        else { const k = ease((t - 0.8) / 0.2); th = tipMax * (1 - Math.min(1, k * 2)); pos = k < 0.5 ? baseAt(th) : arc(start, home, (k - 0.5) * 2); }
        a.root.position.copy(pos);
        a.root.rotation.set(homeRot.x, homeRot.y, homeRot.z - side * th);
        setLevel(t, th);
        const pouring = th > tipMax * 0.6 && t < 0.85;
        const liquid = a.contents || a.lastContents;
        if (pouring && !stream && liquid) {
          const col = liquid.color != null ? liquid.color : (liquid.layers && liquid.layers[0] && liquid.layers[0].color);
          stream = new THREE.Mesh(new THREE.CylinderGeometry(0.2, 0.1, 1, 12), liquidMat(col, liquid.cloudy));
          stream.renderOrder = 1; group.add(stream);           // wide at the lip, drawn thin by gravity
          splash = new THREE.Mesh(new THREE.TorusGeometry(0.35, 0.07, 6, 20),
            new THREE.MeshBasicMaterial({ color: new THREE.Color(col == null ? 0xbcd6e2 : col), transparent: true, opacity: 0.5, depthWrite: false }));
          splash.rotation.x = Math.PI / 2; splash.renderOrder = 2; group.add(splash);
        }
        if (stream) {
          const lip = anchor.clone();
          const land = b.surface != null ? b.root.position.y + b.surface : b.shape ? b.root.position.y + 1 : topY - 2;
          const len = Math.max(0.5, lip.y - land);
          stream.scale.set(1, len, 1); stream.position.set(lip.x, lip.y - len / 2, lip.z);
          stream.visible = pouring;
          if (splash) {                                         // a ring spreading where the stream breaks the surface
            const w = (t * 3) % 1;
            splash.position.set(lip.x, land + 0.05, lip.z);
            splash.scale.setScalar(0.4 + 1.5 * w);
            splash.material.opacity = pouring ? 0.45 * (1 - w) : 0;
            splash.visible = pouring;
          }
        }
      }, done: () => { a.root.position.copy(home); a.root.rotation.copy(homeRot); if (stream) { group.remove(stream); stream.geometry.dispose(); stream.material.dispose(); }
        if (splash) { group.remove(splash); splash.geometry.dispose(); splash.material.dispose(); }
        if (cut) { a.root.remove(cut); cut.geometry.dispose(); cut.material.dispose(); }
        a.pouring = false; a.liquid.forEach(m => { m.visible = true; }); } });
      render();
    }

    /* a reagent going in: a graduated cylinder pours it, a Pasteur pipette drips it, or a weighing boat tips a solid in.
       The tool appears on the bench beside the vessel, travels over it in an arc, delivers, and goes back. The mouth is
       the top of whatever stands on the vessel's neck (a condenser, a funnel): reagents go in there, as at the bench. */
    let tempN = 0;
    const STACKABLE = new Set(['condenser', 'aircondenser', 'adapter', 'claisen', 'funnel', 'buchner', 'hirsch']);
    function mouthOf(id) {
      let cur = items.get(id);
      for (let g = 0; g < 6; g++) {
        const up = onTopOf(cur.spec.id).find(o => STACKABLE.has(o.spec.kind) && !/:(side|neckL|neckR)$/.test(String(o.spec.on)));
        if (!up) break;
        cur = up;
      }
      return cur;
    }
    const canDispense = id => { const it = items.get(id); return !!(it && it.shape && it.root.visible); };
    function dispense(targetId, how, o = {}) {
      if (!canDispense(targetId)) return false;
      const tgt = items.get(targetId), top = mouthOf(targetId), ms = o.ms || 2200, tp = tgt.root.position;
      const home = new THREE.Vector3(Math.max(-44, Math.min(44, tp.x - 8)), 0, Math.max(-22, Math.min(22, tp.z + 8)));
      const side = top.root.position.x >= home.x ? 1 : -1;
      if (how === 'pour') {                                     // a graduated cylinder of the right size, poured in
        const ml = Math.max(0.2, +o.ml || 5), vol = [10, 25, 50, 100, 250, 500].find(v => v >= ml * 1.25) || 1000;
        const spec = { id: '__reagent' + (++tempN), kind: 'cylinder', volume: vol };
        const built = buildVessel(spec);
        built.root.position.copy(home); group.add(built.root);
        items.set(spec.id, Object.assign({ spec, contents: null, liquid: [], home: home.clone(), homeRot: built.root.rotation.clone() }, built));
        setContents(spec.id, { ml: Math.min(ml, vol * 0.9), color: o.color });
        pour(spec.id, top.spec.id, ms);
        anims.push({ t0: performance.now(), ms: ms + 40, step: () => {}, done: () => { group.remove(built.root); disposeTree(built.root); items.delete(spec.id); } });
        return true;
      }
      const mouthY = top.shape ? top.root.position.y + top.shape.height : boundsOf([top.root]).max.y;
      const land = tgt.surface != null ? tp.y + tgt.surface : tp.y + (tgt.shape.plate != null ? tgt.shape.plate : 0.3);
      const over = new THREE.Vector3(top.root.position.x, mouthY + 1.2, top.root.position.z);
      const arc = (p, q, u) => { const c = p.clone().lerp(q, 0.5); c.y = Math.max(p.y, q.y) + 5;
        return p.clone().multiplyScalar((1 - u) * (1 - u)).add(c.multiplyScalar(2 * u * (1 - u))).add(q.clone().multiplyScalar(u * u)); };
      const ease = x => x < 0.5 ? 2 * x * x : 1 - Math.pow(-2 * x + 2, 2) / 2;
      const tool = new THREE.Group(); tool.position.copy(home); group.add(tool);
      const bits = new THREE.InstancedMesh(how === 'solid' ? new THREE.BoxGeometry(0.1, 0.1, 0.1) : new THREE.SphereGeometry(0.11, 8, 6),
        how === 'solid' ? new THREE.MeshStandardMaterial({ color: new THREE.Color(o.color || 0xf4f2ea), roughness: 0.8 })
          : new THREE.MeshPhysicalMaterial({ color: new THREE.Color(o.color || 0xdce8ee), roughness: 0.05, transparent: true, opacity: 0.85 }), how === 'solid' ? 18 : 3);
      bits.renderOrder = 3; bits.visible = false; group.add(bits);
      let bulb = null, mound = null;
      if (how === 'drops') {                                    // a Pasteur pipette: glass stem, fine tip, rubber bulb; the root is the tip
        const tip = glassMesh(new THREE.CylinderGeometry(0.13, 0.05, 3, 10), '__tool'); tip.position.y = 1.5; tool.add(tip);
        const stem = glassMesh(new THREE.CylinderGeometry(0.28, 0.13, 9, 12), '__tool'); stem.position.y = 7.5; tool.add(stem);
        bulb = new THREE.Mesh(new THREE.SphereGeometry(0.55, 14, 10), new THREE.MeshStandardMaterial({ color: 0xc0392b, roughness: 0.55 }));
        bulb.position.y = 12.6; bulb.scale.set(1, 1.5, 1); tool.add(bulb);
        const fill = new THREE.Mesh(new THREE.CylinderGeometry(0.2, 0.1, 5, 10), liquidMat(o.color, false)); fill.position.y = 4.5; tool.add(fill);
      } else {                                                  // a weighing boat with a heap of solid; the root is its pouring edge
        const boat = new THREE.Mesh(new THREE.BoxGeometry(2.8, 0.25, 2.6), new THREE.MeshStandardMaterial({ color: 0xf2f2ee, roughness: 0.6 }));
        boat.position.set(-side * 1.4, 0.12, 0); tool.add(boat);
        mound = new THREE.Mesh(new THREE.SphereGeometry(0.9, 16, 8, 0, Math.PI * 2, 0, Math.PI / 2), new THREE.MeshStandardMaterial({ color: new THREE.Color(o.color || 0xf4f2ea), roughness: 0.85 }));
        mound.scale.set(1, 0.45, 1); mound.position.set(-side * 1.4, 0.25, 0); tool.add(mound);
      }
      tool.traverse(m => { if (m.isMesh) m.userData = { id: '__tool', part: 'glass' }; });
      const seeds = Array.from({ length: 18 }, (_, i) => ({ dx: (Math.sin(i * 12.9) * 0.5), dz: (Math.cos(i * 7.3) * 0.5), t: (i * 0.618) % 1 }));
      anims.push({ t0: performance.now(), ms, step: t => {
        let pos = over;
        if (t < 0.3) pos = arc(home, over, ease(t / 0.3)); else if (t > 0.8) pos = arc(over, home, ease((t - 0.8) / 0.2));
        tool.position.copy(pos);
        const u = Math.max(0, Math.min(1, (t - 0.32) / 0.46));  // the delivery, 0 → 1
        if (how === 'drops') {
          bulb.scale.set(1, 1.5 - 0.35 * Math.abs(Math.sin(u * Math.PI * 3)), 1);
          bits.visible = t > 0.3 && t < 0.82;
          for (let i = 0; i < 3; i++) { const k = Math.max(0, Math.min(1, (u - i * 0.3) / 0.35));
            M.makeTranslation(over.x, over.y - (over.y - land) * k * k, over.z); M.scale(new THREE.Vector3(k > 0 && k < 1 ? 1 : 0, k > 0 && k < 1 ? 1.3 : 0, k > 0 && k < 1 ? 1 : 0)); bits.setMatrixAt(i, M); }
        } else {
          tool.rotation.z = -side * 1.0 * (t < 0.3 ? 0 : t < 0.42 ? ease((t - 0.3) / 0.12) : t < 0.8 ? 1 : 1 - ease((t - 0.8) / 0.12));
          mound.scale.set(1 - 0.7 * u, 0.45 * (1 - 0.7 * u), 1 - 0.7 * u);
          bits.visible = u > 0 && u < 1;
          seeds.forEach((sd, i) => { const k = (u * 2.2 + sd.t) % 1;
            M.makeTranslation(over.x + sd.dx * 0.4, over.y - (over.y - land) * k * k, over.z + sd.dz * 0.4); bits.setMatrixAt(i, M); });
        }
        bits.instanceMatrix.needsUpdate = true;
      }, done: () => { group.remove(tool); disposeTree(tool); group.remove(bits); bits.geometry.dispose(); bits.material.dispose(); } });
      render();
      return true;
    }

    /* ---------------- camera: orbit, zoom, glide ---------------- */
    const cam = { target: new THREE.Vector3(0, 12, 0), dist: 90, yaw: 0.25, pitch: 0.3 };
    let glide = null;
    function applyCamera() {
      const { target, dist, yaw, pitch } = cam;
      camera.position.set(target.x + dist * Math.cos(pitch) * Math.sin(yaw), target.y + dist * Math.sin(pitch), target.z + dist * Math.cos(pitch) * Math.cos(yaw));
      camera.lookAt(target);
    }
    let speed = 1;                                             // playback speed: camera glides run faster or slower
    function glideTo(target, dist, ms = 700) {
      ms /= speed;
      glide = { t0: performance.now(), ms, from: { t: cam.target.clone(), d: cam.dist }, to: { t: target, d: dist } };
      render();
    }
    function boundsOf(objs) {
      const box = new THREE.Box3();
      for (const o of objs) box.expandByObject(o);
      return box;
    }
    function frame(instant) {
      const box = boundsOf([...items.values()].map(it => it.root).concat(hoses));
      if (box.isEmpty()) { glideTo(new THREE.Vector3(0, 12, 0), 90, instant ? 0 : 600); return; }
      const c = box.getCenter(new THREE.Vector3()), s = box.getSize(new THREE.Vector3());
      const d = Math.max(30, 1.25 * Math.max(s.x, s.y * 1.5, s.z));
      if (instant) { cam.target.copy(c); cam.dist = d; applyCamera(); } else glideTo(c, d);
    }
    function focus(id, o = {}) {
      const list = [].concat(id || []).map(k => items.get(k)).filter(Boolean);   // one piece, or several framed together
      if (!list.length) { frame(false); return; }
      const box = boundsOf(list.map(it => it.root));
      const c = box.getCenter(new THREE.Vector3()), s = box.getSize(new THREE.Vector3());
      if (o.above) { c.y = box.max.y; }                          // leave room over the mouth (a reagent being added)
      glideTo(c, Math.max(o.above ? 38 : 16, 2.2 * Math.max(s.x, s.y, s.z) * (o.above ? 1.4 : 1)));
    }
    let ring = null;
    function highlight(id) {
      if (ring) { group.remove(ring); ring.geometry.dispose(); ring.material.dispose(); ring = null; }
      const it = id && items.get(id);
      if (it) {
        const box = boundsOf([it.root]), s = box.getSize(new THREE.Vector3()), c = box.getCenter(new THREE.Vector3());
        ring = new THREE.Mesh(new THREE.TorusGeometry(0.6 * Math.max(s.x, s.z) + 1.4, 0.18, 8, 56), new THREE.MeshBasicMaterial({ color: 0x3fbf8f }));
        ring.rotation.x = Math.PI / 2; ring.position.set(c.x, Math.max(0.2, box.min.y + 0.2), c.z);
        group.add(ring);
      }
      render();
    }

    /* ---------------- pointer: orbit, zoom, pick ---------------- */
    const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
    /* three.js ray casting ignores visibility, so a hidden piece (a bath taken away) must be skipped by hand */
    const shown = o => { for (let n = o; n && n !== group; n = n.parent) if (!n.visible) return false; return true; };
    function pickAt(ev) {
      const r = renderer.domElement.getBoundingClientRect();
      ndc.set(((ev.clientX - r.left) / r.width) * 2 - 1, -((ev.clientY - r.top) / r.height) * 2 + 1);
      ray.setFromCamera(ndc, camera);
      const hits = ray.intersectObjects(group.children, true);
      /* pick what the eye sees: clear glass is looked through, so the first liquid along the ray wins even behind
         the glass of another vessel; failing that, the first glass; anything opaque stops the ray */
      let glass = null;
      for (const h of hits) {
        const u = h.object.userData;
        if (!u || !u.id || !shown(h.object)) continue;
        if (u.part === 'liquid') return u;
        if (u.part === 'glass') { if (!glass) glass = u; continue; }
        return glass || u;                                     // an opaque body (or a hose)
      }
      return glass;
    }
    let drag = null, lastHover = '';
    const el = renderer.domElement;
    el.addEventListener('pointerdown', e => {
      drag = { x: e.clientX, y: e.clientY, yaw: cam.yaw, pitch: cam.pitch, moved: false, pan: e.button === 2 || e.shiftKey, t: cam.target.clone() };
      el.setPointerCapture(e.pointerId); el.style.cursor = 'grabbing'; glide = null;
    });
    el.addEventListener('contextmenu', e => e.preventDefault());
    el.addEventListener('pointermove', e => {
      if (drag) {
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
        if (drag.pan) {                                         // shift-drag or right-drag pans
          const k = cam.dist / 600, right = new THREE.Vector3(Math.cos(cam.yaw), 0, -Math.sin(cam.yaw));
          cam.target.copy(drag.t).addScaledVector(right, -dx * k).add(new THREE.Vector3(0, dy * k, 0));
        } else {
          cam.yaw = drag.yaw - dx * 0.008;
          cam.pitch = Math.max(-0.05, Math.min(1.3, drag.pitch + dy * 0.006));
        }
        applyCamera(); render(); return;
      }
      const u = pickAt(e);
      const k = u ? u.id + ':' + u.part : '';
      el.style.cursor = u ? 'pointer' : 'grab';
      if (k !== lastHover) { lastHover = k; if (opts.onHover) opts.onHover(u ? u.id : null, u ? u.part : null, e); }
      else if (u && opts.onHover) opts.onHover(u.id, u.part, e, true);
    });
    el.addEventListener('pointerup', e => {
      if (drag && !drag.moved) { const u = pickAt(e); if (opts.onPick) opts.onPick(u ? u.id : null, u ? u.part : null); }
      drag = null; el.style.cursor = 'grab';
    });
    el.addEventListener('pointercancel', () => { drag = null; });
    el.addEventListener('pointerleave', () => { if (lastHover && opts.onHover) { lastHover = ''; opts.onHover(null, null); } });
    el.addEventListener('dblclick', e => { const u = pickAt(e); focus(u ? u.id : null); });
    el.addEventListener('wheel', e => {
      e.preventDefault(); glide = null;
      cam.dist = Math.max(6, Math.min(220, cam.dist * (1 + Math.sign(e.deltaY) * 0.11)));
      applyCamera(); render();
    }, { passive: false });

    /* ---------------- the loop: only draws while something moves ---------------- */
    let dirty = true, running = true;
    const t0 = performance.now();
    function render() { dirty = true; }
    const M = new THREE.Matrix4(), _p = new THREE.Vector3(), _q = new THREE.Quaternion(), _s = new THREE.Vector3();
    function loop() {
      if (!running) return;
      requestAnimationFrame(loop);
      const now = performance.now(), t = (now - t0) / 1000;
      let moving = false;
      if (glide) {
        const k = Math.min(1, (now - glide.t0) / Math.max(1, glide.ms)), e = k < 0.5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
        cam.target.copy(glide.from.t).lerp(glide.to.t, e); cam.dist = glide.from.d + (glide.to.d - glide.from.d) * e;
        applyCamera(); moving = true;
        if (k >= 1) glide = null;
      }
      if (anims.length) {
        moving = true;
        anims = anims.filter(a => { const k = Math.min(1, (now - a.t0) / a.ms); a.step(k); if (k >= 1) { if (a.done) a.done(); return false; } return true; });
      }
      for (const it of items.values()) {
        if (it.bubbles) {
          moving = true;
          const { mesh, seeds, floor, top } = it.bubbles, H = Math.max(0.2, top - floor);
          seeds.forEach((sd, i) => {
            const k = (t * sd.s + sd.t) % 1;
            const rise = Math.pow(k, 0.78);                     // quick off the bottom, then steadier as drag builds
            const sway = Math.sin(t * 5 + sd.ph) * sd.w * H * 0.25;
            const grow = sd.sc * (0.7 + 0.5 * k);               // it swells as the pressure above it drops
            const burst = k > 0.9 ? Math.max(0, (1 - k) / 0.1) : 1;   // and bursts at the surface
            const g = grow * burst;
            _p.set(sd.x + sway, floor + rise * H, sd.z + Math.cos(t * 4.3 + sd.ph) * sd.w * H * 0.2);
            M.compose(_p, _q, _s.set(g, g, g)); mesh.setMatrixAt(i, M);
          });
          mesh.instanceMatrix.needsUpdate = true;
        }
        if (it.stirbar) { moving = true; it.stirbar.mesh.rotation.y = t * 9 + 0.12 * Math.sin(t * 2.3); }   // it slips slightly in the liquid
        if (it.fumes) {
          moving = true;
          /* fumes rise from the neck (r0 0) or, for a steam bath, from all round its rim (r0 = the rim radius) */
          const { mesh, y0, seeds } = it.fumes, r0 = it.fumes.r0 || 0, rise = it.fumes.rise || 12;
          seeds.forEach((sd, i) => {
            const k = (t * 0.25 + sd.t) % 1;
            const a = (r0 ? sd.a : sd.a + t * 0.7) + Math.sin(t * 0.6 + sd.t * 9) * 0.35;   // the column wanders
            const r = r0 + 1.2 * k + Math.sin(t * 0.9 + sd.t * 12) * 0.35 * k;
            const g = (1 + (it.fumes.grow || 2) * k) * (k > 0.72 ? Math.max(0, (1 - k) / 0.28) : 1);   // thins away at the top
            _p.set(Math.cos(a) * r, y0 + rise * Math.pow(k, 0.85), Math.sin(a) * r);
            M.compose(_p, _q, _s.set(g, g, g)); mesh.setMatrixAt(i, M);
          });
          mesh.instanceMatrix.needsUpdate = true;
        }
        if (it.steam && it.heat > 0) { moving = true; }
        if (it.drops) { moving = true; const k = (t * 1.6) % 1; it.drops.mesh.position.set(0, -0.2 - 3.2 * k * k, 0); it.drops.mesh.visible = k < 0.9; }
        if (it.flame && it.heat > 0) {
          moving = true;
          const f = 1 + 0.05 * Math.sin(t * 23) + 0.03 * Math.sin(t * 37);
          const h = 1 + 0.07 * Math.sin(t * 17 + 1.3) + 0.04 * Math.sin(t * 31);
          it.flame[0].scale.set(f, (0.5 + 0.7 * it.heat) * h, f);
          it.flame[1].scale.set(2 - f, h, 2 - f);               // the inner cone steadies as the outer one flares
          it.flame[0].position.x = Math.sin(t * 9) * 0.04; it.flame[1].position.x = Math.sin(t * 9) * 0.03;
        }
        if (it.spinner && it.heat > 0) { moving = true; it.spinner.rotation.y = t * 12 * it.heat; }
      }
      for (const d of drips) {
        moving = true;
        for (let i = 0; i < d.n; i++) {
          const k = (t * (0.62 + 0.17 * i) + i * 0.41) % 1;
          const hang = k < 0.34 ? k / 0.34 : 1;                 // it gathers at the tip before it lets go
          const fall = k < 0.34 ? 0 : Math.pow((k - 0.34) / 0.66, 2);
          const g = 0.45 + 0.55 * hang;
          _p.set(d.origin.x, d.origin.y - d.fall * fall, d.origin.z);
          M.compose(_p, _q, _s.set(g, g * (1 + 1.1 * fall), g)); d.mesh.setMatrixAt(i, M);   // stretched by the fall
        }
        d.mesh.instanceMatrix.needsUpdate = true;
        if (d.ring) d.ring.scale.setScalar(1 + 0.06 * Math.sin(t * 3));
      }
      if (dirty || moving) { renderer.render(scene, camera); dirty = false; }
    }
    function resize() {
      const w = container.clientWidth || 640, h = container.clientHeight || 420;
      renderer.setSize(w, h, false);
      camera.aspect = w / h; camera.updateProjectionMatrix(); render();
    }
    const ro = typeof ResizeObserver === 'function' ? new ResizeObserver(resize) : null;
    if (ro) ro.observe(container);
    resize(); applyCamera(); loop();

    function dispose() {
      running = false;
      if (ro) ro.disconnect();
      for (const it of items.values()) disposeTree(it.root);
      renderer.dispose();
      if (el.parentNode) el.parentNode.removeChild(el);
    }
    return {
      setSpeed: v => { speed = Math.max(0.1, +v || 1); }, dispense, canDispense,
      setBench, setContents, setHeat, setFlame, setReading, setVisible, pour, focus, highlight, resize, dispose, render,
      contentsOf: id => (items.get(id) || {}).contents || null,
      isVisible: id => { const it = items.get(id); return !!(it && it.root.visible); },
      specOf: id => (items.get(id) || {}).spec || null,
      _debug: { items, hoses, scene, camera, cam, pickAt: (clientX, clientY) => pickAt({ clientX, clientY }), volumeOf: id => { const it = items.get(id); return it && it.shape ? volumeBelow(it.shape.profile, it.shape.neck) : null; } },
    };
  }
  return { create, SHAPES, _volumeBelow: volumeBelow, _heightFor: heightFor };
})();
