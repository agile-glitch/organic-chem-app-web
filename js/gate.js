/* Gate — a password in front of the Simulate and Possibilities pages.
   Each locked page is covered by a box asking for the password; what is typed shows as * (one per character) until Enter
   or "Open" checks it. The right password unlocks both pages until the browser tab is closed (sessionStorage, which the
   home screen's tabs share, so it is asked once per visit).
   The site has no server, so this is a light gate, not security: only a salted SHA-256 of the password is in this file,
   but anyone who edits the page in their browser can get past it. To change the password, put the hash of
   SALT + the new password in HASH:   python3 -c "import hashlib;print(hashlib.sha256(b'agilechem-gate:NEWPASSWORD').hexdigest())"   */
(() => {
  'use strict';
  const PAGES = ['simulate', 'possib'];
  const SALT = 'agilechem-gate:';
  const HASH = '04fc2e1763f1a738801e959a6842af0f4494b7f49a9a8bd4ebdae45b2dcc553f';
  const KEY = 'agilechem.gate';

  /* SHA-256 of a string (UTF-8), hex. Written out so it also works where crypto.subtle is missing (a page opened from a file). */
  function sha256(text) {
    const K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
      0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
      0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
      0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
      0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
      0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2];
    const H = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    const bytes = Array.from(new TextEncoder().encode(text)), bits = bytes.length * 8;
    bytes.push(0x80);
    while (bytes.length % 64 !== 56) bytes.push(0);
    for (let i = 7; i >= 0; i--) bytes.push(i >= 4 ? 0 : (bits >>> (i * 8)) & 0xff);
    const rot = (x, n) => (x >>> n) | (x << (32 - n)), W = new Array(64);
    for (let o = 0; o < bytes.length; o += 64) {
      for (let i = 0; i < 16; i++) W[i] = (bytes[o + 4 * i] << 24) | (bytes[o + 4 * i + 1] << 16) | (bytes[o + 4 * i + 2] << 8) | bytes[o + 4 * i + 3];
      for (let i = 16; i < 64; i++) {
        const s0 = rot(W[i - 15], 7) ^ rot(W[i - 15], 18) ^ (W[i - 15] >>> 3), s1 = rot(W[i - 2], 17) ^ rot(W[i - 2], 19) ^ (W[i - 2] >>> 10);
        W[i] = (W[i - 16] + s0 + W[i - 7] + s1) | 0;
      }
      let [a, b, c, d, e, f, g, h] = H;
      for (let i = 0; i < 64; i++) {
        const t1 = (h + (rot(e, 6) ^ rot(e, 11) ^ rot(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + W[i]) | 0;
        const t2 = ((rot(a, 2) ^ rot(a, 13) ^ rot(a, 22)) + ((a & b) ^ (a & c) ^ (b & c))) | 0;
        h = g; g = f; f = e; e = (d + t1) | 0; d = c; c = b; b = a; a = (t1 + t2) | 0;
      }
      [a, b, c, d, e, f, g, h].forEach((v, i) => { H[i] = (H[i] + v) | 0; });
    }
    return H.map(v => (v >>> 0).toString(16).padStart(8, '0')).join('');
  }

  const unlocked = () => { try { return sessionStorage.getItem(KEY) === HASH; } catch (e) { return false; } };
  const remember = () => { try { sessionStorage.setItem(KEY, HASH); } catch (e) { /* storage blocked: asked again next time */ } };
  const boxes = [];

  /* a text box that shows * for each character and keeps what was typed */
  function maskedInput(input) {
    let real = '';
    const show = (caret) => { input.value = '*'.repeat(real.length); input.setSelectionRange(caret, caret); };
    input.addEventListener('beforeinput', e => {
      const s = input.selectionStart, t = input.selectionEnd, data = e.data || '';
      e.preventDefault();
      if (e.inputType.startsWith('insert')) {
        const add = e.inputType === 'insertFromPaste' && e.dataTransfer ? e.dataTransfer.getData('text') : data;
        if (e.inputType === 'insertLineBreak') return;
        real = real.slice(0, s) + add + real.slice(t); show(s + add.length);
      } else if (e.inputType === 'deleteContentBackward') {
        const from = s === t ? Math.max(0, s - 1) : s; real = real.slice(0, from) + real.slice(t); show(from);
      } else if (e.inputType === 'deleteContentForward') {
        const to = s === t ? Math.min(real.length, t + 1) : t; real = real.slice(0, s) + real.slice(to); show(s);
      } else if (e.inputType.startsWith('delete')) { real = real.slice(0, s) + real.slice(t); show(s); }
    });
    input.addEventListener('paste', e => {          // browsers without beforeinput paste data
      if (!e.clipboardData) return;
      e.preventDefault();
      const add = e.clipboardData.getData('text'), s = input.selectionStart, t = input.selectionEnd;
      real = real.slice(0, s) + add + real.slice(t); show(s + add.length);
    });
    return { get value() { return real; }, clear() { real = ''; input.value = ''; } };
  }

  function lock(section) {
    section.classList.add('gated');                     // css/style.css hides everything else in the page, also what its script adds later
    const box = document.createElement('div');
    box.className = 'gate';
    box.innerHTML = `<form class="gate-card" autocomplete="off">
        <div class="gate-icon" aria-hidden="true"><svg viewBox="0 0 24 24" width="28" height="28"><path fill="currentColor" d="M12 2a5 5 0 0 0-5 5v3H6a2 2 0 0 0-2 2v8a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2v-8a2 2 0 0 0-2-2h-1V7a5 5 0 0 0-5-5zm-3 8V7a3 3 0 1 1 6 0v3H9z"/></svg></div>
        <p class="gate-status">Under Maintenance</p>
        <h2>This page needs a password</h2>
        <label class="gate-label" for="gate-${section.id}">Password</label>
        <div class="gate-row"><input id="gate-${section.id}" class="gate-input" type="text" inputmode="text" spellcheck="false" autocapitalize="off" autocorrect="off" autocomplete="off"><button type="submit">Open</button></div>
        <p class="gate-msg" role="alert"></p>
      </form>`;
    section.prepend(box);
    const input = box.querySelector('input'), msg = box.querySelector('.gate-msg'), pw = maskedInput(input);
    const form = box.querySelector('form');
    input.addEventListener('input', () => { msg.textContent = ''; });
    input.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); check(); } });
    form.addEventListener('submit', e => { e.preventDefault(); check(); });
    function check() {
      if (sha256(SALT + pw.value) === HASH) { remember(); unlockAll(); return; }
      msg.textContent = pw.value ? 'That password is not right. Try again.' : 'Type the password, then press Enter.';
      pw.clear(); input.focus();
    }
    boxes.push({ section, box, input });
  }

  function unlockAll() {
    boxes.splice(0).forEach(({ section, box }) => { box.remove(); section.classList.remove('gated'); });
    window.dispatchEvent(new Event('resize'));          // pages that size a canvas while hidden (the 3D bench) measure again
  }

  function start() {
    if (unlocked()) return;
    PAGES.forEach(p => { const s = document.getElementById('page-' + p); if (s) lock(s); });
    // focus the box when its page is shown
    const focus = () => { const b = boxes.find(x => x.section.classList.contains('active')); if (b) setTimeout(() => b.input.focus(), 0); };
    window.addEventListener('hashchange', focus);
    document.addEventListener('click', e => { if (e.target.closest && e.target.closest('.tab')) focus(); });
    setTimeout(focus, 0);
  }

  window.Gate = { locked: page => PAGES.includes(page) && !unlocked(), sha256 };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();
})();
