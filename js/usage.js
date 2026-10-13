/* Usage — a random username for each visitor and a count of the visits and of each app opened, sent to Supabase
   (supabase/usage.sql makes the table; the counts are read in the Supabase dashboard, see supabase/README.md). The username (e.g. BlueFalcon42) is made up in the
   browser and kept in localStorage, so the same browser keeps the same name; nothing else about the visitor is sent:
   only the username, what happened ('visit' or 'open') and which app. With no project set in js/usage_config.js,
   nothing leaves the browser.
     Usage.user        the username
     Usage.track(event, app)    event 'visit' (the home screen, once per browser tab) or 'open' (an app)            */
(() => {
  'use strict';
  const KEY = 'agilechem.user';
  const ADJ = ['Amber', 'Azure', 'Bold', 'Brave', 'Bright', 'Calm', 'Clever', 'Coral', 'Crimson', 'Swift', 'Gentle', 'Golden', 'Happy', 'Jade',
    'Keen', 'Lucky', 'Lunar', 'Mellow', 'Misty', 'Noble', 'Olive', 'Quiet', 'Rapid', 'Ruby', 'Sage', 'Silver', 'Solar', 'Sunny', 'Teal', 'Witty'];
  const NOUN = ['Atom', 'Benzene', 'Beaker', 'Carbon', 'Comet', 'Crystal', 'Dolphin', 'Falcon', 'Flask', 'Fox', 'Helium', 'Ion', 'Koala', 'Lynx',
    'Neon', 'Otter', 'Owl', 'Panda', 'Photon', 'Quark', 'Raven', 'Salmon', 'Tiger', 'Turtle', 'Wolf', 'Zinc', 'Proton', 'Ketone', 'Ester', 'Radical'];
  const pick = a => a[Math.floor(Math.random() * a.length)];
  const store = { get(k, s) { try { return s.getItem(k); } catch (e) { return null; } }, set(k, v, s) { try { s.setItem(k, v); } catch (e) { /* storage blocked */ } } };

  let user = store.get(KEY, window.localStorage);
  if (!/^[A-Za-z]{2,20}[0-9]{1,4}$/.test(user || '')) {
    user = pick(ADJ) + pick(NOUN) + (10 + Math.floor(Math.random() * 90));
    store.set(KEY, user, window.localStorage);
  }

  function track(event, app) {
    const C = window.USAGE_CONFIG || {};
    if (!C.url || !C.anonKey) return;
    try {
      fetch(C.url.replace(/\/+$/, '') + '/rest/v1/usage_events', {
        method: 'POST', keepalive: true,
        headers: { apikey: C.anonKey, Authorization: 'Bearer ' + C.anonKey, 'Content-Type': 'application/json', Prefer: 'return=minimal' },
        body: JSON.stringify({ username: user, event, app: app || null }),
      }).catch(() => { /* offline or blocked: the count is simply missed */ });
    } catch (e) { /* no fetch */ }
  }

  window.Usage = { user, track };
})();
