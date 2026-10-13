/* Where usage counts go (js/usage.js; set-up in supabase/README.md). Both values are in your Supabase project under
   Project Settings > API: the Project URL and the "anon" "public" key (the public one, never the service_role key).
   While either is empty, nothing is sent anywhere. */
window.USAGE_CONFIG = {
  url: '',          // e.g. 'https://abcdefghijklmnop.supabase.co'
  anonKey: '',      // the long key that starts with 'eyJ'
};
