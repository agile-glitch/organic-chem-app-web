# Usage counts with Supabase

Each visitor gets a random username (e.g. `BlueFalcon42`), made up in their browser and kept there, shown at the bottom of the home screen.
The site counts, under that name, each visit to the home screen and each app opened, and `usage.html` shows the counts per username and app.
Only the username, the kind of event (`visit` / `open`), which app and the time are stored: no IP address, no browser details, nothing typed.

## Set-up (once, about 5 minutes)

1. Create a free project at <https://supabase.com> (New project; any name and region; keep the database password somewhere safe).
2. Open **SQL Editor > New query**, paste all of [`usage.sql`](usage.sql), change `change-me` on the line marked `<-- CHANGE THIS` to the password
   you want for the usage page, and press **Run**.
3. Open **Project Settings > API** (or **Data API**) and copy two values into [`js/usage_config.js`](../js/usage_config.js):
   * `url`: the **Project URL**, like `https://abcdefghijklmnop.supabase.co`
   * `anonKey`: the **anon / public** key (starts with `eyJ`). Never the `service_role` key: that one can read and delete everything.
4. Commit and publish the site. Open `usage.html` on the site and enter the password from step 2.

The anon key is meant to be public: with the rules in `usage.sql` it can only add usage rows, and read the totals only through
`usage_summary` with the right password. It cannot read, change or delete the individual rows.

## Notes

* The same person on another browser or device, or after clearing site data, gets a new username.
* While `js/usage_config.js` is empty, nothing is sent anywhere.
* To change the usage-page password later, run in the SQL Editor:
  `update public.usage_secret set password_hash = crypt('new password', gen_salt('bf')) where id = 1;`
* To start the counts again: `truncate public.usage_events;`
