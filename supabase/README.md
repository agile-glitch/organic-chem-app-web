# Usage counts with Supabase

Each visitor gets a random username (e.g. `BlueFalcon42`), made up in their browser and kept there, shown at the bottom of the home screen.
The site counts, under that name, each visit to the home screen and each app opened. Only the username, the kind of event (`visit` / `open`),
which app and the time are stored: no IP address, no browser details, nothing typed.
The counts are only visible to you, in the Supabase dashboard; the site has no page that shows them.

## Set-up (once, about 5 minutes)

1. Create a free project at <https://supabase.com> (New project; any name and region; keep the database password somewhere safe).
2. Open **SQL Editor > New query**, paste all of [`usage.sql`](usage.sql) and press **Run**.
3. Open **Project Settings > API** (or **Data API**) and copy two values into [`js/usage_config.js`](../js/usage_config.js):
   * `url`: the **Project URL**, like `https://abcdefghijklmnop.supabase.co`
   * `anonKey`: the **anon / public** key (starts with `eyJ`). Never the `service_role` key: that one can read and delete everything.
4. Commit and publish the site.

The anon key is meant to be public: with the rules in `usage.sql` it can only add usage rows. It cannot read, change or delete anything.

## Reading the counts

Sign in to supabase.com, open the project, then either:

* **Table Editor > `usage_counts`**: one row per username and app, with how many times, first and last seen
  (`usage_events` has every single visit and app opened);
* or **SQL Editor**, for example:
  ```sql
  -- per user and app
  select * from usage_counts order by username, app;
  -- per app, all users
  select app, sum(times) as times, count(*) as users from usage_counts group by app order by times desc;
  -- the last 7 days
  select username, app, count(*) from usage_events where created_at > now() - interval '7 days' group by 1, 2 order by 1, 2;
  ```

## Notes

* The same person on another browser or device, or after clearing site data, gets a new username.
* While `js/usage_config.js` is empty, nothing is sent anywhere.
* To start the counts again: `truncate public.usage_events;`
