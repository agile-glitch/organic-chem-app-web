-- AgileCHEM usage counts: one row per visit to the home screen and per app opened, under the visitor's random username.
-- Run this once in your Supabase project: Dashboard > SQL Editor > New query > paste > Run.
--
-- The public (anon) key in js/usage_config.js can only ADD a row to usage_events: it cannot read, change or delete anything.
-- You read the counts in the Supabase dashboard, signed in to your own account (Table Editor > usage_counts, or the SQL Editor).

create table if not exists public.usage_events (
  id         bigint generated always as identity primary key,
  username   text not null check (username ~ '^[A-Za-z]{2,20}[0-9]{1,4}$'),
  event      text not null check (event in ('visit', 'open')),
  app        text check (app is null or app in ('home', 'draw', 'sketcher', 'mol', 'retro', 'reactions', 'simulate', 'possib', 'valid')),
  created_at timestamptz not null default now()
);
create index if not exists usage_events_user_idx on public.usage_events (username);
create index if not exists usage_events_time_idx on public.usage_events (created_at);

alter table public.usage_events enable row level security;
drop policy if exists "anyone can add a usage event" on public.usage_events;
create policy "anyone can add a usage event" on public.usage_events for insert to anon, authenticated with check (true);
-- no select / update / delete policy: the public key cannot read or change the events

-- counts per username and app, for reading in the dashboard (not reachable with the public key)
create or replace view public.usage_counts with (security_invoker = true) as
  select username, coalesce(app, 'home') as app, count(*) as times, min(created_at) as first_seen, max(created_at) as last_seen
  from public.usage_events group by username, coalesce(app, 'home');
revoke all on public.usage_counts from anon, authenticated;

-- (a project set up with the earlier version of this file also had a password-protected summary function; it is no longer used)
drop function if exists public.usage_summary(text);
drop table if exists public.usage_secret;
