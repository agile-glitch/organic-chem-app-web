-- AgileCHEM usage counts: one row per visit to the home screen and per app opened, under the visitor's random username.
-- Run this once in your Supabase project: Dashboard > SQL Editor > New query > paste > Run.
-- Before running, change the dashboard password on the line marked  <-- CHANGE THIS  (usage.html asks for it).
--
-- What the public (anon) key can do with this setup:
--   * add a row to usage_events (nothing else: it cannot read, change or delete rows);
--   * call usage_summary(password), which returns the counts only when the password is right.

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

-- the dashboard password, kept where the public key cannot read it
create table if not exists public.usage_secret (id int primary key default 1 check (id = 1), password_hash text not null);
alter table public.usage_secret enable row level security;
create extension if not exists pgcrypto;
insert into public.usage_secret (id, password_hash)
values (1, crypt('change-me', gen_salt('bf')))                                   -- <-- CHANGE THIS
on conflict (id) do update set password_hash = excluded.password_hash;

-- counts per username and app, for usage.html; only with the right password
create or replace function public.usage_summary(p_password text)
returns table (username text, app text, opens bigint, first_seen timestamptz, last_seen timestamptz)
language plpgsql security definer set search_path = public, extensions as $$
begin
  if not exists (select 1 from public.usage_secret s where s.password_hash = crypt(p_password, s.password_hash)) then
    raise exception 'wrong password' using errcode = '28P01';
  end if;
  return query
    select e.username, coalesce(e.app, 'home'), count(*), min(e.created_at), max(e.created_at)
    from public.usage_events e
    group by e.username, coalesce(e.app, 'home')
    order by e.username, 2;
end $$;
revoke all on function public.usage_summary(text) from public;
grant execute on function public.usage_summary(text) to anon, authenticated;
