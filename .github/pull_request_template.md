<!--
Write this yourself, for a reviewer who wasn't in the session. Short and
concrete — no chat/session context, no "as discussed".
-->

## Purpose

The intent this implements — the problem and the agreed approach (link it if
there's a written one). Keep the PR to one small, coherent slice.

## Approach & key decisions

How it's built, the design/architecture choices worth a reviewer's attention,
and anything that **deviates from the agreed approach** (and why).

## Outputs and side-effects

What this creates, modifies, or writes (files, buckets, records, infra).

## How it was checked

- [ ] pre-commit (lint/format/hygiene) and the type checker pass
- [ ] fresh-agent correctness review done — and a security review if it touches
      data access, secrets, or infrastructure
- Behaviours / edge cases covered, with evidence (paste the relevant output, not
  a wall of logs).

## Out of scope

What this deliberately does **not** do (and where it'll be handled instead).
