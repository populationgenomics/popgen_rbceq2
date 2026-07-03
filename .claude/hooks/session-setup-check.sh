#!/usr/bin/env bash
# SessionStart hook: a one-time nudge if the contributor hasn't done the
# user-global setup. Non-blocking (SessionStart can't block) — it points at
# docs/contributor-setup.md. The two CPG plugins are also declared in
# settings.json, so Claude Code offers to install those on trust; this catches
# the skills that can't be auto-prompted (the Matt Pocock collection).
set -uo pipefail

input=$(cat)
src=$(printf '%s' "$input" \
  | python3 -c 'import sys, json; print(json.load(sys.stdin).get("source", ""))' \
  2>/dev/null || true)

# Only nudge on a fresh start, not on resume/clear/compact.
[ "$src" = "startup" ] || exit 0

missing=()
[ -d "$HOME/.claude/plugins/marketplaces/cpg" ] || missing+=("cpg-skills plugin (cpg-infra grounding)")
[ -d "$HOME/.claude/skills/grill-me" ] || missing+=("workflow skills (grill-me / grill-with-docs / tdd)")

[ ${#missing[@]} -eq 0 ] && exit 0

msg="CPG contributor setup looks incomplete — missing: ${missing[*]}. See docs/contributor-setup.md for the one-time, user-global setup."

# stdout JSON surfaces the message to the user; stderr is a transcript fallback.
printf '%s' "$msg" >&2
python3 -c 'import json, sys; print(json.dumps({"systemMessage": sys.argv[1]}))' "$msg" 2>/dev/null || true
exit 0
