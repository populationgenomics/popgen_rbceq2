#!/usr/bin/env bash
# PreToolUse hook for Bash commands.
#
# Enforces the team policy that prose alone can't guarantee:
#   - read-only gcloud only; write/delete cloud-storage commands are not allowed
#   - no execution against live systems (test-only) for now
#
# The permission deny-list in settings.json catches these when they start a
# command; this hook also catches them mid-command (e.g. after `cd x && ...`).
# Exit code 2 blocks the tool call and shows the message to the agent.
#
# Policy is deliberately conservative — the team plans to revisit it soon.
set -uo pipefail

input=$(cat)
cmd=$(printf '%s' "$input" \
  | python3 -c 'import sys, json; print(json.load(sys.stdin).get("tool_input", {}).get("command", ""))' \
  2>/dev/null || true)

[ -z "$cmd" ] && exit 0

block() {
  echo "BLOCKED by .claude/hooks/block-destructive-cloud.sh: $1" >&2
  exit 2
}

# Cloud-storage writes / deletes.
if echo "$cmd" | grep -Eq '(gcloud[[:space:]]+storage[[:space:]]+(rm|mv|rewrite)|gcloud[[:space:]]+storage[[:space:]]+(objects|buckets)[[:space:]]+delete|gsutil[[:space:]]+(rm|mv))'; then
  block "write/delete cloud-storage commands are not allowed (read-only gcloud only)."
fi

# Execution against live systems is off for now.
if echo "$cmd" | grep -Eq 'analysis-runner.*--access-level[[:space:]]+(standard|full)'; then
  block "no execution on live systems: use --access-level test (policy under review)."
fi

exit 0
