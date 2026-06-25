# `.claude/` — enforced execution policy

`settings.json` and the PreToolUse hook encode what an agent may run in a repo —
the binding execution policy that `AGENTS.md` norm #11 ("mind what you may
execute") points at. Prose alone isn't binding; the settings are. `settings.json`
also declares the CPG plugins and the SessionStart setup nudge.

The default policy here is conservative on purpose:

- **Allowed**: read-only context-gathering — read-only `gcloud` (`ls`, `cat`,
  `du`) and read-only shell.
- **Ask first**: `gcloud storage cp`, `gsutil` (a human approves each).
- **Denied / blocked**: cloud-storage `rm` / `mv` / `rewrite` / `delete`, and
  `analysis-runner --access-level standard|full` (no live execution for now).

This `.claude/` is the **baseline a new repo copies** (see
[`../docs/repo-setup.md`](../docs/repo-setup.md)) and then tailors — set each
repo's `allow` / `ask` / `deny` lists for what its agents may run. When you open
a repo, Claude Code asks you to trust the project hook.
