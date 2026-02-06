# HARNESS_PATCHLOG.md

Change log for `tools/run.sh` — the universal run harness.

---

## Format

Each entry follows:

```
### PATCH-<NNN> — <date> — <short title>
**Trigger:** Why the change was needed.
**Change:** What was modified.
**Files:** List of touched files.
**Tested:** How it was verified.
```

---

## Log

### PATCH-001 — 2026-02-06 — Initial creation
**Trigger:** Need a reproducible, logged execution harness across all research projects.
**Change:** Created `tools/run.sh` with CWD pinning, venv/conda auto-activation, pip-install guard, stdout/stderr capture, git-state recording, and JSON manifest output.
**Files:** `tools/run.sh`
**Tested:** Manual dry-run review.

---

<!-- Append new entries above this line -->
