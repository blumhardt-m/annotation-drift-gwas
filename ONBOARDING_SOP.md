# ONBOARDING_SOP.md

Session-initialization procedure for the `annotation-drift-gwas` repository. Run at the start of **every** agent session and every human re-entry. It is a startup checklist, not a design document — its purpose is reproducible entry and protection against six-month-later confusion about what is frozen and what is in progress.

Rationale (why this exists): the project's own history produced a session that began in the wrong working tree (an empty scaffold whose git remote was an unrelated repo), against a stale memory index that misreported which files were on disk. This SOP makes that failure impossible to repeat silently: nothing runs until identity, baseline, and locks are confirmed.

Do not perform any implementation, edit, or pipeline run until every step below passes.

## Step 1 — Confirm repo, branch, and HEAD
- `git remote -v` → origin URL contains `annotation-drift-gwas` (NOT `zmet-nid` or any monorepo).
- `git rev-parse --show-toplevel` → this directory is the repo root, not a subdirectory nested in another repo.
- `git rev-parse --abbrev-ref HEAD` → an authorized feature branch (Study 3: `study3-ranked-enrichment`), never `main`/`master` for implementation work.
- `git merge-base --is-ancestor v1.0_submission HEAD` → HEAD descends from the frozen submission baseline.
- If any check fails: STOP. You are in the wrong place. Re-clone or re-checkout before doing anything else. (See the recovery precedent: the authoritative repo is `github.com/blumhardt-m/annotation-drift-gwas` @ tag `v1.0_submission`.)

## Step 2 — Run the gate
- Execute `./pre_auth_check.sh`.
- Exit 0 = authorized to proceed. Any non-zero exit = STOP and escalate; do not weaken, bypass, or comment out the check.
- The gate is the only authorized execution entry point for implementation work. Inspection (reading files, `git status`, `git log`, running `validate_config.py` to look) is always permitted before the gate; but no edit, stage run, or pipeline execution proceeds until `pre_auth_check.sh` exits 0.

## Step 3 — Verify locked files
- Confirm every artifact in `LOCKED_PROTOCOL_FILES.md` is present and, where hash-locked, matches its recorded SHA-256.
- Confirm the locked inventory is consistent with `validate_config.py:FROZEN_FIELDS` and `TASK.md` (governance-inventory consistency; a mismatch is REVIEW_TRIGGER #11). `pre_auth_check.sh` performs the structural cross-check, but a human re-entering the repo should eyeball that the three documents describe the same freeze.
- Do not modify any locked file. Modification is a human act requiring a `DECISIONS.md` entry and deliberate re-freeze.

## Step 4 — Read the contracts
- Read `TASK.md` — the active session scope: task ID, branch, allowed files, prohibited actions, definition of done.
- Read `EXECUTION_BOUNDARY.md` — the scientific-vs-implementation boundary and the "import, don't reimplement" rule.
- Read `AGENTS.md` if this is an agent session — the hard constraints.
- Work only within the files `TASK.md` authorizes.

## Step 5 — Gate on all checks
- Perform no implementation, edit, or run until Steps 1–4 have all passed.
- If any check is ambiguous or a required file is missing, escalate rather than proceed. An unspecified condition is a human decision, never an agent default.

## Step 6 — Stop on any trigger
- Throughout the session, halt immediately on any `REVIEW_TRIGGER.md` condition (integrity or statistical), including the Study 3 degenerate-GSEA case and the ORA consistency-check failure.
- On a trigger: stop, record the condition and evidence in the session log, surface to the human, consult `REVIEW_CHECKLIST.md`, and await instruction. Do not remediate frozen state autonomously; do not re-freeze.

---

Onboarding is complete only when Steps 1–4 pass and the active `TASK.md` scope is understood. Then, and only then, implementation may begin within that scope.
