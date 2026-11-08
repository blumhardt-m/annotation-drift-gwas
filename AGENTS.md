# AGENTS.md

Agent contract for the `annotation-drift-gwas` repository. Applies to Claude Code, Codex, and any autonomous coding agent. Reading this file and `EXECUTION_BOUNDARY.md` is mandatory at session start (`ONBOARDING_SOP.md` / `pre_auth_check.sh` enforce it).

## Roles

- **Human (Michael):** protocol authority. Sole party permitted to freeze, re-freeze, amend frozen parameters, edit governance-runtime files, or resolve escalations.
- **Agent:** implementation executor. Builds and edits only the non-locked implementation surface under an active `TASK.md`, on an authorized branch.

## Hard constraints

1. **Additive-only.** Never modify, overwrite, delete, or move any file in `LOCKED_PROTOCOL_FILES.md`. Create only the files the active `TASK.md` authorizes.
2. **Branch discipline.** Work only on the branch named in `TASK.md` (Study 3: `study3-ranked-enrichment`). Never commit to `main` or move the `v1.0_submission` tag.
3. **Import, don't reimplement.** Frozen computations (nearest-TSS mapping, `build_alias_table()`, flip definition, bootstrap, ORA machinery) must be imported and called. Reimplementing any of them is forbidden — a re-implementation that diverges by one line is silent scientific drift.
4. **Never resolve an underspecified parameter.** If the prereg or config does not fully specify a value the code needs, STOP and escalate. Do not choose a default, infer intent, or "make it work." An unspecified analytic choice is a human decision.
5. **Never re-freeze.** Do not run the `validate_config.py` re-freeze remediation, do not `cp config.yaml prereg/frozen_params.yaml`, do not regenerate any `.sha256` lock. See `EXECUTION_BOUNDARY.md`.
6. **Never downgrade a gate.** If `pre_auth_check.sh` or `validate_config.py` fails, do not bypass, comment out, or weaken the check. Failure means stop.
7. **Escalate on trigger.** On any condition in `REVIEW_TRIGGER.md`, halt and surface the condition to the human. Do not proceed through an invalid state, even if a fix looks obvious.

## Execution mode

- Default: **plan-then-confirm** for any action touching the filesystem outside the authorized implementation surface. Read-only inspection is always permitted.
- The agent may run the pipeline stages it is authorized to build, but must run `pre_auth_check.sh` first and abort on non-zero exit.

## Study 3 session summary (current)

The task is to implement the registered Study 3 ranked-enrichment contrast (stages 19–21; OSF-registered, controlling artifact `29fe3ff1…`, approved 2026-07-09 — DEC-008/DEC-011) strictly per `prereg/study3_prereg_v1.0.md`. The prereg is the scientific contract; this file and `EXECUTION_BOUNDARY.md` are the execution contract. Build to the prereg, additively, on `study3-ranked-enrichment`, importing frozen machinery, escalating on any ambiguity or trigger.

**Registration verified (Gate 0 cleared, 2026-07-17).** The temporary registration-verification hold is released (DEC-011). The controlling Study 3 prereg (SHA-256 `29fe3ff1…`) was byte-matched to the immutable embargoed OSF registration snapshot; registration was initiated 2026-07-08 and approved 2026-07-09, before any execution of stages 19–21. A local freeze is still never, by itself, an external registration, and `pre_auth_check.sh` passing verifies local lock integrity only. Stages 19–21 proceed under the registered sequencing (prereg §8) and the amended §3.4 ORA preflight (DEC-009), on `study3-ranked-enrichment` created from `review/study3-freeze-and-gate0-20260716`; the additive-only boundary and all prohibitions remain in force — see `EXECUTION_BOUNDARY.md` "Registration verification (Gate 0 — cleared)".


## Prospective per-trait pilot authority (DEC-016 candidate; inactive until adoption)

Baseline main: 847fd7cd98f68695b709ae155a1d2573ff7b7ee7. Installed DEC-015 controls prereg/pertrait_pilot_prereg_v1.1.md, SHA-256 b78907970faafa656c91ba55b2a948f405045dd8165fcceac9ee1bf70333d1db.
Order SHA-256: b28fdff1a374f3bfb17aaa3ca68ccd0bfe862fcac7385ae1bfe5c73e4a338043; rider SHA-256: fedfcce238d6ff66170ed45698d145fa31ee6c8b5115a46dc501c7d874907224.
INPUT_MANIFEST_rev2.json SHA-256: dc838286666ee2d35c88282a05a991570230373b627d3189e84457ece61f1288; copy its exact bytes to manifests/pertrait_pilot/INPUT_MANIFEST_rev2.json after adoption.
Rscript absolute path: /home/michaelblumhardt/miniconda3/envs/study3probe/bin/Rscript; R 4.5.3, fgsea 1.36.2, data.table 1.18.4, BiocParallel 1.44.0. Python and all dependency/input hashes are bound by that manifest.
Scientific scope: v29/v45 with Reactome only; ORA primary, preranked fgsea secondary, gene-to-gene C2, conditional chromosome deletion, and only the one rule-triggered extension. GENCODE v37/v45 and GO:BP are deferred by v1.1 §11 and have no implementation, restoration or execution authority here.

Upon exact owner adoption only, the new worktree /home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot on pertrait-pilot uses PILOT_EXECUTION_AUTHORITY.md and the replacement TASK.md as its scoped execution contract. Their pilot scope supersedes historical Study 3 task/branch language for this worktree only. Study 3 remains CLOSED. The owner must explicitly apply or narrowly delegate the exact six repository candidates; no standing protected-write authority is inferred. The sole guard edit supplies its missing config arguments, retaining every check; the additional pilot gate and all binding checks remain mandatory. No other locked file or scientific definition may be edited. Main and all earlier evidence remain intact.
