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

The task is to implement the locally frozen Study 3 ranked-enrichment contrast (stages 19–21; external OSF registration pending verification) strictly per `prereg/study3_prereg_v1.0.md`. The prereg is the scientific contract; this file and `EXECUTION_BOUNDARY.md` are the execution contract. Build to the prereg, additively, on `study3-ranked-enrichment`, importing frozen machinery, escalating on any ambiguity or trigger.

External OSF registration has not been verified from repository evidence. References to the locally frozen preregistration do not establish an immutable external registration and do not authorize confirmatory Study 3 execution. Confirmatory execution is prohibited until the investigator verifies the registration entity, effective timestamp, exact attached artifact, and byte identity — see `EXECUTION_BOUNDARY.md` "Registration hold".
