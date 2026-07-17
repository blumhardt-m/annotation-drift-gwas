# TASK.md

Active session authorization. Rewritten per session by the human. An agent runs only within this scope.

## Task ID
`study3-build-001`

## Objective
Implement the registered Study 3 ranked-enrichment contrast (count-based) strictly per `prereg/study3_prereg_v1.0.md` (SHA in `prereg/study3_prereg.sha256`; controlling artifact `29fe3ff1…` verified against the immutable embargoed OSF registration snapshot, approved 2026-07-09 — DEC-008/DEC-011).

## Branch
`study3-ranked-enrichment` (create from `review/study3-freeze-and-gate0-20260716`, not `main` — the reviewed protocol, checksum lineage, Gate 0 record, and DEC-004 through DEC-011 are not present on `main`; never commit to `main`).

## Allowed files (create/edit)
- `src/19_gene_counts.py` — build per-gene, per-version count vectors over the shared union universe from `comparison.tsv` (prereg §3.2).
- `src/20_gsea_preranked.py` — preranked fgsea per version, frozen parameters (prereg §3.3).
- `src/21_boundary_flip_localization.py` — primary endpoint C, degenerate-case rule, descriptive A, chromosome-block bootstrap (prereg §4).
- `results/study3_ranked_enrichment/**` — new outputs and diagnostics (prereg §6).
- `RUN_LOG.md` — append run records only.

## Prohibited actions
- Modifying anything in `LOCKED_PROTOCOL_FILES.md`.
- Reimplementing frozen machinery (import from `src/04`, `src/12`, `src/16`, `src/05`, `build_alias_table()`).
- Resolving any parameter not fully specified in the prereg — escalate instead.
- Running the `validate_config.py` re-freeze remediation, or any `.sha256` regeneration.
- Proceeding through any `REVIEW_TRIGGER.md` condition, including the degenerate-GSEA case and the ORA consistency-check failure.

## Preconditions
`pre_auth_check.sh` must exit 0 before any stage runs. If it fails, stop and escalate; do not weaken the check.

Gate 0 (registration verification): CLEARED 2026-07-17 (DEC-011). The external OSF registration of the controlling prereg (`29fe3ff1…`) was verified — initiated 2026-07-08, approved 2026-07-09, before any confirmatory execution. `pre_auth_check.sh` must still exit 0 (local lock integrity), and the amended §3.4 ORA consistency preflight (DEC-009) must pass against the frozen parent inputs before the contrast; either failure halts execution. See `EXECUTION_BOUNDARY.md` "Registration verification (Gate 0 — cleared)".

## Definition of done
Stages 19–21 implemented and run to the prereg; all §6 diagnostics emitted; H1 evaluated per §5 (or correctly halted at a trigger); results under `results/study3_ranked_enrichment/`; `RUN_LOG.md` appended with the run record and the Study 3 prereg SHA.
