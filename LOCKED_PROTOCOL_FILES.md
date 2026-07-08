# LOCKED_PROTOCOL_FILES.md

Machine-readable freeze layer for the `annotation-drift-gwas` repository. Any file listed here is **immutable after lock**. Agents (Claude Code, Codex) must never modify, overwrite, delete, move, or re-freeze these. Modification is a human act requiring a `DECISIONS.md` entry and an explicit re-freeze; see `EXECUTION_BOUNDARY.md`.

This file supersedes nothing already in `LOCK_SUMMARY.md` / `prereg/frozen_params.sha256` — it is the agent-facing enumeration of the same freeze those artifacts already encode, plus the Study 3 additions. Where this list and `LOCK_SUMMARY.md` disagree, that disagreement is itself a `REVIEW_TRIGGER` condition (do not proceed).

## Layer 1 — Scientific immutability (parent study, v1.0_submission @ 37272b2)

### 1a. Scientific inputs (content-addressed, immutable)
- `data/processed/mappings/comparison.tsv` — the frozen locus-to-gene table (721,406 rows; 291,292 rsIDs; 4.15% reassignment). The sole required input to the Study 3 count statistic.
- GWAS Catalog raw file, SHA-256 `bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b`, release `v1.0.2-associations_e115_r2026-01-19_full`. When restored to `data/raw/gwas_catalog.tsv`, it is locked by hash.
- Reactome GMT v95, SHA-256 `a34f43fd763a6800a9a94310c416ae511e84724d87a3c255f35d13764d12bebb`.
- GENCODE v29 GTF, SHA-256 `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9`.
- GENCODE v45 GTF, SHA-256 `56a166e4a43c934a1910c0d09dec786ac048862fc65130731d41f3c7f451ae68`.

### 1b. Frozen algorithms (import-only; may be read and called, never edited)
Reimplementing any of these is forbidden — a re-implementation that diverges by one line is silent scientific drift.
- `src/12_ensembl_id_enrichment.py`
- `src/primary_reactome_export.py`
- `scripts/study2_v37_v45/export_primary_reactome.py`
- `src/04_map_variants.py` (nearest-TSS mapping)
- `src/16_flip_rate.py` (flip definition)
- `src/05_bootstrap.py` (bootstrap machinery)

### 1c. Protocol & preregistration
- `prereg/prereg.md`
- `prereg/frozen_params.yaml`
- `prereg/frozen_params.sha256`

  Reconciled 2026-07-08: removed `PROTOCOL_v1.3.md` and `LOCK_SUMMARY.md` — those were
  L2G-audit v1.3 documents, disposed as foreign artifacts (see "Foreign artifacts
  (disposition)" below). The annotation-drift study of record has no separate
  PROTOCOL/LOCK_SUMMARY: its protocol of record is `config.yaml` (live authoritative),
  locked via `prereg/frozen_params.yaml` + `prereg/frozen_params.sha256` and described in
  `prereg/prereg.md` §0.

### 1d. Published outputs (v1.0_submission)
- `results/**` (all committed tables, figures, and study1/study2 outputs as of commit `37272b2`)
- `manuscript/**`
- `annotation_drift_manuscript.pdf`

## Layer 1b — Study 3 freeze (locked at Study 3 registration)

- `prereg/study3_prereg_v1.0.md` — the frozen Study 3 preregistration (SHA recorded in `prereg/study3_prereg.sha256` at freeze).
- `prereg/study3_prereg.sha256`

## Governance runtime (self-locking)

These files define the freeze and are therefore themselves protocol-locked. An agent that could edit these could grant itself arbitrary permission; changes require human action + `DECISIONS.md`:
- `LOCKED_PROTOCOL_FILES.md` (this file)
- `EXECUTION_BOUNDARY.md`
- `AGENTS.md`
- `ONBOARDING_SOP.md`
- `REVIEW_TRIGGER.md`
- `REVIEW_CHECKLIST.md`
- `pre_auth_check.sh`
- `validate_config.py`

These files are mutually constraining: they are meant to describe the same freeze from different angles. A disagreement between them — a file locked here but not enforced by `validate_config.py`, an artifact in `TASK.md`'s allowed set that appears in the locked list, a frozen field absent from one document — is itself a governance failure requiring human review (REVIEW_TRIGGER #11), not something to reconcile silently.

## Explicitly NOT locked (Study 3 implementation surface — additive-only)

Agents may create and edit only these for the Study 3 build:
- `src/19_gene_counts.py`, `src/20_gsea_preranked.py`, `src/21_boundary_flip_localization.py` (new stages)
- `results/study3_ranked_enrichment/**` (new output directory)
- `TASK.md` (session scope; rewritten per session by human authorization)
- `RUN_LOG.md` (append-only; agents append run records, never rewrite history)

## Foreign artifacts (disposition — DISPOSE, do NOT relocate)

The working tree contains files foreign to this study: a **two-generations-obsolete
(v1.3) snapshot** of the L2G Instability Audit's governance/protocol layer, injected in
commit `90f6631` alongside legitimate Study 1/Study 2 work. These are **not locked and
never will be** — they are scheduled for disposal, not freeze.

**Why dispose, not relocate.** The sibling audit at
`~/Study/research/08_GENOMICS/l2g_instability_audit/` runs a *superseded* design: the live
study is protocol **v1.6** (v29→v45, RA-only pilot, kill-criterion gate, Phase 1 already
executed), having walked away from the v1.3 design (v38/v45, multi-trait RA/SCZ/BMI,
perturbation matrix, LD panels, MSigDB, spectral AUC). Relocating v1.3 into that tree would
inject a conflicting, obsolete protocol next to a frozen v1.2 brief — actively harmful, not
restorative. The v1.3 files carry no unique value: v1.6 is preserved in ≥4 places (live
tree, Insync backup, `~/Downloads/archive/pre_2026-04-11/`, Trash). **There is no valid home
to relocate to.** (This supersedes the earlier "verify-then-remove against the sibling repo"
instruction, which assumed the sibling was a valid destination.)

### Class 1 — Dispose outright (no annotation-drift value)
Remove, or quarantine to `~/Downloads/anndrift_v1.3_contamination_removed/`:
- `PROTOCOL_v1.3.md`, `LOCK_SUMMARY.md`, `TASKS_WEEK1.md`, `EXECUTION_MEMO_v1.0.md`, `OUTPUT_SCHEMA.md`
- `configs/` (`ld_reference.yaml`, `pathway_arbiter.yaml`, `perturbation_matrix.yaml`, `traits_to_run.yaml`)
- `results/locus_instability_atlas/`

### Class 2 — Discard audit content, recreate fresh (filename is one this framework expects)
- `DECISIONS.md` — current entries are audit DEC-001…006; recreate empty for annotation-drift governance decisions.
- `RUN_LOG.md` — recreate as this study's append-only run log.

### Class 3 — Extract annotation-drift content, THEN discard the shell (carve-out)
- `DATA_SOURCES.md` — mixed: real annotation-drift inputs (GENCODE v37/v45, Reactome, GWAS
  Catalog) alongside audit-only rows (L2G scores, LD panel, MSigDB, GTEx, ABC). Pull the
  annotation-drift rows into this study's fresh `DATA_SOURCES.md`; drop the audit rows.
- `MANUSCRIPT_POSITIONING.md` — audit-headed but frames *this* study as "paper 1"; that
  positioning may still be valid. Extract any paper-1 positioning worth keeping into this
  study's own doc before discarding the shell.

**Rule:** no foreign artifact is ever locked or relocated into a live L2G tree. Class 1
removed/quarantined; Class 2 recreated fresh; Class 3 content-extracted then discarded.
Disposal is recorded in `DECISIONS.md`.
