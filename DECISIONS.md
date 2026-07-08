# Decision Log — Annotation Drift in GWAS Locus-to-Gene Mapping

Append-only. Do not edit or delete existing entries; new entries at the bottom.
This log belongs to the **annotation-drift study of record** (protocol v1.0). It was
authored fresh on 2026-07-08; it is not the decision log of the L2G Instability Audit
(that study's DEC-001…006 live in its own archive).

---

## DEC-001 — Agent-governance layer adopted (2026-07-08)

Decision: Adopted the agent-governance package (`AGENTS.md`, `EXECUTION_BOUNDARY.md`,
`LOCKED_PROTOCOL_FILES.md`, `ONBOARDING_SOP.md`, `REVIEW_TRIGGER.md`,
`REVIEW_CHECKLIST.md`, `pre_auth_check.sh`, `TASK.md`) plus the `validate_config.py`
FROZEN_FIELDS patch, on branch `governance-adoption`, merged to `main` as baseline
infrastructure (not Study 3 work).
Rationale: The repository had no runtime governance; this installs the anti-drift gate,
freeze inventory, and review triggers before Study 3 begins.
Provenance / scope: see `GOVERNANCE_ADOPTION.md` (one-time migration record; archivable
once `pre_auth_check.sh` passes on `main`).

---

## DEC-002 — GWAS Catalog snapshot recovered and release ID documented (2026-07-08)

Decision: The frozen GWAS Catalog input (`data/raw/gwas_catalog.tsv`, SHA-256
`bfaffa3a…`) was recovered from Trash on 2026-07-07 and verified byte-identical to the
pre-registered pin. Its containing directory identified the release as
`v1.0.2-associations_e115_r2026-01-19_full` (GWAS Catalog v1.0.2, Ensembl e115, release
2026-01-19).
Rationale: The original lock pinned the SHA-256 but left the human-readable release ID as
the placeholder `REPLACE_WITH_RELEASE_ID`; the recovered directory name resolves it.
Nature: provenance restoration of metadata — the analytic input is unchanged (bytes match
the pin). **Pending (human re-freeze, NOT done here):** fill `snapshot` in `config.yaml`
and `prereg/frozen_params.yaml`, add `inputs.gwas_loci.snapshot` to FROZEN_FIELDS, and
regenerate `frozen_params.sha256`. See `validate_config_FROZEN_FIELDS_patch.md`.

---

## DEC-003 — v1.3 L2G-audit contamination disposed (2026-07-08)

Decision: Removed the two-generations-obsolete v1.3 L2G-Instability-Audit governance/
protocol files that had been injected into this repo in commit `90f6631`
(`PROTOCOL_v1.3.md`, `LOCK_SUMMARY.md`, `TASKS_WEEK1.md`, `EXECUTION_MEMO_v1.0.md`,
`OUTPUT_SCHEMA.md`, `DATA_SOURCES.md`, `MANUSCRIPT_POSITIONING.md`, `DECISIONS.md`,
`RUN_LOG.md`, `configs/`, `results/locus_instability_atlas/`).
Rationale: Disposed, not relocated — the sibling L2G audit is a superseded design (v1.6,
v29→v45 RA-only) with no valid home for v1.3; the files carry no unique value. The two
that appeared to hold annotation-drift content (`DATA_SOURCES.md`,
`MANUSCRIPT_POSITIONING.md`) on inspection held a stale, unpinned, partially-wrong (v37,
~8.3%) view of this study that `config.yaml`/`frozen_params.yaml` already supersede
(v29/v45, 4.15%); extracting would have seeded a Trigger-11 inconsistency. Nothing
extracted; `DATA_SOURCES.md`, `DECISIONS.md`, `RUN_LOG.md` re-authored fresh from
authoritative sources.
Safety: audit content preserved in the L2G audit archive (`archived-2026-07-08` bundle,
3 copies), in this repo's git history, and in
`~/Downloads/anndrift_v1.3_contamination_removed/`. See `LOCKED_PROTOCOL_FILES.md`
"Foreign artifacts (disposition)" and `GOVERNANCE_ADOPTION.md`.
