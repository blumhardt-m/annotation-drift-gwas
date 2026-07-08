# GOVERNANCE_ADOPTION.md

> **One-time migration record — NOT part of the standing governance package.**
> Documents the single migration in which `annotation-drift-gwas` adopted the agent-
> governance layer (AGENTS.md, EXECUTION_BOUNDARY.md, LOCKED_PROTOCOL_FILES.md,
> ONBOARDING_SOP.md, REVIEW_TRIGGER.md, REVIEW_CHECKLIST.md, pre_auth_check.sh + the
> validate_config.py patch) and disposed of foreign v1.3 L2G-audit contamination.
> A historical artifact, like a changelog entry — not enforced or recurring.
>
> **End-of-life:** archivable once `pre_auth_check.sh` passes on `main` post-adoption.

Date: 2026-07-08 · Branch: governance-adoption (off main) · Authority: Michael Blumhardt

---

## Corrected adoption sequence
1. Create branch `governance-adoption` off `main` — DONE.
2. §1b frozen-algorithm filename check (src/04, 05, 12, 16, primary_reactome_export,
   scripts/study2_v37_v45/export_primary_reactome) — DONE, all six resolve.
3. Copy governance files into repo root (7 files + ONBOARDING_SOP.md); chmod +x pre_auth_check.sh.
4. Dispose foreign v1.3 artifacts per LOCKED_PROTOCOL_FILES.md "Foreign artifacts (disposition)":
   - Class 1 (dispose outright): PROTOCOL_v1.3.md, LOCK_SUMMARY.md, TASKS_WEEK1.md,
     EXECUTION_MEMO_v1.0.md, OUTPUT_SCHEMA.md, configs/, results/locus_instability_atlas/
     — remove or quarantine to ~/Downloads/anndrift_v1.3_contamination_removed/.
   - Class 2 (recreate fresh): DECISIONS.md, RUN_LOG.md.
   - Class 3 (extract then discard): DATA_SOURCES.md (keep annotation-drift inputs),
     MANUSCRIPT_POSITIONING.md (keep paper-1 positioning if valuable).
5. Reconcile LOCKED_PROTOCOL_FILES.md §1c to the annotation-drift protocol of record
   (config.yaml + prereg/prereg.md + prereg/frozen_params.{yaml,sha256}); drop the
   PROTOCOL_v1.3.md / LOCK_SUMMARY.md references (disposed).
6. Apply both validate_config.py patches (Change 1: add inputs.gwas_loci.snapshot to
   FROZEN_FIELDS; Change 2: hard-fail on frozen_params SHA mismatch).
7. Add fresh DECISIONS.md entry recording (a) governance adoption and (b) v1.3 disposal (2026-07-08).
8. Commit on governance-adoption; merge to main (baseline infrastructure, not Study 3 work).

## Supersedes
This supersedes the earlier "verify-then-remove against the sibling repo" instruction,
which assumed the v1.3 files had a valid home in
`~/Study/research/08_GENOMICS/l2g_instability_audit/`. They do not: that sibling runs a
superseded design — the live L2G audit is protocol v1.6 (v29→v45, RA-only pilot,
kill-criterion gate, Phase 1 already executed), having abandoned the v1.3 design (v38/v45,
multi-trait, perturbation matrix, LD panels, MSigDB, spectral AUC). Relocating v1.3 there
would inject a conflicting obsolete protocol beside a frozen v1.2 brief. v1.6 is preserved
in ≥4 places, so the v1.3 snapshot carries no unique value: disposal, not relocation.

## Out of scope (higher-priority, separate track)
Recovery of the live L2G audit's own governance state is NOT part of this migration:
recover PROTOCOL_v1.6.md from Trash/Downloads-archive into the live l2g_instability_audit/
tree; place that tree under version control (currently untracked in the zmet-nid monorepo);
verify the Insync backup's git remote pushes somewhere durable. That frozen, executed study
— protocol in Trash, only tracked copy an unverified backup — is a more acute integrity
risk than this cleanup, but a separate study and a separate decision.
