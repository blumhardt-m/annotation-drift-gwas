# REVIEW_TRIGGER.md

Conditions that require the agent to HALT and escalate to the human. On any trigger, the agent stops, reports the condition, and does not proceed. Proceeding through a triggered state is a governance violation.

## Integrity triggers (each of these actually occurred or nearly occurred in this project's history)

1. **SHA mismatch on any locked artifact.** Any file in `LOCKED_PROTOCOL_FILES.md` whose current bytes do not match its recorded SHA-256. Includes `frozen_params.yaml` vs `frozen_params.sha256` — which `validate_config.py` currently only *warns* on; `pre_auth_check.sh` promotes it to a hard stop.
2. **Placeholder in a frozen field.** Any frozen config/prereg value still holding an unresolved placeholder (e.g. `REPLACE_WITH_RELEASE_ID`). A locked protocol containing a placeholder is not validly frozen.
3. **Missing locked input.** Any hash-locked data file absent from disk where config declares it present (the parent study's `MEMORY.md` asserted the Catalog was present when it was in Trash).
4. **Repo-identity mismatch.** Working tree's `git remote` / toplevel does not match the expected `annotation-drift-gwas` repo (this session began in a scaffold whose remote was `zmet-nid`). `pre_auth_check.sh` verifies identity before any run.
5. **Hash mismatch on restoration.** A file being restored to a hash-locked path whose SHA-256 does not equal the locked value. Restoration of matching bytes is permitted; non-matching bytes escalate.
6. **`validate_config.py` non-zero exit.** Any frozen-field divergence. The agent must not attempt the re-freeze remediation; it escalates.

## Scientific / statistical triggers (Study 3)

7. **Degenerate GSEA result (prereg §4.1 sub-case ii).** GSEA produces a trivially small significant set (< 10 significant pathways under both annotation versions), so no meaningful boundary-flip statistic can be computed. Per the prereg this is *not-evaluable, not support* — the agent halts and reports; it must not score this as corroboration of H1.
8. **ORA consistency-check failure.** Regenerated ORA ρ or Jaccard diverges from the parent-study expectations (ρ ≈ 0.923, Jaccard ≈ 0.613) beyond floating-point tolerance. "Within tolerance" means agreement to at least 3 significant figures (e.g. 0.922998 vs 0.923 is a *pass*; 0.918 is a *failure*); trivial last-digit differences from platform floating-point are not triggers. A genuine divergence indicates the frozen machinery or inputs are not reproducing the published result; halt before the contrast.
9. **Distribution outside preregistered sanity limits.** Count-vector diagnostics violate the preregistered diagnostic thresholds (prereg §6 sanity limits): fraction of genes with changed count < 0.5% or > 25% (the 4.15% locus-reassignment rate should induce a count change in a bounded minority of genes, not near-zero and not most of them); nonzero-head size deviating from the verified 17,256 (v29) / 17,296 (v45) by more than ±5%; or zero-tail fraction implying a universe size inconsistent with `union_protein_coding_v29_v45`. These are fixed numeric limits, not a judgment of "plausibility." Any violation halts for human sanity-check.
10. **Universe or size-filter mismatch.** The pathway count after the 15/500 filter differs between the GSEA arm and the regenerated ORA arm, indicating the arms are not testing the same pathway set.

11. **Governance-inventory mismatch.** The set of protocol artifacts referenced by `LOCKED_PROTOCOL_FILES.md`, `validate_config.py:FROZEN_FIELDS`, and `TASK.md` disagree with one another — e.g. a file locked in one but absent from another's list, or a frozen field enforced by the guard but not enumerated in the locked inventory. This is silent governance drift: the layers of the runtime no longer describe the same freeze. `pre_auth_check.sh` cross-checks these inventories and halts on disagreement. Reconciliation is a human act with a `DECISIONS.md` entry.

## Escalation action

On any trigger: stop, write the condition and evidence to the session log, surface to the human, and await instruction. Do not remediate frozen state autonomously. See `REVIEW_CHECKLIST.md` for what the human verifies per trigger.
