# REVIEW_CHECKLIST.md

What the human verifies when a `REVIEW_TRIGGER.md` condition fires. Standardizes intervention so review is not ad hoc.

## On any integrity trigger (1–6)
- [ ] Identify which locked artifact / field is implicated and its expected vs actual value.
- [ ] Confirm whether this is provenance restoration (matching bytes) or an input/parameter change (differing bytes/values). Only the former is ever permitted, and only by the documented restore path.
- [ ] If a placeholder or missing `FROZEN_FIELDS` entry is the cause: decide whether to amend (human act) with a `DECISIONS.md` entry, or revert. Never let the agent re-freeze.
- [ ] For SHA mismatch on `frozen_params.yaml`: determine what modified the locked file post-freeze; restore from `v1.0_submission` or re-lock deliberately with a decision record.
- [ ] For repo-identity mismatch: confirm the working tree is the correct clone before anything else.

## On degenerate GSEA (trigger 7)
- [ ] Confirm the significant-set counts under both versions and that sub-case (ii), < 10 under both, is genuinely met (not a bug in stage 20).
- [ ] Verify the count vector and fgsea inputs are well-formed (tie structure, non-empty nonzero head).
- [ ] Record the result as not-evaluable per prereg §4.1; do NOT allow it to be scored as H1 support.
- [ ] Decide whether the arm is reported as inconclusive or whether a prereg-level reconsideration is warranted (a new registration, not an edit).

## On ORA consistency failure (trigger 8)
- [ ] **First, check the gates.** Review `pre_auth_check.sh` and `validate_config.py` output from this session — if either failed or was bypassed earlier, resolve that before investigating scientific causes. A guard failure upstream can masquerade as a consistency failure downstream.
- [ ] Compare regenerated ρ/Jaccard to 0.923 / 0.613; quantify divergence (fail = disagreement beyond 3 significant figures).
- [ ] Check inputs: is `comparison.tsv` the frozen file (hash), are TSS tables regenerated from the correct GENCODE hashes, is the Reactome GMT v95 by hash.
- [ ] Do not proceed to the contrast until the parent result reproduces.

## On unexpected distribution / universe mismatch (triggers 9–10)
- [ ] Review count-vector diagnostics against the 4.15% reassignment expectation.
- [ ] Confirm the GSEA and ORA arms share the same union universe and identical 15/500-filtered pathway set.
- [ ] Sanity-check hub-gene tail (upper-tail quantiles) for domination before interpreting any ranking result.

## Sign-off
- [ ] Resolution recorded in `DECISIONS.md` (for any frozen-state change) and/or `RUN_LOG.md` (for run-level notes).
- [ ] Only after sign-off does the agent resume, under a possibly-revised `TASK.md`.
