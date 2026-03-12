# Decision Log

This file is append-only. Do not edit or delete existing entries.
New entries are appended at the bottom.

---

## Entry Format

```
Date:           YYYY-MM-DD
Decision ID:    DEC-NNN
Decision:       [one sentence]
Rationale:      [one to three sentences]
Affected files: [comma-separated list]
```

---

## DEC-001 — Protocol Lock

Date:           2026-03-11
Decision ID:    DEC-001
Decision:       Protocol L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3 is declared locked as of 2026-03-11 and must not be modified.
Rationale:      Locking the protocol before execution prevents post-hoc adjustment of scientific criteria, preserves the integrity of the instability audit as a pre-registered study, and ensures that any future change is traceable through a new protocol version number.
Affected files: PROTOCOL_v1.3.md, LOCK_SUMMARY.md

---

## DEC-002 — Spectral Analysis Demoted to Phase C

Date:           2026-03-11
Decision ID:    DEC-002
Decision:       Topology-derived spectral predictors (e.g., graph-spectral features of LD structure) are deferred to Phase C and evaluated only if Phase B baseline AUC falls below 0.75.
Rationale:      Spectral features require additional preprocessing infrastructure and carry a risk of overfitting in small locus panels. Baseline structural predictors (gene density, LD entropy, transcript overlap, etc.) are sufficient to determine whether the instability regime is predictive or descriptive. Spectral analysis is computationally expensive and should not be run speculatively.
Affected files: PROTOCOL_v1.3.md, LOCK_SUMMARY.md, configs/perturbation_matrix.yaml

---

## DEC-003 — Baseline-Shift Perturbation Design

Date:           2026-03-11
Decision ID:    DEC-003
Decision:       The perturbation design follows a baseline-shift structure: baseline pipeline, then +annotation perturbation, then +evidence perturbation, then +distance perturbation, with a minimum of 4 runs per trait.
Rationale:      A baseline-shift design isolates the marginal contribution of each perturbation dimension, enabling clean decomposition of instability sources. Running fewer than 4 perturbations per trait would prevent meaningful decomposition and would not meet the minimum evidence threshold for the instability heterogeneity test.
Affected files: PROTOCOL_v1.3.md, configs/perturbation_matrix.yaml, TASKS_WEEK1.md

---

## DEC-004 — Reactome as Primary Pathway Arbiter

Date:           2026-03-11
Decision ID:    DEC-004
Decision:       Reactome is designated the primary pathway arbiter for all downstream impact evaluation; MSigDB C2 is designated sensitivity-only.
Rationale:      Reactome provides manually curated, hierarchically organized pathways with explicit versioning, making it the most reproducible and auditable choice for primary evaluation. MSigDB C2 includes heterogeneous gene set collections of variable curation quality and is appropriate only for sensitivity confirmation. Consistency with the predecessor annotation drift manuscript also motivates this choice.
Affected files: PROTOCOL_v1.3.md, configs/pathway_arbiter.yaml, OUTPUT_SCHEMA.md

---

## DEC-005 — LD Reference Lock

Date:           2026-03-11
Decision ID:    DEC-005
Decision:       The LD reference panel is locked to 1000 Genomes Phase 3 EUR and must not vary across any perturbation run.
Rationale:      Varying the LD reference would confound instability attributable to annotation or evidence-source changes with instability attributable to LD structure differences. Fixing a single reference panel eliminates this confound and ensures that all perturbation comparisons are interpretable as annotation or evidence effects rather than population effects.
Affected files: PROTOCOL_v1.3.md, configs/ld_reference.yaml, LOCK_SUMMARY.md

---

## DEC-006 — Trait Replication Rule

Date:           2026-03-11
Decision ID:    DEC-006
Decision:       The instability pattern must replicate in at least 2 of 3 traits, where replication requires both instability heterogeneity (mean Q4 ≥ 3× mean Q1) and material downstream pathway impact (primary and secondary rules both satisfied); failure triggers project termination.
Rationale:      A single-trait finding would be insufficient to support a general claim about L2G instability across genetic architectures. Requiring replication in two architecturally distinct traits provides cross-architecture generalizability. The dual requirement (heterogeneity + pathway impact) prevents scenarios where instability is present but downstream-inconsequential, which the protocol explicitly designates as a termination trigger.
Affected files: PROTOCOL_v1.3.md, LOCK_SUMMARY.md, TASKS_WEEK1.md
