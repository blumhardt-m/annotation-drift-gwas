# Execution Memo v1.0

Project: L2G Instability Audit Under Annotation and Evidence-Source Drift
Protocol: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Memo date: 2026-03-11
Operator: Michael Blumhardt

---

## Week-1 Objectives

1. Confirm trait selections and verify that GWAS summary statistics and L2G scores are available for all three trait categories.
2. Confirm GENCODE v38 and v45 GTF files are present and checksummed.
3. Verify LD reference panel (1000 Genomes Phase 3 EUR) is accessible and locked to a fixed version.
4. Establish the baseline pipeline and produce a single baseline L2G run for one trait.
5. Configure the full perturbation matrix (4 runs minimum per trait) and confirm all input paths resolve.
6. Produce the first locus_instability_atlas.tsv stub with correct schema and at least one trait's baseline row populated.
7. Initialize RUN_LOG.md with at least one entry.

---

## In Scope

- Dataset verification and path confirmation for all required inputs
- Baseline pipeline setup (annotation v38, distance-only evidence, nearest-gene distance rule)
- Perturbation matrix configuration (configs/perturbation_matrix.yaml)
- Locus instability atlas generation (Phase A: Existence)
- Top-1 change rate computation across perturbation pairs
- Run logging per RUN_LOG.md schema
- Reactome pathway evaluation against baseline

---

## Out of Scope (Week 1)

- Phase B predictive modeling (AUC computation, feature engineering)
- Phase C spectral analysis
- MSigDB C2 sensitivity analysis
- Cross-trait replication assessment (requires all three traits to be complete)
- Manuscript writing
- Any modification to the locked protocol

---

## Deliverables

| Deliverable | Format | Location |
|-------------|--------|----------|
| Verified input manifest | text checklist | DATA_SOURCES.md (status column updated) |
| Baseline run (1 trait) | pipeline output + log entry | results/locus_instability_atlas/ + RUN_LOG.md |
| Perturbation matrix config | YAML | configs/perturbation_matrix.yaml |
| Atlas stub (1 trait) | TSV | results/locus_instability_atlas/locus_instability_atlas.tsv |
| Top-1 change rate table (1 trait, baseline vs annotation perturb) | TSV | results/locus_instability_atlas/ |
| Updated RUN_LOG.md | Markdown | RUN_LOG.md |

---

## Locked Decision Rules Active This Week

The following protocol rules govern execution this week and are not negotiable:

- LD reference must not vary across any run. Lock it before the first pipeline execution.
- Reactome is the primary pathway arbiter. Do not substitute GO or MSigDB for primary evaluation.
- Reactome, MSigDB, and STRING are forbidden from influencing L2G assignment.
- The baseline-shift design is: baseline → +annotation → +evidence → +distance. Do not reorder.
- Minimum 4 runs per trait before any instability metric is computed.
- Any protocol deviation must be recorded in RUN_LOG.md under "protocol deviations" and in DECISIONS.md.
