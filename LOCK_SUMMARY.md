# Lock Summary

Protocol ID: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Lock Date: 2026-03-11
Locked by: Michael Blumhardt (Research Director)

---

## Locked Components

| Component | Locked Value |
|-----------|-------------|
| Protocol version | v1.3 |
| Trait panel | Autoimmune + Psychiatric + Metabolic (one each) |
| Trait replication rule | ≥2 of 3 traits; instability heterogeneity + pathway impact both required |
| Annotation versions | GENCODE v38 and GENCODE v45 |
| Minimum pipeline runs per trait | 4 |
| LD reference | 1000 Genomes Phase 3 EUR (constant across all perturbations) |
| Primary pathway arbiter | Reactome |
| Sensitivity pathway database | MSigDB C2 |
| Circularity constraint | Reactome, MSigDB, STRING forbidden from L2G assignment |
| Instability metrics | top-1 change rate, top-k instability, rank entropy |
| Material impact threshold (primary) | Top 10% unstable loci → ≥30% Jaccard degradation |
| Material impact threshold (secondary) | Jaccard similarity < 0.65 |
| Instability heterogeneity rule | mean(Q4) ≥ 3 × mean(Q1) |
| Phase B baseline predictors | gene density, LD entropy, transcript overlap complexity, nearest-gene distance asymmetry, regulatory link multiplicity, regulatory evidence entropy |
| Phase B AUC thresholds | Predictive: >0.60; Descriptive: ≤0.60; Phase C gate: <0.75 |
| Spectral retention threshold | ΔAUC ≥ 0.05 |
| Project termination triggers | Instability pathway-inconsequential OR uniform stability |
| Baseline-shift perturbation design | baseline → +annotation → +evidence → +distance |

---

## What Remains Open

The following have not been locked and must be confirmed before execution:

- Final trait selection within each category (RA vs IBD; SCZ vs MDD; BMI vs T2D) — must verify GWAS summary stat availability and L2G score availability
- Evidence-source file paths (eQTL, ABC scores) — require dataset verification
- Exact GENCODE v38 GTF source URL and checksum
- MSigDB C2 version to use for sensitivity analysis
- Number of top-k assignments for top-k instability metric
- Phase C spectral feature set (intentionally deferred)

---

## Rule for Scientific Changes

Any modification to scientific content, decision thresholds, outcome definitions, trait selection criteria, or perturbation design requires issuance of a new protocol version (v1.4 or higher).

Changes to operational details (file paths, execution order, logging format) do not require a new protocol version but must be recorded in DECISIONS.md with a rationale.

The locked protocol document (PROTOCOL_v1.3.md) must never be edited after lock date.
