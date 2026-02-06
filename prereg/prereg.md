# Preregistration: Annotation Drift in GWAS Locus-to-Gene Mapping

**Study Title:** Annotation Version Sensitivity of GWAS Locus-to-Gene Mapping and Pathway Enrichment
**Registration Date:** 2026-01-24
**Author:** Michael Blumhardt
**Status:** Protocol locked; full analysis pending
**Protocol Version:** 1.0

---

## 0. Protocol of Record

`config.yaml` (project root) is the authoritative, machine-readable protocol of record.
`prereg/frozen_params.yaml` is a human-readable mirror that must match `config.yaml`
at all frozen fields. Any divergence constitutes a protocol amendment and must be
documented in Section 7 below before the pipeline is re-run.

`validate_config.py` (project root) enforces this at runtime: every `make` target
calls it before executing, and the pipeline halts on any violation.

---

## 1. Hypotheses

### Primary Hypothesis (H1)
For common-variant GWAS loci, a non-trivial fraction of locus-to-gene assignments
change when using different GENCODE annotation versions (v29 vs v45, both GRCh38),
and this instability may propagate to downstream pathway enrichment results.

### Quantitative Predictions
- **Instability rate (I):** We predict I > 5% under the primary mapping rule
- **Pathway impact:** We predict Jaccard overlap of significant pathways < 0.90

### Null Hypothesis (H0)
Annotation version choice has negligible effect: I < 5% AND pathway Jaccard > 0.90

---

## 2. Methods

### 2.1 Data Sources (Frozen)

| Source | Version | Build | Rationale |
|--------|---------|-------|-----------|
| GWAS Catalog | Current snapshot | GRCh38 | Standard source |
| GENCODE | v29 (July 2018) | GRCh38.p12 | First mature GRCh38 release |
| GENCODE | v45 (Oct 2023) | GRCh38.p14 | Current release |

**Critical design choice:** Both annotation versions use GRCh38 to isolate
annotation drift from liftover/coordinate artifacts.

### 2.2 Inclusion Criteria (Frozen)
- European ancestry
- Genome-wide significant (p < 5x10^-8)
- Autosomal only (chr1-22)
- Sentinel SNP per locus; duplicates collapsed by locus_id

### 2.3 Mapping Rule (Frozen)
**Primary:** Nearest TSS within symmetric window
- TSS definition: strand-aware (+ strand = min(start), - strand = max(end))
- All transcripts considered; gene-level TSS = nearest transcript TSS to locus
- Ties broken by lexicographic gene_id (smallest Ensembl stable ID)
- Gene filter: protein_coding only
- Gene ID: version-stripped Ensembl stable ID (suffix after '.' removed)
- Window grid: 10, 25, 50, 100, 250, 500 kb (fixed; no post-hoc additions)

### 2.4 Instability Metric (Primary Endpoint)

I = (# loci whose assigned gene_id differs) / (# loci mapped in both versions)

Confidence intervals: Wilson score 95% CI (primary); block bootstrap stratified
by chromosome (1,000 replicates, seed = 42) as supplementary.

### 2.5 Pathway Impact (Co-Primary Endpoint)
- Over-representation analysis (Fisher's exact test) on Reactome gene sets
- Gene-set universe: union of all protein_coding genes in v29 or v45
- Pathway size filter: 15-500 genes
- Benjamini-Hochberg FDR < 0.05
- Cross-version comparison: unweighted Jaccard over significant pathways

### 2.6 Sensitivity Analyses (Pre-Committed)

| Analysis | Label | Specification |
|----------|-------|---------------|
| Window grid | Primary | Full grid (10-500 kb); kill if I range > 20% |
| Biotype expansion | Exploratory | protein_coding + lncRNA + pseudogene |
| Mapping rule | Exploratory | Overlap-first (gene body overlap before nearest-TSS fallback) |

No other sensitivity analyses will be run without a protocol amendment.

---

## 3. Stop Rules (Hard)

| Condition | Criterion | Action |
|-----------|-----------|--------|
| Negative result | I < 5% at all W AND Jaccard > 0.90 | Publish as negative / validation result |
| Instability excess | I range > 20% across window grid | Reframe scope or terminate |
| Pipeline failure | < 90% of loci mapped in either version | Debug before proceeding |

---

## 4. Frozen Fields (Enforced by Validator)

The following config.yaml fields are locked by `validate_config.py` and will
cause a hard pipeline failure if modified without a formal amendment:

1. `map.mapping_rule.type` — Mapping rule
2. `map.windows_kb` — Window grid
3. `map.window_shape` — Window shape
4. `map.tie_break.rule` — Tie-break rule
5. `analysis.primary_instability.denominator` — Instability denominator
6. `enrichment.gene_universe` — Gene-set universe
7. `enrichment.pathway_db.sha256` — Pathway DB SHA-256
8. `analysis.uncertainty.bootstrap.seed` — Bootstrap RNG seed
9. `annotations.gencode` — GENCODE annotation definitions (tags, builds, files, hashes)
10. `inputs.gwas_loci.raw_file` — GWAS input file
11. `inputs.gwas_loci.sha256` — GWAS input SHA-256

---

## 5. Deliverables

1. Quantitative bound: I estimate with CI for GWAS locus-to-gene mapping
2. Pathway propagation test: does gene-list drift translate to enrichment drift?
3. Reproducible pipeline with pinned versions and config validation

---

## 6. Pre-Registered Figures and Tables

| Output | Description |
|--------|-------------|
| I vs W curve | Primary instability across the window grid |
| Delta-distance histogram | Distribution of distance-to-gene changes across versions |
| Enrichment Jaccard by W | Pathway overlap summary across window sizes |
| Instability by chromosome | Per-chromosome breakdown |
| Sensitivity table | I under alternative biotype/mapping-rule specifications |

---

## 7. Amendments

Any deviations from this protocol will be documented here with justification
before the affected analysis step is re-run.

| Date | Amendment | Justification |
|------|-----------|---------------|
| 2026-01-24 | v19 -> v29 | v19 is GRCh37; v29 and v45 are both GRCh38, isolating annotation drift from liftover artifacts |
| 2026-02-06 | Protocol hardening | Reconciled config.yaml with Locked Analyses doc; added validate_config.py with 11 frozen fields; designated config.yaml as protocol of record |
