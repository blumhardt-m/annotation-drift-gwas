**Publication status — October 30, 2026:** Unpublished manuscript. No preprint has been posted.

# Annotation Drift in GWAS Pathway Interpretation

Reproducible analysis pipeline and manuscript for the study:

**Annotation version as a hidden degree of freedom in GWAS pathway interpretation: threshold amplification of modest locus-level drift**

## Status

Study 3 has been superseded without execution (DEC-013). Corrections to the parent study are recorded in DEC-012. The successor per-trait annotation-version pilot remains a nonmember candidate pending a portfolio membership decision; no pilot execution is authorized. The manuscript is unpublished; no preprint has been posted.

## Summary

Genome annotation resources such as GENCODE are updated regularly, but the choice of annotation version is rarely reported in GWAS pathway analyses.

This project quantifies how annotation evolution propagates through the standard GWAS interpretation pipeline:

GWAS loci → locus-to-gene mapping → pathway enrichment.

Using a controlled GRCh38-on-GRCh38 comparison between **GENCODE v29 (2018)** and **GENCODE v45 (2023)**, we measure how annotation drift affects downstream biological interpretation.

## Key Results

Dataset  
721,406 GWAS loci from the NHGRI-EBI GWAS Catalog.

Primary findings:

| Metric | Result |
|------|------|
| Locus-to-gene instability | **4.15%** |
| Spearman correlation of pathway p-values | **ρ = 0.923** |
| Jaccard overlap of significant pathways | **0.61** |
| Significant pathways | v29: 74 • v45: 47 |
| Shared significant pathways | **46** |

Interpretation:

Annotation updates perturb locus-to-gene assignments modestly but can produce larger differences in **thresholded pathway discoveries** due to **threshold amplification**, while the global evidence structure remains stable.

## Repository Structure
