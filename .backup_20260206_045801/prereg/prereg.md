# Preregistration: Annotation Drift in GWAS Locus-to-Gene Mapping

**Registration Date:** 2026-01-24  
**Author:** Michael Blumhardt  
**Status:** Pre-analysis

---

## 1. Hypotheses

**H1:** A non-trivial fraction (>5%) of GWAS locus→gene assignments change between GENCODE v29 and v45 (both GRCh38).

**H0:** I < 5% AND pathway Jaccard > 0.90 (negligible drift)

## 2. Data Sources (Frozen)

| Source | Version | Build |
|--------|---------|-------|
| GWAS Catalog | 2026-01 snapshot | GRCh38 |
| GENCODE | v29 (July 2018) | GRCh38 |
| GENCODE | v45 (Oct 2023) | GRCh38 |

## 3. Mapping Rule (Frozen)

- **Primary:** Nearest TSS within ±500 kb
- **TSS:** Most 5' position across transcripts per gene
- **Ties:** Smallest Ensembl stable ID
- **Filter:** Protein-coding only

## 4. Metrics

**Instability rate:** I = (# genes differ) / (# valid in both)

**Pathway impact:** Jaccard overlap of significant pathways (FDR < 0.05)

## 5. Stop Rules

| Condition | Action |
|-----------|--------|
| I < 5% AND Jaccard > 0.90 | Negative result |
| I range > 20% across windows | Reframe or kill |

## 6. Sensitivity Analyses (Bounded)

1. Window: ±250 kb, ±500 kb, ±1 Mb
2. Biotype: protein-coding vs +lncRNA
3. Versions: v44 vs v45 (adjacent)

## 7. Amendments

| Date | Change | Rationale |
|------|--------|-----------|
| 2026-01-24 | v19→v29 | Both v29/v45 are GRCh38 |
