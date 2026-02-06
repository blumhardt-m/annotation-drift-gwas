# Annotation Drift in GWAS Locus-to-Gene Mapping

**Status:** Pilot complete, proceeding to full analysis  
**Started:** 2026-01-24  
**Target:** Bioinformatics or NAR Genomics & Bioinformatics

## One-Line Summary

Quantifying how often GWAS locus→gene assignments change between GENCODE versions, and whether this propagates to pathway enrichment differences.

## The Question

> If the same GWAS loci are processed through two annotation versions under identical mapping rules, how often does the assigned gene change?

## Key Design

**GRCh38-on-GRCh38** (v29 vs v45) to isolate annotation drift from liftover artifacts.

## Pilot Results (Synthetic)

| Metric | Value |
|--------|-------|
| Instability rate | 8.34% |
| 95% CI | [6.63%, 10.05%] |
| Window stability | 0.4% range |

**Verdict:** Non-trivial effect. Proceed.

## Quick Start

```bash
# 1. Setup environment
python3 -m venv .venv && source .venv/bin/activate
pip install numpy pyyaml

# 2. Download data (~500 MB)
./src/01_download_data.sh

# 3. Parse to TSS tables
make parse

# 4. Run analysis
make map analyze
```

## Stop Rules

| Condition | Action |
|-----------|--------|
| I < 5% AND Jaccard > 0.90 | Publish as negative |
| I range > 20% across windows | Reframe or kill |
