# Data Sources Registry

Protocol: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Last updated: 2026-03-11

---

## Required — Already In Hand

| Source Name | Purpose | Expected Path | Status | Verification Note |
|-------------|---------|---------------|--------|-------------------|
| GENCODE v37 GTF | Annotation version for predecessor study reference | data/raw/gencode.v37.annotation.gtf.gz | present | Used in annotation drift manuscript; SHA-256 not yet recorded for audit |
| GENCODE v45 GTF | Annotation version for perturbation runs | data/raw/gencode.v45.annotation.gtf.gz | present | Present on disk; SHA-256 not yet recorded |
| Reactome pathways GMT | Primary pathway arbiter for downstream impact evaluation | data/raw/ReactomePathways.gmt | present | File present; version string must be extracted and recorded |
| Reactome pathways TXT | Reactome pathway-gene membership (alternative format) | data/raw/ReactomePathways.txt | present | Present on disk |
| Reactome pathway relations | Hierarchical pathway structure | data/raw/ReactomePathwaysRelation.txt | present | Present on disk |
| GWAS Catalog associations | Trait locus input for all pipeline runs | data/raw/gwas_catalog.tsv | present | Version and download date must be confirmed |
| GWAS Catalog full download | Extended associations for trait filtering | data/raw/gwas-catalog-download-associations-alt-full.tsv | present | File present; download date not recorded |
| Eligible GWAS traits list | Preprocessed eligible trait set | data/processed/gwas_eligible.tsv | present | Output of src/03_parse_gwas.py; regenerable |

---

## Required — External / Public

| Source Name | Purpose | Expected Path | Status | Verification Note |
|-------------|---------|---------------|--------|-------------------|
| GENCODE v38 GTF | Primary annotation version for baseline pipeline | data/raw/gencode.v38.annotation.gtf.gz | TODO_VERIFY | Must be downloaded and checksummed before first run; see https://www.gencodegenes.org/human/release_38.html |
| L2G scores — autoimmune trait | Locus-to-gene priority scores for autoimmune trait (RA or IBD) | data/raw/l2g_autoimmune_TODO_VERIFY.tsv | TODO_VERIFY | Source: Open Targets Genetics; trait ID TBD pending DEC-003 trait selection |
| L2G scores — psychiatric trait | Locus-to-gene priority scores for psychiatric trait (SCZ or MDD) | data/raw/l2g_psychiatric_TODO_VERIFY.tsv | TODO_VERIFY | Source: Open Targets Genetics; trait ID TBD |
| L2G scores — metabolic trait | Locus-to-gene priority scores for metabolic trait (BMI or T2D) | data/raw/l2g_metabolic_TODO_VERIFY.tsv | TODO_VERIFY | Source: Open Targets Genetics; trait ID TBD |
| GWAS summary stats — autoimmune trait | Full summary statistics for locus definition | data/raw/gwas_sumstats_autoimmune_TODO_VERIFY.tsv | TODO_VERIFY | Source depends on trait selection; must use GRCh38 coordinates |
| GWAS summary stats — psychiatric trait | Full summary statistics for locus definition | data/raw/gwas_sumstats_psychiatric_TODO_VERIFY.tsv | TODO_VERIFY | Source depends on trait selection |
| GWAS summary stats — metabolic trait | Full summary statistics for locus definition | data/raw/gwas_sumstats_metabolic_TODO_VERIFY.tsv | TODO_VERIFY | Source depends on trait selection |
| 1000 Genomes Phase 3 EUR LD panel | LD reference panel; must be constant across all runs | data/raw/ld_reference/1kg_phase3_eur_TODO_VERIFY/ | TODO_VERIFY | Version must be fixed before first run per DEC-005; plink .bed/.bim/.fam format expected |
| GTEx eQTL (latest) | Evidence source for eQTL-only and eQTL+ABC perturbation configs | data/raw/gtex_eqtl_TODO_VERIFY.tsv | TODO_VERIFY | Version and tissue scope must be confirmed; GTEx v8 or v10 |
| ABC scores | Evidence source for eQTL+ABC perturbation config | data/raw/abc_scores_TODO_VERIFY.tsv | TODO_VERIFY | Fulco et al. or Nasser et al. dataset; cell type scope TBD |
| MSigDB C2 GMT | Sensitivity pathway database | data/raw/msigdb_c2_TODO_VERIFY.gmt | TODO_VERIFY | Version must match between baseline and all sensitivity runs; download from MSigDB portal |

---

## Optional — Sensitivity Only

| Source Name | Purpose | Expected Path | Status | Verification Note |
|-------------|---------|---------------|--------|-------------------|
| MSigDB Hallmarks GMT | Additional sensitivity check (not primary) | data/raw/msigdb_hallmarks_TODO_VERIFY.gmt | optional | Not required for Phase A or B; only if sensitivity analysis is extended |
| Chromatin interaction data (HiC/PLAC-seq) | Full evidence stack perturbation config | data/raw/hic_TODO_VERIFY/ | optional | Required only for full-evidence-stack configuration; cell type and resolution TBD |
| Variant effect predictions (VEP or CADD) | Full evidence stack perturbation config | data/raw/vep_scores_TODO_VERIFY.tsv | optional | Required only for full-evidence-stack config |
| Regulatory annotation tracks (ENCODE) | Full evidence stack perturbation config | data/raw/encode_regulatory_TODO_VERIFY/ | optional | Required only for full-evidence-stack config |
