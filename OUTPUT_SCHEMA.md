# Output Schema

Protocol: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
File: results/locus_instability_atlas/locus_instability_atlas.tsv
Last updated: 2026-03-11

---

## Composite Key

Uniqueness is enforced on the combination of: `gwas_trait_id` + `locus_id`

Each row represents one locus from one trait. Instability metrics are aggregated across all perturbation runs for that trait.

---

## Schema

| Field Name | Type | Description |
|------------|------|-------------|
| gwas_trait_id | string | Trait identifier (e.g., Open Targets study ID or EFO code). Part of composite key. |
| locus_id | string | Locus identifier. Format: CHR_LEAD_VARIANT (e.g., 6_32575460_A_G) or equivalent stable identifier. Part of composite key. |
| lead_variant_rsid | string | rsID of the lead variant for this locus. May be null if rsID not available. |
| chrom | string | Chromosome (1–22, X). |
| pos_hg38 | integer | GRCh38 position of the lead variant. |
| top1_gene_baseline | string | HGNC symbol of the top-1 L2G gene under the baseline pipeline configuration. |
| top1_change_rate | float | Fraction of perturbation runs (excluding baseline) in which the top-1 gene assignment differs from baseline. Range [0, 1]. |
| top1_change_count | integer | Number of perturbation runs in which the top-1 gene differs from baseline. |
| n_perturbation_runs | integer | Total number of perturbation runs used to compute instability metrics (excluding baseline). Minimum 3 per protocol. |
| top_k | integer | Value of k used for top-k instability metric. |
| topk_instability | float | Fraction of top-k gene sets that differ from baseline top-k set, averaged across perturbation runs. Range [0, 1]. |
| rank_entropy | float | Shannon entropy of the rank distribution of the top-1 gene across perturbation runs. Higher values indicate greater instability. |
| instability_quantile | integer | Quartile bin of top1_change_rate for this trait (1=most stable, 4=most unstable). Computed within-trait. |
| annotation_perturbation_top1_change | integer | 0/1 flag: top-1 gene differs between baseline and annotation-perturbed run. |
| evidence_perturbation_top1_change | integer | 0/1 flag: top-1 gene differs between baseline and evidence-perturbed run. |
| distance_perturbation_top1_change | integer | 0/1 flag: top-1 gene differs between baseline and distance-perturbed run. |
| reactome_jaccard_baseline_vs_annotation | float | Jaccard similarity of Reactome pathway enrichment sets between baseline and annotation-perturbed run for this trait (trait-level, repeated per locus). |
| reactome_jaccard_baseline_vs_evidence | float | Jaccard similarity of Reactome enrichment sets between baseline and evidence-perturbed run. |
| reactome_jaccard_baseline_vs_distance | float | Jaccard similarity of Reactome enrichment sets between baseline and distance-perturbed run. |
| gene_density_1mb | float | Number of protein-coding genes within 1 Mb of the lead variant (GRCh38). Baseline predictor for Phase B. |
| ld_entropy | float | Shannon entropy of the LD r² distribution within the locus credible set. Baseline predictor for Phase B. |
| transcript_overlap_complexity | float | Count of distinct transcript isoforms overlapping the locus window, normalized by window size. Baseline predictor for Phase B. |
| nearest_gene_distance_asymmetry | float | Absolute difference in distance between the nearest upstream and nearest downstream gene. Baseline predictor for Phase B. |
| regulatory_link_multiplicity | integer | Number of distinct regulatory links (eQTL, chromatin, VEP) connecting lead variant to any gene. Baseline predictor for Phase B. |
| regulatory_evidence_entropy | float | Shannon entropy of regulatory evidence source weights at this locus. Baseline predictor for Phase B. |
| phase_b_predicted_unstable | integer | 0/1 flag populated in Phase B: model predicts this locus as unstable. Null until Phase B is run. |
| notes | string | Free-text notes. Use for anomalies, data quality flags, or manual review flags. |

---

## File Format

- Tab-separated values (.tsv)
- Header row required
- UTF-8 encoding
- No trailing whitespace
- Null values represented as empty string (not "NA", "null", or ".")
- One row per (gwas_trait_id, locus_id) pair
- Rows sorted by gwas_trait_id ascending, then chrom ascending, then pos_hg38 ascending
