# Study 3 Preregistration — Ranked-Enrichment Contrast (count-based)

**Parent study:** Annotation version as a hidden degree of freedom in GWAS pathway interpretation (BIOINF-2026-0616; repo `annotation-drift-gwas` @ tag `v1.0_submission`, commit `37272b2`).
**Study 3 ID:** `annotation-drift-gwas / study3_ranked_enrichment`
**Status:** FROZEN v1.0 — frozen 2026-07-08. SHA-256 recorded in `prereg/study3_prereg.sha256`; the mapping input of record is `comparison.tsv` SHA `8f27c1db…` (see §3.1), and `inputs.gwas_loci.snapshot` is pinned and enforced by `validate_config.py`. OSF registration is the sole remaining external step (to be completed by the human before any Study 3 analysis code executes). Post-freeze changes require a documented protocol amendment.
**Registration basis:** All parameters in §3 are fixed before first execution. The single confirmatory endpoint (§4.1) and its directional prediction are declared before observation of any GSEA result. Secondary/descriptive quantities (§4.2) and diagnostics (§6) are reported unconditionally.

---

## 1. Purpose

The parent study established, under a nearest-TSS / Reactome / overrepresentation-analysis (ORA) pipeline on a fixed GRCh38 build, that annotation revision from GENCODE v29 to v45 produces modest locus-level reassignment (4.15%), high preservation of pathway p-value rank structure (Spearman ρ = 0.923), and substantially lower concordance of binary significance classifications (Jaccard = 0.61 at FDR < 0.05). The parent study attributed this dissociation to threshold amplification: a discrete decision rule applied at a continuous statistical boundary converts small upstream perturbations into large differences in reported findings.

Threshold amplification is a claim about the decision rule, not about ORA specifically. It predicts a testable contrast: an enrichment method that integrates a rank-ordered per-gene evidence statistic, rather than thresholding a discrete gene list, should convert the same upstream annotation perturbation into fewer boundary-localized classification flips.

Study 3 tests this prediction by contrasting preranked GSEA (fgsea) against the frozen ORA result, using a trait-agnostic count of genome-wide-significant loci per gene as the evidence statistic. The scientific object is the per-gene **count**; the preranked algorithm consumes the **ranking induced by that count**. This distinction is maintained throughout: the count is the measured quantity, and the rank ordering is a derived input to fgsea, not itself the object of study. The design matches the ORA arm's trait-agnostic evidentiary basis and avoids introducing a trait-dependent p-value aggregation rule.

**Design rationale — why count-based, not p-value-based.** The parent ORA foreground pooled all genome-wide-significant loci trait-agnostically; set membership does not read p-value magnitude, only threshold crossing. Recovered eligibility data (`gwas_eligible.tsv`, GWAS Catalog release `v1.0.2-associations_e115_r2026-01-19_full`, SHA-256 `bfaffa3a…`) shows 31.3% of SNPs carry multiple trait associations with widely varying p-values (e.g. rs1260326: 2,224 associations, 941 distinct p-values across unrelated phenotypes). Any per-locus p-value therefore requires a trait-collapse rule; the only trait-agnostic collapse (min-p across traits) is dominated by pleiotropic, heavily-studied loci — a study-intensity confound. To contrast against the trait-agnostic pooled ORA on the same basis without smuggling in this confound, the gene statistic is defined as the **count of mapped genome-wide-significant loci per gene**, computed directly from the frozen `comparison.tsv`. This keeps Study 3 a strictly additive extension of frozen artifacts and reads no p-value magnitude. This does not eliminate analytic choices — rsID-within-gene deduplication, integer counting, the shared universe, and tie handling are all frozen choices (§3.2) — but it specifically avoids a trait-dependent p-value aggregation rule, which was the dominant and least-defensible degree of freedom the recovered eligibility data exposed.

**Provenance and the viability/scientific-support boundary.** An exploratory viability probe (2026-07-08, quarantined in `scratch/study3_viability/`, not part of the confirmatory run) established two things before this prereg was frozen: (a) the regenerated `comparison.tsv` (SHA-256 `8f27c1db…`, rebuilt from the hash-verified raw GWAS Catalog `bfaffa3a…` through the frozen mapping pipeline) reproduces the parent ORA result within 3 significant figures (ρ = 0.923, Jaccard = 0.613, sig sets 74/47/46), confirming continuity with the published study; and (b) the count-based GSEA arm is computationally non-degenerate (47 significant pathways under v29, 63 under v45 — well above the §4.1 sub-case-(ii) threshold of 10). The probe established *viability* — that the proposed analysis executes on the recovered data and produces analyzable output. It is **not** evidence for H1: the significant-set counts it produced are descriptive and tie-sensitive (§3.2, §4.1), and the threshold-amplification hypothesis is tested only by the confirmatory boundary-flip-localization contrast, which the probe deliberately did not compute. This distinction between demonstrated viability and untested scientific support is maintained throughout.

## 2. Hypothesis

**H1 (boundary-flip localization, confirmatory, sole primary).** Classification flips — pathways significant under exactly one annotation version at FDR < 0.05 — are more sharply localized near the FDR decision boundary under ORA than under count-based preranked GSEA. The boundary-localization statistic C (§4.1) satisfies C_ORA > C_GSEA. This prediction follows because ORA applies a discrete threshold during construction of the foreground gene set, whereas preranked enrichment delays discretization until pathway-level inference; both methods ultimately threshold pathways at BH-adjusted FDR < 0.05, so the contrast is early versus late thresholding, not threshold versus no threshold. The mechanism predicts that annotation drift is amplified into boundary-localized pathway flips specifically where an early, gene-level discretization intervenes between the perturbed evidence and the pathway test.

Support for H1 corroborates threshold amplification as a property of the decision rule. Failure narrows the parent claim to an ORA-specific annotation-sensitivity result (§5).

## 3. Frozen parameters

Enforced by the pipeline configuration-validation layer (`validate_config.py`) at each stage. Inherited constants are imported from existing frozen modules, not restated in code.

### 3.1 Inputs (verified provenance)
- Mapping artifact: `data/processed/mappings/comparison.tsv` — the frozen parent-study locus-to-gene table (721,406 rows; 291,292 unique rsIDs; 29,920 reassignments = 4.15%), SHA-256 `8f27c1dba92c2af3ec4c427b1bf89cdb03d67f68adef09c1c9131111f6ba2bb1` (recorded in `frozen_params.yaml:inputs.mapping_artifact` and `DATA_SOURCES.md`). This is the sole required input for the gene-count statistic.
- GWAS Catalog release: `v1.0.2-associations_e115_r2026-01-19_full`, SHA-256 `bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b`. Restored and hash-verified into the clone `data/raw/` prior to execution; release ID pinned to `inputs.gwas_loci.snapshot` and enforced. Required for provenance completeness and the governance amendment (§7); NOT required for the count statistic itself, which derives entirely from `comparison.tsv`.
- Reactome GMT: v95, SHA-256 `a34f43fd…` per `frozen_params.yaml`.
- Annotations: GENCODE v29 (GRCh38.p12) and v45 (GRCh38.p14), version-locked in `01_download_data.sh`; TSS tables regenerable and provenance-clean.

### 3.2 Per-gene rank statistic (count-based)
- Statistic: for each gene under each annotation version, the **count of distinct genome-wide-significant loci assigned to that gene** under the frozen nearest-TSS mapping (imported from stage `04_map_variants.py`; ±500 kb window, TSS = most-5′ across transcripts, tie-break lexicographically smallest Ensembl ID — all inherited unchanged).
- Locus counting unit: distinct rsID assigned to the gene under that version. Where the mapping table represents a single rsID in multiple rows, the rsID is counted once per gene (deduplicated on rsID within gene, consistent with `merge_loci: false` / sentinel-SNP eligibility in the parent prereg). No trait dimension enters; the count is trait-agnostic.
- Universe (frozen to match ORA arm): the single shared union universe `union_protein_coding_v29_v45` (`config.yaml:enrichment.gene_universe`; the `ensembl_universe` constructed in `primary_reactome_export.py`). The universe is NOT per-version.
- Zero assignment: genes with no assigned genome-wide-significant locus under a version receive count 0 in that version's vector. Each preranked vector comprises a nonzero head (genes with ≥1 assigned significant locus: 17,256 under v29, 17,296 under v45 per verification) and a tied-zero tail (union universe minus head).
- Ties: the preranked vector is integer-valued (locus counts), so ties are frequent and expected. The tie fraction is high — approximately 88% of the union universe shares a value with at least one other gene, dominated by the large tied-zero tail (~10.5% of the universe is zero) and by low-integer plateaus. fgsea evaluates numeric scores; equal counts are genuine ties whose treatment is fixed by the pinned fgsea version, permutation mode, and seed (§3.3). Any secondary ordering (e.g. Ensembl ID) is a list-reproducibility device only and does not determine the enrichment score. Because the tie fraction is this high, the raw per-pathway NES and p-values are tie-order-sensitive; consequently the significant-set *counts* are treated as non-confirmatory (see §4.1), and the confirmatory claim rests on boundary-flip localization, which is computed from the significance *classifications* rather than from tie-sensitive rank magnitudes. Tie structure is reported as a mandatory diagnostic (§6).

### 3.3 Enrichment (GSEA arm)
- Implementation: `fgsea` (preranked). Broad GSEA is not used (no phenotype labels; the score is a locus-derived count).
- Gene sets: the same Reactome GMT (v95) used in the ORA arm, identical symbol→Ensembl harmonization. The harmonization MUST import and use the identical `build_alias_table()` (MARCH/MARCHF, SEPT/SEPTIN renames) from `primary_reactome_export.py` / `src/12_ensembl_id_enrichment.py`; reimplementation is prohibited.
- minSize: 15 (frozen to `config.yaml:enrichment.pathway_size.min` and `MIN_SIZE` in the ORA export script — identical filter to the ORA arm)
- maxSize: 500
- Multiple-testing correction: Benjamini–Hochberg, FDR < 0.05 (identical to ORA).
- Seed: 42 (`config.yaml:repro.determinism.global_seed`), recorded in §6.
- Permutation: `fgseaMultilevel` with `eps = 0` (frozen), identical across annotation versions within a run. Toolchain frozen to the versions established by the 2026-07-08 viability probe: R 4.5.3, fgsea 1.36.2 (Bioconductor). Mode, eps, and versions are recorded as mandatory §6 diagnostics. (Rationale: with the high tie fraction of the count-based ranking, `fgseaMultilevel` is the appropriate estimator; the standard permutation path is not used.)

### 3.4 ORA arm (inherited)
- One-sided Fisher's exact test per Reactome pathway; foreground = trait-agnostic pooled genes mapped from genome-wide-significant loci; background = shared union universe; BH FDR < 0.05. Implementation: frozen ORA machinery in `primary_reactome_export.py` / `src/12_ensembl_id_enrichment.py`.
- ORA reference significant sets (used to compute C_ORA and the descriptive A_ORA) are regenerated by running the frozen `primary_reactome_export.py` against `comparison.tsv` and the TSS tables. Consistency check: regenerated ρ and Jaccard must match the parent-study expectations (ρ ≈ 0.923, Jaccard ≈ 0.613) within floating-point tolerance; divergence halts execution for investigation before the contrast proceeds.

## 4. Endpoints

### 4.1 Primary (confirmatory) — boundary-flip localization
A classification flip is a pathway significant under exactly one annotation version at FDR < 0.05. For each method m ∈ {ORA, GSEA}, the signed distance of each pathway from the FDR boundary is computed on the −log₁₀(FDR) scale using the more-significant of its two version-specific values. The boundary-localization statistic is:

  C_m = (number of flips within ±Δ of the FDR boundary on the −log₁₀(FDR) scale) / (total flips under method m)

with Δ = 0.5 (frozen). The flip *definition* is imported from `16_flip_rate.py` so ORA and GSEA flip sets are constructed identically.

The number of significant pathways per version, and the between-version Jaccard of those sets, are **not** confirmatory quantities under this design. Because the count-based ranking has a high tie fraction (§3.2), those raw counts are tie-order-sensitive. The confirmatory claim rests solely on the *localization* of flips relative to the FDR boundary (C_ORA > C_GSEA) — a comparison of where classification changes occur, which is robust to the tie sensitivity of the underlying magnitudes. Significant-set counts are reported descriptively only.

- Directional prediction: **C_ORA > C_GSEA** (ORA flips concentrate near the cutoff; GSEA flips, if present, are not preferentially boundary-localized).
- Reported alongside C_m: the full flip-count-versus-distance-from-boundary distribution for both methods, so localization is visible independent of Δ. (Artifact name `boundary_flip_localization` — deliberately distinct from the existing spatial `08_concentration_analysis.py`.)
- Uncertainty: flip counts and C_m are bootstrapped chromosome-block, seed 42, 1000 reps (`frozen_params.yaml:analysis.uncertainty.bootstrap`), resampling at the variant level and propagating through mapping → count-statistic → enrichment → flip classification. This matches the frozen prereg convention; the i.i.d. default in `05_bootstrap.py` is not used.
- **Degenerate-case rule (predeclared).** C_GSEA is undefined if GSEA produces zero flips. Two sub-cases, distinguished before observation: (i) GSEA produces a non-trivial significant set (≥ 10 significant pathways under at least one version) but zero or few flips → this is *supportive* of H1 (ranked enrichment resists the thresholded instability ORA exhibits) and is scored as corroboration, with the near-zero flip count reported directly rather than as a ratio C. (ii) GSEA produces a trivially small significant set (< 10 significant pathways under both versions), such that no meaningful flip statistic can be computed → the contrast is *void, not supportive*: the count statistic failed to produce a usable enrichment signal, H1 is scored not-evaluable, and the disposition defers to reporting the count-distribution and enrichment diagnostics without a mechanistic claim. The 10-pathway threshold is frozen. This rule prevents a degenerate GSEA result from being read as spurious support for the mechanism.

### 4.2 Secondary / descriptive — amplification ratio A
For each method: A_m = (1 − J_m) / (1 − ρ_m), where J_m is the between-version Jaccard of significant sets and ρ_m the between-version Spearman correlation of pathway statistics. J_m, ρ_m, A_m reported separately, with chromosome-block bootstrap CIs (seed 42, 1000 reps). A is descriptive only; it is not a confirmatory endpoint and does not gate the disposition. Rationale for demotion: the frozen count statistic is rsID-deduplicated per gene (§3.2), and its distribution across the nonzero head is such that most genes aggregate several distinct significant loci (deduplicated: mean ≈ 16.9 loci/gene, median ≈ 9 across the ~17,300-gene head; see the §6 upper-tail quantiles for the full shape). Because a gene's count changes only when a *distinct* assigned rsID is added or removed, single-locus reassignment is unlikely to move a well-supported gene's rank; perturbation concentrates in sparsely supported genes. It is the *shape* of this distribution — the broad mid-body plus a large tied-zero tail — not any single summary figure, that governs the buffering. (Note: an earlier draft cited "mean 41.8 loci/gene"; that figure is the non-deduplicated row count, 721,406 rows ÷ ~17,300 head genes, and does not describe the frozen rsID-deduplicated statistic. Corrected here to the deduplicated distribution the statistic actually implements.) ρ_GSEA is therefore expected near-degenerate (1 − ρ_GSEA small), making A numerically unstable in exactly the regime most favorable to the mechanism. The mechanistic test therefore rests on boundary-flip localization (§4.1), which is less sensitive to near-degenerate rank concordance (though it carries its own assumptions, notably a non-trivial GSEA flip set — see §4.1 degenerate-case rule).

## 5. Disposition and kill criterion

- **Mechanism corroborated:** C_ORA > C_GSEA (with the GSEA flip distribution not boundary-concentrated relative to ORA). The manuscript may advance the threshold-amplification framing as a property of the decision rule, demonstrated across ORA (thresholded set) and count-based ranked enrichment. Generalization beyond this pipeline (other mappings, databases, p-value-weighted rankings) remains out of scope and is not licensed by this result alone.
- **Mechanism narrowed (kill):** C_ORA ≤ C_GSEA, i.e. ranked enrichment does not reduce boundary-localized flips relative to ORA. The threshold-amplification mechanism claim is withdrawn from its general form; the manuscript is reframed as an ORA-specific annotation-sensitivity result and the higher-ceiling target (PLOS Computational Biology) is not pursued on mechanistic grounds.

No post-hoc adjustment of Δ or the endpoint structure after observation. A (§4.2) is reported regardless but never converts a kill into a pass.

## 6. Mandatory diagnostics (reported unconditionally)

- Nonzero-head size (genes with ≥1 assigned genome-wide-significant locus) under v29 and v45.
- Tied-zero-tail size under v29 and v45, and tail as a fraction of the union universe.
- Full tie structure of the count vectors (distribution of count values; size of each tied block, especially the zero block).
- Upper-tail quantiles of the per-gene count distribution under each version: maximum, 99.9th, 99th, 95th percentiles, and median. This directly addresses whether a small number of hub genes dominate the ranking (the buffering argument rests on this distribution), which reviewers will ask about.
- Fraction of genes whose count changes between v29 and v45. This bridges the 4.15% locus-reassignment rate and downstream pathway behavior — it quantifies how much of the count vector actually moved — without introducing a new hypothesis.
- `fgsea` version and computing environment; enrichment seed and bootstrap seed (both 42).
- Permutation mode used (standard N, or multilevel with eps).
- Reactome pathways retained vs excluded after the 15/500 size filter, under each version (interacts with the parent thin-category finding).
- Consistency check: regenerated ORA ρ and Jaccard vs parent expectations (0.923 / 0.613).
- Provenance record: Catalog release `v1.0.2-associations_e115_r2026-01-19_full`, SHA-256 verification of restored raw file against `bfaffa3a…`.

## 7. Governance prerequisites (must complete before freeze)

These are writes to the clone, required before SHA-locking Study 3:
1. Fill `snapshot: "REPLACE_WITH_RELEASE_ID"` → `"v1.0.2-associations_e115_r2026-01-19_full"` in both `config.yaml` and `prereg/frozen_params.yaml`. Add a `DECISIONS.md` entry recording that the release was recovered retroactively (2026-07-07) from the Trash-directory name and hash-verified against the pre-registered pin, and that it was a placeholder at original freeze. The recovered Catalog file's SHA-256 matches the pre-registered pin `bfaffa3a…` exactly; the governance amendment therefore **restores missing metadata (the human-readable release identifier), it does not change any analytical input** — the bytes the parent study ran on are provably identical.
2. Restore the hash-verified Catalog file (`bfaffa3a…`) from Trash into the clone `data/raw/`; re-verify SHA-256; regenerate `gwas_eligible.tsv` via `03_parse_gwas.py` and confirm it matches the recovered eligible file (closes the derived-artifact provenance loop rather than trusting the loose Trash copy).
3. Add the first `RUN_LOG.md` entry (parent pipeline had zero recorded runs): document the recovery, the release ID, and the Study 3 registration.
4. Regenerate `MEMORY.md` against actual clone disk state; note the prior index described remote/other-workspace contents, not the working tree.
5. Consider (not required for freeze) whether the parent study warrants a provenance note: the published ORA result is reproducible (file hashes correctly) but the snapshot pin was a placeholder and the file survived only in Trash across three workspace deletions.

**Completion status (2026-07-08):** #1 done (snapshot pinned in both files; DEC-004 recorded). #2 done (raw Catalog restored + hash-verified `bfaffa3a…`; `comparison.tsv` regenerated `8f27c1db…` and reproduces the parent ORA result — see the viability probe). #3 done (RUN_LOG entry appended). #4/#5 remain as noted; neither blocks the freeze.

## 8. Analysis sequencing

1. Complete §7 governance prerequisites.
2. Freeze this prereg; compute SHA-256; append to `prereg/` lineage; OSF-register before any Study 3 code executes.
3. On a new branch (`study3-ranked-enrichment`), add stages following the numbered convention (e.g. `src/19_gene_counts.py`, `src/20_gsea_preranked.py`, `src/21_boundary_flip_localization.py`), importing stage-04 mapping, `build_alias_table()`, the 15/500 filter, and the `16_flip_rate.py` flip definition. Additive-only: no modification of any existing artifact, `comparison.tsv`, `results/`, or the manuscript.
4. Build per-gene count vectors (v29, v45) over the shared union universe from `comparison.tsv` (§3.2).
5. Run preranked fgsea under each version with frozen parameters (§3.3); record §6 diagnostics.
6. Regenerate ORA reference (`primary_reactome_export.py`); run the ρ/Jaccard consistency check.
7. Compute C_ORA, C_GSEA and the flip-distance distributions; chromosome-block bootstrap (seed 42, 1000 reps).
8. Compute descriptive A_m components with CIs.
9. Evaluate H1 (C_ORA > C_GSEA); apply §5 disposition.
10. Report all §6 diagnostics unconditionally regardless of disposition.
