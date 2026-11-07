# Per-trait Annotation-Version Pilot — Preregistration v1.0 (candidate for freeze)

**Parent study:** Annotation version as a hidden degree of freedom in GWAS pathway interpretation (BIOINF-2026-0616; repo `annotation-drift-gwas` @ tag `v1.0_submission`, commit `37272b2`; public mirror HEAD `61df5e0`; local clone of record `~/Study/research/08_GENOMICS/annotation-drift-gwas-clean`).
**Supersedes (without execution):** Study 3, `prereg/study3_prereg_v1.0.md` sha256 `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`, OSF `fv4sg` (repo DEC-013).
**Specification source:** DECISION_PACKAGE_annotation_drift_2026-11-07_v1.3 §D (pilot spec v0.3), owner-adopted 2026-11-07 14:27 -08:00 (D-D1). The text of §1–§10 below is that specification verbatim, re-numbered.
**Basis:** GWAS_ANNOTATION_DRIFT_INDEPENDENT_ASSESSMENT_2026-11-07_v1.2 (md sha256 `91fdc8fbc64ff694ccd061813a3fadaeb060c1e12d14d751c5aec13e786aadd1`); VERIFICATION_REPORT_20261107.md (sha256 `329c432bd55050d56ee5b2e300565d81869c4f6dc7fc56da800ca986c53c837e`).
**Status:** CANDIDATE — becomes FROZEN v1.0 when the owner writes `prereg/pertrait_pilot_prereg.sha256` over this file's bytes and records the freeze in `DECISIONS.md`. Post-freeze changes require a documented protocol amendment.
**Registration basis:** every parameter in §2–§4 and every rule in §8 is fixed before any per-trait enrichment is run. One field is deferred by design: the frozen-roster hash (§3), written by the roster script before the first enrichment run and recorded in `RUN_LOG.md`; until then its value is NOT_RECORDED.
**Authorization:** this document authorizes nothing. Execution requires the §10 conditions, including a portfolio live-slot allocation or an explicit owner exception.

## 0. Frozen inputs (identities)

| Input | Identity |
|---|---|
| GWAS Catalog associations snapshot | `data/raw/gwas_catalog.tsv`, sha256 `bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b`, release `v1.0.2-associations_e115_r2026-01-19_full` |
| Eligible associations | `data/processed/gwas_eligible.tsv`, sha256 `8f7b2994d7e62da7845b832eb959a611a9b7f81603cfc019b69d22faa5edd7e2` (746,850 rows; 307,102 distinct SNPS values, 97 compound) |
| Parent mapping table (reference only) | `data/processed/mappings/comparison.tsv`, sha256 `8f27c1dba92c2af3ec4c427b1bf89cdb03d67f68adef09c1c9131111f6ba2bb1` |
| Parent frozen parameters | `prereg/frozen_params.yaml`, sha256 `8407f51ff47b8ba4cdfd0f58b22a23a94f3d221dd295ae828cbee225ec01322f` |
| GENCODE v29 / v37 / v45 GTFs | as pinned in `01_download_data.sh` and `frozen_params.yaml` (content hashes therein) |
| Reactome GMT v95 | sha256 `a34f43fd…` per `frozen_params.yaml`; Ensembl conversion via the existing alias table |
| GO:BP set | the Study 1 input, identity as recorded in `results/study1_GO_replication/` |
| Frozen roster (≥ 100) | `prereg/pertrait_roster_ge100.tsv` — sha256 NOT_RECORDED until written by the roster script (§3) |
| Extension roster (≥ 50) | `prereg/pertrait_roster_ge50.tsv` — written only if §8 rule 4 or 5 fires; sha256 NOT_RECORDED until then |


## 1. Scientific question and estimand

**Object of study (stated exactly):** Catalog-derived, positionally thinned, per-trait association sets — for each `DISEASE/TRAIT` string in the frozen GWAS Catalog snapshot, the genome-wide-significant rsIDs after per-trait deduplication and ±250 kb positional thinning (§3). This is *not* an ancestry-controlled, LD-defined independent-locus GWAS: the associations file carries no per-row ancestry, and positional thinning is a proxy for LD clumping. Every conclusion transports to that object and is written in those terms.

**Question:** for such per-trait sets (sparse regime), how often and by how much does the set of significant pathways change when only the GENCODE release used for locus-to-gene assignment changes, with everything else frozen?

**Primary estimand**, per trait *t* and release pair *(a, b)*: the number of significance flips *F_t* (pathways significant under exactly one release at BH-FDR 0.05) and the significant-set Jaccard *J_t*, summarized over the informative subset of the frozen roster by distribution (§8), with the uninformative count reported alongside.

**Secondary estimands:** per-trait Spearman ρ of −log10 p across releases; per-trait distance-to-boundary distribution of flipped pathways; the C2 comparison (§6); the ranked-method, secondary-pair and secondary-database results.

## 2. Inputs (all frozen already)

- GWAS Catalog snapshot `bfaffa3a…` (release `v1.0.2-associations_e115_r2026-01-19_full`); `gwas_eligible.tsv` (746,850 rows; columns SNPS, CHR_ID, CHR_POS, P-VALUE, DISEASE/TRAIT, STUDY).
- GENCODE v29, v37, v45 GTFs as frozen in `01_download_data.sh` / `frozen_params.yaml`; the v37 TSS table from Study 2.
- Reactome GMT v95 (`a34f43fd…`) converted to Ensembl IDs with the existing alias table; the GO:BP set used in Study 1 as a secondary database.
- Mapping: nearest-TSS, ±500 kb, imported from `src/04_map_variants.py` unchanged.

## 3. Trait roster rule (frozen before any per-trait run)

- Unit: the `DISEASE/TRAIT` string as recorded (no manual merging of synonyms; this is a known coarse choice and is reported as such).
- Per trait: keep rows with P-VALUE < 5e-8 and a strict single rsID (`rs[0-9]+`); compound identifier strings are excluded and counted.
- Deduplication: one row per (trait, rsID), keeping the minimum p across studies of that trait.
- Sentinel thinning within trait: greedy, by ascending p, retaining an rsID only if no retained rsID lies within ±250 kb on the same chromosome (positional proxy for LD; the parent study's frozen `sentinel_snp_only` was never operationalized — this is the operationalization, and it is a declared degree of freedom). Tie-break for equal p: ascending CHR_ID (numeric), then ascending CHR_POS, then lexicographic rsID string; fully deterministic, no RNG.
- Ancestry: the Catalog associations file carries no per-row ancestry; the frozen "European ancestry" criterion cannot be applied from this file. Declared as not applied, with reason.
- Eligibility: traits with ≥ 100 retained sentinel rsIDs after thinning, autosomes only. Roster size *N* is whatever this rule yields and is reported before any enrichment is run (freeze the roster list and its hash).
- No post-hoc roster edits. Traits may be dropped only by the predeclared rules in §8.

## 4. Per-trait enrichment (two methods)

- ORA: foreground = distinct Ensembl genes assigned to the trait's sentinel rsIDs under release *a*; universe = frozen union protein-coding set (19,315); Fisher one-sided, BH across the 15–500-size pathway set, FDR 0.05 — the parent apparatus, imported.
- Ranked: preranked `fgseaMultilevel` on a per-gene statistic defined *within trait* as max over the gene's assigned sentinel rsIDs of −log10 p (trait-specific, so the cross-trait min-p aggregation confound identified in the Study 3 prereg does not arise); genes with no assigned rsID receive 0 (tied tail; reported). minSize 15, maxSize 500, `fgseaMultilevel(eps = 0, sampleSize = 101, nproc = 1)` preceded by `set.seed(42L)` once per trait × release × database call, R 4.5.3 / fgsea 1.36.2 as pinned in Study 3; ties in the per-gene statistic broken by Ensembl ID lexicographic order before the call.
- Annotation pairs: (v29, v45) primary; (v37, v45) secondary. Reactome primary; GO:BP secondary.

## 5. Per-trait metrics

*J_t*, *F_t* (split into a-only / b-only), ρ_t, top-K overlap at K = 25/50 (reported only when ≥ K pathways are tested), and for every flipped pathway the distance of its adjusted p from the cutoff under both releases and the overlap-count change (Δ genes).

## 6. C2 — perturbation null for gene-to-gene reassignment (executable definition)

Purpose: compare, for each trait, the flip count produced by the *observed gene-to-gene reassignments alone* against the distribution produced by the same number of random gene-to-gene reassignments of matched locality under thresholding. C2 is deliberately restricted to gene-to-gene changes; unmapped→mapped and mapped→unmapped changes are a different perturbation type and are handled only in the full sensitivity analysis (§5), where *J_t* and *F_t* include them.

Definitions, all frozen:
- *Reference release:* *a* (the earlier release of the pair). Candidate genes and TSSs come from release *a*'s frozen TSS table only.
- *Observed gene-to-gene set* for trait *t*: rsIDs mapped under both releases whose assigned Ensembl gene differs; size *n₁*. Also reported, but excluded from C2: *n₂⁺* (unmapped under *a*, mapped under *b*) and *n₂⁻* (mapped under *a*, unmapped under *b*).
- *Observed comparator* *F_t^{g2g}*: starting from release *a*'s mapping, apply exactly the observed gene-to-gene reassignments (release *b*'s gene for each of the *n₁* rsIDs; all other rsIDs unchanged, including those unmapped under *b*), recompute the foreground and ORA, and count flips against release *a*'s unperturbed result. This is the only observed quantity C2 compares against; *F_t* from §5 is not used in C2.
- *Eligible rsID:* mapped under *a* with ≥ 1 protein-coding gene in release *a* other than its assigned gene with a TSS within ±500 kb; its *alternative set* is those genes, ordered lexicographically by Ensembl ID.
- *Synthetic replicate r* (r = 1…1000): draw *k₁ = min(n₁, |eligible|)* eligible rsIDs without replacement; replace each one's gene by an alternative drawn uniformly by seeded index over the lexicographic order; recompute the foreground and ORA; count flips against release *a*'s unperturbed result.
- *Shortfall:* if *k₁ < n₁*, the trait's null is TRUNCATED, the shortfall is reported, and the trait is excluded from the roster-level C2 calibration (not from §5 summaries). A trait with *n₁* = 0 has no C2 (no perturbation to match) and is likewise excluded from the calibration, reported as NO-DRIFT.
- *Seeds:* numpy `default_rng(42 + i)`, *i* = 0-based index in the frozen roster file; one generator per trait, consumed in replicate order.
- *Outputs:* per trait, the synthetic distribution of *F*, its 95th percentile, *F_t^{g2g}*, and *E_t* = [*F_t^{g2g}* > 95th percentile]. Roster-level: *x* = fraction of traits with *E_t* = 1 among the C2-available traits (§8 defines availability), and the C2 denominator *d_x*.

Interpretation rule, frozen: *F_t^{g2g}* within the null → consistent with thresholding of randomly relocated evidence; no mechanism claim beyond it is supported. Roster-level *x* well above the nominal 5% → observed reassignments produce more flips than random local relocation of the same count, and a further explanation is required. Consistency with a model is not proof of it.

## 7. Chromosome-deletion sensitivity (not confidence intervals)

The pipeline deduplicates rsIDs within trait and genes within foreground and consumes gene *sets*, so a chromosome-block bootstrap collapses to subset sampling, and no resampling scheme here would validate interval coverage for thresholded counts, medians over a changing informative subset, or C2 indicators. Accordingly no confidence intervals are reported. The uncertainty-adjacent output is a **leave-one-chromosome-out sensitivity analysis** (22 autosomal deletions, deterministic):

- For each chromosome *c*, remove every rsID on *c* from every trait's thinned set simultaneously; recompute mapping under both releases, enrichment, *J_t*, *F_t*, ρ_t for all traits, and the roster-level summaries *m_J*, *m_F*, *τ*, *c* (coverage).
- Under each deletion, C2 availability and *x* are **held at their full-data values**; only the non-C2 quantities (*J_t*, *F_t*, ρ_t, *u*, *c*, *m_J*, *m_F*, *τ*) are recomputed, on the unchanged frozen roster. Report, per trait and per summary: the full-data value, the minimum and maximum across the 22 deletions, and the deletion(s) attaining them; the number of traits whose informative status changes under any deletion; and the number of deletions under which the §8 outcome — recomputed descriptively with the same thresholds and the full-data C2 — would differ from the full-data outcome.
- Simultaneous deletion preserves the shared perturbation across overlapping traits; this is what makes the roster-level sensitivity meaningful as a sensitivity, and is all that is claimed for it.
- C2 is not re-run under deletion and no C2 sensitivity is reported; C2 is computed once on full data (§6).
- The §8 decision is taken once, on full data. The sensitivity table is descriptive and never changes an outcome. A full-data outcome that would differ under ≥ 3 deletions is flagged FRAGILE; both the outcome-change counts and the FRAGILE flag are labelled "conditional on full-data C2" — they measure sensitivity of the non-C2 inputs, not of the complete decision procedure.

## 8. Stopping rule (frozen before execution; outcomes mutually exclusive, evaluated in the stated precedence)

Status of the thresholds: the numbers below (0.60, 0.80, 0.90, 0.10, 0.25, 10, 0.5 *n_i*, ≥ 3) are **predeclared project-continuation criteria** — commitments about what the owner will do on each outcome — not validated statistical error guarantees; no coverage, power or error-rate property is claimed for them. Any "bounding note" produced under rule 3 describes the observed frozen roster and the specified informative and C2-available subsets, not a population of traits or analyses.

Definitions (primary analysis = pair v29/v45, Reactome, ORA, FDR 0.05):
- *Uninformative trait:* significant set empty under both releases. Excluded from *J*, *F*, *τ* and C2 summaries; the count *u* is always reported. A trait with an empty set under exactly one release is informative with *J_t* = 0, *F_t* = size of the non-empty set.
- *Informative coverage* *c* = (roster − *u*) / roster; *n_i* = roster − *u*.
- *Central tendency:* *m_J* = median *J_t*, *m_F* = median *F_t*, over informative traits.
- *Tail* *τ* = fraction of informative traits with *J_t* < 0.5 **or** *F_t* ≥ 3.
- *C2 availability:* *d_x* = number of informative traits with a non-truncated C2 null and *n₁* > 0. C2 is **available** if *d_x* ≥ 10 and *d_x* ≥ 0.5 *n_i*; otherwise unavailable, *x* is undefined and is reported as C2-UNAVAILABLE with *d_x*. When available, *x* = fraction of those *d_x* traits with *E_t* = 1 (so *x* ≥ 0.25 requires ≥ 3 exceedances).

Scope rule: every roster-level statement is a statement about the informative subset and is written as such ("among the *n_i* of *N* roster traits with any significant Reactome pathway at FDR 0.05…"). No bound is reported over all traits. C2 statements are scoped to the *d_x* C2-available traits.

Precedence and outcomes (first matching rule applies; exactly one outcome):
1. **UNDERPOWERED** if *c* < 0.60. Terminal for the primary analysis: report *u*, *c*, and all descriptive outputs; no PASS/KILL. By policy no extension follows this outcome. FDR 0.10 results are reported as secondary descriptives only and produce no outcome.
2. **PASS** if *m_J* < 0.80, or *τ* > 0.25, or (C2 available and *x* ≥ 0.25). Outcome: the narrowed manuscript (assessment §11), carrying the C2 calibration (if available) and the ranked-method contrast as secondary.
3. **KILL** if C2 is available and all of: *m_J* ≥ 0.90; *m_F* ≤ 1; *τ* ≤ 0.10; *x* ≤ 0.10. Outcome: a bounding note over the observed frozen roster's informative subset and C2-available subset, stating that among the informative traits of this roster release choice changed the Reactome ORA significant set by more than two pathways, or by Jaccard < 0.5, for at most ~10% of traits, with the roster, *u*, *d_x*, denominators and C2 calibration; project parks. (Rule 2 is tested first, so KILL never coincides with any PASS condition.)
4. **C2-UNAVAILABLE** if C2 is unavailable and no rule-2 condition holds. Outcome: no KILL is possible without C2; proceed to the predeclared extension (which may raise *d_x*), then re-evaluate rules 1–4; if C2 is still unavailable after the extension, C2-UNAVAILABLE-FINAL: full report of the §5 results with the non-C2 KILL conditions stated as met or not met, no bounding claim, project parks pending an owner decision made outside this preregistration.
5. **INDETERMINATE** otherwise (C2 available, no PASS, not all KILL conditions). Outcome: the single predeclared extension, then re-evaluation (below).

Predeclared extension (fixed now; the only one): **roster extension to traits with ≥ 50 retained sentinels** (same database, same pair, same thresholds), adding the 50–99-sentinel traits to the frozen ≥ 100 roster as a second frozen roster file with its own hash. Rules 1–5 are re-evaluated once on the extended roster with identical thresholds; rule 5 on the extended roster is terminal — INDETERMINATE-FINAL, full report, no claim beyond description, project parks pending an owner decision made outside this preregistration. GO:BP, the v37/v45 pair, the ranked method, and FDR 0.10 are secondary throughout: reported in full, never decisive, never a substitute extension.

Rules are evaluated once per roster on full data by a script that reads only the frozen parameters, after all per-trait runs for that roster complete; the §7 sensitivity table never changes an outcome.

## 9. Deliverables, effort, runtime (estimates, unmeasured)

Deliverables: frozen roster file + hash; per-trait tables (one row per trait × method × pair × database); flipped-pathway table; C2 null tables (gene-to-gene); chromosome-deletion sensitivity tables; the rules 1–5 outcome computed by a script that reads only the frozen rule parameters; `RUN_LOG.md` entries; SHA256SUMS. Build: 2–4 agent sessions reusing stages 03/04/12 and the Study 3 viability scripts' fgsea path. Runtime: per-trait ORA is seconds; fgsea per trait seconds–minutes; the C2 null (N traits × 2 pairs × 1000 synthetic mappings × ORA) is the dominant cost — minutes per trait, hours overall; the chromosome-deletion sensitivity is 22 full recomputations of the per-trait pipeline (without C2), comparable in cost to the base run. All unmeasured.

## 10. Governance provisions

- Preregister as `prereg/pertrait_pilot_prereg_v1.0.md` (from a reviewed spec; v0.3 is the current draft) with SHA sidecar; OSF registration optional (owner's call; no embargo needed if §C is adopted).
- New branch `pertrait-pilot`; additive-only; frozen computations imported, never reimplemented; `pre_auth_check.sh` must pass; agent prompt carries the seven-field header and bounded-autonomy authorization for the whole pilot as one unit of work.
- Execution requires: §A done, §B deposited, §C adopted (or explicitly deferred with Study 3 left BLOCKED), register allocation of a live slot or an explicit owner exception, and the frozen roster hash recorded before the first enrichment run.


