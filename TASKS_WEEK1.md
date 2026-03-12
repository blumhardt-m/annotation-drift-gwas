# Week-1 Task Checklist

Protocol: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Week starting: 2026-03-11
Operator: Michael Blumhardt

Mark each item [ ] incomplete, [x] complete, [!] blocked.

---

## 1. Dataset Verification

- [ ] 1.1 Confirm GENCODE v38 GTF file is present, record checksum in DATA_SOURCES.md
      Completion criteria: File exists at expected path; SHA-256 checksum recorded; file size plausible (>500 MB).

- [ ] 1.2 Confirm GENCODE v45 GTF file is present, record checksum in DATA_SOURCES.md
      Completion criteria: File exists at expected path; SHA-256 checksum recorded.

- [ ] 1.3 Select final traits within each category (autoimmune, psychiatric, metabolic); confirm GWAS summary statistics exist for each
      Completion criteria: Three traits named; three GWAS summary stat files located; paths recorded in configs/traits_to_run.yaml with TODO_VERIFY markers resolved.

- [ ] 1.4 Confirm L2G score availability for selected traits (OpenTargets or equivalent)
      Completion criteria: L2G score files located or access confirmed; paths recorded in DATA_SOURCES.md.

- [ ] 1.5 Confirm eQTL dataset availability (GTEx or equivalent) for eQTL-only and eQTL+ABC evidence configs
      Completion criteria: eQTL file located; version recorded; path in DATA_SOURCES.md.

- [ ] 1.6 Confirm ABC score dataset availability for eQTL+ABC evidence config
      Completion criteria: ABC score file located or access confirmed; path in DATA_SOURCES.md.

- [ ] 1.7 Confirm 1000 Genomes Phase 3 EUR LD reference panel is present and locked
      Completion criteria: LD panel files present; version confirmed; path recorded in configs/ld_reference.yaml with TODO_VERIFY resolved.

- [ ] 1.8 Confirm Reactome GMT file is present
      Completion criteria: File present at data/raw/ReactomePathways.gmt; version recorded in DATA_SOURCES.md.

- [ ] 1.9 Confirm MSigDB C2 GMT file availability
      Completion criteria: File located or access confirmed; path in DATA_SOURCES.md.

---

## 2. Pipeline Baseline Setup

- [ ] 2.1 Execute baseline pipeline for one selected trait (annotation v38, distance-only evidence, nearest-gene distance rule)
      Completion criteria: Pipeline completes without error; output TSV produced at results/locus_instability_atlas/; row count plausible for trait locus count; run logged in RUN_LOG.md.

- [ ] 2.2 Verify baseline output schema matches OUTPUT_SCHEMA.md
      Completion criteria: All required columns present; types correct; composite key (gwas_trait_id + locus_id) is unique in output.

- [ ] 2.3 Record baseline run in RUN_LOG.md
      Completion criteria: RUN_LOG.md entry includes date, operator, script, inputs, outputs, protocol deviations (none expected), and notes.

---

## 3. Perturbation Configuration

- [ ] 3.1 Resolve all TODO_VERIFY markers in configs/perturbation_matrix.yaml
      Completion criteria: All placeholder paths replaced with confirmed values or explicitly flagged as blocked with a note.

- [ ] 3.2 Resolve all TODO_VERIFY markers in configs/traits_to_run.yaml
      Completion criteria: All three trait entries have confirmed gwas_id, summary_stats_path, and l2g_scores_path.

- [ ] 3.3 Resolve all TODO_VERIFY markers in configs/ld_reference.yaml
      Completion criteria: ld_panel_path confirmed; checksum populated.

- [ ] 3.4 Confirm that perturbation matrix covers minimum 4 runs per trait
      Completion criteria: configs/perturbation_matrix.yaml lists baseline + annotation perturb + evidence perturb + distance perturb for each trait.

---

## 4. Instability Atlas Generation

- [ ] 4.1 Run all 4 minimum perturbation configurations for one trait; produce per-locus top-1 change rate
      Completion criteria: Four output TSVs exist; top-1 change rate column populated for each perturbation pair; values are in [0, 1].

- [ ] 4.2 Compute rank entropy per locus for the one-trait pilot
      Completion criteria: rank_entropy column populated in atlas for all loci; distribution is non-degenerate.

- [ ] 4.3 Verify instability heterogeneity test is computable (Q1/Q4 quantile split)
      Completion criteria: Locus count is sufficient for a 4-quantile split; Q1 and Q4 means computed and recorded; ratio logged.

- [ ] 4.4 Run Reactome pathway evaluation for baseline and annotation-perturbed run (one trait)
      Completion criteria: Two enrichment result files produced; Jaccard similarity computed between the two; value recorded in atlas or separate summary table.

- [ ] 4.5 Commit first atlas TSV to results/locus_instability_atlas/
      Completion criteria: File committed to git with informative commit message referencing protocol ID and run date.

---

## 5. Logging

- [ ] 5.1 Create initial RUN_LOG.md entry for the baseline run
      Completion criteria: Entry present; all template fields populated; no blank required fields.

- [ ] 5.2 Confirm DECISIONS.md reflects current state (all 6 initial entries present)
      Completion criteria: DEC-001 through DEC-006 are present and correctly dated.

- [ ] 5.3 Update DATA_SOURCES.md status column for all verified inputs
      Completion criteria: No entry has status "unverified" unless explicitly blocked with a reason.

- [ ] 5.4 Archive this checklist with completion state at end of week
      Completion criteria: Copy of TASKS_WEEK1.md with final checkbox states committed to archive/ or noted in RUN_LOG.md.
