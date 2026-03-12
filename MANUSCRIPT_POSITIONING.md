# Manuscript Positioning

Protocol: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Last updated: 2026-03-11

---

## What Paper 1 Established

The predecessor manuscript (annotation_drift_manuscript.pdf, nearly complete as of 2026-03-11) demonstrates that genome annotation version drift propagates into substantial pathway interpretation differences in GWAS analysis.

Specifically, paper 1 established:

- Locus-to-gene assignments change between GENCODE versions at a non-trivial rate (pilot: ~8.3%, 95% CI [6.6%, 10.1%]) even when the reference genome is held constant (GRCh38-on-GRCh38).
- These assignment changes are not uniformly distributed across the genome; they concentrate in annotation-dense or structurally complex regions.
- Assignment changes propagate into Reactome pathway enrichment: Jaccard similarity between pathway result sets from different annotation versions falls measurably below 1.0.
- The effect is present across multiple mapping window sizes and is not an artifact of a single window choice.
- Annotation drift is a scientifically consequential source of variability that is currently underacknowledged in GWAS interpretation workflows.

Paper 1 operates under a two-annotation-version (v29/v37 vs v44/v45) design using a single evidence source (genomic distance), establishing the existence and downstream consequences of annotation drift.

---

## What This Project Adds

The L2G Instability Audit extends paper 1 in three complementary directions:

1. **Scope of perturbation.** Paper 1 varied only annotation version. This project varies annotation version, evidence-source configuration (distance-only, eQTL-only, eQTL+ABC, full stack), and distance rules simultaneously, enabling a full decomposition of instability sources.

2. **Locus-level resolution.** Paper 1 characterizes aggregate instability rates. This project constructs a locus-level instability atlas, enabling identification of which specific loci drive pathway instability and why.

3. **Predictability question.** Paper 1 demonstrates that instability exists and matters; it does not attempt to predict which loci will be unstable. This project tests whether structural locus features (gene density, LD entropy, transcript overlap complexity, etc.) are sufficient to predict instability, and whether topology-derived features add further predictive power.

4. **Multi-trait generalizability.** Paper 1 is conducted on a single trait or a narrow trait set. This project applies identical perturbation experiments across three architecturally distinct genetic architectures (autoimmune, psychiatric, metabolic), with an explicit cross-trait replication rule.

---

## Why This Is a Follow-On

This project is a follow-on, not a replication or extension of paper 1:

- It takes paper 1's central empirical finding (annotation drift propagates to pathway instability) as established background, not as a hypothesis to be tested again.
- It uses the same methodological framework (GENCODE annotation, GWAS locus mapping, Reactome pathway evaluation) but adds new perturbation dimensions and a predictive modeling phase.
- The locus instability atlas produced here could directly enrich the discussion or supplementary materials of paper 1 if paper 1 is not yet submitted, or form the basis of a standalone methods paper otherwise.
- The decision rules (termination conditions, phase gates) are structured so that a negative result (uniform stability, or instability without pathway consequence) cleanly terminates this project without invalidating paper 1's findings.

---

## Intended Audience

Primary: Computational genomics methods community; GWAS analysts who rely on L2G or equivalent prioritization tools.

Secondary: Pathway analysis method developers; researchers designing GWAS follow-up studies who need to understand annotation sensitivity of their locus prioritization step.

Tertiary: Genomics tool developers at resources such as Open Targets, who maintain L2G pipelines and may use this atlas as an empirical benchmark for pipeline stability evaluation.

The paper resulting from this project is positioned as a methods/resource paper rather than a biological discovery paper. The scientific contribution is the instability atlas and the characterization of the stability regime (predictive vs. descriptive).
