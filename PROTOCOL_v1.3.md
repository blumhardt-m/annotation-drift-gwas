LOCKED PROTOCOL DOCUMENT v1.3
Project: L2G Instability Audit Under Annotation and Evidence-Source Drift
Program: 08_GENOMICS / annotation_drift_gwas
Protocol ID: L2G_INSTABILITY_AUDIT_PROTOCOL_v1.3
Lock Date: 2026-03-11
Status: FINAL — EXECUTION AUTHORIZED

Strategic Objective
Establish the empirical ground truth of instability in locus-to-gene prioritization pipelines under controlled perturbations of annotation versions and evidence-source configurations.

The project answers three sequential questions:

Existence — how much instability exists in L2G assignments?
Structure — does instability concentrate in specific loci or genomic contexts?
Predictability — can structural locus features predict instability?

Topology-derived predictors (e.g., spectral features) are evaluated only after baseline predictors are tested.

Core Hypothesis
L2G prioritization outputs occupy heterogeneous stability regimes.

Some loci remain stable across pipeline perturbations.

Some loci are intrinsically unstable due to locus architecture or annotation sensitivity.

A minority of unstable loci may disproportionately alter downstream pathway conclusions.

Trait Panel
Instability will be evaluated across three GWAS traits representing distinct genetic architectures.

Autoimmune example: rheumatoid arthritis or inflammatory bowel disease
Psychiatric example: schizophrenia or major depression
Metabolic example: BMI or type-2 diabetes

All perturbation experiments must be applied identically across all traits.

Trait Replication Rule
The instability pattern must replicate in at least 2 of the 3 traits.

Replication requires both instability heterogeneity and material downstream pathway impact.

Failure to satisfy both conditions in two traits triggers project termination.

Pipeline Perturbation Design

Perturbation dimensions:

Annotation version
GENCODE v38
GENCODE v45

Evidence-source configuration examples
distance-only
eQTL-only
eQTL + ABC
full evidence stack

Distance rule examples
nearest gene
±250 kb window
±500 kb window

Baseline-shift design:

baseline pipeline
+ annotation perturbation
+ evidence perturbation
+ distance perturbation

Minimum runs per trait: 4

LD reference control
LD reference panel must remain constant across all perturbations.

Example reference:
1000 Genomes Phase 3 EUR

Pathway evaluation

Primary arbiter: Reactome
Sensitivity database: MSigDB C2

Circularity constraint

Evidence sources allowed in L2G:

genomic distance
eQTL
chromatin interaction
regulatory annotations
variant effect predictions

Forbidden in L2G assignment:

Reactome
MSigDB
STRING

Instability Metrics

top-1 assignment change rate
top-k assignment instability
rank entropy

Material Downstream Impact

Primary rule:
Top 10% of unstable loci account for ≥30% of pathway Jaccard degradation.

Secondary rule:
Jaccard similarity < 0.65.

Instability heterogeneity

mean(top1_change_rate_Q4) ≥ 3 × mean(top1_change_rate_Q1)

Phase B predictive modeling

Baseline predictors include:

gene density
LD entropy
transcript overlap complexity
nearest-gene distance asymmetry
regulatory link multiplicity
regulatory evidence entropy

Decision thresholds:

Predictive success
AUC > 0.60

Descriptive success
AUC ≤ 0.60

Proceed to Phase C only if baseline AUC < 0.75.

Spectral retention threshold

ΔAUC ≥ 0.05

If bootstrap CI width exceeds threshold, spectral results are considered underpowered exploratory analysis and reported in supplementary material.

Project outcomes

Predictive instability regime → predictive methods paper
Descriptive instability regime → descriptive instability audit
Instability present but pathway-inconsequential → terminate project
Uniform stability → terminate project

Authorization
Research Director
Michael Blumhardt
