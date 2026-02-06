#!/bin/bash
# deploy.sh — One-shot deploy of all updated files into the annotation_drift_gwas project.
# Run from project root:
#   bash deploy.sh
#
# This script writes files directly — no dependency on ~/Downloads.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_ROOT"

echo "Deploying to: $PROJECT_ROOT"
echo ""

# ── 1. Backup existing files ─────────────────────────────────────────────
BACKUP_DIR=".backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"
for f in config.yaml Makefile validate_config.py prereg/prereg.md src/02_parse_gencode.py; do
    if [ -f "$f" ]; then
        mkdir -p "$BACKUP_DIR/$(dirname "$f")"
        cp "$f" "$BACKUP_DIR/$f"
    fi
done
echo "Backups saved to $BACKUP_DIR/"
echo ""

# ── 2. config.yaml (with real SHA-256 hashes) ────────────────────────────
cat > config.yaml << 'CONFIGEOF'
# config.yaml
# Annotation Drift in GWAS Locus-to-Gene Mapping
# Protocol version: 1.0
# Frozen prior to full analysis run
project:
  name: "annotation_drift_gwas_l2g"
  genome_build: "GRCh38"
  protocol_version: "1.0"

paths:
  raw_dir: "data/raw"
  processed_dir: "data/processed"
  results_dir: "results"
  tables_dir: "results/tables"
  figures_dir: "results/figures"
  manuscript_dir: "results/manuscript"
  prereg_dir: "prereg"
  docs_dir: "docs"
  src_dir: "src"

inputs:
  gwas_loci:
    source: "NHGRI-EBI GWAS Catalog"
    snapshot: "REPLACE_WITH_RELEASE_ID"
    format: "snp_list"
    raw_file: "data/raw/gwas_catalog.tsv"
    eligible_file: "data/processed/gwas_eligible.tsv"
    sha256: "bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b"
    columns:
      chr: "chr"
      pos: "pos"
      locus_id: "locus_id"
    inclusion:
      ancestry: "European"
      p_threshold: 5.0e-8
      autosomal_only: true
      locus_definition: "sentinel_snp_only"
    qc:
      allowed_chromosomes: ["1","2","3","4","5","6","7","8","9","10",
                             "11","12","13","14","15","16","17","18",
                             "19","20","21","22"]
      drop_missing: true
      drop_duplicates_by: ["locus_id"]
      merge_loci: false

annotations:
  gencode:
    - tag: "v29"
      build: "GRCh38.p12"
      gtf: "data/raw/gencode.v29.annotation.gtf.gz"
      sha256: "e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9"
    - tag: "v45"
      build: "GRCh38.p14"
      gtf: "data/raw/gencode.v45.annotation.gtf.gz"
      sha256: "56a166e4a43c934a1910c0d09dec786ac048862fc65130731d41f3c7f451ae68"
  gene_def:
    id_field: "gene_id"
    version_strip: true          # ENSG000001234.5 -> ENSG000001234; strip suffix after '.'
    biotypes: ["protein_coding"]
  transcript_def:
    mode: "all_transcripts"
    tss_rule: "strand_aware"

parse:
  tss_table:
    out_template: "data/processed/tss_tables/tss_{tag}.tsv"
    allow_multiple_tss_per_gene: true
    collapse_rule: "nearest_transcript"

map:
  mapping_rule:
    type: "nearest_tss"
    locus_representation: "point"
  windows_kb: [10, 25, 50, 100, 250, 500]
  window_shape: "symmetric"
  allow_unmapped: true
  multi_map:
    enabled: false
  tie_break:
    rule: "lexicographic_gene_id"
  outputs:
    per_version_template: "data/processed/mappings/l2g_{tag}_W{W}kb.tsv"
    comparison_template: "data/processed/mappings/l2g_compare_W{W}kb.tsv"

analysis:
  primary_instability:
    symbol: "I"
    definition: "assigned_gene_id_differs"
    denominator: "mapped_in_both"
  window_stability:
    metric: "max_minus_min"
  uncertainty:
    ci_method: "wilson"
    bootstrap:
      enabled: true
      method: "block"
      block_unit: "chromosome"
      n_reps: 1000
      seed: 42
  outputs:
    tables:
      - "results/tables/instability_by_window.tsv"
      - "results/tables/instability_overall.tsv"
      - "results/tables/instability_by_chromosome.tsv"
      - "results/tables/instability_by_distance.tsv"
      - "results/tables/instability_by_margin.tsv"
      - "results/tables/bootstrap.txt"
      - "results/tables/sensitivity.tsv"
      - "results/tables/difference_classification.tsv"
      - "results/tables/difference_classification_summary.txt"
    figures:
      - "results/figures/I_vs_W.png"
      - "results/figures/delta_distance_hist.png"

enrichment:
  enabled: true
  gene_universe: "union_protein_coding_v29_v45"
  pathway_db:
    name: "Reactome"
    version: "PENDING_DOWNLOAD"
    file: "data/raw/reactome.gmt"
    sha256: "PENDING_DOWNLOAD"
  method:
    type: "ORA"
    test: "fisher"
    multiple_testing: "BH"
    fdr_alpha: 0.05
  pathway_size:
    min: 15
    max: 500
  comparison:
    metric: "jaccard"
    weighting: "unweighted"
    significance_basis: "fdr"
  outputs:
    tables:
      - "results/tables/enrichment_jaccard_by_window.tsv"

sensitivity:
  window_grid:
    enabled: true
    kill_threshold_range: 0.20
  biotype_expansion:
    enabled: true
    biotypes: ["protein_coding","lncRNA","pseudogene"]
    label: "exploratory"
  mapping_rule_alt:
    enabled: true
    rule: "overlap_first"
    label: "exploratory"

repro:
  environment:
    python_version: "3.13.11"
    requirements_file: "requirements.txt"
    requirements_sha256: "7ab802fb0caeb7b1b41aa0a01e5d351653688bc8afbf90cf8f109ed68362855e"
  determinism:
    global_seed: 42
  manifest_out: "results/MANIFEST.sha256"
  locked_config: "prereg/frozen_params.yaml"
CONFIGEOF
echo "  config.yaml"

# ── 3. Makefile ───────────────────────────────────────────────────────────
cat > Makefile << 'MAKEEOF'
# Makefile — Annotation Drift in GWAS Locus-to-Gene Mapping
# Every pipeline target calls `check` first to enforce frozen config fields.

.PHONY: all check freeze download parse map analyze enrich sensitivity clean help

PYTHON  ?= python3
CONFIG  ?= config.yaml
LOCKED  ?= prereg/frozen_params.yaml

# ── Guard ────────────────────────────────────────────────────────────────
check:
	@$(PYTHON) validate_config.py $(CONFIG) $(LOCKED)

# ── Freeze (run once) ────────────────────────────────────────────────────
freeze:
	cp $(CONFIG) $(LOCKED)
	sha256sum $(LOCKED) > $(LOCKED:.yaml=.sha256)
	@echo "Frozen: $(LOCKED)"
	@echo "Hash:   $$(cat $(LOCKED:.yaml=.sha256))"

# ── Pipeline ─────────────────────────────────────────────────────────────
download:
	$(MAKE) check
	@echo "=== Downloading data ==="
	chmod +x src/01_download_data.sh
	./src/01_download_data.sh

parse:
	$(MAKE) check
	@echo "=== Parsing GENCODE v29 ==="
	$(PYTHON) src/02_parse_gencode.py \
		--gtf data/raw/gencode.v29.annotation.gtf.gz \
		--output data/processed/tss_tables/tss_v29.tsv
	@echo "=== Parsing GENCODE v45 ==="
	$(PYTHON) src/02_parse_gencode.py \
		--gtf data/raw/gencode.v45.annotation.gtf.gz \
		--output data/processed/tss_tables/tss_v45.tsv
	@echo "=== Parsing GWAS Catalog ==="
	$(PYTHON) src/03_parse_gwas.py

map:
	$(MAKE) check
	@echo "=== Mapping variants ==="
	$(PYTHON) src/04_map_variants.py

analyze: bootstrap sensitivity classify
	@echo "=== Analysis complete ==="

bootstrap:
	$(MAKE) check
	@echo "=== Bootstrap CI ==="
	$(PYTHON) src/05_bootstrap.py

sensitivity:
	$(MAKE) check
	@echo "=== Sensitivity analysis ==="
	$(PYTHON) src/06_sensitivity.py

classify:
	$(MAKE) check
	@echo "=== Classifying differences ==="
	$(PYTHON) src/07_classify_differences.py

concentration:
	$(MAKE) check
	@echo "=== Concentration analysis ==="
	$(PYTHON) src/08_concentration_analysis.py

# ── Full pipeline ────────────────────────────────────────────────────────
all: download parse map analyze concentration
	@echo "=== Pipeline complete ==="

# ── Utilities ────────────────────────────────────────────────────────────
clean:
	rm -f data/processed/tss_tables/*.tsv
	rm -f data/processed/mappings/*.tsv
	rm -f results/tables/*.tsv results/tables/*.txt
	rm -f results/figures/*.png
	@echo "Cleaned generated outputs."

help:
	@echo "Annotation Drift GWAS — Available targets:"
	@echo ""
	@echo "  make freeze       Snapshot config.yaml -> prereg/frozen_params.yaml"
	@echo "  make check        Validate config against frozen snapshot"
	@echo ""
	@echo "  make download     Download GWAS Catalog and GENCODE GTFs"
	@echo "  make parse        Parse GTFs to TSS tables and filter GWAS"
	@echo "  make map          Run variant-to-gene mapping"
	@echo "  make analyze      Bootstrap + sensitivity + classification"
	@echo "  make concentration  Concentration analysis"
	@echo "  make all          Full pipeline (download -> concentration)"
	@echo ""
	@echo "  make clean        Remove generated files"
	@echo "  make help         This message"
MAKEEOF
echo "  Makefile"

# ── 4. validate_config.py ────────────────────────────────────────────────
cat > validate_config.py << 'VALEOF'
#!/usr/bin/env python3
"""
validate_config.py — Frozen-field guard for annotation_drift_gwas_l2g.

Compares a live config.yaml against a locked snapshot and fails hard if any
protocol-critical field has changed.  Intended to run as the first step of
`make parse`, `make map`, and `make analyze`.

Usage:
    python validate_config.py config.yaml prereg/frozen_params.yaml
    python validate_config.py config.yaml                          # uses prereg/frozen_params.yaml default

Setup (run once after freezing):
    cp config.yaml prereg/frozen_params.yaml
    sha256sum prereg/frozen_params.yaml > prereg/frozen_params.sha256

Exit codes:
    0  All frozen fields match.
    1  At least one frozen field differs (details printed to stderr).
    2  IO / parse error.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Any

import yaml

# ── Frozen fields ──────────────────────────────────────────────────────────
# Each entry: (dotted key path into the YAML, human label for error messages).
# Add new fields here if the protocol is amended.
FROZEN_FIELDS: list[tuple[str, str]] = [
    # Core mapping
    ("map.mapping_rule.type",                   "Mapping rule"),
    ("map.windows_kb",                          "Window grid"),
    ("map.window_shape",                        "Window shape"),
    ("map.tie_break.rule",                      "Tie-break rule"),
    # Instability definition
    ("analysis.primary_instability.denominator", "Instability denominator"),
    # Enrichment
    ("enrichment.gene_universe",                "Gene-set universe"),
    ("enrichment.pathway_db.sha256",            "Pathway DB SHA-256"),
    # Reproducibility
    ("analysis.uncertainty.bootstrap.seed",     "Bootstrap RNG seed"),
    # Input data
    ("annotations.gencode",                     "GENCODE annotation definitions"),
    ("inputs.gwas_loci.raw_file",               "GWAS input file"),
    ("inputs.gwas_loci.sha256",                 "GWAS input SHA-256"),
]


def resolve(cfg: dict, dotted_key: str) -> Any:
    """Walk a nested dict by dotted path.  Raises KeyError if missing."""
    node = cfg
    for part in dotted_key.split("."):
        node = node[part]
    return node


def load_yaml(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def validate(live_path: Path, locked_path: Path) -> list[str]:
    """Return a list of violation messages (empty = pass)."""
    live = load_yaml(live_path)
    locked = load_yaml(locked_path)
    violations: list[str] = []

    for dotted_key, label in FROZEN_FIELDS:
        try:
            locked_val = resolve(locked, dotted_key)
        except KeyError:
            violations.append(
                f"  {label} ({dotted_key}): key missing from locked snapshot"
            )
            continue

        try:
            live_val = resolve(live, dotted_key)
        except KeyError:
            violations.append(
                f"  {label} ({dotted_key}): key missing from live config"
            )
            continue

        if live_val != locked_val:
            violations.append(
                f"  {label} ({dotted_key}):\n"
                f"    locked : {locked_val!r}\n"
                f"    live   : {live_val!r}"
            )

    return violations


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python validate_config.py <config.yaml> [prereg/frozen_params.yaml]",
              file=sys.stderr)
        sys.exit(2)

    live_path = Path(sys.argv[1])
    locked_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("prereg/frozen_params.yaml")

    if not live_path.exists():
        print(f"ERROR: live config not found: {live_path}", file=sys.stderr)
        sys.exit(2)
    if not locked_path.exists():
        print(f"ERROR: locked snapshot not found: {locked_path}", file=sys.stderr)
        print("       Run `cp config.yaml prereg/frozen_params.yaml` after freezing.",
              file=sys.stderr)
        sys.exit(2)

    # Integrity check: warn if the locked file itself was modified
    lock_hash_file = locked_path.with_suffix(".sha256")
    if lock_hash_file.exists():
        expected = lock_hash_file.read_text().strip().split()[0]
        actual = file_sha256(locked_path)
        if actual != expected:
            print(
                f"WARNING: locked snapshot has been modified since freezing.\n"
                f"  expected SHA-256: {expected}\n"
                f"  actual SHA-256:   {actual}\n",
                file=sys.stderr,
            )

    violations = validate(live_path, locked_path)

    if violations:
        print(
            "PROTOCOL VIOLATION: frozen config fields have changed.\n"
            "The following fields differ between the live config and the locked snapshot:\n",
            file=sys.stderr,
        )
        for v in violations:
            print(v, file=sys.stderr)
        print(
            "\nTo proceed, either:\n"
            "  (a) revert config.yaml to match the locked snapshot, or\n"
            "  (b) document a protocol amendment and re-freeze:\n"
            "      cp config.yaml prereg/frozen_params.yaml\n"
            "      sha256sum prereg/frozen_params.yaml > prereg/frozen_params.sha256\n",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print(f"validate_config: all {len(FROZEN_FIELDS)} frozen fields OK.")
        sys.exit(0)


if __name__ == "__main__":
    main()
VALEOF
echo "  validate_config.py"

# ── 5. prereg/prereg.md ──────────────────────────────────────────────────
mkdir -p prereg
cat > prereg/prereg.md << 'PREREGEOF'
# Preregistration: Annotation Drift in GWAS Locus-to-Gene Mapping

**Study Title:** Annotation Version Sensitivity of GWAS Locus-to-Gene Mapping and Pathway Enrichment
**Registration Date:** 2026-01-24
**Author:** Michael Blumhardt
**Status:** Protocol locked; full analysis pending
**Protocol Version:** 1.0

---

## 0. Protocol of Record

`config.yaml` (project root) is the authoritative, machine-readable protocol of record.
`prereg/frozen_params.yaml` is a human-readable mirror that must match `config.yaml`
at all frozen fields. Any divergence constitutes a protocol amendment and must be
documented in Section 7 below before the pipeline is re-run.

`validate_config.py` (project root) enforces this at runtime: every `make` target
calls it before executing, and the pipeline halts on any violation.

---

## 1. Hypotheses

### Primary Hypothesis (H1)
For common-variant GWAS loci, a non-trivial fraction of locus-to-gene assignments
change when using different GENCODE annotation versions (v29 vs v45, both GRCh38),
and this instability may propagate to downstream pathway enrichment results.

### Quantitative Predictions
- **Instability rate (I):** We predict I > 5% under the primary mapping rule
- **Pathway impact:** We predict Jaccard overlap of significant pathways < 0.90

### Null Hypothesis (H0)
Annotation version choice has negligible effect: I < 5% AND pathway Jaccard > 0.90

---

## 2. Methods

### 2.1 Data Sources (Frozen)

| Source | Version | Build | Rationale |
|--------|---------|-------|-----------|
| GWAS Catalog | Current snapshot | GRCh38 | Standard source |
| GENCODE | v29 (July 2018) | GRCh38.p12 | First mature GRCh38 release |
| GENCODE | v45 (Oct 2023) | GRCh38.p14 | Current release |

**Critical design choice:** Both annotation versions use GRCh38 to isolate
annotation drift from liftover/coordinate artifacts.

### 2.2 Inclusion Criteria (Frozen)
- European ancestry
- Genome-wide significant (p < 5x10^-8)
- Autosomal only (chr1-22)
- Sentinel SNP per locus; duplicates collapsed by locus_id

### 2.3 Mapping Rule (Frozen)
**Primary:** Nearest TSS within symmetric window
- TSS definition: strand-aware (+ strand = min(start), - strand = max(end))
- All transcripts considered; gene-level TSS = nearest transcript TSS to locus
- Ties broken by lexicographic gene_id (smallest Ensembl stable ID)
- Gene filter: protein_coding only
- Gene ID: version-stripped Ensembl stable ID (suffix after '.' removed)
- Window grid: 10, 25, 50, 100, 250, 500 kb (fixed; no post-hoc additions)

### 2.4 Instability Metric (Primary Endpoint)

I = (# loci whose assigned gene_id differs) / (# loci mapped in both versions)

Confidence intervals: Wilson score 95% CI (primary); block bootstrap stratified
by chromosome (1,000 replicates, seed = 42) as supplementary.

### 2.5 Pathway Impact (Co-Primary Endpoint)
- Over-representation analysis (Fisher's exact test) on Reactome gene sets
- Gene-set universe: union of all protein_coding genes in v29 or v45
- Pathway size filter: 15-500 genes
- Benjamini-Hochberg FDR < 0.05
- Cross-version comparison: unweighted Jaccard over significant pathways

### 2.6 Sensitivity Analyses (Pre-Committed)

| Analysis | Label | Specification |
|----------|-------|---------------|
| Window grid | Primary | Full grid (10-500 kb); kill if I range > 20% |
| Biotype expansion | Exploratory | protein_coding + lncRNA + pseudogene |
| Mapping rule | Exploratory | Overlap-first (gene body overlap before nearest-TSS fallback) |

No other sensitivity analyses will be run without a protocol amendment.

---

## 3. Stop Rules (Hard)

| Condition | Criterion | Action |
|-----------|-----------|--------|
| Negative result | I < 5% at all W AND Jaccard > 0.90 | Publish as negative / validation result |
| Instability excess | I range > 20% across window grid | Reframe scope or terminate |
| Pipeline failure | < 90% of loci mapped in either version | Debug before proceeding |

---

## 4. Frozen Fields (Enforced by Validator)

The following config.yaml fields are locked by `validate_config.py` and will
cause a hard pipeline failure if modified without a formal amendment:

1. `map.mapping_rule.type` — Mapping rule
2. `map.windows_kb` — Window grid
3. `map.window_shape` — Window shape
4. `map.tie_break.rule` — Tie-break rule
5. `analysis.primary_instability.denominator` — Instability denominator
6. `enrichment.gene_universe` — Gene-set universe
7. `enrichment.pathway_db.sha256` — Pathway DB SHA-256
8. `analysis.uncertainty.bootstrap.seed` — Bootstrap RNG seed
9. `annotations.gencode` — GENCODE annotation definitions (tags, builds, files, hashes)
10. `inputs.gwas_loci.raw_file` — GWAS input file
11. `inputs.gwas_loci.sha256` — GWAS input SHA-256

---

## 5. Deliverables

1. Quantitative bound: I estimate with CI for GWAS locus-to-gene mapping
2. Pathway propagation test: does gene-list drift translate to enrichment drift?
3. Reproducible pipeline with pinned versions and config validation

---

## 6. Pre-Registered Figures and Tables

| Output | Description |
|--------|-------------|
| I vs W curve | Primary instability across the window grid |
| Delta-distance histogram | Distribution of distance-to-gene changes across versions |
| Enrichment Jaccard by W | Pathway overlap summary across window sizes |
| Instability by chromosome | Per-chromosome breakdown |
| Sensitivity table | I under alternative biotype/mapping-rule specifications |

---

## 7. Amendments

Any deviations from this protocol will be documented here with justification
before the affected analysis step is re-run.

| Date | Amendment | Justification |
|------|-----------|---------------|
| 2026-01-24 | v19 -> v29 | v19 is GRCh37; v29 and v45 are both GRCh38, isolating annotation drift from liftover artifacts |
| 2026-02-06 | Protocol hardening | Reconciled config.yaml with Locked Analyses doc; added validate_config.py with 11 frozen fields; designated config.yaml as protocol of record |
PREREGEOF
echo "  prereg/prereg.md"

# ── 6. src/02_parse_gencode.py ────────────────────────────────────────────
cat > src/02_parse_gencode.py << 'PARSEEOF'
#!/usr/bin/env python3
"""Parse GENCODE GTF to TSS table.

Reads a GENCODE comprehensive GTF and extracts one TSS per gene
(most-upstream transcript TSS on + strand, most-downstream on - strand).

Protocol-critical behaviour (see Locked Analyses items 4, 5, 7, 8, 10):
  - Autosomes only: chr1-22 (X/Y/MT excluded)
  - Biotype filter: protein_coding (default; overridable via --biotype)
  - Gene ID: Ensembl stable ID, version-stripped (ENSG000001234.5 -> ENSG000001234)
  - TSS: strand-aware, collapsed per gene via min(+) / max(-)
"""

import gzip
import csv
import argparse
from collections import defaultdict

AUTOSOMES = frozenset(str(i) for i in range(1, 23))


def parse_attributes(attr_string):
    """Parse GTF attribute column into a dict."""
    attrs = {}
    for attr in attr_string.strip().split(';'):
        attr = attr.strip()
        if not attr or ' ' not in attr:
            continue
        key, value = attr.split(' ', 1)
        attrs[key] = value.strip('"')
    return attrs


def parse_gtf_to_tss(gtf_path, output_path, biotype_filter="protein_coding"):
    """Extract one TSS per gene from a GENCODE GTF.

    Parameters
    ----------
    gtf_path : str
        Path to GENCODE comprehensive GTF (plain or .gz).
    output_path : str
        Output TSV path.
    biotype_filter : str
        Gene biotype to retain (default: protein_coding).
    """
    genes = defaultdict(lambda: {
        "tss_positions": [], "symbol": None, "strand": None, "chrom": None
    })

    opener = gzip.open if gtf_path.endswith('.gz') else open
    print(f"Parsing {gtf_path}...")

    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            chrom, _, _, start, end, _, strand, _, attributes = fields
            attrs = parse_attributes(attributes)

            # Version-stripped Ensembl stable ID (protocol item 8)
            gene_id = attrs.get('gene_id', '').split('.')[0]
            gene_type = attrs.get('gene_type', '')

            if biotype_filter and gene_type != biotype_filter:
                continue

            # Autosomal only: chr1-22 (protocol items 4-5)
            chrom_clean = chrom.replace('chr', '')
            if chrom_clean not in AUTOSOMES:
                continue

            # Strand-aware TSS (protocol item 10)
            tss = int(start) if strand == '+' else int(end)
            genes[gene_id]["tss_positions"].append(tss)
            genes[gene_id]["symbol"] = attrs.get('gene_name', '')
            genes[gene_id]["strand"] = strand
            genes[gene_id]["chrom"] = chrom_clean

    output_rows = []
    for gene_id, data in genes.items():
        if not data["tss_positions"]:
            continue
        # Gene-level collapse: most-upstream (+) or most-downstream (-)
        tss = (min(data["tss_positions"]) if data["strand"] == '+'
               else max(data["tss_positions"]))
        output_rows.append({
            'ensembl_id': gene_id,
            'gene_symbol': data["symbol"],
            'chromosome': data["chrom"],
            'tss': tss,
            'strand': data["strand"],
        })

    output_rows.sort(key=lambda r: (int(r['chromosome']), r['tss']))

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=['ensembl_id', 'gene_symbol', 'chromosome', 'tss', 'strand'],
            delimiter='\t',
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"  Wrote {len(output_rows)} genes to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse GENCODE GTF to TSS table (autosomal, protein_coding)."
    )
    parser.add_argument('--gtf', required=True, help="Path to GENCODE GTF (.gz)")
    parser.add_argument('--output', required=True, help="Output TSV path")
    parser.add_argument('--biotype', default='protein_coding',
                        help="Gene biotype filter (default: protein_coding)")
    args = parser.parse_args()
    parse_gtf_to_tss(args.gtf, args.output, args.biotype)
PARSEEOF
echo "  src/02_parse_gencode.py"

# ── 7. Summary ────────────────────────────────────────────────────────────
echo ""
echo "All files deployed. Verify:"
echo ""
echo "  grep AUTOSOMES src/02_parse_gencode.py"
echo "  grep sha256 config.yaml | head -5"
echo ""
echo "Then freeze:"
echo ""
echo "  make freeze"
echo "  make check"
