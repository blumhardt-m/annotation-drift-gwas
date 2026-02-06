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

# ── Paths (derived from config.yaml) ──────────────────────────────────────
GWAS_RAW     := data/raw/gwas_catalog.tsv
GWAS_ELIGIBLE:= data/processed/gwas_eligible.tsv
TSS_V29      := data/processed/tss_tables/tss_v29.tsv
TSS_V45      := data/processed/tss_tables/tss_v45.tsv
COMPARISON   := data/processed/mappings/comparison.tsv
BOOTSTRAP_OUT:= results/tables/bootstrap.txt
SENSITIVITY  := results/tables/sensitivity.tsv
CLASSIFY_OUT := results/tables/difference_classification.tsv
CLASSIFY_SUM := results/tables/difference_classification_summary.txt
RESULTS_DIR  := results/tables

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
		--output $(TSS_V29)
	@echo "=== Parsing GENCODE v45 ==="
	$(PYTHON) src/02_parse_gencode.py \
		--gtf data/raw/gencode.v45.annotation.gtf.gz \
		--output $(TSS_V45)
	@echo "=== Parsing GWAS Catalog ==="
	$(PYTHON) src/03_parse_gwas.py \
		--input $(GWAS_RAW) \
		--output $(GWAS_ELIGIBLE)

map:
	$(MAKE) check
	@echo "=== Mapping variants ==="
	$(PYTHON) src/04_map_variants.py \
		--variants $(GWAS_ELIGIBLE) \
		--tss-v29 $(TSS_V29) \
		--tss-v45 $(TSS_V45) \
		--window-kb 500 \
		--output $(COMPARISON)

analyze: check bootstrap sensitivity classify
	@echo "=== Analysis complete ==="

bootstrap:
	$(MAKE) check
	@echo "=== Bootstrap CI ==="
	$(PYTHON) src/05_bootstrap.py \
		--comparison $(COMPARISON) \
		--output $(BOOTSTRAP_OUT)

sensitivity:
	$(MAKE) check
	@echo "=== Sensitivity analysis ==="
	$(PYTHON) src/06_sensitivity.py \
		--variants $(GWAS_ELIGIBLE) \
		--tss-v29 $(TSS_V29) \
		--tss-v45 $(TSS_V45) \
		--output $(SENSITIVITY)

classify:
	$(MAKE) check
	@echo "=== Classifying differences ==="
	$(PYTHON) src/07_classify_differences.py \
		--comparison $(COMPARISON) \
		--tss-v29 $(TSS_V29) \
		--tss-v45 $(TSS_V45) \
		--output $(CLASSIFY_OUT) \
		--summary $(CLASSIFY_SUM)

concentration:
	$(MAKE) check
	@echo "=== Concentration analysis ==="
	$(PYTHON) src/08_concentration_analysis.py \
		--variants $(GWAS_ELIGIBLE) \
		--tss-v29 $(TSS_V29) \
		--tss-v45 $(TSS_V45) \
		--output-dir $(RESULTS_DIR)

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
