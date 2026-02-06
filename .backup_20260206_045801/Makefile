.PHONY: download parse map analyze clean help

help:
	@echo "Targets: download, parse, map, analyze, clean"
	@echo "Run: make download && make parse && make map && make analyze"

download:
	./src/01_download_data.sh

parse:
	python src/02_parse_gencode.py --gtf data/raw/gencode.v29.annotation.gtf.gz --output data/processed/tss_tables/tss_v29.tsv
	python src/02_parse_gencode.py --gtf data/raw/gencode.v45.annotation.gtf.gz --output data/processed/tss_tables/tss_v45.tsv
	python src/03_parse_gwas.py --input data/raw/gwas_catalog.tsv --output data/processed/gwas_eligible.tsv

map:
	python src/04_map_variants.py --variants data/processed/gwas_eligible.tsv \
		--tss-v29 data/processed/tss_tables/tss_v29.tsv \
		--tss-v45 data/processed/tss_tables/tss_v45.tsv \
		--output data/processed/mappings/comparison.tsv

analyze:
	python src/05_bootstrap.py --comparison data/processed/mappings/comparison.tsv --output results/tables/bootstrap.txt
	python src/06_sensitivity.py --variants data/processed/gwas_eligible.tsv \
		--tss-v29 data/processed/tss_tables/tss_v29.tsv \
		--tss-v45 data/processed/tss_tables/tss_v45.tsv \
		--output results/tables/sensitivity.tsv

clean:
	rm -f data/processed/**/*.tsv results/tables/*
