#!/bin/bash
set -e
DATA_DIR="data/raw"
mkdir -p "$DATA_DIR"

echo "Downloading data..."

echo "[1/4] GWAS Catalog..."
wget -O "$DATA_DIR/gwas_catalog.tsv" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv"

echo "[2/4] GENCODE v29..."
wget -O "$DATA_DIR/gencode.v29.annotation.gtf.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"

echo "[3/4] GENCODE v45..."
wget -O "$DATA_DIR/gencode.v45.annotation.gtf.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"

echo "[4/4] GENCODE v44..."
wget -O "$DATA_DIR/gencode.v44.annotation.gtf.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

echo "Done!"
ls -lh "$DATA_DIR"
