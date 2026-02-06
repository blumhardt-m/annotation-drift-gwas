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
