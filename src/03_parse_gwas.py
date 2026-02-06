#!/usr/bin/env python3
"""Parse GWAS Catalog to eligible associations."""

import csv
import argparse

def parse_gwas_catalog(input_path, output_path, max_pvalue=5e-8):
    eligible = []
    stats = {'total': 0, 'no_rsid': 0, 'no_chrom': 0, 'non_auto': 0, 'no_pos': 0, 'pval_fail': 0}
    
    print(f"Parsing {input_path}...")
    with open(input_path, 'r', encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            stats['total'] += 1
            snp = row.get('SNPS', '').strip()
            chrom = row.get('CHR_ID', '').strip()
            pos = row.get('CHR_POS', '').strip()
            pval = row.get('P-VALUE', '').strip()
            
            if not snp.startswith('rs'): stats['no_rsid'] += 1; continue
            if not chrom: stats['no_chrom'] += 1; continue
            if chrom not in [str(i) for i in range(1,23)]: stats['non_auto'] += 1; continue
            if not pos.isdigit(): stats['no_pos'] += 1; continue
            try:
                if float(pval) > max_pvalue: stats['pval_fail'] += 1; continue
            except: stats['pval_fail'] += 1; continue
            
            eligible.append({'SNPS': snp, 'CHR_ID': chrom, 'CHR_POS': pos, 'P-VALUE': pval,
                           'DISEASE/TRAIT': row.get('DISEASE/TRAIT',''), 'STUDY': row.get('STUDY ACCESSION','')})
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['SNPS','CHR_ID','CHR_POS','P-VALUE','DISEASE/TRAIT','STUDY'], delimiter='\t')
        writer.writeheader()
        writer.writerows(eligible)
    
    print(f"  Total: {stats['total']:,} → Eligible: {len(eligible):,}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-pvalue', type=float, default=5e-8)
    args = parser.parse_args()
    parse_gwas_catalog(args.input, args.output, args.max_pvalue)
