#!/usr/bin/env python3
"""Map GWAS variants to nearest gene by TSS."""

import csv
import argparse
from bisect import bisect_left
from collections import defaultdict

def load_tss_table(path):
    """Load TSS table, return dict: chrom -> sorted list of (pos, gene_id, symbol)"""
    tss_by_chrom = defaultdict(list)
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            tss_by_chrom[row['chromosome']].append(
                (int(row['tss']), row['ensembl_id'], row['gene_symbol'])
            )
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort(key=lambda x: x[0])
    return tss_by_chrom

def find_nearest_gene(chrom, pos, tss_by_chrom, window_bp):
    """Find nearest gene within window. Returns (gene_id, symbol, distance) or None."""
    if chrom not in tss_by_chrom:
        return None
    
    genes = tss_by_chrom[chrom]
    positions = [g[0] for g in genes]
    idx = bisect_left(positions, pos)
    
    candidates = []
    for i in range(max(0, idx-1), min(len(genes), idx+2)):
        dist = abs(genes[i][0] - pos)
        if dist <= window_bp:
            candidates.append((dist, genes[i][1], genes[i][2]))
    
    if not candidates:
        return None
    
    candidates.sort(key=lambda x: (x[0], x[1]))  # by distance, then by ID
    return (candidates[0][1], candidates[0][2], candidates[0][0])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--window-kb', type=int, default=500)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    window_bp = args.window_kb * 1000
    
    print("Loading TSS tables...")
    tss_v29 = load_tss_table(args.tss_v29)
    tss_v45 = load_tss_table(args.tss_v45)
    
    print("Mapping variants...")
    results = []
    n_both, n_differ = 0, 0
    
    with open(args.variants, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            chrom = row['CHR_ID']
            pos = int(row['CHR_POS'])
            snp = row['SNPS']
            
            g29 = find_nearest_gene(chrom, pos, tss_v29, window_bp)
            g45 = find_nearest_gene(chrom, pos, tss_v45, window_bp)
            
            if g29 and g45:
                n_both += 1
                differs = g29[0] != g45[0]
                if differs:
                    n_differ += 1
                
                results.append({
                    'snp': snp, 'chrom': chrom, 'pos': pos,
                    'gene_v29': g29[0], 'symbol_v29': g29[1], 'dist_v29': g29[2],
                    'gene_v45': g45[0], 'symbol_v45': g45[1], 'dist_v45': g45[2],
                    'differs': 1 if differs else 0
                })
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['snp','chrom','pos','gene_v29','symbol_v29','dist_v29',
                                               'gene_v45','symbol_v45','dist_v45','differs'], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    rate = n_differ / n_both if n_both > 0 else 0
    print(f"\nResults:")
    print(f"  Mapped in both: {n_both:,}")
    print(f"  Gene differs:   {n_differ:,}")
    print(f"  Instability:    {rate:.2%}")
    print(f"\nOutput: {args.output}")

if __name__ == "__main__":
    main()
