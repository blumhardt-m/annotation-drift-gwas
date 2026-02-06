#!/usr/bin/env python3
"""
Classify gene differences: naming drift vs true mapping drift.

Categories:
1. SAME_GENE_DIFF_SYMBOL: Same Ensembl ID, different symbol (naming drift)
2. SAME_LOCUS_OVERLAP: Different IDs but genes overlap genomically (model change)
3. TRUE_DIFFERENT: Genuinely different genomic genes (real mapping drift)

Usage:
    python 07_classify_differences.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --output results/tables/difference_classification.tsv
"""

import csv
import argparse
from collections import defaultdict


def load_tss_with_details(path):
    """Load TSS table into dict: gene_id -> {symbol, chrom, tss, strand}"""
    genes = {}
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            genes[row['ensembl_id']] = {
                'symbol': row['gene_symbol'],
                'chrom': row['chromosome'],
                'tss': int(row['tss']),
                'strand': row.get('strand', '.')
            }
    return genes


def classify_difference(gene_v29, symbol_v29, gene_v45, symbol_v45, 
                        tss_v29_lookup, tss_v45_lookup, overlap_threshold=50000):
    """
    Classify a gene difference into categories.
    
    Returns: (category, details)
    """
    # Case 1: Same gene ID (shouldn't happen if differs=1, but check)
    if gene_v29 == gene_v45:
        if symbol_v29 != symbol_v45:
            return "SAME_ID_DIFF_SYMBOL", f"{symbol_v29}->{symbol_v45}"
        return "NO_DIFFERENCE", ""
    
    # Different gene IDs - check if they're close/overlapping
    info_v29 = tss_v29_lookup.get(gene_v29)
    info_v45 = tss_v45_lookup.get(gene_v45)
    
    if not info_v29 or not info_v45:
        return "LOOKUP_FAILED", f"v29:{gene_v29 in tss_v29_lookup} v45:{gene_v45 in tss_v45_lookup}"
    
    # Same chromosome?
    if info_v29['chrom'] != info_v45['chrom']:
        return "DIFFERENT_CHROM", f"{info_v29['chrom']} vs {info_v45['chrom']}"
    
    # How far apart are the TSSs?
    tss_distance = abs(info_v29['tss'] - info_v45['tss'])
    
    # Case 2: Very close TSSs suggest overlapping/related genes
    if tss_distance < overlap_threshold:
        return "NEARBY_GENES", f"TSS_dist={tss_distance}bp"
    
    # Case 3: Truly different genes
    return "TRUE_DIFFERENT", f"TSS_dist={tss_distance}bp"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--summary', default=None, help='Summary output file')
    args = parser.parse_args()
    
    print("Loading TSS tables...")
    tss_v29 = load_tss_with_details(args.tss_v29)
    tss_v45 = load_tss_with_details(args.tss_v45)
    
    print(f"  v29: {len(tss_v29)} genes")
    print(f"  v45: {len(tss_v45)} genes")
    
    print("Classifying differences...")
    
    categories = defaultdict(int)
    results = []
    
    with open(args.comparison, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['differs'] != '1':
                continue
            
            category, details = classify_difference(
                row['gene_v29'], row['symbol_v29'],
                row['gene_v45'], row['symbol_v45'],
                tss_v29, tss_v45
            )
            
            categories[category] += 1
            
            results.append({
                'snp': row['snp'],
                'chrom': row['chrom'],
                'pos': row['pos'],
                'gene_v29': row['gene_v29'],
                'symbol_v29': row['symbol_v29'],
                'gene_v45': row['gene_v45'],
                'symbol_v45': row['symbol_v45'],
                'category': category,
                'details': details
            })
    
    # Write detailed results
    with open(args.output, 'w', newline='') as f:
        fieldnames = ['snp', 'chrom', 'pos', 'gene_v29', 'symbol_v29', 
                      'gene_v45', 'symbol_v45', 'category', 'details']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    # Print summary
    total = sum(categories.values())
    print(f"\n{'='*60}")
    print(f"DIFFERENCE CLASSIFICATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total differences classified: {total:,}")
    print(f"\nBreakdown:")
    print(f"{'-'*60}")
    
    for cat in sorted(categories.keys()):
        count = categories[cat]
        pct = 100 * count / total if total > 0 else 0
        print(f"  {cat:25s}: {count:8,} ({pct:5.1f}%)")
    
    print(f"{'-'*60}")
    
    # Key interpretation
    naming_drift = categories.get('SAME_ID_DIFF_SYMBOL', 0)
    nearby = categories.get('NEARBY_GENES', 0)
    true_diff = categories.get('TRUE_DIFFERENT', 0)
    diff_chrom = categories.get('DIFFERENT_CHROM', 0)
    
    real_mapping_drift = true_diff + diff_chrom
    ambiguous = nearby
    
    print(f"\nINTERPRETATION:")
    print(f"  Naming drift only:     {naming_drift:,} ({100*naming_drift/total:.1f}%)")
    print(f"  Nearby genes (<50kb):  {nearby:,} ({100*nearby/total:.1f}%) - ambiguous")
    print(f"  True different genes:  {real_mapping_drift:,} ({100*real_mapping_drift/total:.1f}%) - REAL DRIFT")
    
    # Write summary file
    summary_path = args.summary or args.output.replace('.tsv', '_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("Difference Classification Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total differences: {total:,}\n\n")
        for cat in sorted(categories.keys()):
            count = categories[cat]
            pct = 100 * count / total if total > 0 else 0
            f.write(f"{cat}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        f.write(f"Naming drift: {naming_drift:,} ({100*naming_drift/total:.1f}%)\n")
        f.write(f"Nearby genes: {nearby:,} ({100*nearby/total:.1f}%)\n")
        f.write(f"True different: {real_mapping_drift:,} ({100*real_mapping_drift/total:.1f}%)\n")
    
    print(f"\nOutput: {args.output}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
