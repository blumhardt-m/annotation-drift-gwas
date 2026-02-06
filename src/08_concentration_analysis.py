#!/usr/bin/env python3
"""
Concentration analysis: where does instability concentrate?

Analyses:
1. Instability vs distance to assigned TSS
2. Instability vs margin (gap between 1st and 2nd nearest gene)
3. Instability by chromosome

Usage:
    python 08_concentration_analysis.py \
        --variants data/processed/gwas_eligible.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --output-dir results/tables
"""

import csv
import argparse
from bisect import bisect_left
from collections import defaultdict
import os


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


def find_nearest_two_genes(chrom, pos, tss_by_chrom, window_bp=500000):
    """Find nearest and second-nearest gene. Returns [(gene_id, dist), ...] sorted by dist."""
    if chrom not in tss_by_chrom:
        return []
    
    genes = tss_by_chrom[chrom]
    positions = [g[0] for g in genes]
    idx = bisect_left(positions, pos)
    
    candidates = []
    # Check genes in a window around the position
    for i in range(max(0, idx-5), min(len(genes), idx+6)):
        dist = abs(genes[i][0] - pos)
        if dist <= window_bp:
            candidates.append((genes[i][1], dist))
    
    # Sort by distance, then by gene_id for determinism
    candidates.sort(key=lambda x: (x[1], x[0]))
    return candidates[:2]  # Return top 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("Loading TSS tables...")
    tss_v29 = load_tss_table(args.tss_v29)
    tss_v45 = load_tss_table(args.tss_v45)
    
    # Distance bins
    dist_bins = [(0, 10000), (10000, 50000), (50000, 100000), 
                 (100000, 250000), (250000, 500000)]
    dist_stats = {b: {'total': 0, 'differ': 0} for b in dist_bins}
    
    # Margin bins (distance from 1st to 2nd nearest)
    margin_bins = [(0, 1000), (1000, 5000), (5000, 10000), 
                   (10000, 50000), (50000, 100000), (100000, float('inf'))]
    margin_stats = {b: {'total': 0, 'differ': 0} for b in margin_bins}
    
    # Chromosome stats
    chrom_stats = defaultdict(lambda: {'total': 0, 'differ': 0})
    
    print("Analyzing variants...")
    n_processed = 0
    
    with open(args.variants, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row['CHR_ID']
            pos = int(row['CHR_POS'])
            
            # Get nearest genes in both versions
            top2_v29 = find_nearest_two_genes(chrom, pos, tss_v29)
            top2_v45 = find_nearest_two_genes(chrom, pos, tss_v45)
            
            if not top2_v29 or not top2_v45:
                continue
            
            nearest_v29, dist_v29 = top2_v29[0]
            nearest_v45, dist_v45 = top2_v45[0]
            differs = 1 if nearest_v29 != nearest_v45 else 0
            
            # Use v45 distance as primary
            dist = dist_v45
            
            # Margin: if there's a second gene, compute gap
            if len(top2_v45) >= 2:
                margin = top2_v45[1][1] - top2_v45[0][1]
            else:
                margin = float('inf')
            
            # Bin by distance
            for (lo, hi) in dist_bins:
                if lo <= dist < hi:
                    dist_stats[(lo, hi)]['total'] += 1
                    dist_stats[(lo, hi)]['differ'] += differs
                    break
            
            # Bin by margin
            for (lo, hi) in margin_bins:
                if lo <= margin < hi:
                    margin_stats[(lo, hi)]['total'] += 1
                    margin_stats[(lo, hi)]['differ'] += differs
                    break
            
            # Chromosome
            chrom_stats[chrom]['total'] += 1
            chrom_stats[chrom]['differ'] += differs
            
            n_processed += 1
            if n_processed % 100000 == 0:
                print(f"  Processed {n_processed:,}...")
    
    print(f"  Total processed: {n_processed:,}")
    
    # Output: Distance analysis
    dist_file = os.path.join(args.output_dir, 'instability_by_distance.tsv')
    with open(dist_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['dist_min', 'dist_max', 'n_total', 'n_differ', 'instability_rate'])
        for (lo, hi) in dist_bins:
            s = dist_stats[(lo, hi)]
            rate = s['differ'] / s['total'] if s['total'] > 0 else 0
            writer.writerow([lo, hi, s['total'], s['differ'], f"{rate:.4f}"])
    
    # Output: Margin analysis
    margin_file = os.path.join(args.output_dir, 'instability_by_margin.tsv')
    with open(margin_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['margin_min', 'margin_max', 'n_total', 'n_differ', 'instability_rate'])
        for (lo, hi) in margin_bins:
            s = margin_stats[(lo, hi)]
            rate = s['differ'] / s['total'] if s['total'] > 0 else 0
            hi_str = 'inf' if hi == float('inf') else hi
            writer.writerow([lo, hi_str, s['total'], s['differ'], f"{rate:.4f}"])
    
    # Output: Chromosome analysis
    chrom_file = os.path.join(args.output_dir, 'instability_by_chromosome.tsv')
    with open(chrom_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['chromosome', 'n_total', 'n_differ', 'instability_rate'])
        for chrom in sorted(chrom_stats.keys(), key=lambda x: int(x) if x.isdigit() else 99):
            s = chrom_stats[chrom]
            rate = s['differ'] / s['total'] if s['total'] > 0 else 0
            writer.writerow([chrom, s['total'], s['differ'], f"{rate:.4f}"])
    
    # Print summary
    print(f"\n{'='*70}")
    print("CONCENTRATION ANALYSIS RESULTS")
    print(f"{'='*70}")
    
    print("\n1. INSTABILITY BY DISTANCE TO TSS")
    print("-" * 50)
    print(f"{'Distance (kb)':<20} {'N':>10} {'Differ':>10} {'Rate':>10}")
    print("-" * 50)
    for (lo, hi) in dist_bins:
        s = dist_stats[(lo, hi)]
        rate = s['differ'] / s['total'] if s['total'] > 0 else 0
        label = f"{lo//1000}-{hi//1000} kb"
        print(f"{label:<20} {s['total']:>10,} {s['differ']:>10,} {rate:>10.2%}")
    
    print("\n2. INSTABILITY BY MARGIN TO 2ND NEAREST GENE")
    print("-" * 50)
    print(f"{'Margin (bp)':<20} {'N':>10} {'Differ':>10} {'Rate':>10}")
    print("-" * 50)
    for (lo, hi) in margin_bins:
        s = margin_stats[(lo, hi)]
        rate = s['differ'] / s['total'] if s['total'] > 0 else 0
        hi_str = '∞' if hi == float('inf') else f"{hi//1000}k"
        lo_str = f"{lo}" if lo < 1000 else f"{lo//1000}k"
        label = f"{lo_str}-{hi_str}"
        print(f"{label:<20} {s['total']:>10,} {s['differ']:>10,} {rate:>10.2%}")
    
    print("\n3. TOP 5 CHROMOSOMES BY INSTABILITY RATE")
    print("-" * 50)
    chrom_rates = [(c, s['differ']/s['total'] if s['total'] > 0 else 0, s['total']) 
                   for c, s in chrom_stats.items()]
    chrom_rates.sort(key=lambda x: -x[1])
    for chrom, rate, n in chrom_rates[:5]:
        print(f"  chr{chrom}: {rate:.2%} (n={n:,})")
    
    print(f"\nOutputs saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
