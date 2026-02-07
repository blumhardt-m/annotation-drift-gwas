#!/usr/bin/env python3
"""
Pathway enrichment analysis: ORA (Fisher exact) on v29 vs v45 gene sets,
then Jaccard similarity of significant pathway sets.

Usage:
    python 09_enrichment.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --fdr-alpha 0.05 \
        --min-size 15 --max-size 500 \
        --output results/tables/enrichment_jaccard_by_window.tsv
"""

import csv
import argparse
import math
from collections import defaultdict


def load_gmt(path, min_size, max_size):
    """Load GMT file. Returns dict: pathway_name -> set of gene symbols."""
    pathways = {}
    with open(path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = set(parts[2:])
            if min_size <= len(genes) <= max_size:
                pathways[name] = genes
    return pathways


def load_gene_universe(tss_v29_path, tss_v45_path):
    """Union of protein-coding gene symbols from both versions."""
    symbols = set()
    for path in [tss_v29_path, tss_v45_path]:
        with open(path, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                symbols.add(row['gene_symbol'])
    return symbols


def load_mapped_genes(comparison_path):
    """Extract v29 and v45 gene symbol sets from comparison table."""
    v29_genes = set()
    v45_genes = set()
    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            v29_genes.add(row['symbol_v29'])
            v45_genes.add(row['symbol_v45'])
    return v29_genes, v45_genes


def fisher_exact_pvalue(a, b, c, d):
    """One-tailed Fisher exact test p-value (overrepresentation).
    Contingency table:
        a = in_pathway & in_geneset
        b = in_pathway & not_in_geneset
        c = not_in_pathway & in_geneset
        d = not_in_pathway & not_in_geneset
    """
    n = a + b + c + d
    # Hypergeometric: P(X >= a) where X ~ Hypergeom(n, a+b, a+c)
    # Use log-space for numerical stability
    def log_choose(n, k):
        if k < 0 or k > n:
            return float('-inf')
        return log_factorial(n) - log_factorial(k) - log_factorial(n - k)

    def log_factorial(x):
        if x <= 1:
            return 0.0
        return sum(math.log(i) for i in range(2, x + 1))

    K = a + b      # pathway size in universe
    n_draws = a + c  # gene set size
    N = n          # universe size

    # P(X >= a)
    log_denom = log_choose(N, n_draws)
    p_val = 0.0
    for i in range(a, min(K, n_draws) + 1):
        log_p = log_choose(K, i) + log_choose(N - K, n_draws - i) - log_denom
        p_val += math.exp(log_p)

    return min(p_val, 1.0)


def bh_correction(pvalues):
    """Benjamini-Hochberg FDR correction. Returns list of adjusted p-values."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    cummin = 1.0
    for rank_from_end, (orig_idx, pval) in enumerate(reversed(indexed)):
        rank = n - rank_from_end
        adj = pval * n / rank
        cummin = min(cummin, adj)
        adjusted[orig_idx] = min(cummin, 1.0)
    return adjusted


def run_ora(gene_set, universe, pathways):
    """Run ORA for a gene set against pathways. Returns list of (name, pval, overlap_count)."""
    gene_set_in_universe = gene_set & universe
    n_geneset = len(gene_set_in_universe)
    N = len(universe)
    results = []

    for name, pw_genes in pathways.items():
        pw_in_universe = pw_genes & universe
        K = len(pw_in_universe)
        if K == 0:
            continue
        overlap = gene_set_in_universe & pw_in_universe
        a = len(overlap)
        b = K - a
        c = n_geneset - a
        d = N - K - c
        pval = fisher_exact_pvalue(a, b, c, d)
        results.append((name, pval, a))

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--fdr-alpha', type=float, default=0.05)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    print("Loading GMT...")
    pathways = load_gmt(args.gmt, args.min_size, args.max_size)
    print(f"  {len(pathways)} pathways after size filter [{args.min_size}, {args.max_size}]")

    print("Loading gene universe (union v29 + v45)...")
    universe = load_gene_universe(args.tss_v29, args.tss_v45)
    print(f"  {len(universe)} genes in universe")

    print("Loading mapped gene sets...")
    v29_genes, v45_genes = load_mapped_genes(args.comparison)
    print(f"  v29: {len(v29_genes)} unique symbols, v45: {len(v45_genes)} unique symbols")

    # ORA for v29
    print("Running ORA for v29 gene set...")
    ora_v29 = run_ora(v29_genes, universe, pathways)
    names_v29 = [r[0] for r in ora_v29]
    pvals_v29 = [r[1] for r in ora_v29]
    fdr_v29 = bh_correction(pvals_v29)

    # ORA for v45
    print("Running ORA for v45 gene set...")
    ora_v45 = run_ora(v45_genes, universe, pathways)
    names_v45 = [r[0] for r in ora_v45]
    pvals_v45 = [r[1] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)

    # Significant pathway sets
    sig_v29 = {name for name, fdr in zip(names_v29, fdr_v29) if fdr < args.fdr_alpha}
    sig_v45 = {name for name, fdr in zip(names_v45, fdr_v45) if fdr < args.fdr_alpha}

    # Jaccard
    intersection = sig_v29 & sig_v45
    union = sig_v29 | sig_v45
    jaccard = len(intersection) / len(union) if len(union) > 0 else 1.0

    print(f"\n{'='*60}")
    print("ENRICHMENT RESULTS")
    print(f"{'='*60}")
    print(f"  Significant pathways (v29):  {len(sig_v29)}")
    print(f"  Significant pathways (v45):  {len(sig_v45)}")
    print(f"  Intersection:                {len(intersection)}")
    print(f"  Union:                       {len(union)}")
    print(f"  Jaccard similarity:          {jaccard:.4f}")
    print(f"{'='*60}")

    # v29-only and v45-only pathways
    v29_only = sig_v29 - sig_v45
    v45_only = sig_v45 - sig_v29

    if v29_only:
        print(f"\n  v29-only pathways ({len(v29_only)}):")
        for p in sorted(v29_only)[:10]:
            print(f"    - {p}")
        if len(v29_only) > 10:
            print(f"    ... and {len(v29_only) - 10} more")

    if v45_only:
        print(f"\n  v45-only pathways ({len(v45_only)}):")
        for p in sorted(v45_only)[:10]:
            print(f"    - {p}")
        if len(v45_only) > 10:
            print(f"    ... and {len(v45_only) - 10} more")

    # Write output
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'pathway', 'pval_v29', 'fdr_v29', 'pval_v45', 'fdr_v45',
            'sig_v29', 'sig_v45', 'overlap_v29', 'overlap_v45'
        ])

        # Build lookup for v45 results
        v45_lookup = {}
        for i, (name, pval, overlap) in enumerate(ora_v45):
            v45_lookup[name] = (pval, fdr_v45[i], overlap)

        for i, (name, pval, overlap) in enumerate(ora_v29):
            fdr29 = fdr_v29[i]
            if name in v45_lookup:
                pval45, fdr45, overlap45 = v45_lookup[name]
            else:
                pval45, fdr45, overlap45 = 1.0, 1.0, 0
            writer.writerow([
                name,
                f"{pval:.6e}", f"{fdr29:.6e}",
                f"{pval45:.6e}", f"{fdr45:.6e}",
                1 if fdr29 < args.fdr_alpha else 0,
                1 if fdr45 < args.fdr_alpha else 0,
                overlap, overlap45
            ])

    print(f"\nOutput: {args.output}")

    # Summary line for stop-rule check
    print(f"\nSTOP-RULE CHECK: I < 5% AND Jaccard = {jaccard:.4f}", end="")
    if jaccard > 0.90:
        print(" > 0.90 → Publish as negative result")
    else:
        print(f" <= 0.90 → Annotation drift has pathway-level impact")


if __name__ == "__main__":
    main()
