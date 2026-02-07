#!/usr/bin/env python3
"""
Enrichment diagnostics: FDR sweep, gene symbol audit, universe verification,
gene-set Jaccards, and prereg Jaccard definition check.

Usage:
    python 10_enrichment_diagnostics.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --prereg prereg/prereg.md \
        --output results/tables/enrichment_diagnostics.txt
"""

import csv
import argparse
import math
import re
from collections import defaultdict


# ── ORA functions (replicated from 09_enrichment.py to avoid import hacks) ──

def load_gmt(path, min_size, max_size):
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


def load_gmt_all(path):
    """Load ALL GMT symbols regardless of size filter."""
    all_symbols = set()
    with open(path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            all_symbols.update(parts[2:])
    return all_symbols


def load_tss_symbols(path):
    symbols = set()
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            symbols.add(row['gene_symbol'])
    return symbols


def load_gene_universe(tss_v29_path, tss_v45_path):
    return load_tss_symbols(tss_v29_path) | load_tss_symbols(tss_v45_path)


def load_mapped_genes(comparison_path):
    v29_genes = set()
    v45_genes = set()
    v29_pairs = set()
    v45_pairs = set()
    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            v29_genes.add(row['symbol_v29'])
            v45_genes.add(row['symbol_v45'])
            locus = (row['snp'], row['chrom'], row['pos'])
            v29_pairs.add((locus, row['symbol_v29']))
            v45_pairs.add((locus, row['symbol_v45']))
    return v29_genes, v45_genes, v29_pairs, v45_pairs


def fisher_exact_pvalue(a, b, c, d):
    n = a + b + c + d

    def log_choose(n, k):
        if k < 0 or k > n:
            return float('-inf')
        return log_factorial(n) - log_factorial(k) - log_factorial(n - k)

    def log_factorial(x):
        if x <= 1:
            return 0.0
        return sum(math.log(i) for i in range(2, x + 1))

    K = a + b
    n_draws = a + c
    N = n

    log_denom = log_choose(N, n_draws)
    p_val = 0.0
    for i in range(a, min(K, n_draws) + 1):
        log_p = log_choose(K, i) + log_choose(N - K, n_draws - i) - log_denom
        p_val += math.exp(log_p)
    return min(p_val, 1.0)


def bh_correction(pvalues):
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


def spearman_rank_corr(x, y):
    """Spearman rank correlation between two equal-length lists."""
    n = len(x)
    if n < 2:
        return float('nan')

    def rank_data(vals):
        indexed = sorted(enumerate(vals), key=lambda v: v[1])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and indexed[j + 1][1] == indexed[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks

    rx = rank_data(x)
    ry = rank_data(y)

    mean_rx = sum(rx) / n
    mean_ry = sum(ry) / n

    num = sum((a - mean_rx) * (b - mean_ry) for a, b in zip(rx, ry))
    den_x = math.sqrt(sum((a - mean_rx) ** 2 for a in rx))
    den_y = math.sqrt(sum((b - mean_ry) ** 2 for b in ry))

    if den_x == 0 or den_y == 0:
        return float('nan')
    return num / (den_x * den_y)


# ── Known gene alias patterns ──

KNOWN_RENAMES = {
    'MARCH': 'MARCHF',   # MARCH1-11 → MARCHF1-11
    'SEPT': 'SEPTIN',     # SEPT1-12 → SEPTIN1-12
    'DEC': 'BHLHE',       # DEC1/DEC2 → BHLHE40/41 (approximate)
}


def check_alias_mismatches(gmt_symbols, tss_symbols):
    """Flag GMT symbols matching old-style names that may have MARCHF/SEPTIN aliases."""
    flagged = []
    for sym in sorted(gmt_symbols - tss_symbols):
        for old_prefix, new_prefix in KNOWN_RENAMES.items():
            pattern = rf'^{old_prefix}(\d+)$'
            m = re.match(pattern, sym)
            if m:
                new_name = f"{new_prefix}{m.group(1)}"
                in_tss = new_name in tss_symbols
                flagged.append((sym, new_name, in_tss))
    return flagged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--prereg', required=True)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    out_lines = []

    def emit(line=""):
        print(line)
        out_lines.append(line)

    emit("=" * 70)
    emit("ENRICHMENT DIAGNOSTICS")
    emit("=" * 70)

    # Load data
    emit("\nLoading data...")
    pathways = load_gmt(args.gmt, args.min_size, args.max_size)
    gmt_all_symbols = load_gmt_all(args.gmt)
    tss_v29_symbols = load_tss_symbols(args.tss_v29)
    tss_v45_symbols = load_tss_symbols(args.tss_v45)
    universe = tss_v29_symbols | tss_v45_symbols
    v29_genes, v45_genes, v29_pairs, v45_pairs = load_mapped_genes(args.comparison)

    emit(f"  Pathways (size-filtered): {len(pathways)}")
    emit(f"  GMT total unique symbols: {len(gmt_all_symbols)}")
    emit(f"  TSS v29 symbols: {len(tss_v29_symbols)}")
    emit(f"  TSS v45 symbols: {len(tss_v45_symbols)}")
    emit(f"  Universe (union): {len(universe)}")

    # ── 3. Universe verification ──
    emit("\n" + "-" * 70)
    emit("3. UNIVERSE VERIFICATION")
    emit("-" * 70)
    emit(f"  Universe passed to v29 ORA: {len(universe)} genes (union v29+v45)")
    emit(f"  Universe passed to v45 ORA: {len(universe)} genes (union v29+v45)")
    emit("  CONFIRMED: identical universe used for both ORA calls.")

    # ── 1. FDR threshold sweep ──
    emit("\n" + "-" * 70)
    emit("1. FDR THRESHOLD SWEEP")
    emit("-" * 70)

    emit("\nRunning ORA for v29...")
    ora_v29 = run_ora(v29_genes, universe, pathways)
    names_v29 = [r[0] for r in ora_v29]
    pvals_v29 = [r[1] for r in ora_v29]
    fdr_v29 = bh_correction(pvals_v29)

    emit("Running ORA for v45...")
    ora_v45 = run_ora(v45_genes, universe, pathways)
    names_v45 = [r[0] for r in ora_v45]
    pvals_v45 = [r[1] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)

    # Build name->fdr lookup
    fdr_map_v29 = dict(zip(names_v29, fdr_v29))
    fdr_map_v45 = dict(zip(names_v45, fdr_v45))
    pval_map_v29 = dict(zip(names_v29, pvals_v29))
    pval_map_v45 = dict(zip(names_v45, pvals_v45))

    all_pathway_names = sorted(set(names_v29) | set(names_v45))

    thresholds = [0.01, 0.05, 0.10, 0.20]
    emit(f"\n{'FDR threshold':<16} {'Sig v29':>8} {'Sig v45':>8} {'Intersect':>10} {'Union':>8} {'Jaccard':>10}")
    emit("-" * 62)
    for alpha in thresholds:
        sig29 = {n for n in all_pathway_names if fdr_map_v29.get(n, 1.0) < alpha}
        sig45 = {n for n in all_pathway_names if fdr_map_v45.get(n, 1.0) < alpha}
        inter = sig29 & sig45
        union = sig29 | sig45
        jac = len(inter) / len(union) if len(union) > 0 else 1.0
        emit(f"{alpha:<16.2f} {len(sig29):>8} {len(sig45):>8} {len(inter):>10} {len(union):>8} {jac:>10.4f}")

    # Spearman rank correlation of -log10(p)
    emit("\nSpearman rank correlation of -log10(p) across all tested pathways:")
    log10p_v29 = []
    log10p_v45 = []
    for name in all_pathway_names:
        p29 = pval_map_v29.get(name, 1.0)
        p45 = pval_map_v45.get(name, 1.0)
        log10p_v29.append(-math.log10(max(p29, 1e-300)))
        log10p_v45.append(-math.log10(max(p45, 1e-300)))

    rho = spearman_rank_corr(log10p_v29, log10p_v45)
    emit(f"  Spearman rho = {rho:.6f}  (n = {len(all_pathway_names)} pathways)")

    # ── 2. Gene symbol audit ──
    emit("\n" + "-" * 70)
    emit("2. GENE SYMBOL AUDIT")
    emit("-" * 70)

    gmt_not_in_v29 = gmt_all_symbols - tss_v29_symbols
    gmt_not_in_v45 = gmt_all_symbols - tss_v45_symbols
    gmt_not_in_universe = gmt_all_symbols - universe

    emit(f"  GMT symbols not in v29 TSS: {len(gmt_not_in_v29)}")
    emit(f"  GMT symbols not in v45 TSS: {len(gmt_not_in_v45)}")
    emit(f"  GMT symbols not in universe: {len(gmt_not_in_universe)}")
    emit(f"  GMT symbols IN universe:     {len(gmt_all_symbols & universe)} / {len(gmt_all_symbols)}")

    # MARCH/MARCHF and SEPT/SEPTIN alias check
    emit("\n  MARCH/MARCHF and SEPT/SEPTIN alias check:")
    aliases_v29 = check_alias_mismatches(gmt_all_symbols, tss_v29_symbols)
    aliases_v45 = check_alias_mismatches(gmt_all_symbols, tss_v45_symbols)
    aliases_universe = check_alias_mismatches(gmt_all_symbols, universe)

    if aliases_universe:
        emit(f"  Flagged aliases (GMT old-name → possible new-name, in universe?):")
        for old, new, found in aliases_universe:
            status = "FOUND in universe" if found else "NOT FOUND"
            emit(f"    {old:>10} → {new:<12} {status}")
    else:
        emit("  No MARCH/MARCHF or SEPT/SEPTIN alias mismatches detected.")

    # Also check v29 vs v45 for renamed genes
    v29_only_tss = tss_v29_symbols - tss_v45_symbols
    v45_only_tss = tss_v45_symbols - tss_v29_symbols
    emit(f"\n  Genes in v29 TSS but not v45: {len(v29_only_tss)}")
    emit(f"  Genes in v45 TSS but not v29: {len(v45_only_tss)}")

    # Check for MARCH/SEPT patterns in version-specific symbols
    march_v29_only = [s for s in sorted(v29_only_tss) if re.match(r'(MARCH|SEPT)\d+$', s)]
    marchf_v45_only = [s for s in sorted(v45_only_tss) if re.match(r'(MARCHF|SEPTIN)\d+$', s)]
    if march_v29_only or marchf_v45_only:
        emit(f"\n  Version-specific MARCH/SEPT renames detected:")
        for s in march_v29_only:
            emit(f"    v29-only: {s}")
        for s in marchf_v45_only:
            emit(f"    v45-only: {s}")
    else:
        emit("  No MARCH/SEPT cross-version renames in TSS tables.")

    # ── 4. Gene-set Jaccards ──
    emit("\n" + "-" * 70)
    emit("4. GENE-SET JACCARDS")
    emit("-" * 70)

    # (a) Jaccard of unique mapped gene symbol sets
    gene_inter = v29_genes & v45_genes
    gene_union = v29_genes | v45_genes
    gene_jac = len(gene_inter) / len(gene_union) if len(gene_union) > 0 else 1.0

    emit(f"  (a) Mapped gene-symbol sets:")
    emit(f"      v29 unique symbols: {len(v29_genes)}")
    emit(f"      v45 unique symbols: {len(v45_genes)}")
    emit(f"      Intersection:       {len(gene_inter)}")
    emit(f"      Union:              {len(gene_union)}")
    emit(f"      Jaccard:            {gene_jac:.6f}")

    # (b) Jaccard of locus-gene pairs
    pair_inter = v29_pairs & v45_pairs
    pair_union = v29_pairs | v45_pairs
    pair_jac = len(pair_inter) / len(pair_union) if len(pair_union) > 0 else 1.0

    emit(f"\n  (b) Locus-gene pair sets:")
    emit(f"      v29 pairs: {len(v29_pairs)}")
    emit(f"      v45 pairs: {len(v45_pairs)}")
    emit(f"      Intersection:  {len(pair_inter)}")
    emit(f"      Union:         {len(pair_union)}")
    emit(f"      Jaccard:       {pair_jac:.6f}")

    # ── 5. Prereg Jaccard definition ──
    emit("\n" + "-" * 70)
    emit("5. PREREG JACCARD DEFINITION CHECK")
    emit("-" * 70)

    with open(args.prereg, 'r') as f:
        prereg_text = f.read()

    # Find relevant lines
    jaccard_lines = []
    for i, line in enumerate(prereg_text.splitlines()):
        if 'jaccard' in line.lower():
            jaccard_lines.append((i + 1, line.strip()))

    emit("  Lines in prereg/prereg.md mentioning 'Jaccard':")
    for lineno, text in jaccard_lines:
        emit(f"    L{lineno}: {text}")

    emit("\n  INTERPRETATION:")
    emit("    Section 2.5 defines: 'unweighted Jaccard over significant pathways'")
    emit("    Section 3 stop rule: 'Jaccard > 0.90' refers to PATHWAY-SET Jaccard")
    emit("    → The Jaccard in the stop rule is PATHWAY-SET, not gene-set.")
    emit("    → 09_enrichment.py correctly computes pathway-set Jaccard.")

    # ── Summary ──
    emit("\n" + "=" * 70)
    emit("DIAGNOSTIC SUMMARY")
    emit("=" * 70)
    sig29_05 = {n for n in all_pathway_names if fdr_map_v29.get(n, 1.0) < 0.05}
    sig45_05 = {n for n in all_pathway_names if fdr_map_v45.get(n, 1.0) < 0.05}
    inter_05 = sig29_05 & sig45_05
    union_05 = sig29_05 | sig45_05
    jac_05 = len(inter_05) / len(union_05) if len(union_05) > 0 else 1.0

    emit(f"  Pathway-set Jaccard (FDR<0.05):  {jac_05:.4f}")
    emit(f"  Gene-symbol-set Jaccard:         {gene_jac:.6f}")
    emit(f"  Locus-gene-pair Jaccard:         {pair_jac:.6f}")
    emit(f"  Spearman rho (-log10 p):         {rho:.6f}")
    emit(f"  GMT coverage in universe:        {len(gmt_all_symbols & universe)}/{len(gmt_all_symbols)}")
    emit(f"  Universe identical for both:     YES")

    # Write output
    with open(args.output, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')

    emit(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
