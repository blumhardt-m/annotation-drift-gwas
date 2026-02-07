#!/usr/bin/env python3
"""
Ensembl-ID-based enrichment: use stable identifiers to separate true mapping
drift from symbol-namespace drift (MARCH/MARCHF, SEPT/SEPTIN).

Usage:
    python 12_ensembl_id_enrichment.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --output results/tables/enrichment_ensembl_id.txt
"""

import csv
import argparse
import math
import re
from collections import defaultdict


# ── HGNC alias table ────────────────────────────────────────────────────

def build_alias_table():
    aliases = {}
    for i in range(1, 12):
        aliases[f"MARCH{i}"] = f"MARCHF{i}"
        aliases[f"MARCHF{i}"] = f"MARCH{i}"
    for i in range(1, 15):
        aliases[f"SEPT{i}"] = f"SEPTIN{i}"
        aliases[f"SEPTIN{i}"] = f"SEPT{i}"
    return aliases


# ── Data loading ────────────────────────────────────────────────────────

def load_tss_tables(v29_path, v45_path):
    """Load both TSS tables. Returns:
    - symbol_to_ensembl: dict mapping symbol -> ensembl_id (union, alias-aware)
    - ensembl_universe: set of all ensembl_ids
    """
    alias_table = build_alias_table()
    symbol_to_ensembl = {}
    ensembl_universe = set()

    for path in [v29_path, v45_path]:
        with open(path, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                eid = row['ensembl_id']
                sym = row['gene_symbol']
                ensembl_universe.add(eid)
                symbol_to_ensembl[sym] = eid

    # For alias pairs, ensure both old and new symbols map to the same Ensembl ID.
    # Prefer the mapping that already exists; fill in the missing direction.
    for old_sym, new_sym in alias_table.items():
        if old_sym in symbol_to_ensembl and new_sym not in symbol_to_ensembl:
            symbol_to_ensembl[new_sym] = symbol_to_ensembl[old_sym]
        elif new_sym in symbol_to_ensembl and old_sym not in symbol_to_ensembl:
            symbol_to_ensembl[old_sym] = symbol_to_ensembl[new_sym]

    return symbol_to_ensembl, ensembl_universe


def load_gmt_as_ensembl(gmt_path, symbol_to_ensembl, min_size, max_size):
    """Load GMT, converting symbols to Ensembl IDs. Returns (pathways, stats)."""
    pathways = {}
    n_mapped = 0
    n_unmapped = 0
    unmapped_symbols = set()

    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            symbols = parts[2:]
            ensembl_genes = set()
            for sym in symbols:
                eid = symbol_to_ensembl.get(sym)
                if eid:
                    ensembl_genes.add(eid)
                    n_mapped += 1
                else:
                    n_unmapped += 1
                    unmapped_symbols.add(sym)
            if min_size <= len(ensembl_genes) <= max_size:
                pathways[name] = ensembl_genes

    stats = {
        'n_mapped': n_mapped,
        'n_unmapped': n_unmapped,
        'n_unmapped_unique': len(unmapped_symbols),
        'n_pathways': len(pathways),
    }
    return pathways, stats, unmapped_symbols


def load_comparison_ensembl(comparison_path):
    """Load comparison table. Returns v29/v45 Ensembl gene sets and locus-level data."""
    v29_genes = set()
    v45_genes = set()
    loci = []

    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            g29 = row['gene_v29']
            g45 = row['gene_v45']
            s29 = row['symbol_v29']
            s45 = row['symbol_v45']
            v29_genes.add(g29)
            v45_genes.add(g45)
            loci.append({
                'snp': row['snp'],
                'chrom': row['chrom'],
                'gene_v29': g29,
                'gene_v45': g45,
                'symbol_v29': s29,
                'symbol_v45': s45,
                'differs_ensembl': int(row['differs']),
            })

    return v29_genes, v45_genes, loci


# ── ORA machinery ───────────────────────────────────────────────────────

def fisher_exact_pvalue(a, b, c, d):
    n = a + b + c + d

    def log_choose(nn, k):
        if k < 0 or k > nn:
            return float('-inf')
        return log_factorial(nn) - log_factorial(k) - log_factorial(nn - k)

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


def fdr_sweep(ora_v29, ora_v45, thresholds):
    names_v29 = [r[0] for r in ora_v29]
    pvals_v29 = [r[1] for r in ora_v29]
    fdr_v29 = bh_correction(pvals_v29)

    names_v45 = [r[0] for r in ora_v45]
    pvals_v45 = [r[1] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)

    fdr_map_v29 = dict(zip(names_v29, fdr_v29))
    fdr_map_v45 = dict(zip(names_v45, fdr_v45))
    pval_map_v29 = dict(zip(names_v29, pvals_v29))
    pval_map_v45 = dict(zip(names_v45, pvals_v45))

    all_names = sorted(set(names_v29) | set(names_v45))

    results = []
    for alpha in thresholds:
        sig29 = {n for n in all_names if fdr_map_v29.get(n, 1.0) < alpha}
        sig45 = {n for n in all_names if fdr_map_v45.get(n, 1.0) < alpha}
        inter = sig29 & sig45
        union = sig29 | sig45
        jac = len(inter) / len(union) if len(union) > 0 else 1.0
        results.append((alpha, len(sig29), len(sig45), len(inter), len(union), jac))

    log10p_v29 = [-math.log10(max(pval_map_v29.get(n, 1.0), 1e-300)) for n in all_names]
    log10p_v45 = [-math.log10(max(pval_map_v45.get(n, 1.0), 1e-300)) for n in all_names]
    rho = spearman_rank_corr(log10p_v29, log10p_v45)

    return results, rho


# ── Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    out_lines = []

    def emit(line=""):
        print(line)
        out_lines.append(line)

    emit("=" * 70)
    emit("ENSEMBL-ID-BASED ENRICHMENT ANALYSIS")
    emit("=" * 70)

    # ── Load and prepare data ──
    emit("\nLoading TSS tables and building symbol → Ensembl mapping...")
    symbol_to_ensembl, ensembl_universe = load_tss_tables(args.tss_v29, args.tss_v45)
    emit(f"  Symbol→Ensembl mappings: {len(symbol_to_ensembl)}")
    emit(f"  Ensembl universe:        {len(ensembl_universe)}")

    emit("\nConverting GMT to Ensembl IDs...")
    pathways_eid, gmt_stats, unmapped_syms = load_gmt_as_ensembl(
        args.gmt, symbol_to_ensembl, args.min_size, args.max_size
    )
    emit(f"  GMT symbols mapped:      {gmt_stats['n_mapped']}")
    emit(f"  GMT symbols unmapped:    {gmt_stats['n_unmapped']} ({gmt_stats['n_unmapped_unique']} unique)")
    emit(f"  Pathways after filter:   {gmt_stats['n_pathways']}")

    emit("\nLoading comparison (Ensembl IDs)...")
    v29_genes, v45_genes, loci = load_comparison_ensembl(args.comparison)
    emit(f"  v29 unique Ensembl IDs:  {len(v29_genes)}")
    emit(f"  v45 unique Ensembl IDs:  {len(v45_genes)}")
    emit(f"  Total loci:              {len(loci)}")

    # ── 6. Decompose differing loci ──
    emit(f"\n{'='*70}")
    emit("6. LOCUS-LEVEL DRIFT DECOMPOSITION")
    emit(f"{'='*70}")

    n_total = len(loci)
    n_ensembl_differs = sum(1 for l in loci if l['differs_ensembl'] == 1)
    n_symbol_differs = sum(1 for l in loci if l['symbol_v29'] != l['symbol_v45'])
    n_symbol_only = sum(
        1 for l in loci
        if l['symbol_v29'] != l['symbol_v45'] and l['gene_v29'] == l['gene_v45']
    )
    n_true_drift = n_ensembl_differs  # same as differs_ensembl=1
    n_ensembl_same_symbol_same = sum(
        1 for l in loci
        if l['gene_v29'] == l['gene_v45'] and l['symbol_v29'] == l['symbol_v45']
    )

    emit(f"\n  Total loci mapped in both:                    {n_total:>10,}")
    emit(f"  Ensembl ID differs (true mapping drift):      {n_ensembl_differs:>10,} ({100*n_ensembl_differs/n_total:.2f}%)")
    emit(f"  Symbol differs:                               {n_symbol_differs:>10,} ({100*n_symbol_differs/n_total:.2f}%)")
    emit(f"  Symbol differs BUT same Ensembl ID:           {n_symbol_only:>10,}")
    emit(f"  True mapping drift (different Ensembl ID):    {n_true_drift:>10,}")

    emit(f"\n  Of the {n_symbol_differs:,} loci with symbol differences:")
    if n_symbol_differs > 0:
        emit(f"    Symbol-only (same gene, renamed):           {n_symbol_only:>10,} ({100*n_symbol_only/n_symbol_differs:.1f}%)")
        emit(f"    True mapping drift (different gene):        {n_symbol_differs - n_symbol_only:>10,} ({100*(n_symbol_differs - n_symbol_only)/n_symbol_differs:.1f}%)")

    # Show some examples of symbol-only differences
    if n_symbol_only > 0:
        emit(f"\n  Examples of symbol-only differences (same Ensembl ID):")
        count = 0
        for l in loci:
            if l['symbol_v29'] != l['symbol_v45'] and l['gene_v29'] == l['gene_v45']:
                emit(f"    {l['gene_v29']}: {l['symbol_v29']} → {l['symbol_v45']}  (locus {l['snp']})")
                count += 1
                if count >= 10:
                    emit(f"    ... and {n_symbol_only - 10} more")
                    break

    # ── Ensembl-ID-based ORA ──
    emit(f"\n{'='*70}")
    emit("ENSEMBL-ID-BASED ORA")
    emit(f"{'='*70}")

    emit("\nRunning ORA (v29 Ensembl IDs)...")
    ora_v29 = run_ora(v29_genes, ensembl_universe, pathways_eid)
    emit("Running ORA (v45 Ensembl IDs)...")
    ora_v45 = run_ora(v45_genes, ensembl_universe, pathways_eid)

    thresholds = [0.01, 0.05, 0.10, 0.20]
    eid_results, eid_rho = fdr_sweep(ora_v29, ora_v45, thresholds)

    # Symbol-based results (from diagnostics) for delta
    SYMBOL_JACCARDS = {0.01: 0.6347, 0.05: 0.6927, 0.10: 0.7170, 0.20: 0.8123}
    SYMBOL_RHO = 0.778094

    emit(f"\n{'FDR':<8} {'EID Jaccard':>12} {'Sym Jaccard':>12} {'Delta':>8}")
    emit("-" * 42)
    for alpha, n29, n45, nint, nunion, jac in eid_results:
        sym_jac = SYMBOL_JACCARDS.get(alpha, float('nan'))
        delta = jac - sym_jac
        emit(f"{alpha:<8.2f} {jac:>12.4f} {sym_jac:>12.4f} {delta:>+8.4f}")

    emit(f"\nSpearman rho (-log10 p):")
    emit(f"  Ensembl-ID-based: {eid_rho:.6f}")
    emit(f"  Symbol-based:     {SYMBOL_RHO:.6f}")
    emit(f"  Delta:            {eid_rho - SYMBOL_RHO:+.6f}")

    # Detailed FDR=0.05
    emit(f"\nDetailed at FDR < 0.05:")
    for alpha, n29, n45, nint, nunion, jac in eid_results:
        if alpha == 0.05:
            emit(f"  v29 sig: {n29}, v45 sig: {n45}, intersect: {nint}, union: {nunion}")
            emit(f"  Jaccard: {jac:.4f}")

    # ── Summary ──
    jac_05_eid = [r[5] for r in eid_results if r[0] == 0.05][0]
    jac_05_sym = SYMBOL_JACCARDS[0.05]

    emit(f"\n{'='*70}")
    emit("SUMMARY")
    emit(f"{'='*70}")
    emit(f"  Locus-level:")
    emit(f"    Instability (Ensembl ID): {n_ensembl_differs:,}/{n_total:,} = {100*n_ensembl_differs/n_total:.2f}%")
    emit(f"    Symbol-only renames:      {n_symbol_only:,} loci (not true drift)")
    emit(f"")
    emit(f"  Pathway-level (FDR < 0.05):")
    emit(f"    Symbol-based Jaccard:     {jac_05_sym:.4f}")
    emit(f"    Ensembl-ID-based Jaccard: {jac_05_eid:.4f}")
    emit(f"    Delta:                    {jac_05_eid - jac_05_sym:+.4f}")
    emit(f"")
    emit(f"  Spearman rho:")
    emit(f"    Symbol-based:             {SYMBOL_RHO:.4f}")
    emit(f"    Ensembl-ID-based:         {eid_rho:.4f}")
    emit(f"")

    if jac_05_eid > 0.90:
        emit(f"  CONCLUSION: Ensembl-ID Jaccard = {jac_05_eid:.4f} > 0.90")
        emit(f"    The symbol-based divergence was driven by namespace drift.")
        emit(f"    Combined with I < 5%, this meets the negative-result stop rule.")
    else:
        emit(f"  CONCLUSION: Ensembl-ID Jaccard = {jac_05_eid:.4f} < 0.90")
        emit(f"    Annotation drift has genuine pathway-level impact even with")
        emit(f"    stable identifiers. The effect is not a symbol-namespace artifact.")

    with open(args.output, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')

    emit(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
