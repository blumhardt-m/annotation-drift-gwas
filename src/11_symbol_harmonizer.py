#!/usr/bin/env python3
"""
Symbol harmonizer: fix MARCH/SEPT alias mismatches in Reactome GMT,
then rerun ORA to measure corrected pathway-set Jaccard.

Usage:
    python 11_symbol_harmonizer.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --output results/tables/enrichment_harmonized.txt
"""

import csv
import argparse
import math
from collections import defaultdict


# ── 1. HGNC alias lookup ────────────────────────────────────────────────

def build_alias_table():
    """Build bidirectional old ↔ new HGNC symbol mapping.

    Covers the 2020–2023 HGNC renames driven by the Excel auto-format problem:
      MARCH1-11  → MARCHF1-11
      SEPT1-14   → SEPTIN1-14  (SEPT6, SEPT13 may not exist in all builds)
    """
    aliases = {}

    for i in range(1, 12):  # MARCH1–MARCH11
        old = f"MARCH{i}"
        new = f"MARCHF{i}"
        aliases[old] = new
        aliases[new] = old

    for i in range(1, 15):  # SEPT1–SEPT14
        old = f"SEPT{i}"
        new = f"SEPTIN{i}"
        aliases[old] = new
        aliases[new] = old

    return aliases


def harmonize_symbol(symbol, universe, alias_table):
    """If symbol is in universe, return it unchanged.
    If not, check alias table and return the alias if it IS in universe.
    Otherwise return the original (unresolvable)."""
    if symbol in universe:
        return symbol
    alias = alias_table.get(symbol)
    if alias and alias in universe:
        return alias
    return symbol


# ── ORA machinery (replicated from 09_enrichment.py) ────────────────────

def load_gmt_harmonized(path, min_size, max_size, universe, alias_table):
    """Load GMT with symbol harmonization. Returns (pathways, n_substitutions)."""
    pathways = {}
    n_subs = 0
    sub_details = defaultdict(int)
    with open(path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            raw_genes = parts[2:]
            harmonized = set()
            for g in raw_genes:
                h = harmonize_symbol(g, universe, alias_table)
                if h != g:
                    n_subs += 1
                    sub_details[f"{g} → {h}"] += 1
                harmonized.add(h)
            if min_size <= len(harmonized) <= max_size:
                pathways[name] = harmonized
    return pathways, n_subs, sub_details


def load_gmt_raw(path, min_size, max_size):
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


def load_tss_symbols(path):
    symbols = set()
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            symbols.add(row['gene_symbol'])
    return symbols


def load_mapped_genes(comparison_path):
    v29 = set()
    v45 = set()
    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            v29.add(row['symbol_v29'])
            v45.add(row['symbol_v45'])
    return v29, v45


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


# ── Sweep helper ────────────────────────────────────────────────────────

def fdr_sweep(ora_v29, ora_v45, thresholds):
    """Compute pathway-set Jaccard at multiple FDR thresholds."""
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

    # Spearman
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
    emit("SYMBOL HARMONIZATION — ENRICHMENT CORRECTION")
    emit("=" * 70)

    # ── Build alias table ──
    alias_table = build_alias_table()
    emit(f"\n1. ALIAS LOOKUP TABLE ({len(alias_table)} entries)")
    emit("-" * 70)
    seen = set()
    for old, new in sorted(alias_table.items()):
        pair = tuple(sorted([old, new]))
        if pair not in seen:
            seen.add(pair)
            emit(f"  {old:<12} ↔ {new}")

    # ── Load data ──
    emit("\nLoading data...")
    tss_v29 = load_tss_symbols(args.tss_v29)
    tss_v45 = load_tss_symbols(args.tss_v45)
    universe = tss_v29 | tss_v45
    v29_genes, v45_genes = load_mapped_genes(args.comparison)
    emit(f"  Universe: {len(universe)} genes")

    # ── Load GMT raw (for baseline) ──
    pathways_raw = load_gmt_raw(args.gmt, args.min_size, args.max_size)

    # ── Load GMT harmonized ──
    pathways_harm, n_subs, sub_details = load_gmt_harmonized(
        args.gmt, args.min_size, args.max_size, universe, alias_table
    )

    emit(f"\n2. HARMONIZATION REPORT")
    emit("-" * 70)
    emit(f"  Total symbol substitutions: {n_subs}")
    emit(f"  Unique substitution types:  {len(sub_details)}")
    emit(f"  Pathways (raw):             {len(pathways_raw)}")
    emit(f"  Pathways (harmonized):      {len(pathways_harm)}")
    emit(f"\n  Substitution details:")
    for sub, count in sorted(sub_details.items(), key=lambda x: -x[1]):
        emit(f"    {sub:<30} (in {count} pathways)")

    # ── Run ORA: raw baseline ──
    emit(f"\n3. ORA COMPARISON: RAW vs HARMONIZED")
    emit("-" * 70)

    thresholds = [0.01, 0.05, 0.10, 0.20]

    emit("\nRunning raw ORA...")
    ora_v29_raw = run_ora(v29_genes, universe, pathways_raw)
    ora_v45_raw = run_ora(v45_genes, universe, pathways_raw)
    raw_results, raw_rho = fdr_sweep(ora_v29_raw, ora_v45_raw, thresholds)

    emit("Running harmonized ORA...")
    ora_v29_harm = run_ora(v29_genes, universe, pathways_harm)
    ora_v45_harm = run_ora(v45_genes, universe, pathways_harm)
    harm_results, harm_rho = fdr_sweep(ora_v29_harm, ora_v45_harm, thresholds)

    # ── Report table ──
    emit(f"\n{'FDR':<8} {'Raw Jaccard':>12} {'Harm Jaccard':>13} {'Delta':>8} {'Direction':>10}")
    emit("-" * 53)
    for (alpha, _, _, _, _, jac_raw), (_, _, _, _, _, jac_harm) in zip(raw_results, harm_results):
        delta = jac_harm - jac_raw
        direction = "+" if delta > 0 else ("-" if delta < 0 else "=")
        emit(f"{alpha:<8.2f} {jac_raw:>12.4f} {jac_harm:>13.4f} {delta:>+8.4f} {direction:>10}")

    emit(f"\nSpearman rho (-log10 p):")
    emit(f"  Raw:        {raw_rho:.6f}")
    emit(f"  Harmonized: {harm_rho:.6f}")
    emit(f"  Delta:      {harm_rho - raw_rho:+.6f}")

    # ── Detailed FDR=0.05 comparison ──
    emit(f"\n4. DETAILED FDR < 0.05 COMPARISON")
    emit("-" * 70)

    names_raw_v29 = [r[0] for r in ora_v29_raw]
    fdr_raw_v29 = bh_correction([r[1] for r in ora_v29_raw])
    names_raw_v45 = [r[0] for r in ora_v45_raw]
    fdr_raw_v45 = bh_correction([r[1] for r in ora_v45_raw])

    names_harm_v29 = [r[0] for r in ora_v29_harm]
    fdr_harm_v29 = bh_correction([r[1] for r in ora_v29_harm])
    names_harm_v45 = [r[0] for r in ora_v45_harm]
    fdr_harm_v45 = bh_correction([r[1] for r in ora_v45_harm])

    sig_raw_v29 = {n for n, f in zip(names_raw_v29, fdr_raw_v29) if f < 0.05}
    sig_raw_v45 = {n for n, f in zip(names_raw_v45, fdr_raw_v45) if f < 0.05}
    sig_harm_v29 = {n for n, f in zip(names_harm_v29, fdr_harm_v29) if f < 0.05}
    sig_harm_v45 = {n for n, f in zip(names_harm_v45, fdr_harm_v45) if f < 0.05}

    raw_union = sig_raw_v29 | sig_raw_v45
    harm_union = sig_harm_v29 | sig_harm_v45
    raw_inter = sig_raw_v29 & sig_raw_v45
    harm_inter = sig_harm_v29 & sig_harm_v45

    emit(f"  Raw:  v29={len(sig_raw_v29)}, v45={len(sig_raw_v45)}, "
         f"inter={len(raw_inter)}, union={len(raw_union)}")
    emit(f"  Harm: v29={len(sig_harm_v29)}, v45={len(sig_harm_v45)}, "
         f"inter={len(harm_inter)}, union={len(harm_union)}")

    # Pathways that changed status
    gained = harm_inter - raw_inter
    lost = raw_inter - harm_inter
    if gained:
        emit(f"\n  Pathways that became concordant after harmonization ({len(gained)}):")
        for p in sorted(gained)[:15]:
            emit(f"    + {p}")
        if len(gained) > 15:
            emit(f"    ... and {len(gained) - 15} more")
    if lost:
        emit(f"\n  Pathways that became discordant after harmonization ({len(lost)}):")
        for p in sorted(lost)[:15]:
            emit(f"    - {p}")
        if len(lost) > 15:
            emit(f"    ... and {len(lost) - 15} more")

    # ── Summary ──
    emit(f"\n{'='*70}")
    emit("SUMMARY")
    emit(f"{'='*70}")
    jac_raw_05 = raw_results[1][5]
    jac_harm_05 = harm_results[1][5]
    emit(f"  Pathway-set Jaccard (FDR<0.05):")
    emit(f"    Before harmonization: {jac_raw_05:.4f}")
    emit(f"    After harmonization:  {jac_harm_05:.4f}")
    emit(f"    Delta:                {jac_harm_05 - jac_raw_05:+.4f}")
    emit(f"  Spearman rho: {raw_rho:.4f} → {harm_rho:.4f} ({harm_rho - raw_rho:+.4f})")
    emit(f"  Symbol substitutions:   {n_subs} across {len(sub_details)} rename patterns")
    residual = jac_harm_05
    emit(f"\n  CONCLUSION: After fixing MARCH/SEPT aliases, pathway Jaccard = {residual:.4f}.")
    if residual > 0.90:
        emit("  → Corrected Jaccard > 0.90: alias mismatch was the primary driver.")
        emit("  → Combined with I < 5%, this satisfies the negative-result stop rule.")
    else:
        emit(f"  → Corrected Jaccard = {residual:.4f} < 0.90: annotation drift has")
        emit("    genuine pathway-level impact beyond the alias artifact.")

    # Write output
    with open(args.output, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')

    emit(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
