#!/usr/bin/env python3
"""
Export primary v29-v45 Reactome enrichment results as per-version TSVs.

Reuses the exact ORA machinery from src/12_ensembl_id_enrichment.py:
  - Ensembl-ID-based ORA (Fisher's exact, one-sided)
  - BH FDR correction
  - Gene universe = union of protein-coding genes from both TSS tables
  - Pathway size filter: 15–500 genes
  - GMT symbol→Ensembl mapping with HGNC alias awareness

Outputs:
  results/primary_reactome/reactome_enrichment_v29.tsv
  results/primary_reactome/reactome_enrichment_v45.tsv
"""

import csv
import math
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMPARISON = os.path.join(BASE, "data", "processed", "mappings", "comparison.tsv")
TSS_V29 = os.path.join(BASE, "data", "processed", "tss_tables", "tss_v29.tsv")
TSS_V45 = os.path.join(BASE, "data", "processed", "tss_tables", "tss_v45.tsv")
GMT = os.path.join(BASE, "data", "raw", "reactome.gmt")
OUT_DIR = os.path.join(BASE, "results", "primary_reactome")

MIN_SIZE = 15
MAX_SIZE = 500


# ── Alias table (identical to 12_ensembl_id_enrichment.py) ─────────────

def build_alias_table():
    aliases = {}
    for i in range(1, 12):
        aliases[f"MARCH{i}"] = f"MARCHF{i}"
        aliases[f"MARCHF{i}"] = f"MARCH{i}"
    for i in range(1, 15):
        aliases[f"SEPT{i}"] = f"SEPTIN{i}"
        aliases[f"SEPTIN{i}"] = f"SEPT{i}"
    return aliases


# ── Data loading ───────────────────────────────────────────────────────

def load_tss_tables():
    alias_table = build_alias_table()
    symbol_to_ensembl = {}
    ensembl_universe = set()

    for path in [TSS_V29, TSS_V45]:
        with open(path, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                eid = row['ensembl_id']
                sym = row['gene_symbol']
                ensembl_universe.add(eid)
                symbol_to_ensembl[sym] = eid

    for old_sym, new_sym in alias_table.items():
        if old_sym in symbol_to_ensembl and new_sym not in symbol_to_ensembl:
            symbol_to_ensembl[new_sym] = symbol_to_ensembl[old_sym]
        elif new_sym in symbol_to_ensembl and old_sym not in symbol_to_ensembl:
            symbol_to_ensembl[old_sym] = symbol_to_ensembl[new_sym]

    return symbol_to_ensembl, ensembl_universe


def load_gmt_as_ensembl(symbol_to_ensembl):
    """Load GMT with pathway_id preserved, converting symbols to Ensembl IDs."""
    pathways = {}
    with open(GMT, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            pw_id = parts[1]
            symbols = parts[2:]
            ensembl_genes = set()
            for sym in symbols:
                eid = symbol_to_ensembl.get(sym)
                if eid:
                    ensembl_genes.add(eid)
            if MIN_SIZE <= len(ensembl_genes) <= MAX_SIZE:
                pathways[name] = {'id': pw_id, 'genes': ensembl_genes}
    return pathways


def load_comparison():
    v29_genes = set()
    v45_genes = set()
    n = 0
    with open(COMPARISON, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            v29_genes.add(row['gene_v29'])
            v45_genes.add(row['gene_v45'])
            n += 1
    return v29_genes, v45_genes, n


# ── ORA machinery (identical to 12_ensembl_id_enrichment.py) ──────────

def fisher_exact_pvalue(a, b, c, d):
    n = a + b + c + d

    def log_factorial(x):
        if x <= 1:
            return 0.0
        return sum(math.log(i) for i in range(2, x + 1))

    def log_choose(nn, k):
        if k < 0 or k > nn:
            return float('-inf')
        return log_factorial(nn) - log_factorial(k) - log_factorial(nn - k)

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
    for name, pw_data in pathways.items():
        pw_genes = pw_data['genes']
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
        results.append({
            'pathway_id': pw_data['id'],
            'pathway_name': name,
            'p_value': pval,
            'gene_count': a,
        })
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


# ── Output ─────────────────────────────────────────────────────────────

def write_enrichment_tsv(ora_results, fdr_values, out_path):
    fieldnames = ['pathway_id', 'pathway_name', 'p_value', 'fdr', 'gene_count']
    rows = []
    for res, fdr_val in zip(ora_results, fdr_values):
        rows.append({
            'pathway_id': res['pathway_id'],
            'pathway_name': res['pathway_name'],
            'p_value': f"{res['p_value']:.6e}",
            'fdr': f"{fdr_val:.6e}",
            'gene_count': res['gene_count'],
        })
    rows.sort(key=lambda r: float(r['p_value']))
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


# ── Main ───────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 65)
    print("PRIMARY REACTOME ENRICHMENT EXPORT (v29 vs v45)")
    print("=" * 65)

    # Load data
    print("\nLoading TSS tables...")
    symbol_to_ensembl, ensembl_universe = load_tss_tables()
    print(f"  Ensembl universe: {len(ensembl_universe)} genes")

    print("Loading GMT...")
    pathways = load_gmt_as_ensembl(symbol_to_ensembl)
    print(f"  Pathways (size {MIN_SIZE}-{MAX_SIZE}): {len(pathways)}")

    print("Loading comparison...")
    v29_genes, v45_genes, n_loci = load_comparison()
    print(f"  Loci: {n_loci:,}")
    print(f"  v29 unique genes: {len(v29_genes):,}")
    print(f"  v45 unique genes: {len(v45_genes):,}")

    # Run ORA
    print("\nRunning ORA (v29)...")
    ora_v29 = run_ora(v29_genes, ensembl_universe, pathways)
    print("Running ORA (v45)...")
    ora_v45 = run_ora(v45_genes, ensembl_universe, pathways)

    # BH correction
    fdr_v29 = bh_correction([r['p_value'] for r in ora_v29])
    fdr_v45 = bh_correction([r['p_value'] for r in ora_v45])

    # Write TSVs
    out_v29 = os.path.join(OUT_DIR, "reactome_enrichment_v29.tsv")
    out_v45 = os.path.join(OUT_DIR, "reactome_enrichment_v45.tsv")

    n_rows_v29 = write_enrichment_tsv(ora_v29, fdr_v29, out_v29)
    n_rows_v45 = write_enrichment_tsv(ora_v45, fdr_v45, out_v45)

    # Compute Spearman ρ on -log10(p)
    pval_map_v29 = {r['pathway_name']: r['p_value'] for r in ora_v29}
    pval_map_v45 = {r['pathway_name']: r['p_value'] for r in ora_v45}
    all_names = sorted(set(pval_map_v29.keys()) | set(pval_map_v45.keys()))

    log10p_v29 = [-math.log10(max(pval_map_v29.get(n, 1.0), 1e-300)) for n in all_names]
    log10p_v45 = [-math.log10(max(pval_map_v45.get(n, 1.0), 1e-300)) for n in all_names]
    rho = spearman_rank_corr(log10p_v29, log10p_v45)

    # Jaccard at FDR < 0.05
    fdr_map_v29 = {r['pathway_name']: fdr for r, fdr in zip(ora_v29, fdr_v29)}
    fdr_map_v45 = {r['pathway_name']: fdr for r, fdr in zip(ora_v45, fdr_v45)}
    sig_v29 = {n for n in all_names if fdr_map_v29.get(n, 1.0) < 0.05}
    sig_v45 = {n for n in all_names if fdr_map_v45.get(n, 1.0) < 0.05}
    inter = sig_v29 & sig_v45
    union = sig_v29 | sig_v45
    jaccard = len(inter) / len(union) if union else 1.0

    # Report
    print(f"\n{'=' * 65}")
    print("RESULTS")
    print(f"{'=' * 65}")
    print(f"  reactome_enrichment_v29.tsv : {n_rows_v29} rows")
    print(f"  reactome_enrichment_v45.tsv : {n_rows_v45} rows")
    print(f"  Spearman ρ (-log10 p):        {rho:.6f}  (expected ≈ 0.923)")
    print(f"  Jaccard (FDR < 0.05):         {jaccard:.4f}  (expected ≈ 0.613)")
    print(f"  Sig v29: {len(sig_v29)}, Sig v45: {len(sig_v45)}, "
          f"Shared: {len(inter)}, Union: {len(union)}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
