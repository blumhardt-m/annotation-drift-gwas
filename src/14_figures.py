#!/usr/bin/env python3
"""
Manuscript figures for Annotation Drift in GWAS Locus-to-Gene Mapping.

Produces:
  fig1_pathway_scatter.pdf  — -log10(p) scatter, Ensembl-ID ORA
  fig2_instability_windows.pdf — I vs window size with CI
  fig3_drift_decomposition.pdf — locus breakdown: same/rename/drift
  fig4_concordance_summary.pdf — multi-panel Jaccard robustness

Usage:
    python 14_figures.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --sensitivity results/tables/sensitivity.tsv \
        --outdir results/figures
"""

import csv
import argparse
import math
import os
import random
import subprocess
import sys
from collections import defaultdict
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ═══════════════════════════════════════════════════════════════════════
# Publication style
# ═══════════════════════════════════════════════════════════════════════

def set_pub_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })


# ═══════════════════════════════════════════════════════════════════════
# ORA machinery
# ═══════════════════════════════════════════════════════════════════════

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
        results.append((name, pval))
    return results


def compute_ora_full(v29_genes, v45_genes, universe, pathways):
    ora_v29 = run_ora(v29_genes, universe, pathways)
    ora_v45 = run_ora(v45_genes, universe, pathways)

    names_v29 = [r[0] for r in ora_v29]
    pvals_v29 = [r[1] for r in ora_v29]
    fdr_v29 = bh_correction(pvals_v29)

    names_v45 = [r[0] for r in ora_v45]
    pvals_v45 = [r[1] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)

    pval_map_v29 = dict(zip(names_v29, pvals_v29))
    pval_map_v45 = dict(zip(names_v45, pvals_v45))
    fdr_map_v29 = dict(zip(names_v29, fdr_v29))
    fdr_map_v45 = dict(zip(names_v45, fdr_v45))
    all_names = sorted(set(names_v29) | set(names_v45))

    return all_names, pval_map_v29, pval_map_v45, fdr_map_v29, fdr_map_v45


# ═══════════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════════

def build_alias_table():
    aliases = {}
    for i in range(1, 12):
        aliases[f"MARCH{i}"] = f"MARCHF{i}"
        aliases[f"MARCHF{i}"] = f"MARCH{i}"
    for i in range(1, 15):
        aliases[f"SEPT{i}"] = f"SEPTIN{i}"
        aliases[f"SEPTIN{i}"] = f"SEPT{i}"
    return aliases


def load_tss_tables_ensembl(v29_path, v45_path):
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
    for old_sym, new_sym in alias_table.items():
        if old_sym in symbol_to_ensembl and new_sym not in symbol_to_ensembl:
            symbol_to_ensembl[new_sym] = symbol_to_ensembl[old_sym]
        elif new_sym in symbol_to_ensembl and old_sym not in symbol_to_ensembl:
            symbol_to_ensembl[old_sym] = symbol_to_ensembl[new_sym]
    return symbol_to_ensembl, ensembl_universe


def load_tss_symbols(path):
    symbols = set()
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            symbols.add(row['gene_symbol'])
    return symbols


def load_gmt_symbols(path, min_size, max_size):
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


def load_gmt_ensembl(gmt_path, symbol_to_ensembl, min_size, max_size):
    pathways = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            ensembl_genes = set()
            for sym in parts[2:]:
                eid = symbol_to_ensembl.get(sym)
                if eid:
                    ensembl_genes.add(eid)
            if min_size <= len(ensembl_genes) <= max_size:
                pathways[name] = ensembl_genes
    return pathways


def load_comparison(path):
    """Returns symbol gene sets, ensembl gene sets, and locus-level data."""
    sym_v29 = set()
    sym_v45 = set()
    eid_v29 = set()
    eid_v45 = set()
    loci = []
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            sym_v29.add(row['symbol_v29'])
            sym_v45.add(row['symbol_v45'])
            eid_v29.add(row['gene_v29'])
            eid_v45.add(row['gene_v45'])
            loci.append(row)
    return sym_v29, sym_v45, eid_v29, eid_v45, loci


def load_sensitivity(path):
    rows = []
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append({
                'window_kb': int(row['window_kb']),
                'n': int(row['n']),
                'n_differ': int(row['n_differ']),
                'rate': float(row['rate']),
            })
    return rows


# ═══════════════════════════════════════════════════════════════════════
# Analysis helpers
# ═══════════════════════════════════════════════════════════════════════

def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0, 0, 0
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    spread = z * math.sqrt((p * (1 - p) + z ** 2 / (4 * n)) / n) / denom
    return p, max(0, centre - spread), min(1, centre + spread)


def size_matched_jaccard(sig_larger, sig_smaller, n_reps, rng):
    target_size = len(sig_smaller)
    larger_list = list(sig_larger)
    jaccards = []
    for _ in range(n_reps):
        sample = set(rng.sample(larger_list, target_size))
        inter = sample & sig_smaller
        union = sample | sig_smaller
        jac = len(inter) / len(union) if len(union) > 0 else 1.0
        jaccards.append(jac)
    jaccards.sort()
    mean_jac = sum(jaccards) / len(jaccards)
    lo = jaccards[int(0.025 * n_reps)]
    hi = jaccards[int(0.975 * n_reps)]
    return mean_jac, lo, hi


def top_k_overlap(pval_v29, pval_v45, all_names, k):
    ranked_v29 = sorted(all_names, key=lambda n: pval_v29.get(n, 1.0))[:k]
    ranked_v45 = sorted(all_names, key=lambda n: pval_v45.get(n, 1.0))[:k]
    inter = set(ranked_v29) & set(ranked_v45)
    return len(inter) / k


# ═══════════════════════════════════════════════════════════════════════
# Manifest
# ═══════════════════════════════════════════════════════════════════════

def write_manifest(outdir,
                   eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45,
                   sym_names, sym_pv29, sym_pv45, sym_fdr29, sym_fdr45,
                   eid_universe_size, sym_universe_size,
                   eid_pathways_count, sym_pathways_count,
                   loci_count):
    """Write FIGURE_MANIFEST.txt with ORA verification and environment snapshot."""
    lines = []

    def emit(line=""):
        lines.append(line)

    emit("=" * 70)
    emit("FIGURE_MANIFEST.txt")
    emit(f"Generated: {datetime.now().isoformat()}")
    emit(f"Python: {sys.version}")
    emit("=" * 70)

    emit("")
    emit("Figures recompute ORA from frozen inputs (comparison.tsv, TSS tables,")
    emit("reactome.gmt). This is a reproducible re-run, not a read-only rendering")
    emit("of pre-saved results.")

    # ── Ensembl-ID ORA counts (used by Figure 1) ──
    emit("")
    emit("-" * 70)
    emit("ENSEMBL-ID ORA COUNTS (Figure 1 basis)")
    emit("-" * 70)

    eid_sig29 = {n for n in eid_names if eid_fdr29.get(n, 1.0) < 0.05}
    eid_sig45 = {n for n in eid_names if eid_fdr45.get(n, 1.0) < 0.05}
    eid_inter = eid_sig29 & eid_sig45
    eid_union = eid_sig29 | eid_sig45
    eid_v29_only = eid_sig29 - eid_sig45
    eid_v45_only = eid_sig45 - eid_sig29
    eid_both_ns = len(eid_names) - len(eid_union)

    emit(f"  Pathways tested:    {len(eid_names)}")
    emit(f"  Universe size:      {eid_universe_size}")
    emit(f"  Pathways loaded:    {eid_pathways_count}")
    emit(f"  Sig v29 (FDR<0.05): {len(eid_sig29)}")
    emit(f"  Sig v45 (FDR<0.05): {len(eid_sig45)}")
    emit(f"  Intersection:       {len(eid_inter)}")
    emit(f"  Union:              {len(eid_union)}")
    emit(f"  v29-only:           {len(eid_v29_only)}")
    emit(f"  v45-only:           {len(eid_v45_only)}")
    emit(f"  Both NS:            {eid_both_ns}")
    emit(f"  Jaccard:            {len(eid_inter)/len(eid_union):.4f}" if len(eid_union) > 0 else "  Jaccard:            N/A")

    # ── Cross-check against enrichment_ensembl_id.txt ──
    emit("")
    emit("CROSS-CHECK vs results/tables/enrichment_ensembl_id.txt:")
    ref_path = "results/tables/enrichment_ensembl_id.txt"
    ref_sig29 = ref_sig45 = ref_inter = ref_union = None
    if os.path.exists(ref_path):
        with open(ref_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'v29 sig:' in line and 'v45 sig:' in line:
                    # Parse: "  v29 sig: 74, v45 sig: 47, intersect: 46, union: 75"
                    parts = line.replace(',', '').split()
                    for i, tok in enumerate(parts):
                        if tok == 'sig:' and i > 0 and 'v29' in parts[i - 1]:
                            ref_sig29 = int(parts[i + 1])
                        if tok == 'sig:' and i > 0 and 'v45' in parts[i - 1]:
                            ref_sig45 = int(parts[i + 1])
                        if tok == 'intersect:':
                            ref_inter = int(parts[i + 1])
                        if tok == 'union:':
                            ref_union = int(parts[i + 1])

    if ref_sig29 is not None:
        emit(f"  Reference: sig_v29={ref_sig29}, sig_v45={ref_sig45}, inter={ref_inter}, union={ref_union}")
        emit(f"  Figure:    sig_v29={len(eid_sig29)}, sig_v45={len(eid_sig45)}, inter={len(eid_inter)}, union={len(eid_union)}")
        match = (len(eid_sig29) == ref_sig29 and len(eid_sig45) == ref_sig45 and
                 len(eid_inter) == ref_inter and len(eid_union) == ref_union)
        if match:
            emit("  STATUS: MATCH — Figure 1 uses identical ORA logic as src/12.")
        else:
            emit("  STATUS: MISMATCH — Figure 1 ORA differs from src/12. This is a BUG.")
    else:
        emit(f"  WARNING: Could not parse reference from {ref_path}")

    # ── Symbol-based ORA counts (used by Figure 4) ──
    emit("")
    emit("-" * 70)
    emit("SYMBOL-BASED ORA COUNTS (Figure 4 basis)")
    emit("-" * 70)

    sym_sig29 = {n for n in sym_names if sym_fdr29.get(n, 1.0) < 0.05}
    sym_sig45 = {n for n in sym_names if sym_fdr45.get(n, 1.0) < 0.05}
    sym_inter = sym_sig29 & sym_sig45
    sym_union = sym_sig29 | sym_sig45

    emit(f"  Pathways tested:    {len(sym_names)}")
    emit(f"  Universe size:      {sym_universe_size}")
    emit(f"  Sig v29 (FDR<0.05): {len(sym_sig29)}")
    emit(f"  Sig v45 (FDR<0.05): {len(sym_sig45)}")
    emit(f"  Intersection:       {len(sym_inter)}")
    emit(f"  Union:              {len(sym_union)}")
    emit(f"  Jaccard:            {len(sym_inter)/len(sym_union):.4f}" if len(sym_union) > 0 else "  Jaccard:            N/A")

    # ── Drift decomposition (Figure 3) ──
    emit("")
    emit("-" * 70)
    emit(f"DRIFT DECOMPOSITION (Figure 3, N = {loci_count:,} loci)")
    emit("-" * 70)

    # ── Environment ──
    emit("")
    emit("-" * 70)
    emit("ENVIRONMENT")
    emit("-" * 70)
    try:
        freeze = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'],
                                          text=True, stderr=subprocess.DEVNULL)
        for pkg_line in freeze.strip().split('\n'):
            emit(f"  {pkg_line}")
    except Exception as e:
        emit(f"  pip freeze failed: {e}")

    manifest_path = os.path.join(outdir, 'FIGURE_MANIFEST.txt')
    with open(manifest_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"  Wrote {manifest_path}")
    return manifest_path


# ═══════════════════════════════════════════════════════════════════════
# Figure 1: Pathway p-value scatter (Ensembl-ID ORA)
# ═══════════════════════════════════════════════════════════════════════

def figure1(all_names, pval_v29, pval_v45, fdr_v29, fdr_v45, outdir):
    print("  Generating Figure 1: pathway scatter...")
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    # Find FDR=0.05 p-value threshold for each version
    # The threshold is the largest raw p-value with FDR < 0.05
    thresh_v29 = 0
    for n in all_names:
        if fdr_v29.get(n, 1.0) < 0.05:
            p = pval_v29.get(n, 1.0)
            if p > thresh_v29:
                thresh_v29 = p
    thresh_v45 = 0
    for n in all_names:
        if fdr_v45.get(n, 1.0) < 0.05:
            p = pval_v45.get(n, 1.0)
            if p > thresh_v45:
                thresh_v45 = p

    log_thresh_v29 = -math.log10(max(thresh_v29, 1e-300))
    log_thresh_v45 = -math.log10(max(thresh_v45, 1e-300))

    # Classify and collect points
    xs_grey, ys_grey = [], []
    xs_red, ys_red = [], []
    xs_blue, ys_blue = [], []

    n_agree_sig = 0
    n_agree_ns = 0
    n_v29_only = 0
    n_v45_only = 0

    for name in all_names:
        p29 = pval_v29.get(name, 1.0)
        p45 = pval_v45.get(name, 1.0)
        f29 = fdr_v29.get(name, 1.0)
        f45 = fdr_v45.get(name, 1.0)
        lp29 = -math.log10(max(p29, 1e-300))
        lp45 = -math.log10(max(p45, 1e-300))

        sig29 = f29 < 0.05
        sig45 = f45 < 0.05

        if sig29 and sig45:
            xs_grey.append(lp29)
            ys_grey.append(lp45)
            n_agree_sig += 1
        elif not sig29 and not sig45:
            xs_grey.append(lp29)
            ys_grey.append(lp45)
            n_agree_ns += 1
        elif sig29 and not sig45:
            xs_red.append(lp29)
            ys_red.append(lp45)
            n_v29_only += 1
        else:
            xs_blue.append(lp29)
            ys_blue.append(lp45)
            n_v45_only += 1

    ax.scatter(xs_grey, ys_grey, c='#AAAAAA', s=6, alpha=0.4, linewidths=0, rasterized=True, zorder=1)
    ax.scatter(xs_red, ys_red, c='#D62728', s=10, alpha=0.7, linewidths=0, rasterized=True, zorder=2, label=f'v29 only (n={n_v29_only})')
    ax.scatter(xs_blue, ys_blue, c='#1F77B4', s=10, alpha=0.7, linewidths=0, rasterized=True, zorder=2, label=f'v45 only (n={n_v45_only})')

    # Threshold lines
    ax.axvline(log_thresh_v29, color='#D62728', ls='--', lw=0.6, alpha=0.6)
    ax.axhline(log_thresh_v45, color='#1F77B4', ls='--', lw=0.6, alpha=0.6)

    # Identity line
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], color='black', ls=':', lw=0.5, alpha=0.4, zorder=0)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect('equal')

    # Quadrant counts
    fontdict = {'fontsize': 6, 'color': '#555555'}
    ax.text(0.03, 0.97, f'Both NS\nn={n_agree_ns}', transform=ax.transAxes, va='top', ha='left', **fontdict)
    ax.text(0.97, 0.03, f'Both sig\nn={n_agree_sig}', transform=ax.transAxes, va='bottom', ha='right', **fontdict)
    ax.text(0.97, 0.97, f'v29 only\nn={n_v29_only}', transform=ax.transAxes, va='top', ha='right',
            fontsize=6, color='#D62728')
    ax.text(0.03, 0.03, f'v45 only\nn={n_v45_only}', transform=ax.transAxes, va='bottom', ha='left',
            fontsize=6, color='#1F77B4')

    ax.set_xlabel('$-\\log_{10}(p)$, GENCODE v29')
    ax.set_ylabel('$-\\log_{10}(p)$, GENCODE v45')
    ax.set_title('Pathway enrichment: v29 vs v45\n(Ensembl-ID-based ORA)')

    path = os.path.join(outdir, 'fig1_pathway_scatter.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(f"    Saved {path}")


# ═══════════════════════════════════════════════════════════════════════
# Figure 2: Instability vs window size
# ═══════════════════════════════════════════════════════════════════════

def figure2(sensitivity_rows, outdir):
    print("  Generating Figure 2: instability vs window...")
    fig, ax = plt.subplots(figsize=(3.5, 2.5))

    windows = [r['window_kb'] for r in sensitivity_rows]
    rates = [r['rate'] * 100 for r in sensitivity_rows]

    # Wilson CIs
    lowers, uppers = [], []
    for r in sensitivity_rows:
        _, lo, hi = wilson_ci(r['n_differ'], r['n'])
        lowers.append(lo * 100)
        uppers.append(hi * 100)

    yerr_lo = [r - l for r, l in zip(rates, lowers)]
    yerr_hi = [u - r for r, u in zip(rates, uppers)]

    ax.bar(range(len(windows)), rates, width=0.5, color='#4C72B0', edgecolor='white', linewidth=0.5, zorder=2)
    ax.errorbar(range(len(windows)), rates, yerr=[yerr_lo, yerr_hi],
                fmt='none', ecolor='black', elinewidth=0.8, capsize=3, capthick=0.8, zorder=3)

    ax.set_xticks(range(len(windows)))
    ax.set_xticklabels([f'{w}' for w in windows])
    ax.set_xlabel('Window size (kb)')
    ax.set_ylabel('Instability rate (%)')
    ax.set_title('Locus-to-gene instability by window size')

    # 5% reference line
    ax.axhline(5.0, color='#999999', ls='--', lw=0.6, zorder=1)
    ax.text(len(windows) - 0.5, 5.1, '5% threshold', fontsize=6, color='#999999', ha='right', va='bottom')

    ax.set_ylim(0, max(rates) * 1.3)

    path = os.path.join(outdir, 'fig2_instability_windows.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(f"    Saved {path}")


# ═══════════════════════════════════════════════════════════════════════
# Figure 3: Drift decomposition
# ═══════════════════════════════════════════════════════════════════════

def figure3(loci, outdir):
    print("  Generating Figure 3: drift decomposition...")

    n_total = len(loci)
    n_same_eid_same_sym = sum(
        1 for l in loci if l['gene_v29'] == l['gene_v45'] and l['symbol_v29'] == l['symbol_v45']
    )
    n_same_eid_diff_sym = sum(
        1 for l in loci if l['gene_v29'] == l['gene_v45'] and l['symbol_v29'] != l['symbol_v45']
    )
    n_diff_eid = sum(1 for l in loci if l['gene_v29'] != l['gene_v45'])

    categories = ['Same gene,\nsame symbol', 'Same gene,\ndifferent symbol\n(renames)', 'Different gene\n(true drift)']
    counts = [n_same_eid_same_sym, n_same_eid_diff_sym, n_diff_eid]
    colors = ['#88BDE6', '#F4A460', '#D62728']
    pcts = [100 * c / n_total for c in counts]

    fig, ax = plt.subplots(figsize=(4, 2.8))

    bars = ax.bar(range(3), counts, color=colors, edgecolor='white', linewidth=0.5, width=0.6)

    for i, (bar, count, pct) in enumerate(zip(bars, counts, pcts)):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + n_total * 0.01,
                f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=7)

    ax.set_xticks(range(3))
    ax.set_xticklabels(categories, fontsize=7)
    ax.set_ylabel('Number of loci')
    ax.set_title(f'Drift decomposition (N = {n_total:,} loci)')
    ax.set_ylim(0, max(counts) * 1.15)

    # Format y-axis with comma separator
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))

    path = os.path.join(outdir, 'fig3_drift_decomposition.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(f"    Saved {path}")


# ═══════════════════════════════════════════════════════════════════════
# Figure 4: Concordance summary (multi-panel)
# ═══════════════════════════════════════════════════════════════════════

def figure4(sym_names, sym_pv29, sym_pv45, sym_fdr29, sym_fdr45,
            eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45,
            outdir):
    print("  Generating Figure 4: concordance summary...")

    rng = random.Random(42)
    thresholds = [0.01, 0.05, 0.10, 0.20]
    k_values = [25, 50, 100, 200]

    # Panel (a): Jaccard vs FDR
    sym_jacs = []
    eid_jacs = []
    for alpha in thresholds:
        sig29 = {n for n in sym_names if sym_fdr29.get(n, 1.0) < alpha}
        sig45 = {n for n in sym_names if sym_fdr45.get(n, 1.0) < alpha}
        union = sig29 | sig45
        sym_jacs.append(len(sig29 & sig45) / len(union) if len(union) > 0 else 1.0)

        sig29 = {n for n in eid_names if eid_fdr29.get(n, 1.0) < alpha}
        sig45 = {n for n in eid_names if eid_fdr45.get(n, 1.0) < alpha}
        union = sig29 | sig45
        eid_jacs.append(len(sig29 & sig45) / len(union) if len(union) > 0 else 1.0)

    # Panel (b): Top-K overlap
    sym_topk = [top_k_overlap(sym_pv29, sym_pv45, sym_names, k) for k in k_values]
    eid_topk = [top_k_overlap(eid_pv29, eid_pv45, eid_names, k) for k in k_values]

    # Panel (c): Size-matched Jaccard at FDR=0.05
    alpha = 0.05
    sym_sig29 = {n for n in sym_names if sym_fdr29.get(n, 1.0) < alpha}
    sym_sig45 = {n for n in sym_names if sym_fdr45.get(n, 1.0) < alpha}
    eid_sig29 = {n for n in eid_names if eid_fdr29.get(n, 1.0) < alpha}
    eid_sig45 = {n for n in eid_names if eid_fdr45.get(n, 1.0) < alpha}

    # Raw Jaccards
    sym_raw = len(sym_sig29 & sym_sig45) / len(sym_sig29 | sym_sig45) if len(sym_sig29 | sym_sig45) > 0 else 1.0
    eid_raw = len(eid_sig29 & eid_sig45) / len(eid_sig29 | eid_sig45) if len(eid_sig29 | eid_sig45) > 0 else 1.0

    # Size-matched
    if len(sym_sig29) > len(sym_sig45):
        sym_sm_mean, sym_sm_lo, sym_sm_hi = size_matched_jaccard(sym_sig29, sym_sig45, 1000, rng)
    elif len(sym_sig45) > len(sym_sig29):
        sym_sm_mean, sym_sm_lo, sym_sm_hi = size_matched_jaccard(sym_sig45, sym_sig29, 1000, rng)
    else:
        sym_sm_mean, sym_sm_lo, sym_sm_hi = sym_raw, sym_raw, sym_raw

    if len(eid_sig29) > len(eid_sig45):
        eid_sm_mean, eid_sm_lo, eid_sm_hi = size_matched_jaccard(eid_sig29, eid_sig45, 1000, rng)
    elif len(eid_sig45) > len(eid_sig29):
        eid_sm_mean, eid_sm_lo, eid_sm_hi = size_matched_jaccard(eid_sig45, eid_sig29, 1000, rng)
    else:
        eid_sm_mean, eid_sm_lo, eid_sm_hi = eid_raw, eid_raw, eid_raw

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.5))

    col_sym = '#4C72B0'
    col_eid = '#DD8452'

    # (a) Jaccard vs FDR
    ax = axes[0]
    ax.plot(thresholds, sym_jacs, 'o-', color=col_sym, ms=4, lw=1.2, label='Symbol')
    ax.plot(thresholds, eid_jacs, 's-', color=col_eid, ms=4, lw=1.2, label='Ensembl ID')
    ax.axhline(0.90, color='#999999', ls='--', lw=0.6)
    ax.text(0.20, 0.91, '0.90', fontsize=6, color='#999999', ha='right')
    ax.set_xlabel('FDR threshold')
    ax.set_ylabel('Pathway-set Jaccard')
    ax.set_title('(a) Jaccard vs FDR threshold')
    ax.set_ylim(0, 1.0)
    ax.set_xticks(thresholds)
    ax.legend(loc='lower right', frameon=False)

    # (b) Top-K overlap
    ax = axes[1]
    ax.plot(k_values, sym_topk, 'o-', color=col_sym, ms=4, lw=1.2, label='Symbol')
    ax.plot(k_values, eid_topk, 's-', color=col_eid, ms=4, lw=1.2, label='Ensembl ID')
    ax.axhline(0.90, color='#999999', ls='--', lw=0.6)
    ax.set_xlabel('K (top pathways)')
    ax.set_ylabel('Overlap (intersection / K)')
    ax.set_title('(b) Top-K overlap')
    ax.set_ylim(0, 1.0)
    ax.set_xticks(k_values)
    ax.legend(loc='lower right', frameon=False)

    # (c) Size-matched Jaccard at FDR=0.05
    ax = axes[2]
    x = [0, 1]
    labels = ['Symbol', 'Ensembl ID']
    raw_vals = [sym_raw, eid_raw]
    sm_vals = [sym_sm_mean, eid_sm_mean]
    sm_los = [sym_sm_lo, eid_sm_lo]
    sm_his = [sym_sm_hi, eid_sm_hi]

    bar_w = 0.3
    bars_raw = ax.bar([xi - bar_w / 2 for xi in x], raw_vals, width=bar_w,
                       color=[col_sym, col_eid], alpha=0.4, edgecolor='white', label='Raw')
    bars_sm = ax.bar([xi + bar_w / 2 for xi in x], sm_vals, width=bar_w,
                      color=[col_sym, col_eid], edgecolor='white', label='Size-matched')
    # CI error bars on size-matched
    yerr_lo = [sm - lo for sm, lo in zip(sm_vals, sm_los)]
    yerr_hi = [hi - sm for sm, hi in zip(sm_vals, sm_his)]
    ax.errorbar([xi + bar_w / 2 for xi in x], sm_vals, yerr=[yerr_lo, yerr_hi],
                fmt='none', ecolor='black', elinewidth=0.8, capsize=3, capthick=0.8)

    ax.axhline(0.90, color='#999999', ls='--', lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Jaccard')
    ax.set_title('(c) Size-matched Jaccard\n(FDR < 0.05)')
    ax.set_ylim(0, 1.0)

    # Custom legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc='grey', alpha=0.4, ec='white'),
        plt.Rectangle((0, 0), 1, 1, fc='grey', ec='white'),
    ]
    ax.legend(handles, ['Raw', 'Size-matched'], loc='upper right', frameon=False)

    plt.tight_layout()
    path = os.path.join(outdir, 'fig4_concordance_summary.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(f"    Saved {path}")


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--sensitivity', required=True)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    set_pub_style()

    print("Loading data...")

    # Comparison data
    sym_v29, sym_v45, eid_v29, eid_v45, loci = load_comparison(args.comparison)

    # Symbol-based ORA
    sym_universe = load_tss_symbols(args.tss_v29) | load_tss_symbols(args.tss_v45)
    sym_pathways = load_gmt_symbols(args.gmt, args.min_size, args.max_size)
    print("Running symbol-based ORA...")
    sym_names, sym_pv29, sym_pv45, sym_fdr29, sym_fdr45 = compute_ora_full(
        sym_v29, sym_v45, sym_universe, sym_pathways
    )

    # Ensembl-ID-based ORA
    sym2ens, eid_universe = load_tss_tables_ensembl(args.tss_v29, args.tss_v45)
    eid_pathways = load_gmt_ensembl(args.gmt, sym2ens, args.min_size, args.max_size)
    print("Running Ensembl-ID-based ORA...")
    eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45 = compute_ora_full(
        eid_v29, eid_v45, eid_universe, eid_pathways
    )

    # Sensitivity data
    sensitivity = load_sensitivity(args.sensitivity)

    # Write manifest before generating figures
    print("\nWriting manifest...")
    write_manifest(
        args.outdir,
        eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45,
        sym_names, sym_pv29, sym_pv45, sym_fdr29, sym_fdr45,
        eid_universe_size=len(eid_universe),
        sym_universe_size=len(sym_universe),
        eid_pathways_count=len(eid_pathways),
        sym_pathways_count=len(sym_pathways),
        loci_count=len(loci),
    )

    # Generate figures
    print("\nGenerating figures...")
    figure1(eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45, args.outdir)
    figure2(sensitivity, args.outdir)
    figure3(loci, args.outdir)
    figure4(sym_names, sym_pv29, sym_pv45, sym_fdr29, sym_fdr45,
            eid_names, eid_pv29, eid_pv45, eid_fdr29, eid_fdr45,
            args.outdir)

    print("\nAll figures saved to", args.outdir)


if __name__ == "__main__":
    main()
