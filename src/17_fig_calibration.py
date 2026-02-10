#!/usr/bin/env python3
"""
Figure 5: Calibration panel — baselines, saturation diagnostic, and flip rate.

Panel A: Spearman rho across ORA regimes
Panel B: Saturation diagnostic (Jaccard vs foreground size)
Panel C: Category-level pathway counts with flip rate

Usage:
    python 17_fig_calibration.py \
        --baseline-fixed results/tables/baseline_fixed.tsv \
        --baseline-random results/tables/baseline_random.tsv \
        --diff-enrichment results/tables/differential_enrichment.tsv \
        --diff-random results/tables/differential_baseline_random.tsv \
        --saturation results/tables/saturation_curve.tsv \
        --category-summary results/tables/category_summary.tsv \
        --flip-rate results/tables/flip_rate.tsv \
        --enrichment-ensid results/tables/enrichment_ensembl_id.txt \
        --outdir results/figures
"""

import argparse
import csv
import datetime
import hashlib
import math
import os
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── Publication style (from src/14_figures.py) ─────────────────────────

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


# ── Data loaders ───────────────────────────────────────────────────────

def load_metric_value(path):
    """Load metric→value TSV into dict."""
    d = {}
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            d[row['metric']] = row.get('value', row.get('mean', ''))
    return d


def load_metric_summary(path):
    """Load metric summary TSV (with mean, sd, ci_lo, ci_hi)."""
    d = {}
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            d[row['metric']] = row
    return d


def parse_real_rho(path):
    """Extract Ensembl-ID Spearman rho from enrichment_ensembl_id.txt."""
    with open(path, 'r') as f:
        for line in f:
            if 'Ensembl-ID-based:' in line and 'Spearman' not in line:
                # Line: "  Ensembl-ID-based: 0.922857"
                parts = line.strip().split()
                return float(parts[-1])
    return float('nan')


def load_saturation(path):
    """Load saturation curve TSV."""
    rows = []
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            m = int(row['m'])
            mj = row['mean_jaccard']
            sj = row['sd_jaccard']
            mean_j = float(mj) if mj != 'nan' else float('nan')
            sd_j = float(sj) if sj != 'nan' else float('nan')
            rows.append((m, mean_j, sd_j))
    return rows


def load_category_summary(path):
    """Load category summary TSV."""
    rows = []
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append({
                'category': row['category'],
                'n_v29': int(row['n_v29']),
                'n_v45': int(row['n_v45']),
                'status': row['status'],
            })
    return rows


# ── Colors (consistent with src/14_figures.py palette) ─────────────────

COL_V29 = '#4C72B0'
COL_V45 = '#DD8452'
COL_GREY = '#AAAAAA'
COL_RED = '#D62728'


# ── Provenance helpers ────────────────────────────────────────────────

def git_commit_hash():
    """Return short git commit hash, or 'unknown' if not in a repo."""
    try:
        out = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            stderr=subprocess.DEVNULL, text=True)
        return out.strip()
    except Exception:
        return 'unknown'


def file_sha256(path):
    """Return SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


# ── Panel A: Spearman rho across regimes ───────────────────────────────

def panel_a(ax, rho_fixed, rho_real, rho_diff, rho_random_mean,
            rho_random_ci_lo, rho_random_ci_hi):
    labels = ['Fixed\n(v45 vs v45)', 'Real\n(v29 vs v45)',
              'Differential\n(drift genes)', 'Matched\nRandom']
    values = [rho_fixed, rho_real, rho_diff, rho_random_mean]
    colors = [COL_GREY, COL_V29, COL_V45, COL_RED]

    # Error bars: only matched random has CI
    yerr_lo = [0, 0, 0, rho_random_mean - rho_random_ci_lo]
    yerr_hi = [0, 0, 0, rho_random_ci_hi - rho_random_mean]

    for i, (val, col) in enumerate(zip(values, colors)):
        ax.plot(i, val, 'o', color=col, ms=7, zorder=3)
        if yerr_lo[i] > 0 or yerr_hi[i] > 0:
            ax.errorbar(i, val, yerr=[[yerr_lo[i]], [yerr_hi[i]]],
                        fmt='none', ecolor=col, elinewidth=1.0,
                        capsize=4, capthick=0.8, zorder=2)

    ax.axhline(0, color='#999999', ls='--', lw=0.6, zorder=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel('Spearman $\\rho$')
    ax.set_title('(a) Pathway rank correlation')
    ax.set_ylim(-0.15, 1.1)

    # Annotate values
    for i, val in enumerate(values):
        ax.text(i, val + 0.06, f'{val:.3f}', ha='center', va='bottom',
                fontsize=6, color=colors[i])


# ── Panel B: Saturation diagnostic ────────────────────────────────────

def panel_b(ax, sat_data):
    ms = [r[0] for r in sat_data]
    means = [r[1] for r in sat_data]
    sds = [r[2] for r in sat_data]

    # Filter NaN for plotting
    plot_m = []
    plot_mean = []
    plot_sd = []
    for m, mean_j, sd_j in sat_data:
        if not math.isnan(mean_j):
            plot_m.append(m)
            plot_mean.append(mean_j)
            plot_sd.append(sd_j if not math.isnan(sd_j) else 0)

    if plot_m:
        ax.errorbar(plot_m, plot_mean,
                     yerr=plot_sd, fmt='o-', color=COL_V29, ms=4, lw=1.0,
                     elinewidth=0.8, capsize=3, capthick=0.8, zorder=3)

    # Mark NaN points as triangles at y=0
    nan_m = [r[0] for r in sat_data if math.isnan(r[1])]
    if nan_m:
        ax.plot(nan_m, [0] * len(nan_m), 'v', color=COL_GREY, ms=5,
                zorder=2, label='No sig. pathways')

    # Reference lines
    ax.axvline(267, color=COL_V45, ls='--', lw=0.8, alpha=0.7)
    ax.text(267, 0.85, '$|G_{\\Delta}|$=267', fontsize=6, color=COL_V45,
            ha='left', va='top', rotation=90)

    ax.axvline(17000, color=COL_RED, ls='--', lw=0.8, alpha=0.7)
    ax.text(17000, 0.85, 'Full\nforeground', fontsize=5.5, color=COL_RED,
            ha='right', va='top')

    # Annotation about saturation
    ax.annotate('Full-foreground Jaccard ≈ 0.62\n(saturation artifact:\n89% universe coverage)',
                xy=(15000, 0.05), fontsize=5.5, color='#555555',
                ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', fc='#F5F5F5',
                          ec='#CCCCCC', lw=0.5))

    ax.set_xscale('log')
    ax.set_xlabel('Foreground size $m$')
    ax.set_ylabel('Jaccard (FDR < 0.05)')
    ax.set_title('(b) Saturation diagnostic')
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlim(50, 25000)
    if nan_m:
        ax.legend(loc='upper left', frameon=False)


# ── Panel C: Category flip rate ───────────────────────────────────────

def panel_c(ax, cat_data, flip_rate):
    # Sort by v29 count descending
    cat_data_sorted = sorted(cat_data, key=lambda r: r['n_v29'], reverse=True)

    categories = [r['category'] for r in cat_data_sorted]
    n_v29 = [r['n_v29'] for r in cat_data_sorted]
    n_v45 = [r['n_v45'] for r in cat_data_sorted]
    statuses = [r['status'] for r in cat_data_sorted]

    n = len(categories)
    y_pos = range(n)
    bar_h = 0.35

    bars_v29 = ax.barh([y - bar_h / 2 for y in y_pos], n_v29,
                        height=bar_h, color=COL_V29, edgecolor='white',
                        linewidth=0.3, label='v29', zorder=2)
    bars_v45 = ax.barh([y + bar_h / 2 for y in y_pos], n_v45,
                        height=bar_h, color=COL_V45, edgecolor='white',
                        linewidth=0.3, label='v45', zorder=2)

    # Hatch v29-only categories
    for i, (bar, status) in enumerate(zip(bars_v29, statuses)):
        if status == 'v29_only':
            bar.set_hatch('///')
            bar.set_edgecolor(COL_RED)
            bar.set_linewidth(0.6)
            # Mark with star
            ax.text(n_v29[i] + 0.3, i - bar_h / 2, '*',
                    fontsize=9, color=COL_RED, va='center', fontweight='bold')

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(categories)
    ax.invert_yaxis()
    ax.set_xlabel('Significant pathways (FDR < 0.05)')
    ax.set_title('(c) Pathway categories')
    ax.legend(loc='lower right', frameon=False)

    # Flip rate annotation
    ax.text(0.97, 0.03,
            f'Flip rate = {flip_rate:.1%}\n(3 categories v29-only)',
            transform=ax.transAxes, fontsize=6.5, color='#555555',
            ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', fc='#F5F5F5',
                      ec='#CCCCCC', lw=0.5))


# ── Main ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-fixed', required=True)
    parser.add_argument('--baseline-random', required=True)
    parser.add_argument('--diff-enrichment', required=True)
    parser.add_argument('--diff-random', required=True)
    parser.add_argument('--saturation', required=True)
    parser.add_argument('--category-summary', required=True)
    parser.add_argument('--flip-rate', required=True)
    parser.add_argument('--enrichment-ensid', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    set_pub_style()

    # Load data
    fixed = load_metric_value(args.baseline_fixed)
    diff_enrich = load_metric_value(args.diff_enrichment)
    diff_random = load_metric_summary(args.diff_random)
    sat_data = load_saturation(args.saturation)
    cat_data = load_category_summary(args.category_summary)
    flip_data = load_metric_value(args.flip_rate)

    rho_fixed = float(fixed['spearman_rho'])
    rho_real = parse_real_rho(args.enrichment_ensid)
    rho_diff = float(diff_enrich['spearman_rho'])

    rho_rand_row = diff_random['spearman_rho']
    rho_random_mean = float(rho_rand_row['mean'])
    rho_random_ci_lo = float(rho_rand_row['ci_lo'])
    rho_random_ci_hi = float(rho_rand_row['ci_hi'])

    flip_rate_val = float(flip_data['flip_rate'])

    print(f"Panel A values: Fixed={rho_fixed:.3f}, Real={rho_real:.3f}, "
          f"Diff={rho_diff:.3f}, Random={rho_random_mean:.3f} "
          f"[{rho_random_ci_lo:.3f}, {rho_random_ci_hi:.3f}]")
    print(f"Panel C: flip_rate={flip_rate_val:.4f}, {len(cat_data)} categories")

    # Create figure
    fig = plt.figure(figsize=(10, 3.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.4], wspace=0.4)

    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    panel_a(ax_a, rho_fixed, rho_real, rho_diff,
            rho_random_mean, rho_random_ci_lo, rho_random_ci_hi)
    panel_b(ax_b, sat_data)
    panel_c(ax_c, cat_data, flip_rate_val)

    path = os.path.join(args.outdir, 'fig5_calibration.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")

    # Update FIGURE_MANIFEST.txt
    manifest_path = os.path.join(args.outdir, 'FIGURE_MANIFEST.txt')

    # Provenance
    timestamp = datetime.datetime.now().isoformat()
    commit = git_commit_hash()
    script_cmd = ' '.join(sys.argv)
    input_files = [
        args.baseline_fixed, args.baseline_random,
        args.diff_enrichment, args.diff_random,
        args.saturation, args.category_summary,
        args.flip_rate, args.enrichment_ensid,
    ]
    hash_lines = '\n'.join(
        f"  {os.path.basename(p)}: {file_sha256(p)}" for p in input_files)

    entry = (
        "\n"
        + "-" * 70 + "\n"
        "FIGURE 5: CALIBRATION (fig5_calibration.pdf)\n"
        + "-" * 70 + "\n"
        f"Generated: {timestamp}\n"
        f"Git commit: {commit}\n"
        f"Script: {script_cmd}\n"
        f"\n"
        f"Input file SHA-256:\n"
        f"{hash_lines}\n"
        f"\n"
        f"Panel A — Spearman rho across regimes:\n"
        f"  Fixed (v45 vs v45):   {rho_fixed:.6f}\n"
        f"  Real (v29 vs v45):    {rho_real:.6f}\n"
        f"  Differential:         {rho_diff:.6f}\n"
        f"  Matched Random:       {rho_random_mean:.6f} "
        f"[{rho_random_ci_lo:.6f}, {rho_random_ci_hi:.6f}]\n"
        f"\n"
        f"Panel B — Saturation diagnostic:\n"
        f"  Random gene sets produce Jaccard ~ 0 at all foreground sizes.\n"
        f"  Full-foreground Jaccard ~ 0.62 is a saturation artifact (89% coverage).\n"
        f"\n"
        f"Panel C — Category flip rate:\n"
        f"  {len(cat_data)} categories, flip rate = {flip_rate_val:.4f}\n"
        f"  v29-only: Cell-Cell communication, Circadian clock, "
        f"Gene expression (Transcription)\n"
        f"  Shared category rank correlation: "
        f"{float(flip_data.get('spearman_rho_shared', 'nan')):.4f}\n"
    )

    if os.path.exists(manifest_path):
        with open(manifest_path, 'a') as f:
            f.write(entry)
        print(f"Appended to {manifest_path}")
    else:
        with open(manifest_path, 'w') as f:
            f.write("FIGURE_MANIFEST.txt\n" + entry)
        print(f"Created {manifest_path}")


if __name__ == "__main__":
    main()
