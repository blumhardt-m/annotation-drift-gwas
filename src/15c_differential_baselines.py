#!/usr/bin/env python3
"""
Differential-gene baselines and saturation diagnostic (EXPLORATORY).

Replaces the uninformative randomized-GMT null (src/15_baselines.py Task 2)
with a differential-gene ORA that operates in a meaningful regime.

Step 1: Compute G_delta_v29 (v29-only) and G_delta_v45 (v45-only) from
        comparison.tsv Ensembl IDs.
Step 2: ORA on differential gene sets, Jaccard and Spearman rho.
Step 3: Matched random null — sample random gene sets of same size,
        repeat ORA, report Jaccard/rho distribution.
Step 4: Saturation curve — Jaccard vs foreground size to explain why the
        original full-set null was uninformative.

Usage:
    python 15c_differential_baselines.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --outdir results/tables
"""

import argparse
import csv
import importlib.util
import math
import os
import random
import sys


# ── Import ORA machinery from src/12 ──────────────────────────────────

def _import_src12():
    src_dir = os.path.dirname(os.path.abspath(__file__))
    path_12 = os.path.join(src_dir, "12_ensembl_id_enrichment.py")
    spec = importlib.util.spec_from_file_location("src12", path_12)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


src12 = _import_src12()


# ── Shared helpers ─────────────────────────────────────────────────────

def compute_fdr(ora_results):
    """Given list of (name, pval, overlap), return (fdr_map, pval_map)."""
    names = [r[0] for r in ora_results]
    pvals = [r[1] for r in ora_results]
    fdrs = src12.bh_correction(pvals)
    return dict(zip(names, fdrs)), dict(zip(names, pvals))


def jaccard_of_significant(fdr_a, fdr_b, all_names, alpha=0.05):
    sig_a = {n for n in all_names if fdr_a.get(n, 1.0) < alpha}
    sig_b = {n for n in all_names if fdr_b.get(n, 1.0) < alpha}
    inter = sig_a & sig_b
    union = sig_a | sig_b
    return len(inter) / len(union) if union else float('nan')


def spearman_of_pvals(pval_a, pval_b, all_names):
    vals_a = [pval_a.get(n, 1.0) for n in all_names]
    vals_b = [pval_b.get(n, 1.0) for n in all_names]
    return src12.spearman_rank_corr(vals_a, vals_b)


def summary_stats(values):
    n = len(values)
    # Filter out NaN for summary
    finite = [v for v in values if not math.isnan(v)]
    if not finite:
        return float('nan'), float('nan'), float('nan'), float('nan')
    nf = len(finite)
    mean = sum(finite) / nf
    sd = math.sqrt(sum((x - mean) ** 2 for x in finite) / (nf - 1)) if nf > 1 else 0.0
    sorted_vals = sorted(finite)
    ci_lo = sorted_vals[int(0.025 * nf)]
    ci_hi = sorted_vals[min(int(0.975 * nf), nf - 1)]
    return mean, sd, ci_lo, ci_hi


# ── Step 1: Differential gene sets ────────────────────────────────────

def load_differential_sets(comparison_path):
    """Return (G_delta_v29, G_delta_v45, universe) as Ensembl ID sets."""
    all_v29 = set()
    all_v45 = set()
    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            all_v29.add(row['gene_v29'])
            all_v45.add(row['gene_v45'])
    universe = all_v29 | all_v45
    g_delta_v29 = all_v29 - all_v45  # v29-only
    g_delta_v45 = all_v45 - all_v29  # v45-only
    return g_delta_v29, g_delta_v45, universe


# ── Step 2: ORA on differential gene sets ─────────────────────────────

def run_differential_ora(g_delta_v29, g_delta_v45, universe, pathways, outpath):
    print("\nStep 2: ORA on differential gene sets...")

    ora_v29 = src12.run_ora(g_delta_v29, universe, pathways)
    ora_v45 = src12.run_ora(g_delta_v45, universe, pathways)

    fdr_v29, pval_v29 = compute_fdr(ora_v29)
    fdr_v45, pval_v45 = compute_fdr(ora_v45)
    all_names = sorted(set(fdr_v29) | set(fdr_v45))

    sig_v29 = {n for n in all_names if fdr_v29.get(n, 1.0) < 0.05}
    sig_v45 = {n for n in all_names if fdr_v45.get(n, 1.0) < 0.05}
    jac = jaccard_of_significant(fdr_v29, fdr_v45, all_names)
    rho = spearman_of_pvals(pval_v29, pval_v45, all_names)

    print(f"  Sig pathways (FDR<0.05): v29-only={len(sig_v29)}, v45-only={len(sig_v45)}")
    print(f"  Jaccard = {jac:.6f}")
    print(f"  Spearman rho = {rho:.6f}")

    with open(outpath, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"n_sig_v29_delta\t{len(sig_v29)}\n")
        f.write(f"n_sig_v45_delta\t{len(sig_v45)}\n")
        f.write(f"jaccard_fdr005\t{jac:.6f}\n")
        f.write(f"spearman_rho\t{rho:.6f}\n")

    print(f"  Saved {outpath}")
    return jac, rho


# ── Step 3: Matched random null ───────────────────────────────────────

def run_matched_random_null(size_v29, size_v45, universe, pathways,
                            n_perms, outpath, seed=42):
    print(f"\nStep 3: Matched random null ({n_perms} permutations, "
          f"sizes {size_v29}/{size_v45})...")
    rng = random.Random(seed)
    universe_list = sorted(universe)  # sorted for reproducibility

    jaccards = []
    rhos = []

    for i in range(n_perms):
        if (i + 1) % 10 == 0:
            print(f"  Permutation {i + 1}/{n_perms}...")

        fake_v29 = set(rng.sample(universe_list, size_v29))
        fake_v45 = set(rng.sample(universe_list, size_v45))

        ora_v29 = src12.run_ora(fake_v29, universe, pathways)
        ora_v45 = src12.run_ora(fake_v45, universe, pathways)

        fdr_v29, pval_v29 = compute_fdr(ora_v29)
        fdr_v45, pval_v45 = compute_fdr(ora_v45)
        all_names = sorted(set(fdr_v29) | set(fdr_v45))

        jac = jaccard_of_significant(fdr_v29, fdr_v45, all_names)
        rho = spearman_of_pvals(pval_v29, pval_v45, all_names)

        jaccards.append(jac)
        rhos.append(rho)

    jac_mean, jac_sd, jac_lo, jac_hi = summary_stats(jaccards)
    rho_mean, rho_sd, rho_lo, rho_hi = summary_stats(rhos)

    with open(outpath, 'w') as f:
        f.write("metric\tmean\tsd\tci_lo\tci_hi\tn_permutations\n")
        f.write(f"jaccard_fdr005\t{jac_mean:.6f}\t{jac_sd:.6f}\t"
                f"{jac_lo:.6f}\t{jac_hi:.6f}\t{n_perms}\n")
        f.write(f"spearman_rho\t{rho_mean:.6f}\t{rho_sd:.6f}\t"
                f"{rho_lo:.6f}\t{rho_hi:.6f}\t{n_perms}\n")

    print(f"  Jaccard: mean={jac_mean:.4f}, SD={jac_sd:.4f}, "
          f"95% CI=[{jac_lo:.4f}, {jac_hi:.4f}]")
    print(f"  Spearman rho: mean={rho_mean:.4f}, SD={rho_sd:.4f}")
    print(f"  Saved {outpath}")


# ── Step 4: Saturation diagnostic ─────────────────────────────────────

def run_saturation_curve(universe, pathways, outpath, seed=42):
    print("\nStep 4: Saturation diagnostic...")
    rng = random.Random(seed)
    universe_list = sorted(universe)
    N = len(universe_list)

    m_values = [100, 250, 500, 1000, 2000, 5000, 10000, 15000]
    # Skip values larger than universe
    m_values = [m for m in m_values if m <= N]
    n_reps = 20

    rows = []
    for m in m_values:
        jaccards = []
        for _ in range(n_reps):
            set_a = set(rng.sample(universe_list, m))
            set_b = set(rng.sample(universe_list, m))

            ora_a = src12.run_ora(set_a, universe, pathways)
            ora_b = src12.run_ora(set_b, universe, pathways)

            fdr_a, pval_a = compute_fdr(ora_a)
            fdr_b, pval_b = compute_fdr(ora_b)
            all_names = sorted(set(fdr_a) | set(fdr_b))

            jac = jaccard_of_significant(fdr_a, fdr_b, all_names)
            jaccards.append(jac)

        finite = [j for j in jaccards if not math.isnan(j)]
        if finite:
            mean_j = sum(finite) / len(finite)
            sd_j = (math.sqrt(sum((x - mean_j) ** 2 for x in finite) / (len(finite) - 1))
                    if len(finite) > 1 else 0.0)
        else:
            mean_j = float('nan')
            sd_j = float('nan')

        rows.append((m, mean_j, sd_j))
        print(f"  m={m:>6d}: mean Jaccard = {mean_j:.4f} (SD={sd_j:.4f})")

    with open(outpath, 'w') as f:
        f.write("m\tmean_jaccard\tsd_jaccard\n")
        for m, mean_j, sd_j in rows:
            f.write(f"{m}\t{mean_j:.6f}\t{sd_j:.6f}\n")

    print(f"  Saved {outpath}")


# ── Main ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Differential-gene baselines and saturation diagnostic")
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--n-permutations', type=int, default=100)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Load shared data ──
    print("Loading data...")
    symbol_to_ensembl, ensembl_universe = src12.load_tss_tables(
        args.tss_v29, args.tss_v45)
    pathways, gmt_stats, _ = src12.load_gmt_as_ensembl(
        args.gmt, symbol_to_ensembl, args.min_size, args.max_size)

    print(f"  Ensembl universe: {len(ensembl_universe)} genes")
    print(f"  Pathways: {gmt_stats['n_pathways']} (size {args.min_size}-{args.max_size})")

    # ── Step 1: Differential gene sets ──
    print("\nStep 1: Computing differential gene sets...")
    g_delta_v29, g_delta_v45, mapping_universe = load_differential_sets(args.comparison)

    # Use the TSS-derived Ensembl universe (consistent with src/12)
    universe = ensembl_universe

    print(f"  |G_delta_v29| (v29-only): {len(g_delta_v29)} "
          f"({100 * len(g_delta_v29) / len(universe):.2f}% of universe)")
    print(f"  |G_delta_v45| (v45-only): {len(g_delta_v45)} "
          f"({100 * len(g_delta_v45) / len(universe):.2f}% of universe)")
    print(f"  |Universe|: {len(universe)}")

    # ── Step 2 ──
    run_differential_ora(
        g_delta_v29, g_delta_v45, universe, pathways,
        os.path.join(args.outdir, "differential_enrichment.tsv"))

    # ── Step 3 ──
    run_matched_random_null(
        len(g_delta_v29), len(g_delta_v45), universe, pathways,
        args.n_permutations,
        os.path.join(args.outdir, "differential_baseline_random.tsv"),
        args.seed)

    # ── Step 4 ──
    run_saturation_curve(
        universe, pathways,
        os.path.join(args.outdir, "saturation_curve.tsv"),
        args.seed)

    print("\nAll differential baselines complete.")


if __name__ == "__main__":
    main()
