#!/usr/bin/env python3
"""
Baseline calibrations for the enrichment analysis (EXPLORATORY).

Task 1 — Fixed-annotation baseline (algorithmic floor):
  Run ORA with v45 gene set on both sides → Jaccard ≈ 1.0, ρ ≈ 1.0.

Task 2 — Randomized-annotation negative control:
  Permute gene→pathway memberships in reactome.gmt, rerun enrichment
  under each permutation → Jaccard ≈ 0.

Usage:
    python 15_baselines.py \
        --comparison data/processed/mappings/comparison.tsv \
        --tss-v29 data/processed/tss_tables/tss_v29.tsv \
        --tss-v45 data/processed/tss_tables/tss_v45.tsv \
        --gmt data/raw/reactome.gmt \
        --outdir results/tables
"""

import argparse
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
    return len(inter) / len(union) if union else 1.0


def spearman_of_pvals(pval_a, pval_b, all_names):
    vals_a = [pval_a.get(n, 1.0) for n in all_names]
    vals_b = [pval_b.get(n, 1.0) for n in all_names]
    return src12.spearman_rank_corr(vals_a, vals_b)


def top_k_overlap(pval_a, pval_b, all_names, k):
    ranked_a = sorted(all_names, key=lambda n: pval_a.get(n, 1.0))[:k]
    ranked_b = sorted(all_names, key=lambda n: pval_b.get(n, 1.0))[:k]
    return len(set(ranked_a) & set(ranked_b)) / k


# ── Task 1: Fixed-annotation baseline ─────────────────────────────────

def run_fixed_baseline(v45_genes, universe, pathways, outpath):
    print("Running fixed-annotation baseline (v45 vs v45)...")

    ora_a = src12.run_ora(v45_genes, universe, pathways)
    ora_b = src12.run_ora(v45_genes, universe, pathways)

    fdr_a, pval_a = compute_fdr(ora_a)
    fdr_b, pval_b = compute_fdr(ora_b)
    all_names = sorted(set(fdr_a) | set(fdr_b))

    jac = jaccard_of_significant(fdr_a, fdr_b, all_names)
    rho = spearman_of_pvals(pval_a, pval_b, all_names)

    rows = [
        ("jaccard_fdr005", f"{jac:.6f}"),
        ("spearman_rho", f"{rho:.6f}"),
    ]
    for k in [10, 20, 50, 100]:
        overlap = top_k_overlap(pval_a, pval_b, all_names, k)
        rows.append((f"top_k_overlap_K{k}", f"{overlap:.6f}"))

    # Validation gate
    if jac < 0.999:
        print(f"  FATAL: Fixed baseline Jaccard = {jac:.6f} != 1.0 — algorithmic bug.")
        sys.exit(1)

    with open(outpath, 'w') as f:
        f.write("metric\tvalue\n")
        for metric, value in rows:
            f.write(f"{metric}\t{value}\n")

    print(f"  Jaccard = {jac:.6f}, rho = {rho:.6f}")
    print(f"  Saved {outpath}")


# ── Task 2: Randomized-annotation negative control ────────────────────

def permute_pathways(pathways, all_genes_list, rng):
    """Return pathways with gene members replaced by random samples of the same size."""
    permuted = {}
    for name, genes in pathways.items():
        permuted[name] = set(rng.sample(all_genes_list, len(genes)))
    return permuted


def run_random_baseline(v29_genes, v45_genes, universe, pathways,
                        n_perms, outpath, seed=42):
    print(f"Running randomized-annotation baseline ({n_perms} permutations)...")
    rng = random.Random(seed)

    # Gene pool: union of all genes appearing in any pathway (intersected with universe)
    all_gmt_genes = set()
    for pw_genes in pathways.values():
        all_gmt_genes |= (pw_genes & universe)
    all_gmt_genes_list = sorted(all_gmt_genes)  # sorted for reproducibility

    jaccards = []
    rhos = []

    for i in range(n_perms):
        if (i + 1) % 10 == 0:
            print(f"  Permutation {i + 1}/{n_perms}...")

        perm_pw = permute_pathways(pathways, all_gmt_genes_list, rng)

        ora_v29 = src12.run_ora(v29_genes, universe, perm_pw)
        ora_v45 = src12.run_ora(v45_genes, universe, perm_pw)

        fdr_v29, pval_v29 = compute_fdr(ora_v29)
        fdr_v45, pval_v45 = compute_fdr(ora_v45)
        all_names = sorted(set(fdr_v29) | set(fdr_v45))

        jac = jaccard_of_significant(fdr_v29, fdr_v45, all_names)
        rho = spearman_of_pvals(pval_v29, pval_v45, all_names)

        jaccards.append(jac)
        rhos.append(rho)

    def summary(values):
        n = len(values)
        mean = sum(values) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in values) / (n - 1)) if n > 1 else 0.0
        sorted_vals = sorted(values)
        ci_lo = sorted_vals[int(0.025 * n)]
        ci_hi = sorted_vals[min(int(0.975 * n), n - 1)]
        return mean, sd, ci_lo, ci_hi

    jac_mean, jac_sd, jac_lo, jac_hi = summary(jaccards)
    rho_mean, rho_sd, rho_lo, rho_hi = summary(rhos)

    if jac_mean > 0.3:
        print(f"  WARNING: Mean random Jaccard = {jac_mean:.4f} > 0.3 — "
              "permutation may not be shuffling correctly.")

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


# ── Main ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fixed and randomized annotation baselines")
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

    print("Loading data...")
    symbol_to_ensembl, ensembl_universe = src12.load_tss_tables(
        args.tss_v29, args.tss_v45)
    pathways, gmt_stats, _ = src12.load_gmt_as_ensembl(
        args.gmt, symbol_to_ensembl, args.min_size, args.max_size)
    v29_genes, v45_genes, _ = src12.load_comparison_ensembl(args.comparison)

    print(f"  Universe: {len(ensembl_universe)} genes")
    print(f"  Pathways: {gmt_stats['n_pathways']} (size {args.min_size}-{args.max_size})")
    print(f"  v29 genes: {len(v29_genes)}, v45 genes: {len(v45_genes)}")

    # Task 1
    run_fixed_baseline(
        v45_genes, ensembl_universe, pathways,
        os.path.join(args.outdir, "baseline_fixed.tsv"))

    # Task 2
    run_random_baseline(
        v29_genes, v45_genes, ensembl_universe, pathways,
        args.n_permutations,
        os.path.join(args.outdir, "baseline_random.tsv"),
        args.seed)

    print("\nAll baselines complete.")


if __name__ == "__main__":
    main()
