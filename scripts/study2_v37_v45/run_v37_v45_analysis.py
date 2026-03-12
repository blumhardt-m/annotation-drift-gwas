#!/usr/bin/env python3
"""
Study 2: GENCODE v37 vs v45 annotation drift analysis.

PURPOSE
-------
Replicates the primary v29-v45 pipeline with GENCODE v37 substituted for v29.
Tests whether annotation drift is specific to the v29-v45 pair or persists
across a more recent version interval (v37-v45, spanning 4 years vs 7 years).

INPUTS (project-relative)
-------------------------
  data/raw/gencode.v37.annotation.gtf.gz   — GENCODE v37 GTF
  data/raw/gencode.v45.annotation.gtf.gz   — GENCODE v45 GTF (via tss_v45.tsv)
  data/processed/tss_tables/tss_v45.tsv    — v45 protein-coding TSS table
  data/processed/gwas_eligible.tsv         — GWAS catalog eligible loci
  data/raw/reactome.gmt                    — Reactome pathway definitions

  (auto-generated on first run)
  data/processed/tss_tables/tss_v37.tsv    — v37 protein-coding TSS table

OUTPUTS (written to results/study2_v37_v45/)
---------------------------------------------
  locus_mapping_v37.tsv        — per-locus gene assignment for v37
  locus_mapping_v45.tsv        — per-locus gene assignment for v45
  instability_summary.tsv      — instability rate, CI, decomposition
  reactome_enrichment_v37.tsv  — per-pathway enrichment for v37
  reactome_enrichment_v45.tsv  — per-pathway enrichment for v45
  rank_scatter.png             — scatter of -log10(p) v37 vs v45
  study2_summary.tsv           — side-by-side comparison vs primary analysis

METHOD
------
  - Parse GENCODE v37 GTF to TSS table (protein-coding, autosomes)
  - Map all GWAS loci to nearest gene within 500kb in both v37 and v45
  - Bootstrap 95% CI on instability rate (10,000 resamples)
  - Decompose reassignments: ID-only vs ID+symbol vs symbol-only
  - Ensembl-ID-based Reactome ORA (Fisher's exact, BH FDR, alpha=0.05)
  - Compute Spearman rho and Jaccard overlap

USAGE
-----
  cd <project_root>
  python scripts/study2_v37_v45/run_v37_v45_analysis.py
"""

import csv
import gzip
import math
import os
import sys
from bisect import bisect_left
from collections import defaultdict

import numpy as np


# ── Path resolution ──────────────────────────────────────────────────────
# Resolve project root: two levels up from this script's location
BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

DATA_RAW  = os.path.join(BASE, "data", "raw")
DATA_PROC = os.path.join(BASE, "data", "processed")
OUT_DIR   = os.path.join(BASE, "results", "study2_v37_v45")

GTF_V37 = os.path.join(DATA_RAW, "gencode.v37.annotation.gtf.gz")
TSS_V37 = os.path.join(DATA_PROC, "tss_tables", "tss_v37.tsv")
TSS_V45 = os.path.join(DATA_PROC, "tss_tables", "tss_v45.tsv")
GWAS    = os.path.join(DATA_PROC, "gwas_eligible.tsv")
GMT     = os.path.join(DATA_RAW, "reactome.gmt")

# Primary analysis reference values (v29-v45) for comparison
PRIMARY = {
    "instability_rate": 0.0415,
    "ci_low": 0.0410,
    "ci_high": 0.0419,
    "spearman_rho": 0.922857,
    "jaccard_005": 0.6133,
    "sig_v29": 74,
    "sig_v45": 47,
    "shared": 46,
}

AUTOSOMES   = frozenset(str(i) for i in range(1, 23))
WINDOW_BP   = 500_000
MIN_PW_SIZE = 15
MAX_PW_SIZE = 500
N_BOOTSTRAP = 10_000
SEED        = 42


# ── GTF parsing (replicates src/02_parse_gencode.py) ─────────────────────

def parse_attributes(attr_string):
    """Parse GTF attribute string into key-value dict."""
    attrs = {}
    for attr in attr_string.strip().split(';'):
        attr = attr.strip()
        if not attr or ' ' not in attr:
            continue
        key, value = attr.split(' ', 1)
        attrs[key] = value.strip('"')
    return attrs


def parse_gtf_to_tss(gtf_path, output_path):
    """Extract one TSS per protein-coding gene from GENCODE GTF.

    Filters to autosomes and protein_coding gene_type.
    For each gene, takes the most 5' TSS across all transcripts.
    """
    genes = defaultdict(lambda: {
        "tss_positions": [], "symbol": None, "strand": None, "chrom": None
    })

    opener = gzip.open if gtf_path.endswith('.gz') else open
    print(f"  Parsing {os.path.basename(gtf_path)}...")

    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            chrom, _, _, start, end, _, strand, _, attributes = fields
            attrs = parse_attributes(attributes)

            gene_id = attrs.get('gene_id', '').split('.')[0]
            gene_type = attrs.get('gene_type', '')

            if gene_type != 'protein_coding':
                continue

            chrom_clean = chrom.replace('chr', '')
            if chrom_clean not in AUTOSOMES:
                continue

            tss = int(start) if strand == '+' else int(end)
            genes[gene_id]["tss_positions"].append(tss)
            genes[gene_id]["symbol"] = attrs.get('gene_name', '')
            genes[gene_id]["strand"] = strand
            genes[gene_id]["chrom"] = chrom_clean

    output_rows = []
    for gene_id, data in genes.items():
        if not data["tss_positions"]:
            continue
        tss = (min(data["tss_positions"]) if data["strand"] == '+'
               else max(data["tss_positions"]))
        output_rows.append({
            'ensembl_id': gene_id,
            'gene_symbol': data["symbol"],
            'chromosome': data["chrom"],
            'tss': tss,
            'strand': data["strand"],
        })

    output_rows.sort(key=lambda r: (int(r['chromosome']), r['tss']))

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=['ensembl_id', 'gene_symbol', 'chromosome',
                        'tss', 'strand'],
            delimiter='\t',
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"    → {len(output_rows)} protein-coding genes")
    return len(output_rows)


# ── TSS loading & variant mapping (replicates src/04_map_variants.py) ────

def load_tss_table(path):
    """Load TSS table into per-chromosome sorted lists."""
    tss_by_chrom = defaultdict(list)
    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            tss_by_chrom[row['chromosome']].append(
                (int(row['tss']), row['ensembl_id'], row['gene_symbol'])
            )
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort(key=lambda x: x[0])
    return tss_by_chrom


def find_nearest_gene(chrom, pos, tss_by_chrom):
    """Find nearest protein-coding gene TSS within WINDOW_BP."""
    if chrom not in tss_by_chrom:
        return None
    genes = tss_by_chrom[chrom]
    positions = [g[0] for g in genes]
    idx = bisect_left(positions, pos)

    candidates = []
    for i in range(max(0, idx - 1), min(len(genes), idx + 2)):
        dist = abs(genes[i][0] - pos)
        if dist <= WINDOW_BP:
            candidates.append((dist, genes[i][1], genes[i][2]))

    if not candidates:
        return None
    candidates.sort(key=lambda x: (x[0], x[1]))
    return (candidates[0][1], candidates[0][2], candidates[0][0])


def map_variants(tss_v37, tss_v45):
    """Map all GWAS loci to nearest gene in v37 and v45."""
    print("  Mapping variants...")
    results = []
    n_both = n_differ = 0

    with open(GWAS, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            chrom = row['CHR_ID']
            pos = int(row['CHR_POS'])
            snp = row['SNPS']

            g37 = find_nearest_gene(chrom, pos, tss_v37)
            g45 = find_nearest_gene(chrom, pos, tss_v45)

            if g37 and g45:
                n_both += 1
                differs = g37[0] != g45[0]
                if differs:
                    n_differ += 1

                results.append({
                    'snp': snp, 'chrom': chrom, 'pos': pos,
                    'gene_v37': g37[0], 'symbol_v37': g37[1], 'dist_v37': g37[2],
                    'gene_v45': g45[0], 'symbol_v45': g45[1], 'dist_v45': g45[2],
                    'differs': 1 if differs else 0
                })

    return results, n_both, n_differ


# ── Bootstrap CI (replicates src/05_bootstrap.py) ────────────────────────

def bootstrap_ci(differs_arr):
    """Compute 95% bootstrap CI for the instability rate."""
    np.random.seed(SEED)
    n = len(differs_arr)
    point_est = differs_arr.mean()

    boot_rates = np.array([
        np.random.choice(differs_arr, size=n, replace=True).mean()
        for _ in range(N_BOOTSTRAP)
    ])

    ci_low = np.percentile(boot_rates, 2.5)
    ci_high = np.percentile(boot_rates, 97.5)
    return point_est, ci_low, ci_high


# ── Gene reassignment decomposition ─────────────────────────────────────

def decompose_reassignments(loci):
    """Classify differing loci by Ensembl ID and symbol changes."""
    n_total = len(loci)
    n_ensembl_differs = 0
    n_symbol_differs = 0
    n_id_only = 0
    n_symbol_only = 0
    n_both_differ = 0

    for l in loci:
        eid_diff = l['gene_v37'] != l['gene_v45']
        sym_diff = l['symbol_v37'] != l['symbol_v45']

        if eid_diff:
            n_ensembl_differs += 1
        if sym_diff:
            n_symbol_differs += 1
        if eid_diff and not sym_diff:
            n_id_only += 1
        if sym_diff and not eid_diff:
            n_symbol_only += 1
        if eid_diff and sym_diff:
            n_both_differ += 1

    return {
        'n_total': n_total,
        'n_ensembl_differs': n_ensembl_differs,
        'n_symbol_differs': n_symbol_differs,
        'n_id_only': n_id_only,
        'n_symbol_only': n_symbol_only,
        'n_both_differ': n_both_differ,
    }


# ── Enrichment (replicates src/12_ensembl_id_enrichment.py) ──────────────

def build_alias_table():
    """HGNC alias table for MARCH/MARCHF, SEPT/SEPTIN renames."""
    aliases = {}
    for i in range(1, 12):
        aliases[f"MARCH{i}"] = f"MARCHF{i}"
        aliases[f"MARCHF{i}"] = f"MARCH{i}"
    for i in range(1, 15):
        aliases[f"SEPT{i}"] = f"SEPTIN{i}"
        aliases[f"SEPTIN{i}"] = f"SEPT{i}"
    return aliases


def load_tss_for_enrichment(v37_path, v45_path):
    """Build symbol→Ensembl mapping and Ensembl universe from TSS tables."""
    alias_table = build_alias_table()
    symbol_to_ensembl = {}
    ensembl_universe = set()

    for path in [v37_path, v45_path]:
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
    """Load Reactome GMT, converting symbols to Ensembl IDs."""
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
            if MIN_PW_SIZE <= len(ensembl_genes) <= MAX_PW_SIZE:
                pathways[name] = {'id': pw_id, 'genes': ensembl_genes}

    return pathways


def fisher_exact_pvalue(a, b, c, d):
    """One-tailed Fisher exact test p-value (overrepresentation)."""
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
    """Benjamini-Hochberg FDR correction."""
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
    """Over-representation analysis for Reactome pathways."""
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
            'pathway_name': name,
            'pathway_id': pw_data['id'],
            'p_value': pval,
            'gene_count': a,
        })
    return results


def spearman_rank_corr(x, y):
    """Compute Spearman rank correlation between two lists."""
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


# ── Figure: rank scatter ─────────────────────────────────────────────────

def make_rank_scatter(pval_v37, pval_v45, all_names, rho, jaccard, out_path):
    """Generate scatter plot of -log10(p) v37 vs v45."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib not available, skipping scatter plot")
        return

    x = [-math.log10(max(pval_v37.get(n, 1.0), 1e-300)) for n in all_names]
    y = [-math.log10(max(pval_v45.get(n, 1.0), 1e-300)) for n in all_names]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x, y, s=8, alpha=0.4, color='#2c7bb6', edgecolors='none')

    lim = max(max(x), max(y)) * 1.05
    ax.plot([0, lim], [0, lim], 'k--', alpha=0.3, linewidth=0.8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel('-log10(p) GENCODE v37', fontsize=12)
    ax.set_ylabel('-log10(p) GENCODE v45', fontsize=12)
    ax.set_title('Study 2: Reactome pathway enrichment\nGENCODE v37 vs v45',
                 fontsize=13, fontweight='bold')
    ax.text(0.05, 0.92, f'Spearman rho = {rho:.4f}\nJaccard = {jaccard:.4f}',
            transform=ax.transAxes, fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"    → {out_path}")


# ── Output writers ───────────────────────────────────────────────────────

def write_locus_mapping(loci, out_path, version_tag):
    """Write locus mapping file for one version."""
    fieldnames = ['locus_id', 'rsid', 'ensembl_id', 'gene_symbol']
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for l in loci:
            gene_key = f'gene_{version_tag}'
            sym_key = f'symbol_{version_tag}'
            writer.writerow({
                'locus_id': f"{l['chrom']}:{l['pos']}",
                'rsid': l['snp'],
                'ensembl_id': l[gene_key],
                'gene_symbol': l[sym_key],
            })


def write_enrichment_tsv(ora_results, fdr_values, out_path):
    """Write per-pathway enrichment results as TSV."""
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


# ══════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 70)
    print("STUDY 2: GENCODE v37 vs v45 ANNOTATION DRIFT ANALYSIS")
    print("=" * 70)

    # ── Step 1: Parse v37 GTF ─────────────────────────────────────────────
    print("\n[1/7] PARSE GENCODE v37 GTF")
    if os.path.exists(TSS_V37):
        print(f"  TSS table already exists: {TSS_V37}")
        n_v37 = sum(1 for _ in open(TSS_V37)) - 1
    else:
        n_v37 = parse_gtf_to_tss(GTF_V37, TSS_V37)

    n_v45 = sum(1 for _ in open(TSS_V45)) - 1

    print(f"\n  v37 protein-coding genes: {n_v37:>6,}")
    print(f"  v45 protein-coding genes: {n_v45:>6,}")
    print(f"  v45 matches primary:      {'YES' if n_v45 == 19114 else 'CHECK'}")

    # ── Step 2: Map variants ──────────────────────────────────────────────
    print("\n[2/7] MAP GWAS VARIANTS")
    print("  Loading TSS tables...")
    tss_v37 = load_tss_table(TSS_V37)
    tss_v45 = load_tss_table(TSS_V45)

    loci, n_both, n_differ = map_variants(tss_v37, tss_v45)
    rate = n_differ / n_both if n_both > 0 else 0

    print(f"    Mapped in both: {n_both:,}")
    print(f"    Gene differs:   {n_differ:,}")
    print(f"    Instability:    {rate:.4f} ({rate*100:.2f}%)")

    write_locus_mapping(loci, os.path.join(OUT_DIR, "locus_mapping_v37.tsv"), "v37")
    write_locus_mapping(loci, os.path.join(OUT_DIR, "locus_mapping_v45.tsv"), "v45")
    print(f"    → locus_mapping_v37.tsv, locus_mapping_v45.tsv")

    # ── Step 3: Bootstrap CI ──────────────────────────────────────────────
    print("\n[3/7] BOOTSTRAP CI")
    differs_arr = np.array([l['differs'] for l in loci])
    point_est, ci_low, ci_high = bootstrap_ci(differs_arr)
    print(f"    Instability: {point_est:.4f} [{ci_low:.4f}, {ci_high:.4f}]")

    # ── Step 4: Gene reassignment decomposition ──────────────────────────
    print("\n[4/7] GENE REASSIGNMENT DECOMPOSITION")
    decomp = decompose_reassignments(loci)
    pct_id_only = (100 * decomp['n_id_only'] / decomp['n_ensembl_differs']
                   if decomp['n_ensembl_differs'] > 0 else 0)
    pct_both = (100 * decomp['n_both_differ'] / decomp['n_ensembl_differs']
                if decomp['n_ensembl_differs'] > 0 else 0)
    pct_sym_only = (100 * decomp['n_symbol_only'] / decomp['n_total']
                    if decomp['n_total'] > 0 else 0)

    print(f"    Total loci:              {decomp['n_total']:,}")
    print(f"    Ensembl ID differs:      {decomp['n_ensembl_differs']:,} "
          f"({100*decomp['n_ensembl_differs']/decomp['n_total']:.2f}%)")
    print(f"    Symbol differs:          {decomp['n_symbol_differs']:,} "
          f"({100*decomp['n_symbol_differs']/decomp['n_total']:.2f}%)")
    print(f"    ID-only change:          {decomp['n_id_only']:,} "
          f"({pct_id_only:.1f}% of diffs)")
    print(f"    Symbol-only rename:      {decomp['n_symbol_only']:,} "
          f"({pct_sym_only:.2f}% of total)")
    print(f"    Both ID+symbol change:   {decomp['n_both_differ']:,} "
          f"({pct_both:.1f}% of diffs)")

    inst_path = os.path.join(OUT_DIR, "instability_summary.tsv")
    with open(inst_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'instability_rate', 'ci_low', 'ci_high',
            'identifier_only_pct', 'true_mapping_pct',
            'n_total', 'n_differ', 'n_id_only', 'n_both_differ',
            'n_symbol_only'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerow({
            'instability_rate': f"{point_est:.6f}",
            'ci_low': f"{ci_low:.6f}",
            'ci_high': f"{ci_high:.6f}",
            'identifier_only_pct': f"{pct_id_only:.2f}",
            'true_mapping_pct': f"{pct_both:.2f}",
            'n_total': decomp['n_total'],
            'n_differ': decomp['n_ensembl_differs'],
            'n_id_only': decomp['n_id_only'],
            'n_both_differ': decomp['n_both_differ'],
            'n_symbol_only': decomp['n_symbol_only'],
        })
    print(f"    → instability_summary.tsv")

    # ── Step 5: Enrichment analysis ───────────────────────────────────────
    print("\n[5/7] ENRICHMENT ANALYSIS (Ensembl-ID-based ORA)")
    print("  Building gene universe...")
    symbol_to_ensembl, ensembl_universe = load_tss_for_enrichment(TSS_V37, TSS_V45)
    pathways = load_gmt_as_ensembl(symbol_to_ensembl)

    n_pathways = len(pathways)
    n_universe = len(ensembl_universe)

    print(f"  Universe size: {n_universe:,}")
    print(f"  Pathways:      {n_pathways:,}")

    v37_genes = set(l['gene_v37'] for l in loci)
    v45_genes = set(l['gene_v45'] for l in loci)

    print(f"  v37 unique genes: {len(v37_genes):,}")
    print(f"  v45 unique genes: {len(v45_genes):,}")

    print("  Running ORA (v37)...")
    ora_v37 = run_ora(v37_genes, ensembl_universe, pathways)
    print("  Running ORA (v45)...")
    ora_v45 = run_ora(v45_genes, ensembl_universe, pathways)

    pvals_v37 = [r['p_value'] for r in ora_v37]
    fdr_v37 = bh_correction(pvals_v37)
    pvals_v45 = [r['p_value'] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)

    write_enrichment_tsv(ora_v37, fdr_v37,
                         os.path.join(OUT_DIR, "reactome_enrichment_v37.tsv"))
    write_enrichment_tsv(ora_v45, fdr_v45,
                         os.path.join(OUT_DIR, "reactome_enrichment_v45.tsv"))
    print(f"    → reactome_enrichment_v37.tsv, reactome_enrichment_v45.tsv")

    # ── Step 6: Jaccard & Spearman ────────────────────────────────────────
    print("\n[6/7] JACCARD OVERLAP & SPEARMAN CORRELATION")

    fdr_map_v37 = {r['pathway_name']: fdr for r, fdr in zip(ora_v37, fdr_v37)}
    fdr_map_v45 = {r['pathway_name']: fdr for r, fdr in zip(ora_v45, fdr_v45)}
    pval_map_v37 = {r['pathway_name']: r['p_value'] for r in ora_v37}
    pval_map_v45 = {r['pathway_name']: r['p_value'] for r in ora_v45}

    all_names = sorted(set(r['pathway_name'] for r in ora_v37) |
                       set(r['pathway_name'] for r in ora_v45))

    sig_v37 = {n for n in all_names if fdr_map_v37.get(n, 1.0) < 0.05}
    sig_v45 = {n for n in all_names if fdr_map_v45.get(n, 1.0) < 0.05}
    inter = sig_v37 & sig_v45
    union = sig_v37 | sig_v45
    jaccard = len(inter) / len(union) if len(union) > 0 else 1.0

    print(f"    v37 significant (FDR<0.05): {len(sig_v37)}")
    print(f"    v45 significant (FDR<0.05): {len(sig_v45)}")
    print(f"    Shared:                     {len(inter)}")
    print(f"    Union:                      {len(union)}")
    print(f"    Jaccard:                    {jaccard:.4f}")

    log10p_v37 = [-math.log10(max(pval_map_v37.get(n, 1.0), 1e-300))
                  for n in all_names]
    log10p_v45 = [-math.log10(max(pval_map_v45.get(n, 1.0), 1e-300))
                  for n in all_names]
    rho = spearman_rank_corr(log10p_v37, log10p_v45)
    print(f"    Spearman rho:               {rho:.6f}")

    # ── Step 7: Scatter plot ──────────────────────────────────────────────
    print("\n[7/7] RANK SCATTER PLOT")
    scatter_path = os.path.join(OUT_DIR, "rank_scatter.png")
    make_rank_scatter(pval_map_v37, pval_map_v45, all_names, rho, jaccard,
                      scatter_path)

    # ── Write study2_summary.tsv ──────────────────────────────────────────
    summary_path = os.path.join(OUT_DIR, "study2_summary.tsv")
    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['metric', 'study2_v37_v45', 'primary_v29_v45'])
        writer.writerow(['instability_rate', f"{point_est:.6f}",
                         f"{PRIMARY['instability_rate']:.6f}"])
        writer.writerow(['ci_low', f"{ci_low:.6f}",
                         f"{PRIMARY['ci_low']:.6f}"])
        writer.writerow(['ci_high', f"{ci_high:.6f}",
                         f"{PRIMARY['ci_high']:.6f}"])
        writer.writerow(['spearman_rho', f"{rho:.6f}",
                         f"{PRIMARY['spearman_rho']:.6f}"])
        writer.writerow(['jaccard_fdr005', f"{jaccard:.4f}",
                         f"{PRIMARY['jaccard_005']:.4f}"])
        writer.writerow(['sig_older_version', len(sig_v37),
                         PRIMARY['sig_v29']])
        writer.writerow(['sig_v45', len(sig_v45), PRIMARY['sig_v45']])
        writer.writerow(['shared_significant', len(inter), PRIMARY['shared']])
        writer.writerow(['n_loci_mapped', decomp['n_total'], 721406])
        writer.writerow(['n_ensembl_differs', decomp['n_ensembl_differs'],
                         29926])
        writer.writerow(['n_pathways_tested', n_pathways, 1369])
        writer.writerow(['universe_size', n_universe, 19315])
        writer.writerow(['identifier_only_pct', f"{pct_id_only:.2f}", "N/A"])
        writer.writerow(['true_mapping_pct', f"{pct_both:.2f}", "N/A"])
    print(f"    → study2_summary.tsv")

    # ══════════════════════════════════════════════════════════════════════
    #  FINAL REPORT
    # ══════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("STUDY 2 FINAL RESULTS: GENCODE v37 vs v45")
    print("=" * 70)

    print(f"\n  METRIC                      Study 2 (v37-v45)  Primary (v29-v45)")
    print(f"  {'─'*64}")
    print(f"  Instability rate            {point_est:>8.4f}           "
          f"{PRIMARY['instability_rate']:>8.4f}")
    print(f"  95% CI                 [{ci_low:.4f},{ci_high:.4f}]  "
          f"[{PRIMARY['ci_low']:.4f},{PRIMARY['ci_high']:.4f}]")
    print(f"  Spearman rho                {rho:>8.4f}           "
          f"{PRIMARY['spearman_rho']:>8.4f}")
    print(f"  Jaccard (FDR<0.05)          {jaccard:>8.4f}           "
          f"{PRIMARY['jaccard_005']:>8.4f}")
    print(f"  Sig pathways (older ver)    {len(sig_v37):>8}           "
          f"{PRIMARY['sig_v29']:>8}")
    print(f"  Sig pathways (v45)          {len(sig_v45):>8}           "
          f"{PRIMARY['sig_v45']:>8}")
    print(f"  Shared significant          {len(inter):>8}           "
          f"{PRIMARY['shared']:>8}")

    print(f"\n  Gene reassignment decomposition:")
    print(f"    Total loci mapped:             {decomp['n_total']:>10,}")
    print(f"    Ensembl ID differs:            {decomp['n_ensembl_differs']:>10,} "
          f"({100*decomp['n_ensembl_differs']/decomp['n_total']:.2f}%)")
    print(f"    ID-only (symbol stable):       {decomp['n_id_only']:>10,} "
          f"({pct_id_only:.1f}% of diffs)")
    print(f"    Both ID+symbol change:         {decomp['n_both_differ']:>10,} "
          f"({pct_both:.1f}% of diffs)")
    print(f"    Symbol-only rename (no drift): {decomp['n_symbol_only']:>10,} "
          f"({pct_sym_only:.2f}% of total)")

    print(f"\n  Output files in {OUT_DIR}/:")
    for fname in sorted(os.listdir(OUT_DIR)):
        fpath = os.path.join(OUT_DIR, fname)
        size = os.path.getsize(fpath)
        print(f"    {fname:40s} {size:>10,} bytes")

    print("\n" + "=" * 70)
    print("STUDY 2 COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
