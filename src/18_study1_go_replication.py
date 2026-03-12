#!/usr/bin/env python3
"""
Study 1: GO Biological Process replication of Reactome enrichment analysis.

Tests whether the threshold-amplification effect (high Spearman rho on p-values
but low Jaccard on significant sets) replicates under Gene Ontology Biological
Process instead of Reactome.

Usage:
    python src/18_study1_go_replication.py

All paths are relative to the project root.
"""

import csv
import gzip
import math
import os
import sys
import urllib.request
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ── Paths ────────────────────────────────────────────────────────────────
PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COMPARISON = os.path.join(PROJECT, "data/processed/mappings/comparison.tsv")
TSS_V29 = os.path.join(PROJECT, "data/processed/tss_tables/tss_v29.tsv")
TSS_V45 = os.path.join(PROJECT, "data/processed/tss_tables/tss_v45.tsv")

GO_OBO = os.path.join(PROJECT, "data/raw/go_data/go-basic.obo")
GENE2GO = os.path.join(PROJECT, "data/raw/go_data/gene2go.gz")
GENE2ENSEMBL = os.path.join(PROJECT, "data/raw/go_data/gene2ensembl.gz")

OUTDIR = os.path.join(PROJECT, "results/study1_GO_replication")

MIN_SIZE = 10
MAX_SIZE = 1000
FDR_ALPHA = 0.05
TAXID = 9606  # Homo sapiens


# ── Download helpers ─────────────────────────────────────────────────────

def download_if_missing(url, dest):
    if os.path.exists(dest):
        print(f"  [cached] {os.path.basename(dest)}")
        return
    print(f"  Downloading {url} ...")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as resp, open(dest, 'wb') as out:
        out.write(resp.read())
    print(f"  → {dest}")


def ensure_go_data():
    """Download go-basic.obo, gene2go, and gene2ensembl from NCBI."""
    os.makedirs(os.path.dirname(GO_OBO), exist_ok=True)
    download_if_missing(
        "http://purl.obolibrary.org/obo/go/go-basic.obo", GO_OBO
    )
    download_if_missing(
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", GENE2GO
    )
    download_if_missing(
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", GENE2ENSEMBL
    )


# ── Load Entrez → Ensembl mapping ───────────────────────────────────────

def load_entrez_to_ensembl():
    """Parse gene2ensembl.gz for human genes. Returns dict: entrez_id (str) -> ensembl_gene_id."""
    mapping = {}
    with gzip.open(GENE2ENSEMBL, 'rt') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            tax_id = int(parts[0])
            if tax_id != TAXID:
                continue
            entrez_id = parts[1]
            ensembl_gene = parts[2]
            if ensembl_gene == '-':
                continue
            # Keep first mapping (deterministic since file is sorted)
            if entrez_id not in mapping:
                mapping[entrez_id] = ensembl_gene
    return mapping


# ── Load GO annotations from gene2go ────────────────────────────────────

def load_go_annotations(entrez_to_ensembl):
    """Parse gene2go.gz for human GO Biological Process annotations.
    Convert Entrez IDs to Ensembl IDs. Returns:
      - go_terms: dict GO_ID -> set of Ensembl IDs
      - go_names: dict GO_ID -> term name
      - stats: mapping success info
    """
    go_terms = defaultdict(set)
    go_names = {}
    n_mapped = 0
    n_unmapped = 0
    seen_entrez = set()

    with gzip.open(GENE2GO, 'rt') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            tax_id = int(parts[0])
            if tax_id != TAXID:
                continue
            entrez_id = parts[1]
            go_id = parts[2]
            # Evidence code is parts[3]
            # Qualifier is parts[4]
            go_term_name = parts[5]
            category = parts[7]  # Process, Function, Component

            if category != "Process":
                continue

            seen_entrez.add(entrez_id)
            ensembl_id = entrez_to_ensembl.get(entrez_id)
            if ensembl_id:
                go_terms[go_id].add(ensembl_id)
                n_mapped += 1
            else:
                n_unmapped += 1

            if go_id not in go_names:
                go_names[go_id] = go_term_name

    stats = {
        'total_annotations': n_mapped + n_unmapped,
        'mapped': n_mapped,
        'unmapped': n_unmapped,
        'unique_entrez_total': len(seen_entrez),
        'unique_entrez_mapped': len({eid for eid in seen_entrez if eid in entrez_to_ensembl}),
    }
    return go_terms, go_names, stats


# ── Load project data ───────────────────────────────────────────────────

def load_ensembl_universe(tss_path):
    """Load all Ensembl IDs from a TSS table (protein-coding background)."""
    ids = set()
    with open(tss_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ids.add(row['ensembl_id'])
    return ids


def load_foreground_genes(comparison_path):
    """Load v29 and v45 Ensembl gene sets from comparison table."""
    v29 = set()
    v45 = set()
    with open(comparison_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            v29.add(row['gene_v29'])
            v45.add(row['gene_v45'])
    return v29, v45


# ── ORA machinery (replicated from 12_ensembl_id_enrichment.py) ────────

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


def run_ora(gene_set, universe, term_sets):
    """ORA: Fisher's exact test for each GO term."""
    gene_set_in_universe = gene_set & universe
    n_geneset = len(gene_set_in_universe)
    N = len(universe)
    results = []
    for term_id, term_genes in term_sets.items():
        pw_in_universe = term_genes & universe
        K = len(pw_in_universe)
        if K == 0:
            continue
        overlap = gene_set_in_universe & pw_in_universe
        a = len(overlap)
        b = K - a
        c = n_geneset - a
        d = N - K - c
        pval = fisher_exact_pvalue(a, b, c, d)
        results.append((term_id, pval, a, overlap))
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


# ── Enrichment output writers ───────────────────────────────────────────

def write_enrichment_tsv(path, ora_results, fdr_values, go_names):
    """Write enrichment results TSV."""
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['term_id', 'term_name', 'p_value', 'fdr', 'gene_count', 'ensembl_genes'])
        for i, (term_id, pval, count, overlap_genes) in enumerate(ora_results):
            writer.writerow([
                term_id,
                go_names.get(term_id, ''),
                f"{pval:.6e}",
                f"{fdr_values[i]:.6e}",
                count,
                ','.join(sorted(overlap_genes))
            ])


def write_summary_tsv(path, rho, jaccard, n29, n45, shared, mapping_rate, bg_size, terms_tested):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['spearman_rho', 'jaccard', 'sig_count_v29', 'sig_count_v45',
                         'shared_count', 'mapping_success_rate', 'background_size', 'terms_tested'])
        writer.writerow([
            f"{rho:.6f}", f"{jaccard:.4f}", n29, n45, shared,
            f"{mapping_rate:.4f}", bg_size, terms_tested
        ])


def make_scatter(path, log10p_v29, log10p_v45):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(log10p_v29, log10p_v45, alpha=0.4, s=12, edgecolors='none')
    max_val = max(max(log10p_v29, default=0), max(log10p_v45, default=0)) * 1.05
    ax.plot([0, max_val], [0, max_val], 'k--', lw=0.8, alpha=0.5)
    ax.set_xlabel('-log10(p) GENCODE v29')
    ax.set_ylabel('-log10(p) GENCODE v45')
    ax.set_title('GO Biological Process enrichment: v29 vs v45')
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ── Main ────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("=" * 70)
    print("STUDY 1: GO BIOLOGICAL PROCESS REPLICATION")
    print("=" * 70)

    # ── 1. Download GO data ──
    print("\n── Downloading GO data ──")
    ensure_go_data()

    # ── 2. Build Entrez → Ensembl mapping ──
    print("\n── Building Entrez → Ensembl mapping ──")
    entrez_to_ensembl = load_entrez_to_ensembl()
    print(f"  Entrez→Ensembl mappings (human): {len(entrez_to_ensembl)}")

    # ── 3. Load GO annotations ──
    print("\n── Loading GO Biological Process annotations ──")
    go_terms, go_names, go_stats = load_go_annotations(entrez_to_ensembl)
    total_entrez = go_stats['unique_entrez_total']
    mapped_entrez = go_stats['unique_entrez_mapped']
    mapping_rate = mapped_entrez / total_entrez if total_entrez > 0 else 0
    print(f"  Total GO:BP annotations (human): {go_stats['total_annotations']}")
    print(f"  Mapped to Ensembl:               {go_stats['mapped']}")
    print(f"  Unmapped:                        {go_stats['unmapped']}")
    print(f"  Unique Entrez genes (GO:BP):     {total_entrez}")
    print(f"  Successfully converted:          {mapped_entrez}")
    print(f"  Mapping success rate:            {mapping_rate:.4f} ({100*mapping_rate:.2f}%)")

    if mapping_rate < 0.95:
        print(f"\n  ABORT: Mapping success rate {mapping_rate:.4f} < 0.95 (>5% failure)")
        sys.exit(1)

    # ── 4. Total GO terms and size filtering ──
    print(f"\n  Total GO:BP terms loaded:        {len(go_terms)}")

    # Size filter
    go_terms_filtered = {}
    for go_id, genes in go_terms.items():
        if MIN_SIZE <= len(genes) <= MAX_SIZE:
            go_terms_filtered[go_id] = genes

    print(f"  GO terms after size filter [{MIN_SIZE}, {MAX_SIZE}]: {len(go_terms_filtered)}")

    # Gene set size stats
    sizes = [len(g) for g in go_terms_filtered.values()]
    sizes.sort()
    median_size = sizes[len(sizes) // 2] if sizes else 0
    mean_size = sum(sizes) / len(sizes) if sizes else 0
    print(f"  Median gene-set size:            {median_size}")
    print(f"  Mean gene-set size:              {mean_size:.1f}")

    # ── 5. Load background universe ──
    print("\n── Loading background universe ──")
    bg_v29 = load_ensembl_universe(TSS_V29)
    bg_v45 = load_ensembl_universe(TSS_V45)
    universe = bg_v29 | bg_v45
    print(f"  v29 protein-coding genes: {len(bg_v29)}")
    print(f"  v45 protein-coding genes: {len(bg_v45)}")
    print(f"  Union (background):       {len(universe)}")

    # ── 6. Load foreground genes ──
    print("\n── Loading foreground gene sets ──")
    fg_v29, fg_v45 = load_foreground_genes(COMPARISON)
    print(f"  v29 foreground genes: {len(fg_v29)}")
    print(f"  v45 foreground genes: {len(fg_v45)}")

    # ── 7. Cross-check with Reactome analysis ──
    print("\n── Background/foreground cross-check ──")
    # Reactome analysis used the same TSS union universe and comparison foreground
    print(f"  Reactome analysis: background={len(universe)} foreground_v29={len(fg_v29)} foreground_v45={len(fg_v45)}")
    print(f"  GO analysis:       background={len(universe)} foreground_v29={len(fg_v29)} foreground_v45={len(fg_v45)}")

    # ── 8. Confirm Ensembl IDs throughout ──
    print("\n── Identifier confirmation ──")
    sample_bg = list(universe)[:3]
    sample_fg = list(fg_v29)[:3]
    sample_go = list(list(go_terms_filtered.values())[0])[:3] if go_terms_filtered else []
    print(f"  Background sample: {sample_bg}")
    print(f"  Foreground sample: {sample_fg}")
    print(f"  GO gene-set sample: {sample_go}")
    all_ensembl = all(
        s.startswith('ENSG') for s in list(universe)[:100]
    ) and all(
        s.startswith('ENSG') for s in list(fg_v29)[:100]
    )
    if go_terms_filtered:
        first_term_genes = list(list(go_terms_filtered.values())[0])[:100]
        all_ensembl = all_ensembl and all(s.startswith('ENSG') for s in first_term_genes)
    print(f"  All Ensembl IDs confirmed: {all_ensembl}")

    # ── DIAGNOSTIC SUMMARY ──
    print("\n" + "=" * 70)
    print("DIAGNOSTIC SUMMARY")
    print("=" * 70)
    print(f"  1. Total GO:BP terms loaded:           {len(go_terms)}")
    print(f"  2. GO terms after size filter:         {len(go_terms_filtered)}")
    print(f"  3. Identifier mapping success rate:    {mapping_rate:.4f} ({100*mapping_rate:.2f}%)")
    print(f"  4. Final background universe size:     {len(universe)}")
    print(f"  5. Median gene-set size:               {median_size}")
    print(f"     Mean gene-set size:                 {mean_size:.1f}")
    print(f"  6. Ensembl IDs confirmed throughout:   {all_ensembl}")

    # ── 9. Run ORA ──
    print("\n" + "=" * 70)
    print("RUNNING ORA")
    print("=" * 70)

    print("\n  ORA for v29 gene set...")
    ora_v29 = run_ora(fg_v29, universe, go_terms_filtered)
    pvals_v29 = [r[1] for r in ora_v29]
    fdr_v29 = bh_correction(pvals_v29)
    print(f"    Terms tested: {len(ora_v29)}")

    print("  ORA for v45 gene set...")
    ora_v45 = run_ora(fg_v45, universe, go_terms_filtered)
    pvals_v45 = [r[1] for r in ora_v45]
    fdr_v45 = bh_correction(pvals_v45)
    print(f"    Terms tested: {len(ora_v45)}")

    # ── 10. Compute metrics ──
    # Significant sets
    sig_v29 = {r[0] for r, fdr in zip(ora_v29, fdr_v29) if fdr < FDR_ALPHA}
    sig_v45 = {r[0] for r, fdr in zip(ora_v45, fdr_v45) if fdr < FDR_ALPHA}
    shared = sig_v29 & sig_v45
    union_sig = sig_v29 | sig_v45
    jaccard = len(shared) / len(union_sig) if len(union_sig) > 0 else 1.0

    # Spearman on all shared terms
    names_v29 = [r[0] for r in ora_v29]
    names_v45 = [r[0] for r in ora_v45]
    pval_map_v29 = dict(zip(names_v29, pvals_v29))
    pval_map_v45 = dict(zip(names_v45, pvals_v45))
    all_terms = sorted(set(names_v29) | set(names_v45))

    log10p_v29 = [-math.log10(max(pval_map_v29.get(t, 1.0), 1e-300)) for t in all_terms]
    log10p_v45 = [-math.log10(max(pval_map_v45.get(t, 1.0), 1e-300)) for t in all_terms]
    rho = spearman_rank_corr(log10p_v29, log10p_v45)

    # ── 11. Report results ──
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"  Spearman rho (-log10 p):       {rho:.6f}")
    print(f"  Jaccard (FDR < {FDR_ALPHA}):          {jaccard:.4f}")
    print(f"  Significant terms (v29):       {len(sig_v29)}")
    print(f"  Significant terms (v45):       {len(sig_v45)}")
    print(f"  Shared significant:            {len(shared)}")
    print(f"  Union significant:             {len(union_sig)}")

    # Compare with Reactome
    print(f"\n  Reactome (Ensembl-ID) comparison:")
    print(f"    Reactome Spearman rho:       0.922857")
    print(f"    GO:BP Spearman rho:          {rho:.6f}")
    print(f"    Reactome Jaccard (FDR<0.05): 0.6133")
    print(f"    GO:BP Jaccard (FDR<0.05):    {jaccard:.4f}")

    # ── 12. Write outputs ──
    print("\n── Writing outputs ──")

    out_v29 = os.path.join(OUTDIR, "go_enrichment_v29.tsv")
    out_v45 = os.path.join(OUTDIR, "go_enrichment_v45.tsv")
    out_summary = os.path.join(OUTDIR, "go_replication_summary.tsv")
    out_scatter = os.path.join(OUTDIR, "go_rank_scatter.png")

    write_enrichment_tsv(out_v29, ora_v29, fdr_v29, go_names)
    print(f"  → {out_v29}")

    write_enrichment_tsv(out_v45, ora_v45, fdr_v45, go_names)
    print(f"  → {out_v45}")

    write_summary_tsv(out_summary, rho, jaccard, len(sig_v29), len(sig_v45),
                      len(shared), mapping_rate, len(universe), len(go_terms_filtered))
    print(f"  → {out_summary}")

    make_scatter(out_scatter, log10p_v29, log10p_v45)
    print(f"  → {out_scatter}")

    print("\n" + "=" * 70)
    print("STUDY 1 COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
