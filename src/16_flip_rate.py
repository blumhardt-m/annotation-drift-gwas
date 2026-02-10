#!/usr/bin/env python3
"""
Pathway category flip rate: do the biological categories implied by
enrichment change across annotation versions?

Maps significant pathways (Ensembl-ID ORA, FDR < 0.05) to Reactome
top-level categories and computes a narrative flip rate.

Usage:
    python 16_flip_rate.py \
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
import re
import sys
import urllib.request


# ── Import ORA machinery from src/12 ──────────────────────────────────

def _import_src12():
    src_dir = os.path.dirname(os.path.abspath(__file__))
    path_12 = os.path.join(src_dir, "12_ensembl_id_enrichment.py")
    spec = importlib.util.spec_from_file_location("src12", path_12)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


src12 = _import_src12()


# ── GMT loader with Reactome IDs ──────────────────────────────────────

def load_gmt_with_ids(gmt_path):
    """Return dict: pathway_name -> reactome_stable_id (R-HSA-xxxxx)."""
    name_to_id = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            desc = parts[1]
            if desc.startswith('R-HSA-'):
                name_to_id[name] = desc
    return name_to_id


# ── Reactome hierarchy download and parsing ───────────────────────────

def download_reactome_files(data_dir):
    """Download hierarchy and pathway name files. Returns (relations_path, names_path) or (None, None)."""
    os.makedirs(data_dir, exist_ok=True)
    rel_path = os.path.join(data_dir, "ReactomePathwaysRelation.txt")
    names_path = os.path.join(data_dir, "ReactomePathways.txt")

    urls = {
        rel_path: "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
        names_path: "https://reactome.org/download/current/ReactomePathways.txt",
    }
    for path, url in urls.items():
        if os.path.exists(path):
            continue
        print(f"  Downloading {url}...")
        try:
            urllib.request.urlretrieve(url, path)
        except Exception as e:
            print(f"  Download failed: {e}")
            return None, None
    return rel_path, names_path


def build_hierarchy(rel_path, names_path):
    """Build mapping: reactome_id -> top-level category name (human only)."""
    # Load human pathway names
    id_to_name = {}
    with open(names_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[2] == "Homo sapiens":
                id_to_name[parts[0]] = parts[1]

    # Load parent→child relations (filter to human: R-HSA-)
    children_of = {}  # parent -> set of children
    parent_of = {}    # child -> parent
    with open(rel_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            parent, child = parts[0], parts[1]
            if not parent.startswith('R-HSA-') or not child.startswith('R-HSA-'):
                continue
            children_of.setdefault(parent, set()).add(child)
            parent_of[child] = parent

    # Top-level categories: human pathways that have no parent
    top_level = set()
    for pid in id_to_name:
        if pid not in parent_of:
            top_level.add(pid)

    # For each pathway, walk up to find its top-level ancestor
    def find_top(pid, visited=None):
        if visited is None:
            visited = set()
        if pid in visited:
            return pid  # cycle guard
        visited.add(pid)
        if pid in top_level:
            return pid
        if pid in parent_of:
            return find_top(parent_of[pid], visited)
        return pid  # orphan — treat as its own category

    pathway_to_category = {}
    for pid in id_to_name:
        top = find_top(pid)
        cat_name = id_to_name.get(top, top)
        pathway_to_category[pid] = cat_name

    print(f"  Human pathways: {len(id_to_name)}")
    print(f"  Top-level categories: {len(top_level)}")
    for tid in sorted(top_level):
        print(f"    {tid}: {id_to_name.get(tid, '?')}")

    return pathway_to_category


# ── Keyword fallback ──────────────────────────────────────────────────

KEYWORD_CATEGORIES = [
    ("Immune System", ["immune", "interferon", "interleukin", "cytokine",
                       "complement", "toll-like", "antigen", "MHC", "TCR",
                       "BCR", "innate", "adaptive"]),
    ("Signal Transduction", ["signaling", "signalling", "MAPK", "Wnt",
                             "Notch", "Hedgehog", "GPCR", "receptor tyrosine",
                             "RAS", "PI3K", "AKT", "mTOR"]),
    ("Metabolism", ["metabol", "glycolysis", "TCA", "oxidative phosphorylation",
                    "fatty acid", "amino acid", "lipid", "cholesterol",
                    "nucleotide", "vitamin"]),
    ("Cell Cycle", ["cell cycle", "mitotic", "mitosis", "checkpoint",
                    "DNA replication", "S phase", "M phase", "G1", "G2"]),
    ("Gene Expression", ["transcription", "translation", "mRNA", "splicing",
                         "ribosom", "RNA polymerase", "chromatin"]),
    ("Neuronal System", ["neuronal", "synap", "neurotransmit", "GABA",
                         "glutamate", "serotonin", "dopamine"]),
    ("Developmental Biology", ["developmental", "axon guidance", "Robo",
                               "Slit", "ephrin", "morphogen"]),
    ("Hemostasis", ["hemostasis", "coagulation", "platelet", "fibrin",
                    "thrombin"]),
    ("Disease", ["disease", "infection", "HIV", "influenza", "SARS",
                 "cancer", "oncogen", "tumor"]),
]


def keyword_category(pathway_name):
    name_lower = pathway_name.lower()
    for category, keywords in KEYWORD_CATEGORIES:
        for kw in keywords:
            if kw.lower() in name_lower:
                return category
    return "Other"


def build_keyword_mapping(gmt_names):
    """Fallback: assign categories by keyword matching on pathway names."""
    print("  Using keyword fallback for category assignment.")
    mapping = {}
    for name in gmt_names:
        mapping[name] = keyword_category(name)
    return mapping


# ── Spearman rank correlation ─────────────────────────────────────────

def spearman_rho(x, y):
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


# ── Main ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Pathway category flip rate analysis")
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--gmt', required=True)
    parser.add_argument('--min-size', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=500)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── 1. Run Ensembl-ID ORA ──
    print("Loading data and running ORA...")
    symbol_to_ensembl, ensembl_universe = src12.load_tss_tables(
        args.tss_v29, args.tss_v45)
    pathways, _, _ = src12.load_gmt_as_ensembl(
        args.gmt, symbol_to_ensembl, args.min_size, args.max_size)
    v29_genes, v45_genes, _ = src12.load_comparison_ensembl(args.comparison)

    ora_v29 = src12.run_ora(v29_genes, ensembl_universe, pathways)
    ora_v45 = src12.run_ora(v45_genes, ensembl_universe, pathways)

    names_v29 = [r[0] for r in ora_v29]
    fdr_v29 = dict(zip(names_v29, src12.bh_correction([r[1] for r in ora_v29])))
    names_v45 = [r[0] for r in ora_v45]
    fdr_v45 = dict(zip(names_v45, src12.bh_correction([r[1] for r in ora_v45])))

    all_names = sorted(set(names_v29) | set(names_v45))
    sig_v29 = {n for n in all_names if fdr_v29.get(n, 1.0) < 0.05}
    sig_v45 = {n for n in all_names if fdr_v45.get(n, 1.0) < 0.05}

    print(f"  Sig v29: {len(sig_v29)}, Sig v45: {len(sig_v45)}")

    # ── 2. Build pathway→category mapping ──
    print("\nBuilding pathway category mapping...")
    name_to_reactome_id = load_gmt_with_ids(args.gmt)

    # Try downloading Reactome hierarchy
    data_dir = os.path.join(os.path.dirname(args.gmt))
    rel_path, names_path = download_reactome_files(data_dir)

    use_hierarchy = False
    if rel_path and names_path:
        try:
            reactome_id_to_category = build_hierarchy(rel_path, names_path)
            use_hierarchy = True
        except Exception as e:
            print(f"  Hierarchy parsing failed: {e}")

    # Map pathway names → categories
    pathway_category = {}
    if use_hierarchy:
        n_mapped = 0
        for pw_name in all_names:
            rid = name_to_reactome_id.get(pw_name)
            if rid and rid in reactome_id_to_category:
                pathway_category[pw_name] = reactome_id_to_category[rid]
                n_mapped += 1
            else:
                pathway_category[pw_name] = keyword_category(pw_name)
        print(f"  Hierarchy-mapped: {n_mapped}/{len(all_names)}, "
              f"keyword fallback: {len(all_names) - n_mapped}")
    else:
        pathway_category = build_keyword_mapping(all_names)

    # ── 3. Category-level summary ──
    print("\nComputing category summary...")
    cats_v29 = {}
    for pw in sig_v29:
        cat = pathway_category.get(pw, "Other")
        cats_v29[cat] = cats_v29.get(cat, 0) + 1

    cats_v45 = {}
    for pw in sig_v45:
        cat = pathway_category.get(pw, "Other")
        cats_v45[cat] = cats_v45.get(cat, 0) + 1

    all_cats = sorted(set(cats_v29) | set(cats_v45))

    # Category summary table
    cat_rows = []
    for cat in all_cats:
        c29 = cats_v29.get(cat, 0)
        c45 = cats_v45.get(cat, 0)
        status = "both"
        if c29 > 0 and c45 == 0:
            status = "v29_only"
        elif c29 == 0 and c45 > 0:
            status = "v45_only"
        cat_rows.append((cat, c29, c45, status))

    cat_path = os.path.join(args.outdir, "category_summary.tsv")
    with open(cat_path, 'w') as f:
        f.write("category\tn_v29\tn_v45\tstatus\n")
        for cat, c29, c45, status in cat_rows:
            f.write(f"{cat}\t{c29}\t{c45}\t{status}\n")

    print(f"\n  Category summary:")
    for cat, c29, c45, status in cat_rows:
        flag = " ***" if status != "both" else ""
        print(f"    {cat:<35s} v29={c29:>3d}  v45={c45:>3d}  {status}{flag}")

    # ── 4. Flip rate ──
    n_union = len(all_cats)
    v29_only_cats = [cat for cat, _, _, s in cat_rows if s == "v29_only"]
    v45_only_cats = [cat for cat, _, _, s in cat_rows if s == "v45_only"]
    n_flips = len(v29_only_cats) + len(v45_only_cats)
    flip_rate = n_flips / n_union if n_union > 0 else float('nan')

    # ── 5. Rank correlation of shared categories ──
    shared_cats = [cat for cat in all_cats
                   if cats_v29.get(cat, 0) > 0 and cats_v45.get(cat, 0) > 0]
    if len(shared_cats) >= 2:
        counts_29 = [cats_v29[c] for c in shared_cats]
        counts_45 = [cats_v45[c] for c in shared_cats]
        rho = spearman_rho(counts_29, counts_45)
    else:
        rho = float('nan')

    # Save flip rate results
    flip_path = os.path.join(args.outdir, "flip_rate.tsv")
    with open(flip_path, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"n_categories_union\t{n_union}\n")
        f.write(f"n_categories_v29\t{len(cats_v29)}\n")
        f.write(f"n_categories_v45\t{len(cats_v45)}\n")
        f.write(f"n_flips\t{n_flips}\n")
        f.write(f"flip_rate\t{flip_rate:.6f}\n")
        f.write(f"n_shared_categories\t{len(shared_cats)}\n")
        f.write(f"spearman_rho_shared\t{rho:.6f}\n")
        f.write(f"category_source\t{'hierarchy' if use_hierarchy else 'keyword'}\n")

    print(f"\n  Flip rate summary:")
    print(f"    Categories (union):     {n_union}")
    print(f"    v29-only categories:    {len(v29_only_cats)} {v29_only_cats}")
    print(f"    v45-only categories:    {len(v45_only_cats)} {v45_only_cats}")
    print(f"    Flip rate:              {flip_rate:.4f} ({n_flips}/{n_union})")
    print(f"    Shared categories:      {len(shared_cats)}")
    print(f"    Spearman rho (shared):  {rho:.4f}")
    print(f"    Category source:        {'hierarchy' if use_hierarchy else 'keyword'}")
    print(f"\n  Saved {cat_path}")
    print(f"  Saved {flip_path}")


if __name__ == "__main__":
    main()
