"""Frozen computations plus the prospective per-trait protocol orchestration."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import re
import statistics
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

def imported(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

parser03 = imported("pilot_frozen03", ROOT / "src/03_parse_gwas.py")
mapping04 = imported("pilot_frozen04", ROOT / "src/04_map_variants.py")
ora12 = imported("pilot_frozen12", ROOT / "src/12_ensembl_id_enrichment.py")
diagnostics13 = imported("pilot_frozen13", ROOT / "src/13_jaccard_robustness.py")
# Pure-result memoization: every new contingency table is evaluated by the
# original imported function. No replacement numerical algorithm or rounding.
uncached_fisher = ora12.fisher_exact_pvalue
ora12.fisher_exact_pvalue = lru_cache(maxsize=None)(uncached_fisher)

@dataclass(frozen=True)
class Sentinel:
    trait: str
    rsid: str
    chrom: int
    pos: int
    p: Decimal

    def key(self):
        return self.p, self.chrom, self.pos, self.rsid

def thin_rows(rows):
    """Strict filtering, minimum-p deduplication and frozen greedy thinning."""
    counts = Counter()
    minima = {}
    excluded = []
    for row in rows:
        counts["input_rows"] += 1
        reason = None
        rsid = row["SNPS"]
        if not re.fullmatch(r"rs[0-9]+", rsid):
            reason = "NON_SINGLE_RSID"
        else:
            p = Decimal(row["P-VALUE"])
            if not p.is_finite() or p <= 0:
                raise RuntimeError("Nonpositive/nonfinite p needs prospective disposition")
            if p >= Decimal("5e-8"):
                reason = "P_NOT_STRICTLY_BELOW_5E8"
            elif row["CHR_ID"] not in {str(i) for i in range(1, 23)}:
                reason = "NON_AUTOSOMAL"
        if reason:
            counts[reason] += 1
            excluded.append({"trait": row["DISEASE/TRAIT"], "rsid": rsid, "reason": reason})
            continue
        s = Sentinel(row["DISEASE/TRAIT"], rsid, int(row["CHR_ID"]), int(row["CHR_POS"]), p)
        key = s.trait, s.rsid
        old = minima.get(key)
        if old is not None:
            counts["deduplicated_rows"] += 1
            if s.p == old.p and (s.chrom, s.pos) != (old.chrom, old.pos):
                raise RuntimeError("Equal minimum p has conflicting coordinates")
        if old is None or s.p < old.p:
            minima[key] = s
    grouped = defaultdict(list)
    for s in minima.values():
        grouped[s.trait].append(s)
    sentinels = {}
    for trait in sorted(grouped):
        kept = []
        positions = defaultdict(list)
        for s in sorted(grouped[trait], key=Sentinel.key):
            if any(abs(s.pos - p) <= 250000 for p in positions[s.chrom]):
                counts["THINNED_WITHIN_250KB"] += 1
                excluded.append({"trait": trait, "rsid": s.rsid, "reason": "THINNED_WITHIN_250KB"})
            else:
                kept.append(s)
                positions[s.chrom].append(s.pos)
        sentinels[trait] = kept
    return sentinels, dict(counts), excluded

def roster(sentinels, cutoff=100, primary=None):
    eligible = sorted(t for t, rows in sentinels.items() if len(rows) >= cutoff)
    if primary is None:
        return eligible
    assert cutoff == 50 and all(t in eligible for t in primary)
    return list(primary) + [t for t in eligible if t not in set(primary)]

def mapped(rows, tables):
    output = []
    for s in rows:
        result = mapping04.find_nearest_gene(str(s.chrom), s.pos, tables, 500000)
        output.append(result[0] if result is not None else None)
    return output

def foreground(assignments):
    return {g for g in assignments if g is not None}

def gene_statistics(rows, assignments, universe):
    scores = {g: 0.0 for g in sorted(universe)}
    for s, g in zip(rows, assignments):
        if g is not None:
            scores[g] = max(scores[g], float(-s.p.log10()))
    return dict(sorted(scores.items(), key=lambda pair: (-pair[1], pair[0])))

def enrichment(assignments, universe, pathways):
    return ora12.run_ora(foreground(assignments), universe, pathways)

def significance(ora, alpha=.05):
    qs = ora12.bh_correction([r[1] for r in ora])
    return {r[0] for r, q in zip(ora, qs) if q < alpha}

def compare(a, b, alpha=.05, q_a=None, q_b=None):
    """ORA tuples or fgsea tuples, with fgsea's own BH-adjusted values supplied."""
    if q_a is None:
        q_a = dict(zip([r[0] for r in a], ora12.bh_correction([r[1] for r in a])))
    if q_b is None:
        q_b = dict(zip([r[0] for r in b], ora12.bh_correction([r[1] for r in b])))
    pa = {r[0]: r[1] for r in a}
    pb = {r[0]: r[1] for r in b}
    assert set(pa) == set(pb) == set(q_a) == set(q_b)
    names = sorted(pa)
    sa = {n for n in names if q_a[n] < alpha}
    sb = {n for n in names if q_b[n] < alpha}
    union = sa | sb
    # The imported FDR sweep supplies the parent ORA Jaccard/rank machinery.
    _, rho = ora12.fdr_sweep(a, b, [alpha])
    row = {"informative": bool(union), "J": len(sa & sb) / len(union) if union else None,
           "F": len(sa ^ sb), "a_only": len(sa - sb), "b_only": len(sb - sa),
           "sig_v29": len(sa), "sig_v45": len(sb), "rho": rho}
    for k in (25, 50):
        row[f"top{k}_overlap"] = diagnostics13.top_k_overlap(pa, pb, names, k) if len(names) >= k else None
    oa, ob = {r[0]: r[2] for r in a}, {r[0]: r[2] for r in b}
    flips = [{"pathway": n, "p_v29": pa[n], "p_v45": pb[n], "q_v29": q_a[n], "q_v45": q_b[n],
              "distance_v29": abs(q_a[n] - alpha), "distance_v45": abs(q_b[n] - alpha),
              "overlap_v29": oa[n], "overlap_v45": ob[n], "delta_genes": ob[n] - oa[n],
              "direction": "v29_only" if n in sa else "v45_only"} for n in sorted(sa ^ sb)]
    return row, flips

def c2_plan(rows, a, b, reference_tss):
    changed = [j for j, (x, y) in enumerate(zip(a, b)) if x is not None and y is not None and x != y]
    hybrid = list(a)
    for j in changed:
        hybrid[j] = b[j]
    eligible = []
    for j, (s, gene) in enumerate(zip(rows, a)):
        if gene is None:
            continue
        alternatives = sorted({g for pos, g, _ in reference_tss.get(str(s.chrom), [])
                               if g != gene and abs(s.pos - pos) <= 500000})
        if alternatives:
            eligible.append((j, alternatives))
    eligible.sort(key=lambda pair: (rows[pair[0]].chrom, rows[pair[0]].pos, rows[pair[0]].rsid))
    return {"n1": len(changed), "n2_plus": sum(x is None and y is not None for x, y in zip(a, b)),
            "n2_minus": sum(x is not None and y is None for x, y in zip(a, b)),
            "hybrid": hybrid, "eligible": eligible}

def c2_draw(a, eligible, k1, rng):
    synthetic = list(a)
    selected = rng.choice(len(eligible), size=k1, replace=False)
    for j in sorted(selected):
        position, alternatives = eligible[j]
        synthetic[position] = alternatives[int(rng.integers(len(alternatives)))]
    return synthetic

def c2_run(rows, a, b, reference_tss, universe, pathways, trait_index):
    plan = c2_plan(rows, a, b, reference_tss)
    n1 = plan["n1"]
    k1 = min(n1, len(plan["eligible"]))
    result = {k: plan[k] for k in ("n1", "n2_plus", "n2_minus")}
    result.update(k1=k1, eligible=len(plan["eligible"]), shortfall=n1-k1, seed=42+trait_index)
    if n1 == 0:
        result.update(status="NO-DRIFT", observed_F_g2g=None, percentile95=None, E=None, null_sorted=[])
        return result
    reference = significance(enrichment(a, universe, pathways))
    observed = len(reference ^ significance(enrichment(plan["hybrid"], universe, pathways)))
    rng = np.random.default_rng(42 + trait_index)
    null = []
    for _ in range(1000):
        draw = c2_draw(a, plan["eligible"], k1, rng)
        null.append(len(reference ^ significance(enrichment(draw, universe, pathways))))
    null.sort()
    result.update(status="TRUNCATED" if k1 < n1 else "AVAILABLE",
                  observed_F_g2g=observed, percentile95=null[949], E=int(observed > null[949]), null_sorted=null)
    return result

def decision(summary, extension=False):
    c, j, f, tau = (summary[k] for k in ("c", "m_J", "m_F", "tau"))
    available, x = summary["C2_available"], summary["x"]
    if c < .60:
        return {"outcome": "UNDERPOWERED", "rule": 1, "extension_required": False}
    if j < .80 or tau > .25 or (available and x >= .25):
        return {"outcome": "PASS", "rule": 2, "extension_required": False}
    if available and j >= .90 and f <= 1 and tau <= .10 and x <= .10:
        return {"outcome": "KILL", "rule": 3, "extension_required": False}
    if not available:
        return {"outcome": "C2-UNAVAILABLE-FINAL" if extension else "C2-UNAVAILABLE",
                "rule": 4, "extension_required": not extension}
    return {"outcome": "INDETERMINATE-FINAL" if extension else "INDETERMINATE",
            "rule": 5, "extension_required": not extension}

def summarize(metrics, c2, held=None):
    n = len(metrics)
    if n == 0:
        raise RuntimeError("Empty eligible roster: no specified coverage denominator")
    informative = [i for i, row in enumerate(metrics) if row["informative"]]
    ni = len(informative)
    dx_indices = [i for i in informative if c2[i]["status"] == "AVAILABLE" and c2[i]["n1"] > 0]
    dx = len(dx_indices)
    available = dx >= 10 and dx >= .5 * ni
    count = sum(c2[i]["E"] for i in dx_indices)
    row = {"N": n, "n_i": ni, "u": n-ni, "c": ni/n,
           "m_J": statistics.median(metrics[i]["J"] for i in informative) if ni else None,
           "m_F": statistics.median(metrics[i]["F"] for i in informative) if ni else None,
           "tau": sum(metrics[i]["J"] < .5 or metrics[i]["F"] >= 3 for i in informative)/ni if ni else None,
           "d_x": dx, "C2_available": available, "C2_exceedances": count, "x": count/dx if available else None}
    if held is not None:
        for key in ("d_x", "C2_available", "C2_exceedances", "x"):
            row[key] = held[key]
        row["conditional_on_full_data_C2"] = True
    return row

def delete_chromosome(rows, chrom):
    return [s for s in rows if s.chrom != chrom]
