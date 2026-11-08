"""Synthetic protocol checks; no empirical outcomes are inspected."""
import csv
from decimal import Decimal
import hashlib
import math
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/pertrait_pilot"))
import core as c
import pre_auth_pilot as gate

def row(rs, pos, p="1e-9", trait="Trait", chrom="1"):
    return {"SNPS": rs, "CHR_POS": str(pos), "P-VALUE": p, "DISEASE/TRAIT": trait, "CHR_ID": chrom}

def sentinel(rs, pos=1, chrom=1):
    return c.Sentinel("T", rs, chrom, pos, Decimal("1e-10"))

def test_strict_filter_min_p_thinning_and_all_ties():
    rows = [row("rs9", 1, "2e-9"), row("rs9", 1, "1e-10"),
            row("rs2", 250001), row("rs3", 250002), row("rs1", 1, chrom="2"),
            row("rs2", 1, chrom="2"), row("rs10", 1, chrom="2"),
            row("rs77; rs78", 800000), row("RS1", 900000), row("rs0", 1000000, "5e-8"),
            row("rs8", 1, chrom="X")]
    # Repeated rs2 has a tied minimum at conflicting positions and must refuse.
    with pytest.raises(RuntimeError, match="conflicting"):
        c.thin_rows(rows)
    rows[4] = row("rs11", 1, chrom="2")
    rows[5] = row("rs12", 1, chrom="2")
    result, counts, exclusions = c.thin_rows(rows)
    assert [(s.rsid, s.chrom, s.pos) for s in result["Trait"]] == [("rs9", 1, 1), ("rs3", 1, 250002), ("rs10", 2, 1)]
    assert result["Trait"][0].p == Decimal("1e-10")
    assert counts["NON_SINGLE_RSID"] == 2 and counts["deduplicated_rows"] == 1
    assert counts["P_NOT_STRICTLY_BELOW_5E8"] == 1
    assert c.thin_rows(list(reversed(rows)))[0] == result

def test_roster_codepoint_order_and_extension_indices():
    groups = {"a": [0]*100, "Z": [0]*101, "A": [0]*50, "ä": [0]*99, "x": [0]*49}
    primary = c.roster(groups)
    assert primary == ["Z", "a"]
    assert c.roster(groups, 50, primary) == ["Z", "a", "A", "ä"]

def test_mapping_enrichment_analytical_case_and_memoization():
    universe = {f"G{i:03d}" for i in range(100)}
    genes = sorted(universe)
    pathways = {"synthetic": set(genes[:15])}
    rows = [sentinel(f"rs{i}", i*1000000) for i in range(4)]
    tss = {"1": [(i*1000000, genes[i], genes[i]) for i in range(4)]}
    assignments = c.mapped(rows, tss)
    assert assignments == genes[:4]
    result = c.enrichment(assignments, universe, pathways)
    expected = math.comb(15, 4)/math.comb(100, 4)
    assert result[0][1] == pytest.approx(expected, rel=1e-11)
    assert result[0][2] == 4 and c.significance(result) == {"synthetic"}
    baseline = c.uncached_fisher(4, 11, 0, 85)
    assert c.ora12.fisher_exact_pvalue(4, 11, 0, 85) == baseline
    opposite = c.enrichment(genes[30:34], universe, pathways)
    metrics, flips = c.compare(result, opposite)
    assert metrics["J"] == 0 and metrics["F"] == 1
    assert flips[0]["delta_genes"] == -4
    assert c.mapped([sentinel("rs9", 500001)], {"1": [(0, "A", "A")]}) == [None]
    assert c.mapped([sentinel("rs9", 500000)], {"1": [(0, "A", "A")]}) == ["A"]
    assert c.mapped([sentinel("rs9", 10)], {"1": [(0, "Z", "Z"), (20, "A", "A")]}) == ["A"]

def test_gene_statistics_max_zero_tail_and_id_ties():
    rows = [sentinel("rs1"), sentinel("rs2"), sentinel("rs3")]
    rows[1] = c.Sentinel("T","rs2",1,3,Decimal("1e-400"))
    scores = c.gene_statistics(rows, ["G2","G2","G1"], {"G1","G2","G3","G0"})
    assert list(scores) == ["G2","G1","G0","G3"]
    assert scores == {"G2":400., "G1":10., "G0":0., "G3":0.}

def test_g2g_hybrid_excludes_gains_and_losses():
    rows = [sentinel("rs1"),sentinel("rs2",2),sentinel("rs3",3)]
    plan = c.c2_plan(rows,["A",None,"C"],["B","D",None],{"1":[(0,"A",""),(1,"B",""),(2,"C","")]})
    assert (plan["n1"],plan["n2_plus"],plan["n2_minus"]) == (1,1,1)
    assert plan["hybrid"] == ["B",None,"C"]
    assert plan["eligible"] == [(0,["B","C"]),(2,["A","B"])]

def test_rng_exact_call_sequence_and_sorted_selection():
    class Trace:
        def __init__(self): self.calls=[]
        def choice(self,n,size,replace):
            self.calls.append(("choice",n,size,replace));return np.array([2,0])
        def integers(self,n):
            self.calls.append(("integers",n));return n-1
    rng=Trace()
    output=c.c2_draw(["a","b","c"],[(0,["x","y"]),(1,["z"]),(2,["w"])],2,rng)
    assert rng.calls == [("choice",3,2,False),("integers",2),("integers",1)]
    assert output == ["y","b","w"]

def test_zero_drift_no_alternative_and_truncation():
    rows=[sentinel("rs1")]
    u={"A","B"}
    no=c.c2_run(rows,["A"],["A"],{"1":[(1,"A","")]},u,{"P":{"A"}},0)
    assert no["status"]=="NO-DRIFT" and no["null_sorted"]==[]
    truncated=c.c2_run(rows,["A"],["B"],{"1":[(1,"A","")]},u,{"P":{"A"}},0)
    assert truncated["status"]=="TRUNCATED" and truncated["shortfall"]==1
    assert truncated["null_sorted"]==[0]*1000 and truncated["E"]==0

def test_order_statistic_and_strict_exceedance(monkeypatch):
    # Reference and hybrid are followed by 1000 replicates. Use a fixed sequence
    # with 949 zeros and 51 twos: the exact 950th is 2, not interpolated 0.1.
    calls=iter([set(),{"a","b"}]+[set()]*949+[{"a","b"}]*51)
    monkeypatch.setattr(c,"enrichment",lambda *a: next(calls))
    monkeypatch.setattr(c,"significance",lambda x:x)
    out=c.c2_run([sentinel("rs1")],["A"],["B"],{"1":[(1,"A",""),(2,"B","")]},{"A","B"},{},0)
    assert out["percentile95"]==2 and out["observed_F_g2g"]==2 and out["E"]==0

def test_empty_significant_sets_and_strict_fdr():
    a=[("P",.5,0)]
    b=[("P",.001,2)]
    empty,_=c.compare(a,a)
    one,_=c.compare(a,b)
    boundary,_=c.compare([("P",.05,1)],[("P",.05,1)])
    assert empty["informative"] is False and empty["J"] is None and empty["F"]==0
    assert one["informative"] is True and one["J"]==0 and one["F"]==1
    assert boundary["informative"] is False

def summary(**kw):
    s={"c":.6,"m_J":.9,"m_F":1,"tau":.1,"C2_available":True,"x":.1}
    return s|kw

@pytest.mark.parametrize("updates,outcome,rule",[
 ({"c":.599,"m_J":0,"x":1},"UNDERPOWERED",1),
 ({"m_J":.7999},"PASS",2),
 ({"tau":.25001},"PASS",2),
 ({"x":.25},"PASS",2),
 ({},"KILL",3),
 ({"m_J":.8,"tau":.25,"x":.24999},"INDETERMINATE",5),
 ({"C2_available":False,"x":None},"C2-UNAVAILABLE",4),
 ({"C2_available":False,"x":None,"m_J":.79},"PASS",2),
 ({"m_J":.89999},"INDETERMINATE",5),
 ({"m_F":1.0001},"INDETERMINATE",5),
 ({"tau":.10001},"INDETERMINATE",5),
 ({"x":.10001},"INDETERMINATE",5),
])
def test_decision_precedence_and_exact_boundaries(updates,outcome,rule):
    got=c.decision(summary(**updates))
    assert (got["outcome"],got["rule"])==(outcome,rule)

@pytest.mark.parametrize("available,outcome",[(False,"C2-UNAVAILABLE-FINAL"),(True,"INDETERMINATE-FINAL")])
def test_only_one_extension(available,outcome):
    got=c.decision(summary(m_J=.85,C2_available=available,x=.15 if available else None),extension=True)
    assert got["outcome"]==outcome and got["extension_required"] is False

def test_small_absent_c2_denominators_and_half_coverage():
    metrics=[{"informative":True,"J":1,"F":0}]*20
    c2=[{"status":"AVAILABLE","n1":1,"E":1}]*10+[{"status":"NO-DRIFT","n1":0,"E":None}]*10
    exact=c.summarize(metrics,c2)
    assert exact["C2_available"] and exact["d_x"]==10 and exact["x"]==1
    more=metrics+[{"informative":True,"J":1,"F":0}]
    fewer=c.summarize(more,c2+[{"status":"NO-DRIFT","n1":0,"E":None}])
    assert not fewer["C2_available"] and fewer["x"] is None
    absent=c.summarize(metrics,[{"status":"NO-DRIFT","n1":0,"E":None}]*20)
    assert absent["d_x"]==0 and absent["x"] is None
    all_empty=c.summarize([{"informative":False,"J":None,"F":0}]*20,c2)
    assert all_empty["n_i"]==0 and c.decision(all_empty)["outcome"]=="UNDERPOWERED"

def test_deletion_keeps_roster_and_full_c2_fixed():
    rows=[sentinel("rs1",1,1),sentinel("rs2",2,2)]
    assert c.delete_chromosome(rows,1)==rows[1:]
    metrics=[{"informative":False,"J":None,"F":0}]*20
    c2=[{"status":"AVAILABLE","n1":1,"E":1}]*20
    held={"d_x":20,"C2_available":True,"C2_exceedances":20,"x":1}
    got=c.summarize(metrics,c2,held)
    assert got["N"]==20 and got["u"]==20 and got["d_x"]==20 and got["x"]==1
    assert c.decision(got)["outcome"]=="UNDERPOWERED"

def test_refuse_missing_adoption(tmp_path):
    with pytest.raises(RuntimeError,match="missing pin"):
        gate.pin(tmp_path/"OWNER_ADOPTION_rev2_2026-11-07.txt",gate.ADOPTION)

def test_refuse_wrong_allocation():
    row={"id":"gwas_annotation_drift_pertrait_pilot","path":"08_GENOMICS/annotation-drift-gwas-pertrait-pilot",
         "status":"PARKED","empirical_execution_authorized":True,"provider_calls_authorized":False}
    with pytest.raises(RuntimeError,match="not ACTIVE"): gate.check_allocation(row)

def test_refuse_missing_input_pin(tmp_path):
    with pytest.raises(RuntimeError,match="missing pin"): gate.pin(tmp_path/"reactome.gmt","0"*64)

def test_refuse_altered_lock(tmp_path):
    p=tmp_path/"frozen_params.yaml";p.write_text("changed")
    with pytest.raises(RuntimeError,match="altered pin"):gate.pin(p,hashlib.sha256(b"original").hexdigest())


def test_baseline_symlink_identity_without_following_target(tmp_path):
    p=tmp_path/"run.sh"
    p.symlink_to("../../../00_ADMIN/tools/run.sh")
    assert gate.baseline_bytes(p,b"120000")==b"../../../00_ADMIN/tools/run.sh"


def test_fgsea_adapter_integration_and_determinism(tmp_path):
    import pipeline
    genes=[f"G{i:03d}" for i in range(100)]
    stats={g:float(100-i) for i,g in enumerate(genes)}
    pathways={"head":set(genes[:15]),"tail":set(genes[-15:])}
    first,qa=pipeline.fgsea_run(tmp_path/"one","v29",stats,pathways)
    second,qb=pipeline.fgsea_run(tmp_path/"two","v29",stats,pathways)
    assert first==second and qa==qb
    assert {r[0] for r in first}==set(pathways)
    assert qa["head"]<.05 and qa["tail"]<.05
    assert [(r[0],r[2]) for r in first]==[(r[0],15) for r in first]
