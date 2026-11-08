"""v1.2 regression tests; fixtures are never empirical-run inputs."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
from decimal import Decimal
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts/pertrait_pilot"))
import core as c
import pipeline as p
import pre_auth_pilot as gate

@pytest.fixture
def bound_output(tmp_path,monkeypatch):
    monkeypatch.setattr(p,"OUT",tmp_path)
    (tmp_path/"EXECUTION_BINDING.json").write_text("synthetic binding")
    return tmp_path

def test_saved_v11_missing_table_is_trait_not_estimable():
    fixture=Path(__file__).parent/"fixtures/v11_fgsea_v29_trait00002.tsv"
    assert gate.sha(fixture)=="645beced3227d3eb952f688e20a1c069b3aa25c93b0af1e16a8babaa34d96bc5"
    rows=p.read_table(fixture)
    pathways={r["pathway"]:set() for r in rows}
    a,qa=p.consume_ranked(rows,pathways,{})
    b=[(name,.001,n) for name,_,n in a];qb={name:.001 for name,_,_ in b}
    comparison,flips=p.ranked_comparison(a,b,qa,qb,.05)
    assert comparison["ranked_status"]=="NOT_ESTIMABLE"
    assert comparison["ranked_na_pathways"]=={"v29":2,"v45":0}
    assert all(comparison[k] is None for k in ("J","F","rho","top25_overlap","top50_overlap","sig_v29","sig_v45"))
    assert flips==[] and len(a)==1369
    assert sum(math.isnan(q) for q in qa.values())==2
    assert gate.sha(fixture)=="645beced3227d3eb952f688e20a1c069b3aa25c93b0af1e16a8babaa34d96bc5"

def test_log2err_na_alone_is_estimable():
    rows=[dict(pathway="P",pval="0.01",padj="0.01",NES="1.5",log2err="NA",size="15")]
    a,qa=p.consume_ranked(rows,{"P":set()},{})
    comparison,flips=p.ranked_comparison(a,a,qa,qa,.05)
    assert comparison["ranked_status"]=="ESTIMABLE" and comparison["J"]==1
    assert comparison["ranked_na_pathways"]=={"v29":0,"v45":0} and not flips

@pytest.mark.parametrize("field,value",[("pval","oops"),("padj","nan"),("pval","Inf"),("NES","unknown"),("log2err","missing"),("size","NA")])
def test_unrecognized_values_fail_closed(field,value):
    row=dict(pathway="P",pval="0.01",padj="0.02",NES="1.5",log2err="0.2",size="15")
    row[field]=value
    with pytest.raises(ValueError):p.consume_ranked([row],{"P":set()},{})

def test_positive_fgsea_large_zero_tail_and_small_signal(bound_output):
    genes=[f"G{i:04d}" for i in range(1000)]
    stats={g:float(30-i) if i<30 else 0. for i,g in enumerate(genes)}
    pathways={"positive":set(genes[:15]),"zero_tail":set(genes[300:700])}
    rows,qs=p.fgsea_run(bound_output/"zero_case","v29",stats,pathways)
    assert all(math.isfinite(v) for _,v,_ in rows) and all(math.isfinite(v) for v in qs.values())
    assert qs["positive"]<.05 and qs["zero_tail"]>=.05
    full=p.read_table(bound_output/"zero_case/fgsea_v29_all_values.tsv")
    assert set(full[0]) >= {"pval","padj","ES","NES","log2err","size","leadingEdge"}
    assert (bound_output/"zero_case/fgsea_v29_all_values.rds").is_file()

def test_ranked_failure_preserves_primary_and_deletions_do_not_call_fgsea(bound_output,monkeypatch):
    genes=[f"G{i:03d}" for i in range(100)]
    tss={"1":[(1,genes[0],genes[0])],"2":[(2,genes[1],genes[1])]}
    pathways={"P":set(genes[:15])}
    monkeypatch.setattr(p,"context",lambda:(tss,tss,set(genes),pathways))
    rows=[c.Sentinel("T","rs1",1,1,Decimal("1e-10")),c.Sentinel("T","rs2",2,2,Decimal("1e-10"))]
    primary=p.trait_enrichment((0,"T",rows,None))
    folder=bound_output/"traits/00000/full"
    before={str(x.relative_to(folder)):x.read_bytes() for x in folder.rglob("*") if x.is_file()}
    calls=[]
    def fail(*args):
        calls.append(args[1]);raise p.RankedExecutionError("synthetic ranked process failure")
    monkeypatch.setattr(p,"fgsea_run",fail)
    ranked=p.trait_ranked((0,"T",rows))
    assert ranked["metrics"][0]["ranked_status"]=="FAILED" and calls==["v29","v45"]
    assert before=={str(x.relative_to(folder)):x.read_bytes() for x in folder.rglob("*") if x.is_file()}
    assert p.resume_checkpoint(folder)==p.clean_json(primary)
    null=p.trait_c2((0,"T",rows));assert null["status"]=="NO-DRIFT"
    deleted=p.trait_enrichment((0,"T",rows,1));assert {r["method"] for r in deleted["metrics"]}=={"ORA"}
    assert calls==["v29","v45"]
    summary=c.summarize([p.metric(primary)],[null]);assert c.decision(summary)["outcome"]=="C2-UNAVAILABLE"

def test_stage_gate_refuses_old_binding_before_dispatch(bound_output,monkeypatch):
    (bound_output/"EXECUTION_BINDING.json").write_text(json.dumps({"code":{},"input_manifest":gate.MANIFEST_SHA}))
    monkeypatch.setattr(p,"code_identity",lambda:{})
    def unexpected(*a,**kw):pytest.fail("old binding reached a command")
    monkeypatch.setattr(p.subprocess,"run",unexpected)
    with pytest.raises(RuntimeError,match="v1.2 execution authority binding differs"):p.stage_gate("enrichment")

def test_binding_requires_every_v12_authority_field():
    good={"code":{"a":"b"},"authority":gate.execution_authority()}
    p.verify_execution_binding(good,{"a":"b"})
    for key in good["authority"]:
        bad={"code":good["code"],"authority":dict(good["authority"])};bad["authority"].pop(key)
        with pytest.raises(RuntimeError):p.verify_execution_binding(bad,{"a":"b"})

def test_refuse_missing_rev3_instrument(tmp_path):
    with pytest.raises(RuntimeError,match="missing pin.*OWNER_ADOPTION_rev3"):
        gate.check_v12_authority(tmp_path,tmp_path)

def test_refuse_v11_governance_binding(tmp_path):
    (tmp_path/"authority").mkdir()
    (tmp_path/"authority/OWNER_ADOPTION_rev3_2026-11-08.txt").write_bytes((gate.DEPOSIT/"authority/OWNER_ADOPTION_rev3_2026-11-08.txt").read_bytes())
    for name,h in gate.V12_GOVERNANCE.items():
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((gate.ROOT/name).read_bytes())
    (tmp_path/"prereg/pertrait_pilot_prereg.sha256").write_bytes((gate.ROOT/"prereg/pertrait_pilot_prereg_v1.1.sha256").read_bytes())
    with pytest.raises(RuntimeError,match="altered pin.*pertrait_pilot_prereg.sha256"):
        gate.check_v12_authority(tmp_path,tmp_path)
