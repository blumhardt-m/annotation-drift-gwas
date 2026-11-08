#!/usr/bin/env python3
"""Single-run, checkpointed per-trait execution. No scientific CLI overrides."""
import argparse
import csv
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import core as c
import pre_auth_pilot as gate

ROOT = c.ROOT
OUT = ROOT / "results/pertrait_pilot"
D = gate.DEPOSIT
PYTHON = "/home/michaelblumhardt/miniconda3/bin/python3"
RSCRIPT = "/home/michaelblumhardt/miniconda3/envs/study3probe/bin/Rscript"
FGSEA_SOURCE = Path("/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-clean/scratch/study3_viability/run_fgsea.R")
MODEL = "Codex, based on GPT-6; version unavailable"
ENV = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
       "OPENBLAS_NUM_THREADS": "1", "PATH": str(Path(PYTHON).parent)+os.pathsep+os.environ["PATH"]}

def clean_json(x):
    if isinstance(x, float) and not math.isfinite(x): return None
    if isinstance(x, dict): return {k: clean_json(v) for k,v in x.items()}
    if isinstance(x, (list, tuple)): return [clean_json(v) for v in x]
    return x

def js(path, x):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean_json(x), indent=2, ensure_ascii=False)+"\n")

def table(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        writer=csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in clean_json(row).items()})

def read_table(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter="\t"))

def log(message):
    with (D/"commands.log").open("a") as f: f.write(time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())+" "+message+"\n")

def state(status, **kw):
    previous = json.loads((D/"STATE.json").read_text())
    js(D/"STATE.json",previous|{"model":MODEL,"status":status}|kw)

def code_identity():
    paths=sorted((ROOT/"scripts/pertrait_pilot").glob("*"))
    paths += sorted((ROOT/"tests/pertrait_pilot").glob("test*.py"))
    source_pins=json.loads((ROOT/"manifests/pertrait_pilot/ADDITIONAL_SOURCE_PINS.json").read_text())
    for name, h in source_pins.items(): gate.pin(name,h)
    return {str(p.relative_to(ROOT)):gate.sha(p) for p in paths if p.is_file()}|source_pins

def stage_gate(stage):
    ident=code_identity()
    checkpoint=OUT/"EXECUTION_BINDING.json"
    if checkpoint.exists():
        assert json.loads(checkpoint.read_text())["code"]==ident,"implementation changed after run binding"
    for label,args in [("legacy",["bash","pre_auth_check.sh"]),
                       ("pilot",[PYTHON,"scripts/pertrait_pilot/pre_auth_pilot.py","--stage",stage])]:
        name=f"{time.time_ns()}_{stage}_{label}.txt"
        result=subprocess.run(args,cwd=ROOT,env=ENV,capture_output=True,text=True)
        (D/name).write_text(result.stdout+result.stderr)
        log(json.dumps({"argv":args,"returncode":result.returncode,"log":name}))
        if result.returncode: raise RuntimeError(f"{stage} {label} gate refused; see {name}")
    if not checkpoint.exists():
        js(checkpoint,{"code":ident,"input_manifest":gate.MANIFEST_SHA,
                       "implementation_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()})
    state("RUNNING_"+stage.upper())

def freeze_roster(cutoff, groups, primary=None):
    path=ROOT/f"prereg/pertrait_roster_ge{cutoff}.tsv"
    names=c.roster(groups,cutoff,primary)
    if path.exists():
        assert [r["trait"] for r in read_table(path)]==names
        gate.pin(path,Path(str(path)+".sha256").read_text().split()[0])
        return names
    table(path,[{"trait_index":i,"trait":t,"sentinel_count":len(groups[t])} for i,t in enumerate(names)],
          ["trait_index","trait","sentinel_count"])
    h=gate.sha(path)
    Path(str(path)+".sha256").write_text(h+"  "+path.name+"\n")
    with (ROOT/"RUN_LOG.md").open("a") as f:
        f.write(f"\nMODEL: {MODEL}. THINKING LEVEL: high requested; actual setting unavailable.\n"
                f"Roster >= {cutoff} frozen before enrichment: N={len(names)}, SHA-256 {h}; code-point order, existing indices preserved.\n")
    log(f"Frozen roster {cutoff}, N={len(names)}, sha256={h}")
    print(f"ROSTER_FROZEN cutoff={cutoff} N={len(names)} sha256={h}",flush=True)
    return names

def roster_stage():
    stage_gate("roster")
    path=OUT/"sentinels_all.tsv"
    if path.exists():
        record=json.loads((OUT/"ROSTER_PREPARATION.json").read_text())
        gate.pin(path,record["sentinels_sha256"])
        groups={}
        for r in read_table(path):
            groups.setdefault(r["trait"],[]).append(c.Sentinel(r["trait"],r["rsid"],int(r["chrom"]),int(r["pos"]),Decimal(r["p"])))
        return groups,freeze_roster(100,groups)
    eligible=OUT/"eligible_regenerated.tsv"
    if not eligible.exists():
        c.parser03.parse_gwas_catalog(ROOT/"data/raw/gwas_catalog.tsv",eligible)
    gate.pin(eligible,gate.sha(ROOT/"data/processed/gwas_eligible.tsv"))
    with eligible.open() as f:
        groups,counts,excluded=c.thin_rows(csv.DictReader(f,delimiter="\t"))
    table(path,[{"trait":t,"rsid":s.rsid,"chrom":s.chrom,"pos":s.pos,"p":str(s.p)}
                for t in sorted(groups) for s in groups[t]],["trait","rsid","chrom","pos","p"])
    table(OUT/"excluded_associations.tsv",excluded,["trait","rsid","reason"])
    table(OUT/"trait_eligibility.tsv",[{"trait":t,"sentinel_count":len(groups[t]),"eligible_ge100":len(groups[t])>=100,
           "eligible_ge50":len(groups[t])>=50} for t in sorted(groups)])
    js(OUT/"ROSTER_PREPARATION.json",{"counts":counts,"sentinels_sha256":gate.sha(path),
                                     "ancestry_applied":False,"ancestry_reason":"Catalog lacks per-row ancestry"})
    return groups,freeze_roster(100,groups)

_CONTEXT=None
def context():
    global _CONTEXT
    if _CONTEXT is None:
        v29=ROOT/"data/processed/tss_tables/tss_v29.tsv"
        v45=ROOT/"data/processed/tss_tables/tss_v45.tsv"
        syms,universe=c.ora12.load_tss_tables(v29,v45)
        assert len(universe)==19315
        pathways,_,_=c.ora12.load_gmt_as_ensembl(ROOT/"data/raw/reactome.gmt",syms,15,500)
        _CONTEXT=(c.mapping04.load_tss_table(v29),c.mapping04.load_tss_table(v45),universe,pathways)
    return _CONTEXT

def output_checkpoint(folder, payload):
    files={p.relative_to(folder).as_posix():gate.sha(p) for p in folder.rglob("*") if p.is_file() and p.name!="DONE.json"}
    js(folder/"DONE.json",{"binding":gate.sha(OUT/"EXECUTION_BINDING.json"),"files":files,"payload":payload})
    return payload

def resume_checkpoint(folder):
    path=folder/"DONE.json"
    if not path.exists(): return None
    record=json.loads(path.read_text())
    assert record["binding"]==gate.sha(OUT/"EXECUTION_BINDING.json"),"checkpoint code/input binding differs"
    for name,h in record["files"].items():gate.pin(folder/name,h)
    return record["payload"]

def fgsea_run(folder, tag, stats, pathways):
    folder.mkdir(parents=True,exist_ok=True)
    stats_path=folder/f"count_vector_{tag}.tsv"
    table(stats_path,[{"gene":g,"count":v} for g,v in stats.items()],["gene","count"])
    pathway_path=folder/"pathways_ensembl.tsv"
    if not pathway_path.exists():
        pathway_path.write_text("".join(name+"\t"+",".join(sorted(genes))+"\n" for name,genes in pathways.items()))
    target=folder/f"fgsea_{tag}.tsv"
    assert not target.exists(),"attempted fgsea output already exists without complete checkpoint; preserve and stop"
    args=[RSCRIPT,"--vanilla",str(ROOT/"scripts/pertrait_pilot/run_fgsea.R"),str(folder),str(FGSEA_SOURCE),tag]
    result=subprocess.run(args,cwd=ROOT,env=ENV|{"TMPDIR":str(folder)},capture_output=True,text=True)
    (folder/f"fgsea_{tag}.log").write_text(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError(f"fgsea failed; see {folder}")
    rows=read_table(target)
    assert {r["pathway"] for r in rows}==set(pathways),"fgsea/ORA pathway family differs"
    tuples=[];qs={}
    nonzero={g for g,v in stats.items() if v>0}
    for r in rows:
        p,q=float(r["pval"]),float(r["padj"])
        if not math.isfinite(p) or not math.isfinite(q):
            raise RuntimeError("Nonfinite fgsea output: no substitute rule is authorized")
        name=r["pathway"]
        tuples.append((name,p,len(nonzero & pathways[name])))
        qs[name]=q
    return tuples,qs

def trait_enrichment(job):
    index,trait,rows,deleted=job
    folder=OUT/"traits"/f"{index:05d}"/("full" if deleted is None else f"delete_{deleted:02d}")
    saved=resume_checkpoint(folder)
    if saved is not None:
        assert saved["trait"]==trait
        return saved
    folder.mkdir(parents=True,exist_ok=True)
    a_tss,b_tss,universe,pathways=context()
    retained=rows if deleted is None else c.delete_chromosome(rows,deleted)
    a,b=c.mapped(retained,a_tss),c.mapped(retained,b_tss)
    table(folder/"mapping.tsv",[{"rsid":s.rsid,"chrom":s.chrom,"pos":s.pos,"p":str(s.p),"v29":x,"v45":y}
                                for s,x,y in zip(retained,a,b)],["rsid","chrom","pos","p","v29","v45"])
    oa,ob=c.enrichment(a,universe,pathways),c.enrichment(b,universe,pathways)
    for tag,values in [("v29",oa),("v45",ob)]:
        qs=c.ora12.bh_correction([x[1] for x in values])
        table(folder/f"ora_{tag}.tsv",[{"pathway":r[0],"p":r[1],"q":q,"overlap":r[2]} for r,q in zip(values,qs)])
    ga,qga=fgsea_run(folder,"v29",c.gene_statistics(retained,a,universe),pathways)
    gb,qgb=fgsea_run(folder,"v45",c.gene_statistics(retained,b,universe),pathways)
    metrics=[];flips=[]
    for method,aa,bb,qa,qb in [("ORA",oa,ob,None,None),("fgsea",ga,gb,qga,qgb)]:
        for alpha in (.05,.10):
            row,ff=c.compare(aa,bb,alpha,qa,qb)
            metadata={"trait_index":index,"trait":trait,"method":method,"alpha":alpha,"deleted_chromosome":deleted}
            metrics.append(metadata|row)
            flips.extend(metadata|f for f in ff)
    payload={"trait":trait,"index":index,"deleted_chromosome":deleted,"retained_sentinels":len(retained),"metrics":metrics}
    table(folder/"metrics.tsv",metrics)
    table(folder/"flips.tsv",flips,["trait_index","trait","method","alpha","deleted_chromosome","pathway","p_v29","p_v45","q_v29","q_v45","distance_v29","distance_v45","overlap_v29","overlap_v45","delta_genes","direction"])
    return output_checkpoint(folder,payload)

def trait_c2(job):
    index,trait,rows=job
    folder=OUT/"traits"/f"{index:05d}"/"c2"
    saved=resume_checkpoint(folder)
    if saved is not None:
        assert saved["trait"]==trait
        return saved
    folder.mkdir(parents=True,exist_ok=True)
    a_tss,b_tss,universe,pathways=context()
    result=c.c2_run(rows,c.mapped(rows,a_tss),c.mapped(rows,b_tss),a_tss,universe,pathways,index)
    values=result.pop("null_sorted")
    table(folder/"null_sorted.tsv",[{"order_statistic":i+1,"F":x} for i,x in enumerate(values)],["order_statistic","F"])
    result.update(trait=trait,trait_index=index)
    js(folder/"C2.json",result)
    return output_checkpoint(folder,result)

def parallel(function,jobs,label):
    # Independent traits; each trait owns its one RNG and outputs. No agents.
    results=[]
    pool=ProcessPoolExecutor(max_workers=4)
    try:
        pending={pool.submit(function,job):job for job in jobs}
        for future in as_completed(pending):
            result=future.result();results.append(result)
            print(f"{label} completed {len(results)}/{len(jobs)}",flush=True)
    except BaseException:
        for future in pending: future.cancel()
        for process in pool._processes.values(): process.terminate()
        pool.shutdown(wait=True,cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)
    return results

def metric(full):
    return next(r for r in full["metrics"] if r["method"]=="ORA" and r["alpha"]==.05)

def aggregate(names, extension=False):
    cutoff=50 if extension else 100
    folder=OUT/f"ge{cutoff}";folder.mkdir(exist_ok=True)
    if (folder/"OUTCOME.json").exists():
        return json.loads((folder/"OUTCOME.json").read_text())
    full=[resume_checkpoint(OUT/"traits"/f"{i:05d}"/"full") for i in range(len(names))]
    c2=[resume_checkpoint(OUT/"traits"/f"{i:05d}"/"c2") for i in range(len(names))]
    assert all(full) and all(c2)
    primary=[metric(r) for r in full]
    summary=c.summarize(primary,c2)
    outcome=c.decision(summary,extension)
    deletion_rows=[];all_deletion_metrics=[];ranges=[]
    changed=set()
    for chrom in range(1,23):
        deletion=[resume_checkpoint(OUT/"traits"/f"{i:05d}"/f"delete_{chrom:02d}") for i in range(len(names))]
        assert all(deletion)
        metrics=[metric(r) for r in deletion]
        s=c.summarize(metrics,c2,summary)
        descriptive=c.decision(s,extension)
        deletion_rows.append({"deleted_chromosome":chrom}|s|descriptive|
                             {"outcome_differs":descriptive["outcome"]!=outcome["outcome"]})
        for i,(before,after) in enumerate(zip(primary,metrics)):
            if before["informative"]!=after["informative"]:changed.add(i)
        all_deletion_metrics.extend(m for r in deletion for m in r["metrics"])
    full_metrics=[m for r in full for m in r["metrics"]]
    for row in full_metrics:
        subset=[m for m in all_deletion_metrics if m["trait_index"]==row["trait_index"] and m["method"]==row["method"] and m["alpha"]==row["alpha"]]
        for key in ("J","F","rho","sig_v29","sig_v45","top25_overlap","top50_overlap"):
            finite=[m for m in subset if m[key] is not None and math.isfinite(m[key])]
            lo=min((m[key] for m in finite),default=None);hi=max((m[key] for m in finite),default=None)
            ranges.append({"trait_index":row["trait_index"],"trait":row["trait"],"method":row["method"],"alpha":row["alpha"],"metric":key,
                           "full":row[key],"minimum":lo,"maximum":hi,"min_deletions":[m["deleted_chromosome"] for m in finite if m[key]==lo],
                           "max_deletions":[m["deleted_chromosome"] for m in finite if m[key]==hi],
                           "undefined_deletions":[m["deleted_chromosome"] for m in subset if m[key] is None]})
    summary_ranges=[]
    for key in ("c","u","n_i","m_J","m_F","tau"):
        finite=[r for r in deletion_rows if r[key] is not None]
        lo=min((r[key] for r in finite),default=None);hi=max((r[key] for r in finite),default=None)
        summary_ranges.append({"metric":key,"full":summary[key],"minimum":lo,"maximum":hi,
                               "min_deletions":[r["deleted_chromosome"] for r in finite if r[key]==lo],
                               "max_deletions":[r["deleted_chromosome"] for r in finite if r[key]==hi]})
    changes=sum(r["outcome_differs"] for r in deletion_rows)
    result={"model":MODEL,"cutoff":cutoff,"summary":summary,**outcome,"conditional_deletion_outcome_changes":changes,
            "FRAGILE_conditional_on_full_data_C2":changes>=3,"traits_with_informative_status_change":len(changed)}
    table(folder/"per_trait.tsv",full_metrics)
    table(folder/"C2.tsv",c2)
    table(folder/"deletion_summaries.tsv",deletion_rows)
    table(folder/"deletion_per_trait.tsv",all_deletion_metrics)
    table(folder/"deletion_trait_ranges.tsv",ranges)
    table(folder/"deletion_summary_ranges.tsv",summary_ranges)
    flips=[]
    for i in range(len(names)):flips.extend(read_table(OUT/"traits"/f"{i:05d}"/"full"/"flips.tsv"))
    table(folder/"flipped_pathways.tsv",flips,["trait_index","trait","method","alpha","deleted_chromosome","pathway","p_v29","p_v45","q_v29","q_v45","distance_v29","distance_v45","overlap_v29","overlap_v45","delta_genes","direction"])
    js(folder/"OUTCOME.json",result)
    return result

def run_roster(names,groups,extension):
    stage_gate("enrichment")
    parallel(trait_enrichment,[(i,t,groups[t],None) for i,t in enumerate(names)],"FULL")
    stage_gate("c2")
    parallel(trait_c2,[(i,t,groups[t]) for i,t in enumerate(names)],"C2")
    stage_gate("deletion")
    parallel(trait_enrichment,[(i,t,groups[t],chrom) for chrom in range(1,23) for i,t in enumerate(names)],"DELETION")
    stage_gate("decision")
    return aggregate(names,extension)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--roster-only",action="store_true",help="Preparation checkpoint; does not enrich.")
    args=parser.parse_args()
    groups,primary=roster_stage()
    if args.roster_only:return
    state("EMPIRICAL_RUNNING",empirical_run_started=True)
    result=run_roster(primary,groups,False)
    print("PRIMARY_OUTCOME "+json.dumps(result),flush=True)
    if result["extension_required"]:
        stage_gate("roster")
        extended=freeze_roster(50,groups,primary)
        result=run_roster(extended,groups,True)
    else:
        print("EXTENSION_NOT_TRIGGERED",flush=True)
    js(D/"OUTCOME.json",result)
    state("EMPIRICAL_COMPLETE",empirical_run_complete=True,outcome=result)
    print("EMPIRICAL_COMPLETE "+json.dumps(result),flush=True)

if __name__=="__main__":
    try:main()
    except BaseException as exc:
        record={"model":MODEL,"error":repr(exc),"traceback":traceback.format_exc(),"time":time.time()}
        js(D/f"FAILURE_{time.time_ns()}.json",record)
        state("STOPPED",error=repr(exc))
        raise
