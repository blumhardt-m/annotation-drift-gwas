#!/usr/bin/env python3
"""Fail-closed authorization gate for the adopted per-trait pilot."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_ROOT = Path("/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot")
DEPOSIT = Path("/home/michaelblumhardt/agent_work/annotation_drift_pertrait_pilot_execution_20261108")
V11_DEPOSIT = DEPOSIT.with_name("annotation_drift_pertrait_pilot_execution_20261107")
PACKAGE = Path("/home/michaelblumhardt/agent_work/annotation_drift_pilot_activation_preflight_rev2_20261107")
BASE = "847fd7cd98f68695b709ae155a1d2573ff7b7ee7"
ADOPTION = "2c7ccecb68d00e29c66e7f214a2860316c1e6a8a9392c6275ce30bb1f61f117e"
PACKAGE_SHA = "4103b84fef3d57cacbdb259b3e6ba82b5164222b419326ee6f8a7686222d67b3"
MANIFEST_SHA = "dc838286666ee2d35c88282a05a991570230373b627d3189e84457ece61f1288"

MAIN_HEAD = "33245c08cf8cbc04263f3de2c85f5cade7626f0f"
ADOPTION_V3 = "133b89edec0edde6260ea12b0a6f21ffddf6c533b5a217fdf5f6179fe647b6b7"
PREREG_V12 = "f76176bd9ba0b06c40366717f46bbe5639702e8a259cc44b46875198aeb475cd"
SIDECAR_V12 = "ca196a010d13470442d8806155fb8624eba92566042dd176896db91215fad50d"
BRACKET_SHA = "3b6caf3c738c21898eb70946b79990002b925987a4af6df3cd1bdb1041d228c5"
V12_GOVERNANCE = {
    "DECISIONS.md": "655879f8a1521b92a6989c990eecfc031338fdca60113f90ab8db44823fbc334",
    "PILOT_EXECUTION_AUTHORITY.md": "c6b64d9c1db98c2b05043861e63bce09c27e258a35152ede7f60ce43ce9f8f8e",
    "prereg/pertrait_pilot_prereg.sha256": SIDECAR_V12,
    "prereg/pertrait_pilot_prereg_v1.2.md": PREREG_V12,
    "prereg/pertrait_pilot_prereg_v1.1.sha256": "a92ab999e7269d1b25be8127f5ef63ad80780d2c95f47a52a7c57bc03a641b13",
}

def check_v12_authority(root, deposit):
    pin(deposit / "authority/OWNER_ADOPTION_rev3_2026-11-08.txt", ADOPTION_V3)
    for name, expected in V12_GOVERNANCE.items():
        pin(root / name, expected)

def execution_authority():
    return {"schema": "pertrait-pilot-v1.2", "prereg": PREREG_V12,
            "controlling_sidecar": SIDECAR_V12, "adoption_rev3": ADOPTION_V3,
            "input_manifest": MANIFEST_SHA, "reactivation_bracket": BRACKET_SHA}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def pin(path, expected):
    require(Path(path).is_file(), f"missing pin: {path}")
    require(sha(path) == expected, f"altered pin: {path}")

def baseline_bytes(path, mode):
    if mode == b"120000":
        require(path.is_symlink(), f"missing baseline symlink: {path}")
        return os.fsencode(os.readlink(path))
    require(path.is_file(), f"missing baseline tracked file: {path}")
    return path.read_bytes()

def check_allocation(row):
    require(row.get("id") == "gwas_annotation_drift_pertrait_pilot", "wrong allocation identity")
    require(row.get("path") == "08_GENOMICS/annotation-drift-gwas-pertrait-pilot", "wrong allocation path")
    require(row.get("status") == "ACTIVE", "pilot not ACTIVE")
    require(row.get("empirical_execution_authorized") is True, "empirical permission absent")
    require(row.get("provider_calls_authorized") is False, "provider flag differs")

def command(args):
    r = subprocess.run(args, cwd=ROOT, text=True, capture_output=True)
    require(r.returncode == 0, f"command failed {args}: {r.stdout}{r.stderr}")
    return r.stdout.strip()

def measure_software(manifest):
    py = manifest["software"]["python"]
    require(platform.python_version() == py["version"], "Python version mismatch")
    require(Path(sys.executable).resolve() == Path(py["executable"]).resolve(), "Python executable mismatch")
    packages = {k: importlib.metadata.version(k) for k in py["packages"]}
    require(packages == py["packages"], "Python package version mismatch")
    rs = manifest["software"]["Rscript"]
    pin(rs["path"], rs["sha256"])
    expression = 'cat(paste(as.character(getRversion()),as.character(packageVersion("fgsea")),as.character(packageVersion("data.table")),as.character(packageVersion("BiocParallel")),sep="\\n"))'
    versions = dict(zip(["R", "fgsea", "data.table", "BiocParallel"],
                        command([rs["path"], "--vanilla", "-e", expression]).splitlines()))
    require(versions == rs["expected_versions"], "R package version mismatch")
    return {"python": platform.python_version(), "packages": packages, "R": versions}

# Revision-4 administrative authority. Scientific execution_authority() is unchanged.
ADOPTION_V4 = '01c5624e84a553a6ab1ca7a61e35dc67bde3bd706fbb5ab9e4c4bb0326820843'
SCOPE_APPROVAL_V4 = 'a5a8f858732e6b661adb3688f94a674c5ae9856b3ed0f41c702245060c18b368'
ORIGINAL_BINDING_SHA = 'c2410ed014044c2025b0d9095f48fd42b085ec9468d01080d680d8925f6a2e1d'
ORIGINAL_FINALIZERS_SHA = '16540e70b7edc2fef52f7f4b050f65721054f67dcfb80816bb41ece755de6316'
BRACKET_ZERO_SHA = '90827e70580191ad21f51f3bbb2e21c5e5d1d962801e67463f350cfd0c58f675'
STOP_ARCHIVE = 'stop_records/stop1_pdg036_2026-11-10'
ADMIN_ELIGIBLE = {'scripts/pertrait_pilot/pre_auth_pilot.py', 'scripts/pertrait_pilot/pipeline.py',
                  'deposit:run_completion.py', 'deposit:completion_supervisor.py',
                  'deposit:terminal_transition.py', 'deposit:deposit_results.py'}

def administrative_authority():
    pin(DEPOSIT/'authority/OWNER_ADOPTION_rev4_2026-11-11.txt', ADOPTION_V4)
    pin(DEPOSIT/'authority/OWNER_SCOPE_APPROVAL_rev4_2026-11-11.txt', SCOPE_APPROVAL_V4)
    pin(DEPOSIT/'authority/REV4_PENDING_SCOPE.txt', '4d38e5fd7331e4d7f65a5693498713f4f302b79345237853392ed6ee8d23f105')
    return {ADOPTION_V4, SCOPE_APPROVAL_V4}

def original_binding():
    path=ROOT/'results/pertrait_pilot/EXECUTION_BINDING.json'
    pin(path, ORIGINAL_BINDING_SHA)
    return json.loads(path.read_text())

def amendment_roots():
    binding=original_binding()
    roots=dict(binding['code'])
    roots['deposit:run_completion.py']=binding['launcher_sha256']
    path=DEPOSIT/STOP_ARCHIVE/'FINALIZER_IDENTITIES.json'
    pin(path, ORIGINAL_FINALIZERS_SHA)
    pin(DEPOSIT/'FINALIZER_IDENTITIES.json', ORIGINAL_FINALIZERS_SHA)
    roots.update({'deposit:'+n:h for n,h in json.loads(path.read_text()).items()})
    return roots

def validate_amendment_entries(entries, roots, authorities):
    import re
    require(isinstance(entries,list), 'invalid amendment entries')
    effective=dict(roots)
    for i,entry in enumerate(entries,1):
        name=entry.get('file_path')
        require(entry.get('sequence')==i, 'amendment sequence differs')
        require(name in ADMIN_ELIGIBLE, 'ineligible administrative amendment: '+str(name))
        require(entry.get('old_sha256')==effective.get(name), 'amendment old hash does not chain: '+str(name))
        require(entry.get('authorizing_instrument_sha256') in authorities, 'unauthenticated amendment authority')
        require(isinstance(entry.get('new_sha256'),str) and re.fullmatch('[0-9a-f]{64}',entry['new_sha256']), 'invalid amended hash')
        require(bool(entry.get('timestamp')) and bool(entry.get('purpose')), 'incomplete amendment record')
        effective[name]=entry['new_sha256']
    return effective

def effective_bindings():
    authorities=administrative_authority()
    root=ROOT/'results/pertrait_pilot'
    ledger=root/'EXECUTION_BINDING_AMENDMENTS.json'
    sidecar=root/'EXECUTION_BINDING_AMENDMENTS.sha256'
    require(sidecar.is_file(), 'missing amendment authentication')
    expected=sidecar.read_text().strip().split()
    require(len(expected)==2 and expected[1]==ledger.name, 'invalid amendment sidecar')
    pin(ledger,expected[0])
    state=json.loads((DEPOSIT/'STATE.json').read_text())
    require(state.get('execution_binding_amendments_sha256')==expected[0], 'amendment STATE identity differs')
    require(state.get('execution_binding_amendments_sidecar_sha256')==sha(sidecar), 'amendment sidecar STATE identity differs')
    require(expected[0] in (ROOT/'RUN_LOG.md').read_text(), 'amendment not recorded in RUN_LOG')
    record=json.loads(ledger.read_text())
    require(record.get('schema')=='pertrait-pilot-administrative-amendments-v1', 'amendment schema differs')
    return validate_amendment_entries(record['entries'], amendment_roots(), authorities)

def amended_code(code):
    effective=effective_bindings()
    return {name:effective[name] if name in ADMIN_ELIGIBLE else value for name,value in code.items()}

def control_pin(name):
    effective=effective_bindings()
    key='deposit:'+name
    require(key in effective, 'unbound control script')
    pin(DEPOSIT/name,effective[key])

def pipeline_protected_bytes(data):
    import ast
    tree=ast.parse(data)
    lines=data.splitlines(keepends=True)
    # Mask only the two authorized verification helpers; all other bytes remain exact.
    allowed={'verify_execution_binding','stage_gate'}
    intervals=[]
    for node in tree.body:
        if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)) and node.name in allowed:
            require(not node.decorator_list,'verification helper acquired decorator')
            intervals.append((node.lineno-1,node.end_lineno,node.name))
    require({n for _,_,n in intervals}==allowed,'verification helper missing')
    for start,end,name in reversed(intervals):lines[start:end]=[('# administrative helper: '+name+'\n').encode()]
    return b''.join(lines)

def verify_pipeline_scope():
    name='scripts/pertrait_pilot/pipeline.py'
    baseline=DEPOSIT/'worktree_snapshot'/name
    pin(baseline, original_binding()['code'][name])
    require(pipeline_protected_bytes(baseline.read_bytes())==pipeline_protected_bytes((ROOT/name).read_bytes()),
            'computational or other nonverification pipeline bytes changed')

def verify_bound_files():
    effective=effective_bindings()
    for name,expected in effective.items():
        path=DEPOSIT/name[8:] if name.startswith('deposit:') else ROOT/name
        pin(path,expected)
    verify_pipeline_scope()
    return effective

def admin_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.admin.tmp')
    with temp.open('w') as stream:
        json.dump(value,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(temp,path)

def pilot_row(data):
    import re
    match=re.search(rb'^  - id: gwas_annotation_drift_pertrait_pilot\r?\n',data,re.M)
    require(match is not None,'pilot row absent')
    end=re.search(rb'^  - id:',data[match.end():],re.M)
    return data[match.start():match.end()+end.start() if end else len(data)]

def portfolio_names():
    return ['00_ADMIN/PORTFOLIO_STATUS.yaml','00_ADMIN/tools/validate_portfolio_status.py',
            '00_ADMIN/DISCOVERY_GOVERNANCE_DECISIONS.md','PORTFOLIO_GOVERNANCE.md']

def bracket_bytes(folder,record):
    names=portfolio_names()
    require(set(record['after'])==set(names),'bracket file set differs')
    data={}
    for name in names:
        item=record['after'][name]
        require(item['path']==str(EXPECTED_ROOT.parents[1]/name),'bracket installed path differs')
        pin(folder/name,item['sha256'])
        data[name]=(folder/name).read_bytes()
        require(len(data[name])==item['bytes'],'bracket size differs')
    return data

def check_absorption(previous,current,terminal=False):
    import ast,re,yaml
    reg,val,ledger,gov=portfolio_names()
    require(len(yaml.safe_load(current[reg])['projects'])==16,'member count changed')
    original=(DEPOSIT/'portfolio/brackets/0'/reg).read_bytes()
    require(pilot_row(current[reg])==pilot_row(original),'pilot row changed')
    require(current[gov]==previous[gov],'governance changed')
    require(current[ledger].startswith(previous[ledger]) and len(current[ledger])>len(previous[ledger]),'ledger not append-only with new entries')
    suffix=current[ledger][len(previous[ledger]):].decode()
    headings=[line for line in suffix.splitlines() if re.match(r'^(?:#+\s*)?PDG-\d+\s+[—–-]\s+',line)]
    require(bool(headings),'no appended ledger headings')
    changed=current[val]!=previous[val]
    require(not changed or terminal,'validator changed')
    if changed:
        assignments={}
        for node in ast.parse(current[val]).body:
            if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in {'CURRENT_OWNER_AUTHORIZED_LIVE','CURRENT_LIVE_PROJECT_COUNT','OWNER_ENUMERATED_MEMBERSHIP_COUNT'}:
                assignments[node.targets[0].id]=ast.literal_eval(node.value)
        require(assignments=={'CURRENT_OWNER_AUTHORIZED_LIVE':{'gwas_annotation_drift_pertrait_pilot':'ACTIVE','historical_replay_validation':'BLOCKED_COMPLETABLE_PREREQUISITE','bh_xrb_spectroscopic_outburst_inventory_stage0':'BLOCKED_COMPLETABLE_PREREQUISITE'},'CURRENT_LIVE_PROJECT_COUNT':3,'OWNER_ENUMERATED_MEMBERSHIP_COUNT':16},'terminal validator live mirror changed')
    return headings,changed

def accepted_bracket(terminal=False):
    administrative_authority()
    root=DEPOSIT/'portfolio/brackets'
    pin(DEPOSIT/'portfolio/REACTIVATION_BRACKET.json',BRACKET_SHA)
    activation=json.loads((DEPOSIT/'portfolio/REACTIVATION_BRACKET.json').read_text())
    pin(root/'0/BRACKET.json',BRACKET_ZERO_SHA)
    folders=sorted((p for p in root.iterdir() if p.is_dir() and p.name.isdecimal()),key=lambda p:int(p.name))
    require([p.name for p in folders]==[str(i) for i in range(len(folders))] and bool(folders),'bracket sequence gap')
    previous_hash=None;previous=None
    for index,folder in enumerate(folders):
        path=folder/'BRACKET.json';record=json.loads(path.read_text())
        require(record['sequence']==index and record['predecessor_sha256']==previous_hash,'bracket predecessor chain broken')
        data=bracket_bytes(folder,record)
        if index==0:
            require(record['after']==activation['after'] and record['reactivation_bracket_sha256']==BRACKET_SHA,'bracket zero differs from reactivation')
        else:
            require(record['authorizing_instrument_sha256']==ADOPTION_V4,'bracket authority differs')
            headings,changed=check_absorption(previous,data,terminal=terminal)
            require(headings==record['appended_ledger_headings'],'bracket headings differ')
            require(bool(record.get('terminal_validator_exception',False))==changed,'validator exception record differs')
        previous=data;previous_hash=sha(path)
    return record,data,previous_hash

def validate_captured_portfolio(folder):
    reg,val,_,_=portfolio_names()
    result=subprocess.run([sys.executable,str(folder/val),str(folder/reg)],capture_output=True,text=True,
                          env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
    require(result.returncode==0 and result.stdout.startswith('VALID: 16 entries, 3 live'),
            'captured portfolio validation failed: '+result.stdout+result.stderr)
    return result.stdout

def refresh_portfolio_locked(terminal=False):
    """Caller owns the shared flock continuously; this helper never acquires/releases it."""
    import datetime,difflib,tempfile,shutil
    record,previous,previous_hash=accepted_bracket(terminal=terminal)
    root=DEPOSIT/'portfolio/brackets';research=EXPECTED_ROOT.parents[1]
    current={name:(research/name).read_bytes() for name in portfolio_names()}
    if current==previous:
        validate_captured_portfolio(root/str(record['sequence']))
        return record
    headings,changed=check_absorption(previous,current,terminal=terminal)
    number=record['sequence']+1
    scratch=Path(tempfile.mkdtemp(prefix='.capture-',dir=root))
    for name,data in current.items():
        path=scratch/name;path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('wb') as stream:stream.write(data);stream.flush();os.fsync(stream.fileno())
    stdout=validate_captured_portfolio(scratch)
    require(current=={name:(research/name).read_bytes() for name in portfolio_names()},'portfolio changed while flock held')
    updated={'sequence':number,'predecessor_sha256':previous_hash,'authorizing_instrument_sha256':ADOPTION_V4,
             'captured_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'lock':str(research/'00_ADMIN/.portfolio_governance_write.lock'),'lock_timeout_seconds':60,'lock_held_continuously':True,
             'after':{name:{'path':str(research/name),'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)} for name,data in current.items()},
             'register_diff':''.join(difflib.unified_diff(previous[portfolio_names()[0]].decode().splitlines(True),current[portfolio_names()[0]].decode().splitlines(True),fromfile='predecessor',tofile='accepted')),
             'appended_ledger_headings':headings,'validator_stdout':stdout,'conditions_a_through_e':not changed,'terminal_validator_exception':changed}
    admin_json(scratch/'BRACKET.json',updated)
    os.rename(scratch,root/str(number))
    path=root/str(number)/'BRACKET.json'
    state=json.loads((DEPOSIT/'STATE.json').read_text())
    state['accepted_bracket']={'sequence':number,'path':str(path),'sha256':sha(path),'complete':True}
    state['bracket']=updated['after'];state['bracket_file']=str(path.relative_to(DEPOSIT))
    admin_json(DEPOSIT/'STATE.json',state)
    line=updated['captured_at_utc']+' revision-4 bracket refresh '+str(number)+' sha256 '+sha(path)+'; '+ '; '.join(headings)+'\n'
    # The deposit journal remains append-only even after repository finalization.
    logpath=DEPOSIT/'RUN_LOG.md' if terminal else ROOT/'RUN_LOG.md'
    with logpath.open('a') as stream:stream.write('\n'+line);stream.flush();os.fsync(stream.fileno())
    with (DEPOSIT/'commands.log').open('a') as stream:stream.write(line);stream.flush();os.fsync(stream.fileno())
    return updated

def refresh_portfolio():
    import fcntl
    lockpath=EXPECTED_ROOT.parents[1]/'00_ADMIN/.portfolio_governance_write.lock'
    with lockpath.open('rb') as lock:
        # flock(1) blocks with a 60-second timeout; inherited descriptor retains the lock.
        result=subprocess.run(['flock','--exclusive','--timeout','60',str(lock.fileno())],pass_fds=(lock.fileno(),))
        require(result.returncode==0,'portfolio flock timeout')
        try:return refresh_portfolio_locked()
        finally:fcntl.flock(lock,fcntl.LOCK_UN)


def check(stage):
    import yaml
    require(ROOT == EXPECTED_ROOT, "wrong worktree")
    for ancestor in set(ROOT.parents) | set(DEPOSIT.parents):
        require(not (ancestor / "AGENTS.md").exists(), f"external ancestor AGENTS.md: {ancestor}")
    check_v12_authority(ROOT, DEPOSIT)
    verify_bound_files()
    pin(DEPOSIT / "authority/OWNER_ADOPTION_rev2_2026-11-07.txt", ADOPTION)
    pin(V11_DEPOSIT / "SHA256SUMS", "5310964053138e146daf9281c50f1992e5dc8b79923a6311586b29da45e53483")
    for line in (V11_DEPOSIT / "SHA256SUMS").read_text().splitlines():
        h, name = line.split("  ", 1)
        pin(V11_DEPOSIT / name, h)
    pin(PACKAGE / "SHA256SUMS", PACKAGE_SHA)
    entries = (PACKAGE / "SHA256SUMS").read_text().splitlines()
    require(len(entries) == 132, "package entry count")
    for line in entries:
        h, name = line.split("  ", 1)
        pin(PACKAGE / name, h)
    require(command(["git", "rev-parse", "--show-toplevel"]) == str(ROOT), "wrong git root")
    require(command(["git", "branch", "--show-current"]) == "pertrait-pilot", "wrong branch")
    require(command(["git", "rev-parse", "refs/heads/main"]) == MAIN_HEAD, "main changed")
    command(["git", "merge-base", "--is-ancestor", BASE, "HEAD"])
    bracket = refresh_portfolio()
    ids = json.loads((PACKAGE / "proposals_rev2/repository/IDENTITIES.json").read_text())
    for name, item in ids["candidates"].items():
        pin(ROOT / name, V12_GOVERNANCE.get(name, item["new"]["sha256"]))
    # Compare original tracked content; only adopted candidates and append-only log differ.
    original = subprocess.check_output(["git", "ls-tree", "-r", "-z", BASE], cwd=ROOT).split(b"\0")
    for entry in original:
        if not entry:
            continue
        meta, name_bytes = entry.split(b"\t", 1)
        name = name_bytes.decode()
        if name in ids["candidates"] or name in V12_GOVERNANCE:
            continue
        content = subprocess.check_output(["git", "cat-file", "blob", meta.split()[2]], cwd=ROOT)
        path = ROOT / name
        actual = baseline_bytes(path, meta.split()[0])
        if name == "RUN_LOG.md":
            require(actual.startswith(content), "RUN_LOG baseline prefix altered")
        else:
            require(actual == content, f"altered baseline lock: {name}")
    mp = ROOT / "manifests/pertrait_pilot/INPUT_MANIFEST_rev2.json"
    pin(mp, MANIFEST_SHA)
    manifest = json.loads(mp.read_text())
    for name, item in manifest["inputs"].items():
        expected = ids["candidates"][name]["new"]["sha256"] if name in ids["candidates"] else item["sha256"]
        pin(ROOT / name, V12_GOVERNANCE.get(name, expected))
    runtime = json.loads((ROOT / "manifests/pertrait_pilot/RUNTIME_INPUT_MANIFEST.json").read_text())
    for ver in ("v29", "v45"):
        expected = manifest["inputs"][f"data/processed/tss_tables/tss_{ver}.tsv"]["sha256"]
        require(runtime["tss"][ver]["byte_equal"] is True, "TSS binding UNRESOLVED")
        pin(ROOT / f"results/pertrait_pilot/tss/tss_{ver}.tsv", expected)
    for name in ("pertrait_pilot_prereg", "pertrait_pilot_prereg_v1.1", "pertrait_pilot_prereg_v1.0", "frozen_params", "study3_prereg"):
        text = (ROOT / f"prereg/{name}.sha256").read_text().strip()
        h, path = text.split(maxsplit=1)
        p = ROOT / path.lstrip("*")
        if not p.exists():
            p = ROOT / "prereg" / path.lstrip("*")
        pin(p, h)
    software = measure_software(manifest)
    rosters = {}
    for cutoff in (100, 50):
        path = ROOT / f"prereg/pertrait_roster_ge{cutoff}.tsv"
        if path.exists():
            sidecar = Path(str(path) + ".sha256")
            require(sidecar.is_file(), "missing roster sidecar")
            h = sidecar.read_text().split()[0]
            pin(path, h)
            require(h in (ROOT / "RUN_LOG.md").read_text(), "roster not recorded before enrichment")
            rosters[str(cutoff)] = h
    if stage in ("enrichment", "c2", "deletion", "decision"):
        require("100" in rosters, "primary roster not frozen")
    return {"stage": stage, "status": "PASS", "software": software, "rosters": rosters}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["testing", "roster", "enrichment", "c2", "deletion", "decision"], required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(check(args.stage), sort_keys=True))
    except Exception as exc:
        print(f"PILOT PRE-AUTH REFUSED: {exc}", file=sys.stderr)
        sys.exit(1)
