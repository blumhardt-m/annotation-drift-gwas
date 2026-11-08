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

def check(stage):
    import yaml
    require(ROOT == EXPECTED_ROOT, "wrong worktree")
    for ancestor in set(ROOT.parents) | set(DEPOSIT.parents):
        require(not (ancestor / "AGENTS.md").exists(), f"external ancestor AGENTS.md: {ancestor}")
    check_v12_authority(ROOT, DEPOSIT)
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
    pin(DEPOSIT / "portfolio/REACTIVATION_BRACKET.json", BRACKET_SHA)
    bracket = json.loads((DEPOSIT / "portfolio/REACTIVATION_BRACKET.json").read_text())
    require(bracket["owner_sha256"] == ADOPTION_V3, "reactivation authority differs")
    for item in bracket["after"].values():
        pin(item["path"], item["sha256"])
    installation = {"after": {
        "register": bracket["after"]["00_ADMIN/PORTFOLIO_STATUS.yaml"],
        "validator": bracket["after"]["00_ADMIN/tools/validate_portfolio_status.py"]}}
    register = yaml.safe_load(Path(installation["after"]["register"]["path"]).read_text())
    rows = register["projects"]
    require(len(rows) == 16, "membership count")
    check_allocation(next((r for r in rows if r["id"] == "gwas_annotation_drift_pertrait_pilot"), {}))
    validator = command([sys.executable, installation["after"]["validator"]["path"],
                         installation["after"]["register"]["path"]])
    require(validator.startswith("VALID: 16 entries, 3 live"), "live mirror mismatch")
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
