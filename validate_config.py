#!/usr/bin/env python3
"""
validate_config.py — Frozen-field guard for annotation_drift_gwas_l2g.

Compares a live config.yaml against a locked snapshot and fails hard if any
protocol-critical field has changed.  Intended to run as the first step of
`make parse`, `make map`, and `make analyze`.

Usage:
    python validate_config.py config.yaml prereg/frozen_params.yaml
    python validate_config.py config.yaml                          # uses prereg/frozen_params.yaml default

Setup (run once after freezing):
    cp config.yaml prereg/frozen_params.yaml
    sha256sum prereg/frozen_params.yaml > prereg/frozen_params.sha256

Exit codes:
    0  All frozen fields match.
    1  At least one frozen field differs (details printed to stderr).
    2  IO / parse error.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Any

import yaml

# ── Frozen fields ──────────────────────────────────────────────────────────
# Each entry: (dotted key path into the YAML, human label for error messages).
# Add new fields here if the protocol is amended.
FROZEN_FIELDS: list[tuple[str, str]] = [
    # Core mapping
    ("map.mapping_rule.type",                   "Mapping rule"),
    ("map.windows_kb",                          "Window grid"),
    ("map.window_shape",                        "Window shape"),
    ("map.tie_break.rule",                      "Tie-break rule"),
    # Instability definition
    ("analysis.primary_instability.denominator", "Instability denominator"),
    # Enrichment
    ("enrichment.gene_universe",                "Gene-set universe"),
    ("enrichment.pathway_db.sha256",            "Pathway DB SHA-256"),
    # Reproducibility
    ("analysis.uncertainty.bootstrap.seed",     "Bootstrap RNG seed"),
    # Input data
    ("annotations.gencode",                     "GENCODE annotation definitions"),
    ("inputs.gwas_loci.raw_file",               "GWAS input file"),
    ("inputs.gwas_loci.sha256",                 "GWAS input SHA-256"),
]


def resolve(cfg: dict, dotted_key: str) -> Any:
    """Walk a nested dict by dotted path.  Raises KeyError if missing."""
    node = cfg
    for part in dotted_key.split("."):
        node = node[part]
    return node


def load_yaml(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def validate(live_path: Path, locked_path: Path) -> list[str]:
    """Return a list of violation messages (empty = pass)."""
    live = load_yaml(live_path)
    locked = load_yaml(locked_path)
    violations: list[str] = []

    for dotted_key, label in FROZEN_FIELDS:
        try:
            locked_val = resolve(locked, dotted_key)
        except KeyError:
            violations.append(
                f"  {label} ({dotted_key}): key missing from locked snapshot"
            )
            continue

        try:
            live_val = resolve(live, dotted_key)
        except KeyError:
            violations.append(
                f"  {label} ({dotted_key}): key missing from live config"
            )
            continue

        if live_val != locked_val:
            violations.append(
                f"  {label} ({dotted_key}):\n"
                f"    locked : {locked_val!r}\n"
                f"    live   : {live_val!r}"
            )

    return violations


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python validate_config.py <config.yaml> [prereg/frozen_params.yaml]",
              file=sys.stderr)
        sys.exit(2)

    live_path = Path(sys.argv[1])
    locked_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("prereg/frozen_params.yaml")

    if not live_path.exists():
        print(f"ERROR: live config not found: {live_path}", file=sys.stderr)
        sys.exit(2)
    if not locked_path.exists():
        print(f"ERROR: locked snapshot not found: {locked_path}", file=sys.stderr)
        print("       Run `cp config.yaml prereg/frozen_params.yaml` after freezing.",
              file=sys.stderr)
        sys.exit(2)

    # Integrity check: warn if the locked file itself was modified
    lock_hash_file = locked_path.with_suffix(".sha256")
    if lock_hash_file.exists():
        expected = lock_hash_file.read_text().strip().split()[0]
        actual = file_sha256(locked_path)
        if actual != expected:
            print(
                f"WARNING: locked snapshot has been modified since freezing.\n"
                f"  expected SHA-256: {expected}\n"
                f"  actual SHA-256:   {actual}\n",
                file=sys.stderr,
            )

    violations = validate(live_path, locked_path)

    if violations:
        print(
            "PROTOCOL VIOLATION: frozen config fields have changed.\n"
            "The following fields differ between the live config and the locked snapshot:\n",
            file=sys.stderr,
        )
        for v in violations:
            print(v, file=sys.stderr)
        print(
            "\nTo proceed, either:\n"
            "  (a) revert config.yaml to match the locked snapshot, or\n"
            "  (b) document a protocol amendment and re-freeze:\n"
            "      cp config.yaml prereg/frozen_params.yaml\n"
            "      sha256sum prereg/frozen_params.yaml > prereg/frozen_params.sha256\n",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print(f"validate_config: all {len(FROZEN_FIELDS)} frozen fields OK.")
        sys.exit(0)


if __name__ == "__main__":
    main()
