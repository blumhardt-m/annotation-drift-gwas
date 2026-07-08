#!/usr/bin/env bash
# pre_auth_check.sh
# Automated gate that must pass before any agent runs a Study 3 pipeline stage.
# Wraps validate_config.py (does not replace it) and adds the checks the guard
# structurally cannot perform: repo identity, placeholder scan, data-presence-by-hash,
# and promotion of the frozen-params SHA warning to a hard failure.
#
# Exit 0 = authorized to run. Non-zero = STOP (see REVIEW_TRIGGER.md).
# This script is protocol-locked (LOCKED_PROTOCOL_FILES.md). Agents must not weaken it.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

fail() { echo "PRE-AUTH FAIL [$1]: $2" >&2; exit 1; }
ok()   { echo "  ok: $1"; }

echo "== pre_auth_check =="

# ---- 1. Repo identity (this session began in the wrong repo; never again) ----
EXPECTED_REMOTE_SUBSTR="annotation-drift-gwas"
REMOTE_URL="$(git remote get-url origin 2>/dev/null || echo none)"
case "$REMOTE_URL" in
  *"$EXPECTED_REMOTE_SUBSTR"*) ok "repo identity ($REMOTE_URL)";;
  *) fail "repo-identity" "origin is '$REMOTE_URL', expected to contain '$EXPECTED_REMOTE_SUBSTR' (TRIGGER 4)";;
esac

# ---- 2. Branch discipline ----
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [ "$BRANCH" = "main" ] || [ "$BRANCH" = "master" ]; then
  fail "branch" "on '$BRANCH'; Study 3 work requires an authorized feature branch (see TASK.md)"
fi
ok "branch ($BRANCH)"

# ---- 2b. Scientific-baseline ancestry (repo identity is not enough) ----
# HEAD must descend from the frozen submission tag, or the scientific baseline is wrong
# even when the remote and branch look fine (a branch off an unrelated commit passes
# every other check but is not built on v1.0_submission).
BASELINE_TAG="v1.0_submission"
if git rev-parse -q --verify "refs/tags/${BASELINE_TAG}" >/dev/null; then
  if ! git merge-base --is-ancestor "refs/tags/${BASELINE_TAG}" HEAD; then
    fail "baseline-ancestry" "HEAD does not descend from ${BASELINE_TAG}; scientific baseline mismatch (TRIGGER 4)"
  fi
  ok "descends from ${BASELINE_TAG}"
else
  fail "baseline-ancestry" "frozen tag ${BASELINE_TAG} not found in this repo (TRIGGER 4)"
fi

# ---- 3. Governance + locked files present ----
for f in AGENTS.md EXECUTION_BOUNDARY.md LOCKED_PROTOCOL_FILES.md REVIEW_TRIGGER.md \
         TASK.md validate_config.py config.yaml prereg/frozen_params.yaml prereg/frozen_params.sha256; do
  [ -f "$f" ] || fail "missing-file" "required governance/protocol file absent: $f"
done
ok "governance + protocol files present"

# ---- 4. Frozen-params SHA: PROMOTE validate_config's soft warning to hard failure ----
RECORDED="$(cat prereg/frozen_params.sha256 | awk '{print $1}')"
ACTUAL="$(sha256sum prereg/frozen_params.yaml | awk '{print $1}')"
[ "$RECORDED" = "$ACTUAL" ] || fail "lock-integrity" \
  "frozen_params.yaml bytes ($ACTUAL) != recorded lock ($RECORDED). Locked file was modified post-freeze (TRIGGER 1)"
ok "frozen_params lock integrity"

# ---- 5. Placeholder scan in frozen config (the exact hole this session found) ----
if grep -q "REPLACE_WITH_RELEASE_ID" config.yaml prereg/frozen_params.yaml 2>/dev/null; then
  fail "placeholder" "unresolved placeholder REPLACE_WITH_RELEASE_ID in frozen config (TRIGGER 2)"
fi
ok "no unresolved placeholders"

# ---- 6. Hash-locked data inputs: present AND correct bytes if present ----
# GWAS Catalog (required restored for provenance completeness)
CATALOG="data/raw/gwas_catalog.tsv"
CATALOG_SHA="bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b"
if [ -f "$CATALOG" ]; then
  A="$(sha256sum "$CATALOG" | awk '{print $1}')"
  [ "$A" = "$CATALOG_SHA" ] || fail "input-hash" "$CATALOG hash $A != locked $CATALOG_SHA (TRIGGER 5)"
  ok "GWAS Catalog hash-verified"
else
  fail "missing-input" "$CATALOG absent; restore hash-verified copy before running (TRIGGER 3)"
fi
# Frozen mapping table (sole required Study 3 input)
[ -f "data/processed/mappings/comparison.tsv" ] || fail "missing-input" "comparison.tsv absent (TRIGGER 3)"
ok "comparison.tsv present"

# ---- 7. Governance-inventory cross-check (TRIGGER 11) ----
# Every frozen data path enforced by validate_config's FROZEN_FIELDS-linked config
# must also appear in LOCKED_PROTOCOL_FILES.md. Minimal structural check: the GWAS
# raw file and the snapshot release ID must be referenced in the locked inventory.
for token in "comparison.tsv" "gwas_catalog" "frozen_params.yaml" "v1.0.2-associations_e115_r2026-01-19_full"; do
  grep -q "$token" LOCKED_PROTOCOL_FILES.md 2>/dev/null \
    || fail "gov-inventory" "'$token' enforced/expected but not enumerated in LOCKED_PROTOCOL_FILES.md (TRIGGER 11)"
done
# The snapshot release ID must now be a frozen field in the guard itself.
grep -q "inputs.gwas_loci.snapshot" validate_config.py 2>/dev/null \
  || fail "gov-inventory" "inputs.gwas_loci.snapshot not in validate_config.py FROZEN_FIELDS (TRIGGER 11)"
ok "governance inventories consistent"

# ---- 8. Delegate frozen-field enforcement to the existing guard ----
echo "-- validate_config.py --"
if ! python3 validate_config.py; then
  fail "validate_config" "frozen-field divergence (TRIGGER 6). DO NOT re-freeze; escalate."
fi
ok "validate_config passed"

echo "== pre_auth_check PASSED — authorized to run =="
exit 0
