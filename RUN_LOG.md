# Run Log — Annotation Drift in GWAS Locus-to-Gene Mapping

Append-only. One record per pipeline run; agents append, never rewrite history.
Authored fresh 2026-07-08 (annotation-drift study of record, protocol v1.0).

## Entry template

```
---
Date:             YYYY-MM-DD
Operator:         [name / agent]
Branch:           [git branch]
Stage / target:   [script or make target]
Inputs:           [paths + SHA-256 where hash-locked]
Outputs:          [paths]
pre_auth_check:   [PASS / FAIL — must PASS before any stage run]
Study 3 prereg:   [SHA of prereg/study3_prereg_v1.0.md, when applicable]
Protocol devs:    [NONE, or description]
Notes:            [anomalies, runtime observations]
```

---

## Runs

*(No runs recorded yet. First entry to be added on first authorized pipeline execution,
after `pre_auth_check.sh` exits 0.)*
