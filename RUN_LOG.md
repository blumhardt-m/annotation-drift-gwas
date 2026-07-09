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

---
Date:             2026-07-08
Operator:         Claude (agent) under M. Blumhardt authorization
Branch:           study3-prereg-freeze
Stage / target:   src/02_parse_gencode.py, src/03_parse_gwas.py, src/04_map_variants.py (frozen mapping pipeline, run as-is); frozen ORA reference via src/primary_reactome_export.py
Inputs:           data/raw/gwas_catalog.tsv (SHA-256 bfaffa3a…c3480b, verified) ;
                  data/raw/gencode.v29.annotation.gtf.gz (e77a9568…c2b9, verified) ;
                  data/raw/gencode.v45.annotation.gtf.gz (56a166e4…ae68, verified) ;
                  data/raw/reactome.gmt (a34f43fd…)
Outputs:          data/processed/tss_tables/tss_v29.tsv, tss_v45.tsv ;
                  data/processed/gwas_eligible.tsv (746,850 eligible) ;
                  data/processed/mappings/comparison.tsv (SHA-256 8f27c1db…6ba2bb1;
                  721,406 rows, 291,292 unique rsIDs, 29,920 differ = 4.15%) ;
                  scratch/study3_viability/** (EXPLORATORY viability probe, quarantined)
pre_auth_check:   n/a for the freeze itself; regeneration ran as authorized §7 prerequisite work (not a confirmatory Study 3 run). The Study 3 build (stages 19–21) remains gated on a feature branch per TASK.md.
Study 3 prereg:   29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17 (prereg/study3_prereg_v1.0.md, FROZEN v1.0; SHA corrected per DEC-005 — supersedes the earlier freeze)
Protocol devs:    NONE. Pipeline code run unmodified. Snapshot pinned + comparison.tsv SHA recorded per DEC-004 (governance amendment restoring metadata; no analytical input changed).
Notes:            Verification gate PASSED — regenerated comparison.tsv reproduces the parent
                  ORA result within 3 sig figs (ρ=0.922857≈0.923, Jaccard=0.6133≈0.613,
                  sig v29/v45/shared = 74/47/46). Exploratory viability probe (non-confirmatory)
                  found the count-based GSEA arm non-degenerate (sig v29=47, v45=63; §4.1 sub-case i).
                  R 4.5.3 / fgsea 1.36.2, fgseaMultilevel eps=0. Env recorded in
                  scratch/study3_viability/ENV.md.
