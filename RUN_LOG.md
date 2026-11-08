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
Study 3 prereg:   b8211b77cfc6d5cc50b0839a9a7496c82481933eb55f63ea26d18671b3024e47 (prereg/study3_prereg_v1.0.md, FROZEN v1.0)
Protocol devs:    NONE. Pipeline code run unmodified. Snapshot pinned + comparison.tsv SHA recorded per DEC-004 (governance amendment restoring metadata; no analytical input changed).
Notes:            Verification gate PASSED — regenerated comparison.tsv reproduces the parent
                  ORA result within 3 sig figs (ρ=0.922857≈0.923, Jaccard=0.6133≈0.613,
                  sig v29/v45/shared = 74/47/46). Exploratory viability probe (non-confirmatory)
                  found the count-based GSEA arm non-degenerate (sig v29=47, v45=63; §4.1 sub-case i).
                  R 4.5.3 / fgsea 1.36.2, fgseaMultilevel eps=0. Env recorded in
                  scratch/study3_viability/ENV.md.

---

## 2026-07-17 — Study 3 governance branch remote-publication verification

Remote-publication status was manually verified after commits
`19d16beebaabfe344b671a78059b452730628056`
and
`938d7e9bea4c3f9c5937ea08d7b1066ce30b4907`.

At verification:

- local branch `review/study3-freeze-and-gate0-20260716` had no configured upstream;
- the local governance and housekeeping commits were not present on the remote;
- remote branch inventory contained only `main` at
  `4da4768e130c1a90534a2f57735d5aa3ddfd3902`;
- remote tag inventory contained only the pre-existing
  `v1.0_submission` tag at
  `37272b2ccf2974aba56a26d68c84f26af4eb362f`;
- the prior remote hotfix branch was absent; and
- no pull request existed for
  `review/study3-freeze-and-gate0-20260716`.

Disposition: the Study 3 review branch remained local and unpublished. No
branch or tag was pushed, and no pull request was created or updated for the
review branch.

This entry records remote-publication status only. It does not lift the
standing Study 3 execution hold or authorize implementation or execution of
stages 19–21.


MODEL: Codex, based on GPT-6; version unavailable. THINKING LEVEL: high requested; actual setting unavailable.
Pilot revision 1.1 execution: exact governance installed, inputs staged and software versions matched; prescribed TSS preparation imported parse_gtf_to_tss with default biotype.
{
  "v29": {
    "regenerated_sha256": "5c831ca2b14155690fd8c4f59c3cb821ce516267765c25515877d27e7bcc1a5d",
    "retained_sha256": "5c831ca2b14155690fd8c4f59c3cb821ce516267765c25515877d27e7bcc1a5d",
    "byte_equal": true,
    "verdict": "EQUAL"
  },
  "v45": {
    "regenerated_sha256": "c109894e13f56d5971fd6191e99bc46af7dc1b96b01a5e53008f1bb5449e263e",
    "retained_sha256": "c109894e13f56d5971fd6191e99bc46af7dc1b96b01a5e53008f1bb5449e263e",
    "byte_equal": true,
    "verdict": "EQUAL"
  }
}

Both pilot and legacy authorization gates PASS; 31 synthetic checks PASS. No empirical outcome has been inspected. Commits: {"governance": "064e46728c1ad4fa20b5c7fb25e46135b404d213", "implementation": "99a08d7806575022468f01aec43885764e62ed98", "tests": "8a4c0ddbe214cc5f3078418f4e43740d4510eaba"}

MODEL: Codex, based on GPT-6; version unavailable. THINKING LEVEL: high requested; actual setting unavailable.
Roster >= 100 frozen before enrichment: N=699, SHA-256 9d6b812deed0337fb81d4e32224701425eb861ecb377957fe7ad6d845842079b; code-point order, existing indices preserved.

MODEL: Codex, based on GPT-6; version unavailable. THINKING LEVEL: high requested; actual runtime setting unavailable.
2026-11-08T04:36:16.753107+00:00 STOPPED: fgsea v29 for trait index 2 returned NA pval/padj for 2 of 1369 pathways. No scientific rule evaluated; no C2, deletion, extension or rerun. Attempted outputs retained; scientific code unchanged. Owner disposition required. The adopted terminal deposit/parking pair closes this attempt.

MODEL: Codex, based on GPT-6; version unavailable. THINKING LEVEL: high requested; actual setting unavailable.
v1.2 TSS regeneration: pinned producer imported; both release tables byte-identical to retained pinned inputs: {"v29": {"regenerated_sha256": "5c831ca2b14155690fd8c4f59c3cb821ce516267765c25515877d27e7bcc1a5d", "retained_sha256": "5c831ca2b14155690fd8c4f59c3cb821ce516267765c25515877d27e7bcc1a5d", "byte_equal": true, "regenerated_path": "/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot/results/pertrait_pilot/tss/regenerated/tss_v29.tsv", "input_path": "/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot/data/processed/tss_tables/tss_v29.tsv"}, "v45": {"regenerated_sha256": "c109894e13f56d5971fd6191e99bc46af7dc1b96b01a5e53008f1bb5449e263e", "retained_sha256": "c109894e13f56d5971fd6191e99bc46af7dc1b96b01a5e53008f1bb5449e263e", "byte_equal": true, "regenerated_path": "/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot/results/pertrait_pilot/tss/regenerated/tss_v45.tsv", "input_path": "/home/michaelblumhardt/Study/research/08_GENOMICS/annotation-drift-gwas-pertrait-pilot/data/processed/tss_tables/tss_v45.tsv"}}

MODEL: Codex, based on GPT-6; version unavailable. THINKING LEVEL: high requested; actual setting unavailable.
v1.2 roster >= 100: regenerated bytes equal roster of record; N=699, SHA-256 9d6b812deed0337fb81d4e32224701425eb861ecb377957fe7ad6d845842079b; verified before enrichment.
