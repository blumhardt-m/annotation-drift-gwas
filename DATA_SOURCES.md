# Data Sources — Annotation Drift in GWAS Locus-to-Gene Mapping

Protocol version: 1.0 (annotation-drift study of record).
Authored fresh 2026-07-08 from the study's own frozen truth (`prereg/frozen_params.yaml`
mirrored from `config.yaml`) — **not** carried from any prior document. Every input below
is pinned by the SHA-256 recorded in `frozen_params.yaml`; that file is authoritative.

## Inputs of record (content-addressed)

| Source | Path | SHA-256 | Notes |
|--------|------|---------|-------|
| GWAS Catalog associations | `data/raw/gwas_catalog.tsv` | `bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b` | NHGRI-EBI GWAS Catalog. European ancestry, p < 5e-8, autosomal, sentinel-SNP loci. Release ID recovered as `v1.0.2-associations_e115_r2026-01-19_full` (see DECISIONS). |
| GENCODE v29 GTF | `data/raw/gencode.v29.annotation.gtf.gz` | `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9` | GRCh38.p12. Older annotation version. |
| GENCODE v45 GTF | `data/raw/gencode.v45.annotation.gtf.gz` | `56a166e4a43c934a1910c0d09dec786ac048862fc65130731d41f3c7f451ae68` | GRCh38.p14. Newer annotation version. |
| Reactome GMT v95 | `data/raw/reactome.gmt` | `a34f43fd763a6800a9a94310c416ae511e84724d87a3c255f35d13764d12bebb` | Primary pathway arbiter for ORA (min 15 / max 500, BH FDR < 0.05). |

## Derived artifact of record
- `data/processed/mappings/comparison.tsv` — frozen locus-to-gene comparison table
  (v29 vs v45, nearest-TSS), the sole required input to the Study 3 count statistic.
  **SHA-256 `8f27c1dba92c2af3ec4c427b1bf89cdb03d67f68adef09c1c9131111f6ba2bb1`** (input of
  record; also pinned at `frozen_params.yaml:inputs.mapping_artifact` and in
  `manifests/comparison_tsv_regenerated_20260708.sha256`).
  Provenance: regenerated 2026-07-08 by the frozen mapping pipeline (`src/02→03→04`,
  `--window-kb 500`) from the pinned raw inputs above; reproduces the parent ORA result
  within 3 sig figs (Spearman ρ 0.923, Jaccard 0.613; sig sets 74/47/46). Structure:
  721,406 rows, 291,292 unique rsIDs, 29,920 differ (4.15%). Gitignored (derived);
  content-addressed by the SHA above. Locked in `LOCKED_PROTOCOL_FILES.md` §1a.

## Provenance notes
- The annotation pair of record is **GENCODE v29 vs v45** (both on GRCh38). Any document
  citing a v37/v44 pair reflects a superseded pre-freeze framing and is not this study.
- The GWAS Catalog release ID (`v1.0.2-associations_e115_r2026-01-19_full`) was recovered
  and verified byte-identical to the pinned SHA-256, and is now **pinned into
  `inputs.gwas_loci.snapshot`** in both `config.yaml` and `frozen_params.yaml` (2026-07-08,
  DEC-004). The field is in `validate_config.py:FROZEN_FIELDS`, so enforcement is now
  active; `frozen_params.sha256` was regenerated to match.
- Raw input bytes are re-downloadable; the SHA-256 pins above are the archival record.
