# EXECUTION_BOUNDARY.md

Defines the boundary between **scientific changes** (forbidden to agents) and **implementation changes** (permitted). This is the anti-drift layer. When in doubt, a change is scientific and forbidden.

## Scientific changes — FORBIDDEN to any agent

An agent is prohibited from making, and must escalate rather than perform, any change to:

- Any file in `LOCKED_PROTOCOL_FILES.md`.
- Any value governed by `validate_config.py:FROZEN_FIELDS` (mapping rule, window grid, window shape, tie-break, instability denominator, gene universe, pathway DB SHA, bootstrap seed, GENCODE definitions, GWAS raw-file path, GWAS SHA-256, and — once added — GWAS snapshot release ID).
- Any Study 3 frozen parameter defined in the registered prereg §3 (OSF-registered; controlling artifact `29fe3ff1…` verified and approved 2026-07-09 — DEC-008/DEC-011): the count-statistic definition (count of distinct genome-wide-significant loci per gene, rsID-deduplicated within gene, trait-agnostic), the shared union universe, minSize 15 / maxSize 500, FDR < 0.05, fgsea seed 42, the boundary band Δ, the chromosome-block bootstrap (seed 42, 1000 reps), the primary endpoint (C_ORA > C_GSEA), the demotion of A to descriptive, the degenerate-case rule and its 10-pathway threshold, or the kill criterion.
- Data inputs: no re-download, no substitution, no regeneration of any hash-locked input with a different byte content. Restoring a file whose hash matches a locked SHA is permitted (it is provenance restoration, not input change); producing a file whose hash differs is forbidden.

## Implementation changes — PERMITTED

An agent may, within the scope of the active `TASK.md` and on the authorized branch:

- Create the new Study 3 stage scripts (`src/19_`, `20_`, `21_`) and their outputs under `results/study3_ranked_enrichment/`.
- Fix paths, add logging, add descriptive diagnostic plots/tables that do not alter any frozen quantity.
- Import and call frozen machinery (mapping, alias table, flip definition, bootstrap) — importing is required; reimplementing a frozen computation is forbidden because it risks silent divergence.
- Refactor its own new code.
- Append to `RUN_LOG.md`.

## The re-freeze fence (critical)

`validate_config.py` prints remediation guidance on failure that includes a two-line re-freeze: `cp config.yaml prereg/frozen_params.yaml && sha256sum ... > prereg/frozen_params.sha256`. **Agents are prohibited from executing the re-freeze path.** Re-freezing converts a scientific divergence into a new lock and is the primary silent-drift vector. The only permitted agent response to a `validate_config.py` non-zero exit, or to any frozen-field divergence, is to STOP and escalate per `REVIEW_TRIGGER.md`. Filling the snapshot placeholder, adding a `FROZEN_FIELDS` entry, or any re-freeze is a human act requiring a `DECISIONS.md` entry.

## Provenance restoration vs. input change

Restoring the GWAS Catalog file from backup/Trash into `data/raw/` is permitted *only* after its SHA-256 is verified equal to `bfaffa3a…`. If the hash does not match, restoration is forbidden and the condition escalates. This distinction — same bytes = provenance restoration, different bytes = input change — is the rule that let this study recover its snapshot without altering what the parent analysis ran on.

## Authorized execution entry point

The only authorized entry point for running any pipeline stage is `pre_auth_check.sh` (which itself calls `validate_config.py`). Invoking `validate_config.py` directly is permitted for inspection but is never a substitute for the full gate: the gate performs checks the guard structurally cannot (repo identity, baseline ancestry, placeholder scan, data-presence-by-hash, inventory cross-check). Guard and gate are harmonized in severity (a modified lock file is fatal by either path), but only the gate is complete. An agent is prohibited from running a stage without a passing `pre_auth_check.sh` in the same session.

## Registration verification (Gate 0 — cleared 2026-07-17)

Gate 0 has been verified and the temporary registration-verification hold is released (DEC-011). The controlling Study 3 preregistration (`prereg/study3_prereg_v1.0.md`, SHA-256 `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`, recorded in `prereg/study3_prereg.sha256`) was downloaded from and byte-matched to the immutable embargoed OSF registration snapshot. Registration was initiated 2026-07-08 and approved 2026-07-09 — before any implementation or execution of Study 3 stages 19–21, and before computation of `C_ORA`, `C_GSEA`, or their directional comparison (DEC-008).

The local-freeze-versus-external-registration distinction remains a standing principle: a local freeze is never, by itself, an external registration, and `pre_auth_check.sh` passing verifies local lock integrity only, not registration. What changed on 2026-07-17 is solely that the external registration has now been verified for this controlling artifact; the distinction is retained for all future records.

Registered Study 3 branch: `study3-ranked-enrichment`, to be created from `review/study3-freeze-and-gate0-20260716` (not from `main`), because the reviewed protocol, its checksum lineage, the Gate 0 record, and decisions DEC-004 through DEC-011 exist only on the review branch.

Confirmatory execution proceeds under the registered sequencing in prereg §8, subject to the amended §3.4 ORA consistency preflight (DEC-009): the preflight is run against the frozen parent inputs and halts execution mechanically if any DEC-009 condition fails (v29/v45/shared/union counts 74/47/46/75, Jaccard 46/75, Spearman ρ serialized to `0.922857`, and matching input hashes, pathway filtering, alignment, and deterministic implementation). The additive-only execution boundary and every prohibition in this document remain fully in force: no modification of frozen inputs, frozen parameters, `comparison.tsv`, prior results, or the manuscript during Study 3 execution.
