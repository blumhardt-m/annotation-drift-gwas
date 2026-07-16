# EXECUTION_BOUNDARY.md

Defines the boundary between **scientific changes** (forbidden to agents) and **implementation changes** (permitted). This is the anti-drift layer. When in doubt, a change is scientific and forbidden.

## Scientific changes — FORBIDDEN to any agent

An agent is prohibited from making, and must escalate rather than perform, any change to:

- Any file in `LOCKED_PROTOCOL_FILES.md`.
- Any value governed by `validate_config.py:FROZEN_FIELDS` (mapping rule, window grid, window shape, tie-break, instability denominator, gene universe, pathway DB SHA, bootstrap seed, GENCODE definitions, GWAS raw-file path, GWAS SHA-256, and — once added — GWAS snapshot release ID).
- Any Study 3 frozen parameter defined in the locally frozen prereg §3 (external OSF registration pending verification): the count-statistic definition (count of distinct genome-wide-significant loci per gene, rsID-deduplicated within gene, trait-agnostic), the shared union universe, minSize 15 / maxSize 500, FDR < 0.05, fgsea seed 42, the boundary band Δ, the chromosome-block bootstrap (seed 42, 1000 reps), the primary endpoint (C_ORA > C_GSEA), the demotion of A to descriptive, the degenerate-case rule and its 10-pathway threshold, or the kill criterion.
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

## Registration hold (Gate 0) — confirmatory execution prohibited until OSF verified

References in earlier project records to a "registered preregistration" do not establish that an immutable OSF registration exists. Confirmatory Study 3 execution is prohibited until the exact external registration, timestamp, and attached artifact hash have been verified by the investigator.

This hold stands regardless of Gate 0's eventual outcome. The Study 3 preregistration is a **local freeze** (`prereg/study3_prereg_v1.0.md`, SHA in `prereg/study3_prereg.sha256`); a local freeze is not an external registration, and the two must be kept distinct in every record. `pre_auth_check.sh` passing is necessary but **not** sufficient: it verifies local lock integrity, not external registration.

Gate 0 verification is a manual, investigator-only step. It must determine:

- whether the OSF entity referenced by the investigator (`fv4sg`) is an immutable *registration* rather than an ordinary (mutable) project;
- its effective registration timestamp;
- which prereg artifact it contains — the superseded `b8211b77…` or the controlling `29fe3ff1…`;
- whether the registered bytes exactly match the local artifact;
- whether registration preceded any confirmatory execution.

Dispositions:

- **No immutable registration** → registration is absent; confirmatory execution stays prohibited.
- **Registration contains `b8211b77…`** → a registration exists but covers a *superseded* protocol; it must NOT be represented as registration of the controlling protocol (`29fe3ff1…`), and confirmatory execution stays prohibited.
- **Registration contains `29fe3ff1…`** → potentially clean, subject to timestamp and byte verification before the hold is lifted.
- **Artifact cannot be retrieved or byte-matched** → registration identity remains unverified; confirmatory execution stays prohibited.

Only the investigator may lift this hold, and only after a clean disposition. No agent may infer registration status from Git, a README statement, a draft-registration notification, or an uploaded file alone.
