# Decision Log — Annotation Drift in GWAS Locus-to-Gene Mapping

Append-only. Do not edit or delete existing entries; new entries at the bottom.
This log belongs to the **annotation-drift study of record** (protocol v1.0). It was
authored fresh on 2026-07-08; it is not the decision log of the L2G Instability Audit
(that study's DEC-001…006 live in its own archive).

---

## DEC-001 — Agent-governance layer adopted (2026-07-08)

Decision: Adopted the agent-governance package (`AGENTS.md`, `EXECUTION_BOUNDARY.md`,
`LOCKED_PROTOCOL_FILES.md`, `ONBOARDING_SOP.md`, `REVIEW_TRIGGER.md`,
`REVIEW_CHECKLIST.md`, `pre_auth_check.sh`, `TASK.md`) plus the `validate_config.py`
FROZEN_FIELDS patch, on branch `governance-adoption`, merged to `main` as baseline
infrastructure (not Study 3 work).
Rationale: The repository had no runtime governance; this installs the anti-drift gate,
freeze inventory, and review triggers before Study 3 begins.
Provenance / scope: see `GOVERNANCE_ADOPTION.md` (one-time migration record; archivable
once `pre_auth_check.sh` passes on `main`).

---

## DEC-002 — GWAS Catalog snapshot recovered and release ID documented (2026-07-08)

Decision: The frozen GWAS Catalog input (`data/raw/gwas_catalog.tsv`, SHA-256
`bfaffa3a…`) was recovered from Trash on 2026-07-07 and verified byte-identical to the
pre-registered pin. Its containing directory identified the release as
`v1.0.2-associations_e115_r2026-01-19_full` (GWAS Catalog v1.0.2, Ensembl e115, release
2026-01-19).
Rationale: The original lock pinned the SHA-256 but left the human-readable release ID as
the placeholder `REPLACE_WITH_RELEASE_ID`; the recovered directory name resolves it.
Nature: provenance restoration of metadata — the analytic input is unchanged (bytes match
the pin). **Pending (human re-freeze, NOT done here):** fill `snapshot` in `config.yaml`
and `prereg/frozen_params.yaml`, add `inputs.gwas_loci.snapshot` to FROZEN_FIELDS, and
regenerate `frozen_params.sha256`. See `validate_config_FROZEN_FIELDS_patch.md`.

---

## DEC-003 — v1.3 L2G-audit contamination disposed (2026-07-08)

Decision: Removed the two-generations-obsolete v1.3 L2G-Instability-Audit governance/
protocol files that had been injected into this repo in commit `90f6631`
(`PROTOCOL_v1.3.md`, `LOCK_SUMMARY.md`, `TASKS_WEEK1.md`, `EXECUTION_MEMO_v1.0.md`,
`OUTPUT_SCHEMA.md`, `DATA_SOURCES.md`, `MANUSCRIPT_POSITIONING.md`, `DECISIONS.md`,
`RUN_LOG.md`, `configs/`, `results/locus_instability_atlas/`).
Rationale: Disposed, not relocated — the sibling L2G audit is a superseded design (v1.6,
v29→v45 RA-only) with no valid home for v1.3; the files carry no unique value. The two
that appeared to hold annotation-drift content (`DATA_SOURCES.md`,
`MANUSCRIPT_POSITIONING.md`) on inspection held a stale, unpinned, partially-wrong (v37,
~8.3%) view of this study that `config.yaml`/`frozen_params.yaml` already supersede
(v29/v45, 4.15%); extracting would have seeded a Trigger-11 inconsistency. Nothing
extracted; `DATA_SOURCES.md`, `DECISIONS.md`, `RUN_LOG.md` re-authored fresh from
authoritative sources.
Safety: audit content preserved in the L2G audit archive (`archived-2026-07-08` bundle,
3 copies), in this repo's git history, and in
`~/Downloads/anndrift_v1.3_contamination_removed/`. See `LOCKED_PROTOCOL_FILES.md`
"Foreign artifacts (disposition)" and `GOVERNANCE_ADOPTION.md`.

---

## DEC-004 — Study 3 prereg frozen; snapshot pinned; comparison.tsv is input of record (2026-07-08)

Decision: Froze the Study 3 ranked-enrichment prereg (`prereg/study3_prereg_v1.0.md`,
Status → FROZEN v1.0, SHA-256 `b8211b77…` in `prereg/study3_prereg.sha256`) and completed
its §7 governance prerequisites #1–#3:
(a) Pinned `inputs.gwas_loci.snapshot` = `v1.0.2-associations_e115_r2026-01-19_full` in
    **both** `config.yaml` and `prereg/frozen_params.yaml`. The field was already in
    `validate_config.py:FROZEN_FIELDS`; it was inert only because the value was the
    `REPLACE_WITH_RELEASE_ID` placeholder. Pinning activates enforcement.
    `prereg/frozen_params.sha256` regenerated (now `8407f51f…`).
(b) Recorded `comparison.tsv` SHA-256 `8f27c1dba92c…6ba2bb1` as the mapping input of record
    in `frozen_params.yaml:inputs.mapping_artifact`, `config.yaml` (mirror preserved),
    `DATA_SOURCES.md`, and `manifests/comparison_tsv_regenerated_20260708.sha256`.
Rationale: The snapshot recovery (DEC-002) restored human-readable metadata only — the
recovered Catalog bytes match the pre-registered pin `bfaffa3a…` exactly, and the
regenerated `comparison.tsv` reproduces the published parent ORA result to 3 sig figs
(ρ 0.923, Jaccard 0.613; sig sets 74/47/46), confirmed by the exploratory viability probe
(quarantined in `scratch/study3_viability/`, non-confirmatory). No analytical input changed.
Provenance/scope: The viability probe also showed the count-based GSEA arm is non-degenerate
(sig v29 = 47, v45 = 63 ≥ the §4.1 threshold of 10) — this established *viability*, not
support for H1, which is tested only by the confirmatory boundary-flip contrast (§4.1).
Toolchain frozen at R 4.5.3 / fgsea 1.36.2, `fgseaMultilevel` eps=0 (§3.3).
Pending (NOT done here): OSF registration of the frozen prereg (external human step);
§7 #4 (MEMORY.md regeneration) and #5 (parent-study provenance note) — neither blocks the
freeze. Δ = 0.5 ratified (§4.1), already in the frozen text.

---

## DEC-005 — Study 3 prereg §4.2 buffering figure corrected pre-registration (2026-07-08)

Study 3 prereg §4.2 buffering figure corrected pre-registration. An earlier draft cited a
non-deduplicated row-count mean of 41.8 loci/gene (721,406 rows ÷ ~17,300 head genes); the
frozen statistic is rsID-deduplicated (§3.2), with mean ≈ 16.9 / median ≈ 9 across the
nonzero head. The 41.8 figure did not describe the implemented statistic and was removed
from the frozen text; recorded here for development provenance. Supersedes the b8211b77
freeze, which was never registered or pushed.

New frozen prereg SHA-256: `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`
(`prereg/study3_prereg.sha256`). The RUN_LOG entry was updated to this SHA. DEC-004 is left
as originally written (append-only history recording the b8211b77 freeze as it stood at the
time); this DEC-005 entry carries the supersession and the new SHA.
Scope: text/provenance correction only — no frozen scientific parameter changed;
`validate_config` still passes and the config↔frozen_params mirror is intact.
Cross-document check: grepped all committed docs for "41.8"/"41.8 loci"; the only occurrence
was the §4.2 parenthetical now removed. DATA_SOURCES.md and all other tracked docs contained
no "41.8" reference (none to remove). The literal string "41.8" appears nowhere in the frozen
prereg or any working-tree doc except this provenance record; "b8211b77" survives in DEC-004
(the freeze it recorded), in this supersession record, and in the immutable git history of
commit 23f6b7e.

---

## DEC-006 — Execution-contract wording corrected: local freeze ≠ external registration; Gate 0 hold added (2026-07-16)

Decision (human-authorized locked-file modification): corrected the execution-contract
files to stop describing the Study 3 prereg as a "registered preregistration" and to make
the registration hold explicit. `AGENTS.md`, `EXECUTION_BOUNDARY.md`, and `TASK.md` now
read **"locally frozen preregistration; external OSF registration pending verification."**
`AGENTS.md` and `EXECUTION_BOUNDARY.md` are governance-runtime files listed in
`LOCKED_PROTOCOL_FILES.md`; editing them is a human act, recorded here per the governance
rule. No hash-locked file changed (`pre_auth_check.sh` presence-checks these docs but does
not content-hash them; `prereg/frozen_params.yaml` and its `.sha256` are untouched, so lock
integrity is intact).

Added to `EXECUTION_BOUNDARY.md` a new section **"Registration hold (Gate 0)"**:
> References in earlier project records to a "registered preregistration" do not establish
> that an immutable OSF registration exists. Confirmatory Study 3 execution is prohibited
> until the exact external registration, timestamp, and attached artifact hash have been
> verified by the investigator.

This hold stands regardless of Gate 0's eventual outcome: a local freeze is not an external
registration, and `pre_auth_check.sh` passing (local lock integrity) does not clear it.

Gate 0 verification (manual, investigator-only) must determine, for the OSF entity `fv4sg`:
whether it is an immutable *registration* (not a mutable project); its effective timestamp;
which prereg artifact it contains — superseded `b8211b77…` or controlling `29fe3ff1…`;
whether the registered bytes exactly match the local artifact; and whether registration
preceded any confirmatory execution. Dispositions: absent → prohibited; contains
`b8211b77…` → registration of a *superseded* protocol, must not be represented as
registration of the controlling protocol, prohibited; contains `29fe3ff1…` → potentially
clean, subject to timestamp + byte verification; artifact unretrievable/unmatched →
unverified, prohibited. Only the investigator may lift the hold, and only on a clean
disposition.

Scope: documentation/execution-contract only; no frozen scientific parameter changed.

---

## DEC-007 — July 8 chronology: provenance framing of record (2026-07-16)

Records the defensible provenance framing for the relationship between the exploratory
GSEA work and the confirmatory design, each statement tied to committed evidence:

- The threshold-amplification mechanism underlying H1 existed in commit `37272b2`
  (2026-03-12) — months before the probe. (Parent manuscript / README interpretation.)
- The 2026-07-08 exploratory viability probe observed significant-pathway counts
  v29 = 47, v45 = 63. (DEC-004; RUN_LOG 2026-07-08; commit `23f6b7e` message; prereg §1.
  The probe artifacts themselves are git-ignored under `scratch/` and were never committed;
  they are preserved in a private archive with a hash manifest.)
- The primary confirmatory statistics `C_ORA`, `C_GSEA`, and their directional comparison
  `C_ORA > C_GSEA` were **not** computed. (prereg §1/§4.1; no `src/21_*` stage exists; no
  computed value in `results/` or `manifests/`.)

Git cannot establish the intra-day order of the final prereg wording and the probe run on
2026-07-08, but the scientific mechanism motivating the directional prediction clearly
predates the probe. **Therefore:** the final protocol was frozen after or contemporaneously
with exploratory GSEA work, but the directional hypothesis was grounded in a version-
controlled mechanism documented months earlier and was **not fitted to the uncomputed
primary endpoint**. The final artifact is not represented as wholly outcome-blind or fully
pre-probe.

OSF registration status remains UNVERIFIED and is out of scope for Git-based verification
(see DEC-006 / EXECUTION_BOUNDARY.md "Registration hold (Gate 0)"). Scope: provenance record
only; no scientific parameter changed.

---

## DEC-008 — Owner review and prospective ratification of Study 3 protocol lineage (2026-07-17)

I reviewed commits `23f6b7e`, `11bb5ea`, `4342e2e`, `25c1834`, and `eae1e21`, including their
complete diffs and both frozen Study 3 preregistration artifacts.

The original freeze in `23f6b7e` and the descriptive re-freeze in `11bb5ea` were performed
during an agent-assisted session even though protocol freeze was reserved to the human owner.
That operating procedure did not comply with the project's reserved-action rule. Git authorship
under my identity does not cure that procedural defect, and this decision does not represent the
original freeze procedure as governance-compliant.

I independently reviewed and prospectively ratify the scientific content of the controlling
Study 3 protocol, SHA-256 `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`,
as the protocol governing future implementation and execution. The change from
`b8211b77cfc6d5cc50b0839a9a7496c82481933eb55f63ea26d18671b3024e47` to
`29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17` removed an erroneous
descriptive note in §4.2 and did not change any frozen scientific endpoint, parameter,
estimator, threshold, uncertainty procedure, or disposition rule.

The exact controlling artifact was downloaded from and byte-matched to the immutable embargoed
OSF registration snapshot. The registration was initiated on July 8, 2026 and approved on
July 9, 2026 before implementation or execution of Study 3 stages 19–21. No C_ORA, C_GSEA, or
primary directional comparison was computed before registration.

Commits `23f6b7e`, `11bb5ea`, and `4342e2e` are retained as the protocol and checksum lineage.
Commits `25c1834` and `eae1e21` are retained as audit history. Any language in those commits
stating that OSF registration remained pending is superseded by the verified Gate 0 disposition
recorded here and must not be treated as current execution authority.

This decision authorizes prospective retention and use of the reviewed protocol. It does not
retroactively cure the reserved-action violation.

This decision records successful verification of the registration artifact and sequencing
requirement, but it does not by itself lift the standing Study 3 execution hold or authorize
implementation of stages 19–21.

---

## DEC-009 — Protocol amendment and owner interpretation of the §3.4 ORA consistency check (2026-07-17)

No full-precision parent Spearman value is preserved. The registered §3.4 consistency check
states that regenerated ρ must match the parent expectation "ρ ≈ 0.923" within floating-point
tolerance. The highest-fidelity preserved parent value is ρ = 0.922857, serialized to six decimal
places. The rounded expectation 0.923 differs from 0.922857 by 1.43 × 10⁻⁴ and therefore cannot
serve as a floating-point reference value. As written, the registered check is not operationally
executable.

A search across the review branch, commit `37272b2`, tag `v1.0_submission`, and the committed
`results/`, `manifests/`, `RUN_LOG.md`, and `DECISIONS.md` records found no preserved artifact
containing either a higher-precision Spearman reference or the primary significant-pathway
identifier sets. The identified per-pathway file, `enrichment_jaccard_by_window.tsv`, belongs to
a different symbol-based diagnostic analysis and is not an authoritative membership record for
the primary parent ORA result.

The previously proposed exact-membership criterion therefore cannot be evaluated against a
preserved parent artifact and is withdrawn for that reason. This is an evidentiary limitation,
not a discretionary relaxation after observation of the Study 3 endpoint.

The highest-fidelity preserved aggregate references are:

- Spearman ρ = 0.922857, serialized to six decimal places;
- v29 significant-pathway count = 74;
- v45 significant-pathway count = 47;
- shared significant pathways = 46;
- union = 75;
- Jaccard = 46/75.

The pre-execution consistency check passes only if:

1. the regenerated significant-pathway counts are exactly 74 under v29 and 47 under v45, with 46
   shared pathways and 75 in the union;
2. Jaccard is therefore exactly 46/75;
3. the regenerated Spearman value, serialized using the parent implementation's six-decimal
   format, is exactly `0.922857`; and
4. the frozen input hashes, pathway filtering, pathway alignment, and deterministic
   implementation match the authoritative parent configuration.

Exact historical pathway-set membership is not imposed as an independent criterion because the
full parent significant-pathway identifier sets were not preserved as an authoritative artifact.

Any failure to satisfy these conditions, or any discrepancy in frozen input hashes, pathway
filtering, pathway alignment, or deterministic implementation, halts execution before the
Study 3 contrast.

This entry is a documented pre-execution protocol amendment to an unexecutable registered
consistency check. It changes neither the Study 3 endpoint, directional hypothesis, estimator,
significance threshold, uncertainty procedure, nor disposition rule. The amendment was made
before implementation or execution of stages 19–21 and before computation of `C_ORA`, `C_GSEA`,
or their directional comparison.

---

## DEC-010 — Parent-study data provenance and disclosure disposition (2026-07-17)

The original parent-study configuration preserved the GWAS Catalog source-file SHA-256 but left
the human-readable release identifier as `REPLACE_WITH_RELEASE_ID`.

On July 7, 2026, the only surviving copy of the exact source file was recovered from a system
Trash directory associated with prior deleted workspaces. The recovered file matched the
previously recorded SHA-256 `bfaffa3ac3f742a207f968f45e09d66a73b5d617783936c77cc0e09242c3480b`
exactly, establishing byte identity with the parent analytical input.

The release label `v1.0.2-associations_e115_r2026-01-19_full` was reconstructed from the
recovered file's directory name; it was not independently established from GWAS Catalog metadata.
The file was restored to the canonical project workspace and reverified. The derived mapping
artifact was regenerated as `comparison.tsv`, SHA-256
`8f27c1dba92c2af3ec4c427b1bf89cdb03d67f68adef09c1c9131111f6ba2bb1`, and reproduced the parent ORA
result: ρ = 0.922857, Jaccard = 0.6133, and significant-set counts 74/47 with 46 shared.

The parent-study provenance record therefore contains four distinct deficiencies identified
during the July 2026 review:

1. the human-readable GWAS Catalog release identifier remained a placeholder;
2. after three workspace deletions, the only surviving copy of the pinned source file was
   recovered from system Trash;
3. no contemporaneous parent execution records predating July 8, 2026 were preserved in
   `RUN_LOG.md`; and
4. no committed artifact preserved the primary significant-pathway identifier sets.

The source-file and regenerated mapping-artifact hashes establish input-byte identity, and the
parent ORA result is reproducible at the preserved aggregate level: 74 significant pathways under
v29, 47 under v45, 46 shared pathways, union 75, Jaccard 46/75 = 0.6133, and Spearman ρ =
0.922857. However, historical pathway-level membership cannot be independently verified against a
preserved parent output. These are data-custody, execution-provenance, and
reproducibility-granularity limitations; no input-byte discrepancy or aggregate scientific-result
discrepancy was found.

I authorize inclusion of a concise data-provenance clarification in the Data Availability or
Reproducibility section of the resubmitted manuscript. This decision does not alter any
scientific claim.

> **Data provenance clarification.** The parent analysis preserved cryptographic hashes for its
> GWAS Catalog input and regenerated mapping artifact, but several provenance elements were
> incomplete. The human-readable release identifier remained a placeholder; after three workspace
> deletions, the only surviving copy of the pinned source file was recovered from system Trash; no
> contemporaneous parent execution records predating July 8, 2026 were preserved in `RUN_LOG.md`;
> and the primary significant-pathway membership table was not retained. The recovered source file
> matched the previously recorded SHA-256 exactly, and regeneration of the locus-to-gene mapping
> artifact reproduced the preserved aggregate ORA results: 74 significant pathways under v29, 47
> under v45, 46 shared pathways, union 75, Jaccard 46/75 = 0.6133, and Spearman ρ = 0.922857.
> Input-byte identity and aggregate-result reproduction are therefore verified, but pathway-level
> membership cannot be independently checked against an archived parent output. These limitations
> concern data custody, execution provenance, and reproducibility granularity; no input-byte
> discrepancy or aggregate scientific-result discrepancy was found.

The authorized provenance clarification will be inserted during preparation of the next
submission version, after the Study 3 scientific disposition is frozen and before the manuscript
is submitted. This decision does not authorize an immediate edit to the locked manuscript source.
