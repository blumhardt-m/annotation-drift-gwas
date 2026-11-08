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

---

## DEC-011 — Release of the temporary Study 3 registration hold (2026-07-17)

I reviewed the controlling Study 3 preregistration, the Gate 0 verification recorded in DEC-008,
the §3.4 protocol amendment recorded in DEC-009, and the temporary registration hold established
by DEC-006, DEC-007, `EXECUTION_BOUNDARY.md`, `AGENTS.md`, and `TASK.md`.

The factual basis for the temporary hold has been resolved. The exact controlling preregistration
artifact, SHA-256 `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`, was
downloaded from and byte-matched to the immutable embargoed OSF registration snapshot.
Registration was initiated on July 8, 2026 and approved on July 9, 2026 before implementation or
execution of Study 3 stages 19–21.

I therefore release the temporary Study 3 registration-verification hold.

This decision does not create supplementary scientific authority or add a new implementation or
execution gate. Authority to proceed derives from the sequencing and execution provisions already
frozen in §8 of the registered protocol. The registered protocol's prerequisites have been
completed, subject to the amended §3.4 ORA consistency check in DEC-009.

The execution-contract files shall be updated solely to:

1. replace the obsolete statements that OSF registration is pending or unverified;
2. record the verified controlling artifact hash and registration approval date;
3. identify `study3-ranked-enrichment` as the registered Study 3 branch;
4. specify that the branch is to be created from `review/study3-freeze-and-gate0-20260716`, because
   the reviewed protocol, checksum lineage, Gate 0 record, and DEC-004 through DEC-011 are not
   present on `main`; and
5. preserve the registered additive-only execution boundary and all existing prohibitions against
   modifying frozen inputs, frozen parameters, `comparison.tsv`, prior results, or the manuscript
   during Study 3 execution.

The §3.4 ORA consistency preflight is an authorized and required real-data preflight under the
registered protocol. It is not the confirmatory Study 3 endpoint. It shall be run against the
frozen parent inputs and shall halt execution mechanically if any DEC-009 condition fails.

After the execution-contract reconciliation is reviewed and committed, the registered sequence may
proceed on `study3-ranked-enrichment`: implementation of stages 19–21, the §3.4 preflight, and—only
if that preflight passes—the registered ranked-enrichment, bootstrap, endpoint evaluation,
disposition, and output freeze.

This decision does not modify the registered endpoint, directional hypothesis, estimator,
thresholds, random seeds, uncertainty procedure, stopping rules, or disposition criteria.


---

## 2026-11-07 — Record corrections after independent assessment (DEC-012)

Basis: GWAS_ANNOTATION_DRIFT_INDEPENDENT_ASSESSMENT_2026-11-07_v1.2 (md sha256 91fdc8fbc64ff694ccd061813a3fadaeb060c1e12d14d751c5aec13e786aadd1) and
VERIFICATION_REPORT_20261107.md (sha256 329c432bd55050d56ee5b2e300565d81869c4f6dc7fc56da800ca986c53c837e).
Scope: corrections to the record of the parent study (v1.0_submission, 37272b2) and its deposited results.
No artifact listed in LOCKED_PROTOCOL_FILES.md is modified by this entry.

(a) Preregistered prediction I > 5% (prereg/prereg.md §1) was not met: I = 4.15% (row-weighted), Wilson
    CI [4.10, 4.19]; results/tables/bootstrap.txt records "CI upper bound < 5%: EFFECT IS TRIVIAL". The
    negative-result rule (I < 5% AND J > 0.90) did not trigger because J = 0.613. The manuscript did not
    disclose the failed prediction. Status: disclosed here; the manuscript framing ("modest drift") stands
    as a description, not as a preregistered outcome.
(b) Window sensitivity was run at 250/500/1000 kb (results/tables/sensitivity.tsv); the frozen grid
    map.windows_kb is 10/25/50/100/250/500 kb. No amendment was recorded. Status: deviation recorded.
(c) The reported bootstrap (10,000 i.i.d. resamples of association rows; src/05_bootstrap.py) is not the
    frozen chromosome-block bootstrap (1000 reps, seed 42). Rows are not independent. The Wilson CI, which
    the protocol names primary, was reported. Status: deviation recorded; CI widths from (c) are not
    to be cited as evidence of precision.
(d) Frozen inclusion rules (European ancestry; sentinel SNP per locus with duplicates collapsed by
    locus_id) are not implemented in src/03_parse_gwas.py; the mapping table is association-level
    (721,406 rows for 291,292 distinct SNP identifiers). The manuscript Methods omit ancestry.
    Status: deviation recorded.
(e) Manuscript §2.5 describes per-version background universes (19,019 / 19,114) and argues from them;
    the analysis used the frozen union universe (19,315 Ensembl IDs; 20,258 symbols). Status: reporting
    error recorded; the code behaviour is the one of record.
(f) Fig. 5B annotation "Full-foreground Jaccard ≈ 0.62 (saturation artifact: 89% universe coverage)" and
    manuscript §3.3 rest on results/tables/baseline_random.tsv (mean 0.6197, SD 0.4649, CI [0.000, 1.000],
    n = 100). src/15_baselines.py returns Jaccard = 1.0 when both significant sets are empty and warns
    that the permutation may not be shuffling; results/tables/saturation_curve.tsv is 0/NaN throughout.
    The statistic is degenerate. Status: the "by construction ≈ 0.62" claim is withdrawn from the record;
    the saturation point itself (coverage 89.3% / 89.6%) stands and is recorded under (i).
(g) "Drift is structured rather than random" rests on rho = 0.275 (differential sets) versus a size-matched
    null of 0.150 with 95% CI [0.010, 0.314]. The observed value lies inside the null interval. Status:
    claim withdrawn from the record.
(h) Instability denominators: 29,920 / 721,406 association rows = 4.147%; 10,129 / 291,292 distinct literal
    SNP identifiers = 3.477% (strict single-rsID subset 10,125 / 291,199 = 3.4770%); zero identifiers with
    conflicting assignments (VERIFICATION_REPORT_20261107.md item 6). Status: both denominators are to be
    reported wherever the rate is stated.
(i) Scientific limitation of record: the ORA foreground is 89.3% (v29) / 89.6% (v45) of the universe
    because all genome-wide-significant Catalog associations were pooled across traits. Pathway
    significance in this regime is driven by the few unhit genes per pathway; 47 of 110 symbol-based
    significance flips are one-gene overlap changes. The Discussion's per-trait non-identifiability
    claim is not evidenced by this design. Status: limitation recorded; no per-trait claim is made
    on this evidence.
(j) Top-K reporting: Methods K = 10/20/50/100; Results "K >= 100"; results/tables/jaccard_robustness.txt
    K = 25/50/100/200. Status: minor inconsistency recorded.
Open (not corrected): identity of the bytes uploaded to ScholarOne 2026-03-12 and the source of fig6.tif
    in the trashed submission_preparation_log.txt (NOT_COMPUTED).

## 2026-11-07 — Lineage reconciliation (DEC-012-a)

Shared ancestor `4da4768e130c1a90534a2f57735d5aa3ddfd3902`; main was fast-forwarded to
`review/study3-freeze-and-gate0-20260716` (12 commits, tip
`8e68ee6e78e4646e8e43f093e863da67a5fdcbed`), then merged with fetched origin/main
`61df5e0df3fd5c0ac3af38099fe9fb4fe6263be7`. The single reconciliation merge is
`491c90cee28a59b13dbbaf567b0df0cf11a7c67a`, with those review and public tips as its two
parents; tag `reconciled-20261107-491c90c` identifies that merge. README auto-merged,
preserving both status statements; its merge-time SHA-256 is
`f0ee5a2d166c016459923b0da0a3021b890cfebda3fed9a082b1356143455266`.
Not pushed; publication decision D-A1b remains pending. This is the owner's explicit
local reconciliation instruction executed by Codex; it does not authorize publication.

Fresh post-merge audit: `~/agent_work/annotation_drift_lineage_verify_20261107/VERIFICATION_LINEAGE_20261107.md`.
Report SHA-256 `c1137b69e6564d9ad3157119b53aa1aed3d796847504de23e7e352ac6aa795ec`; commands.log SHA-256 `ea117fcf076262eb877748c192fa75c8db576d24262f98f2acd43c801ffce299`, copied from the audit's SHA256SUMS.

The audit verified both histories reachable, both prereg sidecars, the controlling Study 3
prereg SHA-256 `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17`,
and RUN_LOG.md byte identity (SHA-256
`09946bd7ef47c4e7c4b1b961511b6c4238e086ee0bb2a2ba1b06cbd900966456`). No Stage 19–21
implementation or results/study3_ranked_enrichment/ path was found at main. Repository
endpoints were clean and unchanged, and the four-file portfolio bracket matched.
The audit precedes DEC-012/012-a/013/014 and the separate README Status correction.

Locked-file hashes below are transcribed programmatically from audit item 3: 52 available
committed paths are byte-identical between the review tip and the merge. Five listed paths
have no blob in either commit and remain ABSENT / NOT_COMPUTED; no local-file hash is
substituted for historical evidence. The disposed foreign results path is excluded per the
inventory. This is a tracked-lineage verification, not a claim that every named path was hashed.

| Path | Bytes at merge | SHA-256 at both review and merge / limitation |
|---|---:|---|
| `AGENTS.md` | 3677 | `fea7d8bb6bbf8f2446213e88f70751859f7d5a40e4960de9efe3aacf9986fe81` |
| `EXECUTION_BOUNDARY.md` | 6041 | `fbe566c4c01650ec6c3f0771769c5e773133bd513040b240426d2b7c77decb6f` |
| `LOCKED_PROTOCOL_FILES.md` | 7008 | `c962e085e16397a882e2a0b01022f38a4d0fd40a02479756dbfa010072c2ebc9` |
| `ONBOARDING_SOP.md` | 3919 | `edc44c095641b9bd610541a85475581e46f42238597912612d06b6ae3a092a55` |
| `REVIEW_CHECKLIST.md` | 2738 | `12c0dc67f6535346ebf179eb5f109c8a8aba20fc30644e0d6e5007ea8551167b` |
| `REVIEW_TRIGGER.md` | 4208 | `c7ede8753370ed9cf5256691b28b544e3ba0f24bcd4f84c6dbaaffb4cba6b630` |
| `annotation_drift_manuscript.pdf` | — | ABSENT in both Git trees; comparison NOT_COMPUTED |
| `data/processed/mappings/comparison.tsv` | — | ABSENT in both Git trees; comparison NOT_COMPUTED |
| `data/raw/gencode.v29.annotation.gtf.gz` | — | ABSENT in both Git trees; comparison NOT_COMPUTED |
| `data/raw/gencode.v45.annotation.gtf.gz` | — | ABSENT in both Git trees; comparison NOT_COMPUTED |
| `data/raw/gwas_catalog.tsv` | — | ABSENT in both Git trees; comparison NOT_COMPUTED |
| `data/raw/reactome.gmt` | 1006622 | `a34f43fd763a6800a9a94310c416ae511e84724d87a3c255f35d13764d12bebb` |
| `manuscript/annotation_drift_manuscript.pdf` | 32991 | `c331bd69c62d32dd628b2fac51c6419e469c687b69c2e2415308d417c97e90dd` |
| `manuscript/build_manuscript.py` | 48621 | `b00b708f44f834edf5dbfeadf64367dbf0da552de3376414431e373cf7ae6d22` |
| `pre_auth_check.sh` | 5227 | `9e98540049a28da64f5c5bbb1375d27130578ae4bae24d424104e6e6978d5184` |
| `prereg/frozen_params.sha256` | 92 | `758696c08cb936760e855374103854388912e9151e2d2ccd0a0eb17ad89814f5` |
| `prereg/frozen_params.yaml` | 5034 | `8407f51ff47b8ba4cdfd0f58b22a23a94f3d221dd295ae828cbee225ec01322f` |
| `prereg/prereg.md` | 5974 | `68aed7121f16d92ca6d67e9f22ac1e3541c660b5afe1cea3ad2b2b30b67073c5` |
| `prereg/study3_prereg.sha256` | 88 | `0abe0e566a0e4089c82271ad267095339f8d67f4ce1792eb51ced3968565c05d` |
| `prereg/study3_prereg_v1.0.md` | 22370 | `29fe3ff1bfe017cc76fec92c4b118a6a21c877ecd83ef71ae562de5cc2651d17` |
| `results/figures/.gitkeep` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `results/figures/FIGURE_MANIFEST.txt` | 10587 | `51af7acbe5756b13f1a031559c190ec2b7498de03c4d121db0d895aa0dd7e4af` |
| `results/figures/fig1_pathway_scatter.pdf` | 44411 | `ad0eeffe8832eecf4ececb90578d109e64f807707e790806507b464867fbbd7a` |
| `results/figures/fig2_instability_windows.pdf` | 13180 | `9c9ccc96869e227f510d7b390b6606830460340d8d487221f3c115480e9f7991` |
| `results/figures/fig3_drift_decomposition.pdf` | 15684 | `8081d0fbb49d81ec1e49032c030438e9ca26c5a9e1a4e3a2e58dd83edafd7315` |
| `results/figures/fig4_concordance_summary.pdf` | 19023 | `1396fbcca889c0feb99998f16ac69360a3e51f2d87996a712770cb2a25eb295f` |
| `results/figures/fig5_calibration.pdf` | 32128 | `fc428b1cac787b4c21c8ceb0f0093f41278bae2817cfe21985b744762bd89e46` |
| `results/final_manuscript_update_log.txt` | 5087 | `bd62f730f16c43a8965eba0e755831a15f3acff88c5bc76874504f5b06d045d6` |
| `results/study1_GO_replication/go_rank_scatter.png` | 72461 | `f5cc67b718d0ae1a8214fcea3ac2201304892b195bbe57f97091f62b64f8eacb` |
| `results/study1_GO_replication/go_replication_summary.tsv` | 157 | `3529d7088cc1abe03a1bd7a3e8d43333e1fb9f98289658ac6e98880273f738e5` |
| `results/study2_v37_v45/instability_summary.tsv` | 193 | `a9d1b7004bb30cdb3a57a8eccff887a82ceb5ffd634c1621ff6c548bb7bea539` |
| `results/study2_v37_v45/rank_scatter.png` | 133635 | `2bd2d91deea77988b98035f82458369745e53b8067515751d48f1375869c3bb0` |
| `results/study2_v37_v45/study2_summary.tsv` | 430 | `72ecce675913a3ba2fdb24690d97b694411fa1e6c0d3642a4e54f8c92df89124` |
| `results/tables/.gitkeep` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `results/tables/baseline_fixed.tsv` | 168 | `2bd20f8a6dc7edded8b0f98ffe1fcee4112ed972d3244d50564307690821d73f` |
| `results/tables/baseline_random.tsv` | 150 | `88e33351e320bb5d3379a6c8d38adac954b32a94f15c411c436b6ebe71c5108d` |
| `results/tables/bootstrap.txt` | 231 | `c8e8c16dca7b45acc460021ea76206a7339e1dbfdc38f89a15b3bbb423052a02` |
| `results/tables/category_summary.tsv` | 544 | `2f0f305c7a435068493dbd798b9027a6ebe2b4b96ed305dea8e54509208109a9` |
| `results/tables/difference_classification.tsv` | 2980440 | `11329c05b081108963f214a25fe592aa567264516c5b51b276d555816573a27b` |
| `results/tables/difference_classification_summary.txt` | 255 | `613e7f786aefa060b822f6d6f449949da0e07a16e852b97fab73034651f2032d` |
| `results/tables/differential_baseline_random.tsv` | 130 | `db0d5c0f52a07ffd147c76bc61e0d5c4fe4c4ec7c38925d1c5c902df153d4c01` |
| `results/tables/differential_enrichment.tsv` | 90 | `bd4ed15760814c9657efac2fe792df7a051f45a160da70e3b6c24c628a13bda9` |
| `results/tables/enrichment_diagnostics.txt` | 4842 | `fa7a58168a2b0209004d5ff6d6c82d6f1a80f913781d2b26f5ef01e11c9039d1` |
| `results/tables/enrichment_ensembl_id.txt` | 3355 | `63a70dc46ab1756584caf994c7243961b4f213d09befb9cd999b0190449e2c5f` |
| `results/tables/enrichment_harmonized.txt` | 2638 | `2934a83b6062963824764a7ba4847adc942b8c3f70002a9669338fa59b563761` |
| `results/tables/enrichment_jaccard_by_window.tsv` | 145208 | `355e7bf5da3eb1fe5b9b6bb184cc52917eae40c27edd1d1b3bef1ddd21990718` |
| `results/tables/flip_rate.tsv` | 182 | `cf2077217eb3ef6a9ad80f80acab30956fc8b96c68f920587cd00d2d79c0fb12` |
| `results/tables/jaccard_robustness.txt` | 2466 | `9d4ee22fb0aefdb8e22623b12c30758453fab62b6709b8548827e00a802f53eb` |
| `results/tables/saturation_curve.tsv` | 179 | `365d4257f656f561418b755ecd55f278e2fe62c39353ca32063505833bcb6e5c` |
| `results/tables/sensitivity.tsv` | 144 | `f982a012328e4f9f8b8ebdaaed9681232d7bc04d31c81423cf190abfcb7e61d8` |
| `scripts/study2_v37_v45/export_primary_reactome.py` | 11737 | `e202126a27f23dc5520e96388adbea3ae96eacdf40961043bd7756a7064f6493` |
| `src/04_map_variants.py` | 3466 | `e64c1bfd76e22b6db181b24171788c5760e238b7257181737dc88d1d061056d0` |
| `src/05_bootstrap.py` | 1828 | `804229a524408d2548498b61871394d0777505771e47c2c90c07b057582c2dd9` |
| `src/12_ensembl_id_enrichment.py` | 15264 | `2934db6053aa9013eff6570ea1fc7278b8fed1854824f39de975970308b3f87a` |
| `src/16_flip_rate.py` | 14140 | `7ba8410dda81b4a2435deecd2ae3bf7d9a25bb834aed0d50200667785b1e3d4e` |
| `src/primary_reactome_export.py` | 10240 | `a4d3a52cc614fddc79ee2e6ec1ee302b80356102a1c0c9895531f5e4b6b0fdd1` |
| `validate_config.py` | 6153 | `77f599697f9d4299ae1219cf43777481b81bec12e7f88b98e4db9d52efde8c2b` |

Preview-run deviations (owner-specified record):

1. Output-path deviation: the preview was written under `~/Study/research/08_GENOMICS/annotation_drift_gwas/annotation_drift_lineage_preview_20261107/` instead of the prompt's `~/agent_work/annotation_drift_lineage_preview_20261107/` (commands.log batches 1–3 and report provenance paragraphs). On 2026-11-07 the complete directory was moved to the requested agent_work path; all nine files retained their SHA-256 and byte count. Historical path statements inside the preview were preserved.
2. Failed script-construction attempt: "A functions.exec script-construction attempt failed with a JavaScript SyntaxError before executing commands or writing files" (commands.log line 41). No side effects; retained as an execution event in the transcript.
3. Post-seal wording correction: `verify_preview.py` edited one row of `LOCKED_FILES_CHECK.md` after the first `SHA256SUMS` was written, then resealed the manifest (commands.log lines 244–249). The delivered files match the resealed manifest; the first seal is superseded. The relocation preserved the final resealed bytes.

Audit execution disclosure: the initial fresh audit attempt was stopped after it read an
explicitly excluded Downloads instruction file; its deviation receipt is preserved as
`attempt1_deviation.txt` in the verification directory. It supplied no passing verdict.
A fresh replacement audit produced the report cited above without reading that receipt or
the excluded instructions. The launcher did not expose a plan-mode selector: the replacement
report records actual Default mode and workspace-write sandbox, with read-only behavior
outside its output directory and narrowly scoped escalation for the deliverables. Exact
deployed model/version was unavailable and is labelled NOT_VERIFIABLE in the report.

## 2026-11-07 — Study 3 superseded without execution (DEC-013)   [ADOPTED 2026-11-07 14:27 -08:00, D-C1]
Study 3 (prereg/study3_prereg_v1.0.md, sha256 29fe3ff1…, OSF fv4sg) is superseded by the per-trait
annotation-version pilot specified in DECISION_PACKAGE_annotation_drift_2026-11-07 §D. Basis: assessment
v1.2 §7 and §10. No confirmatory Study 3 execution evidence (stages 19–21, results directory, R/fgsea logs, or an
implementation branch) was found in the inspected scope (VERIFICATION_REPORT_20261107 item 4).
Preserved unchanged: the prereg and sidecar, the OSF registration (embargo disposition per D-A1b), the
viability probe under scratch/study3_viability/, branch review/study3-freeze-and-gate0-20260716 and tag
local-safety/study3-pre-review-20260716. No evaluation of the Study 3 kill rule (C_ORA <= C_GSEA) was found in the inspected evidence; the
parent threshold-amplification claim is therefore neither corroborated nor narrowed by Study 3 and is
carried as "mechanism not established" per DEC-012 (i).

---

## 2026-11-07 — Per-trait pilot preregistration v1.0 FROZEN (DEC-014)
prereg/pertrait_pilot_prereg_v1.0.md sha256 d7e26148c6f540156b11cd294268f3513b613e8f46b01e9bdc5247da41f20dad, sidecar prereg/pertrait_pilot_prereg.sha256 written by Codex under the owner's explicit freeze instruction on 2026-11-07.
Source: DECISION_PACKAGE_annotation_drift_2026-11-07_v1.3 §D (md sha256
29f4ee72018e184e31284a36ee5a8b006a428b51745d7e1d1588948691e7d3e3), owner-adopted D-D1. Deferred field: frozen
roster hash (NOT_RECORDED until the roster script runs). No execution authorized by this entry.
The candidate was copied byte-for-byte without edits; its conditional CANDIDATE status is fulfilled by this sidecar and decision. This is a local freeze, not an external OSF registration.



## 2026-11-07 — Per-trait pilot preregistration amended v1.0 → v1.1 before any execution (DEC-015)

Controlling file after adoption: prereg/pertrait_pilot_prereg_v1.1.md, sha256 b78907970faafa656c91ba55b2a948f405045dd8165fcceac9ee1bf70333d1db; sidecar prereg/pertrait_pilot_prereg.sha256 rewritten by owner over the v1.1 bytes. v1.0 preserved unchanged at prereg/pertrait_pilot_prereg_v1.0.md (sha256 d7e26148…, DEC-014); v1.0 sidecar retained as prereg/pertrait_pilot_prereg_v1.0.sha256 for provenance.

Basis: activation preflight package 2026-11-07 (manifest sha256 bd41cd6d…), owner dispositions 17:28 -08:00. Execution state at amendment: none — no roster, enrichment, C2, sensitivity, or implementation exists under any version.

Amendment contents (prereg §11): v37/v45 pair and GO:BP deferred (both v1.0 provenance statements were text errors: no v37 pin exists in the clone; Study 1 GO inputs were unpinned mutable downloads); v29/v45 TSS tables bound by deterministic regeneration from the pinned GTFs with src/02_parse_gencode.py::parse_gtf_to_tss and byte comparison to the retained tables before enrichment, mismatch leaving the binding unresolved; roster order (code-point; extension appends, indices preserved); C2 draw order, generator scope and exact sampling calls; C2 percentile = 950th order statistic with strict exceedance; R executable pinned by absolute path. Thresholds, precedence, outcomes, extension rule, estimand, mapping, ORA and fgsea parameters unchanged.

Authority: owner adoption of the amendment is recorded here; no implementation, roster, or empirical execution is authorized by this entry. Membership and allocation remain as recorded in PDG-033 (nonmember, OWNER_DECISION_PENDING).
