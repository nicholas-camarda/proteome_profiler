## Context

The current pipeline assumes one LI-COR workbook per analysis group, averages duplicate membrane spots inside each workbook, widens directly to one column per treatment group, and computes descriptive fold-change ratios from those collapsed values. That is acceptable for the original exploratory VEGFRi/Dox workflow, but it is not a valid representation for collaborator studies that have multiple biological samples per group.

The next required use case is a factorial mouse-plasma study with `sex x treatment x biological replicate`. The user does not want a direct male-versus-female comparison. They want treatment-versus-vehicle contrasts within each sex, with analyte-level multiplicity correction and an explicit significant-hit list.

The design therefore needs to do three things at once:
- preserve biological replicate observations as first-class data
- support subgroup-restricted treatment contrasts such as `vehicle vs aldosterone` within `male` and within `female`
- separate exploratory threshold screens from inferential statistical results

## Goals / Non-Goals

**Goals:**
- Introduce a sample-level input model that can represent multiple biological replicate workbooks per treatment arm.
- Keep technical replicate handling limited to the duplicate membrane spots within each workbook.
- Run analyte-level differential testing within configured subgroup strata.
- Produce raw p-values, multiplicity-adjusted p-values, and significant-hit tables with configurable adjustment method and `BH` as the default.
- Make output structure and labels explicit so exploratory fold-change screens are not confused with inferential significance.
- Preserve the current one-workbook-per-group exploratory workflow for analyses that do not have biological replicates.

**Non-Goals:**
- Add direct sex-versus-sex contrasts in this change.
- Implement mixed-effects or hierarchical Bayesian models.
- Solve all possible factorial designs; this change targets stratified two-group treatment comparisons with biological replicates.
- Remove the exploratory threshold-based workflow entirely.

## Decisions

### 1. Require an explicit sample manifest for replicate-aware analyses

Replicate-aware runs SHALL be driven by a manifest that maps each workbook to `sample_id`, treatment group, and optional stratification variables such as sex.

Why this decision:
- Filename parsing is brittle and will silently fail on collaborator data.
- A manifest makes the biological analysis unit explicit.
- The same mechanism can later support additional stratifiers or covariates.

Alternatives considered:
- Derive sample metadata from workbook filenames.
  Rejected because it creates an undocumented contract and is too fragile for collaborator-facing analysis.
- Keep one-workbook-per-group input and ask users to average biological replicates upstream.
  Rejected because it destroys within-group variance and makes inferential statistics invalid.

### 1a. Keep the legacy exploratory path for non-replicate analyses

Analyses that truly have one workbook per group and no biological replicate structure SHALL remain runnable without a manifest. Those runs shall stay in the exploratory workflow and must not claim replicate-aware inferential outputs.

Why this decision:
- The existing VEGFRi/Dox use case is still valid as an exploratory screen.
- Forcing every analysis through the replicate-aware path would add unnecessary configuration burden.
- The real failure case is ambiguous multi-workbook input without sample metadata, not the existence of legitimate single-workbook exploratory analyses.

Alternatives considered:
- Make the manifest mandatory for every run.
  Rejected because it would degrade the simple exploratory workflow without adding value for non-replicate projects.

### 1b. Dispatch analysis mode from config, not from filename heuristics

The run configuration SHALL determine whether the pipeline uses the legacy exploratory path or the replicate-aware inferential path. If a sample manifest is configured, the pipeline enters replicate-aware mode. If no manifest is configured and there is exactly one workbook per group, the pipeline stays in the legacy exploratory mode.

Why this decision:
- The user needs the implementation to be explicit and inspectable.
- Mode dispatch should not depend on brittle file-naming conventions.
- This keeps the old and new workflows understandable and testable.

Alternatives considered:
- Auto-detect the mode entirely from the files present in the data directory.
  Rejected because auto-detection becomes ambiguous once collaborators use mixed naming conventions or partial datasets.

### 2. Introduce a long-format sample-level analysis dataset

After workbook import, the canonical intermediate dataset SHALL contain one row per `analyte x sample_id`, along with treatment and subgroup metadata. The existing duplicate membrane spots within a workbook will still be averaged, but biological replicates will remain separate rows.

Why this decision:
- Statistical testing needs access to per-sample variability.
- Long format is the natural representation for grouped summaries and model fitting.
- It prevents the current premature collapse to one column per group.

Alternatives considered:
- Continue widening to one column per group and store list-columns of sample values.
  Rejected because it complicates downstream modeling and validation.

### 2a. Use a minimal manifest schema with optional extra metadata

The replicate-aware manifest will require `sample_id`, `workbook_path`, and `treatment`. Additional columns such as `sex`, `batch`, or `notes` may be present, and the analysis config will name which metadata column should be used for subgrouping.

Why this decision:
- It keeps the minimum contract small and clear.
- It supports the immediate within-sex use case without hard-coding `sex` into every analysis.
- It leaves room for future covariates without redesigning the manifest format.

## Configuration Contract

The replicate-aware path should be driven by explicit config keys rather than implicit file layout rules.

The user-editable analysis metadata should live in a dedicated config script such as `scripts/config/analysis_config.R`, while path-resolution and validation helpers stay in `scripts/helpers/project_paths.R`. This makes it obvious which file collaborators are meant to edit.

Recommended config shape:

```r
list(
  user = "lauren",
  analysis_slug = "aldo_plasma_cytokine_xl",
  sample_manifest = "manifests/aldo_plasma_samples.csv",
  protocol_preset = "cytokine_xl",
  info_fn = "output/cytoXL array kit - protocol.xlsx",
  subgroup_var = "sex",
  treatment_var = "treatment",
  comparisons = list("vehicle" = c("aldosterone")),
  ref_thresh_to_filter = c(150),
  p_adjust_method = "BH",
  alpha = 0.05
)
```

Recommended minimal manifest:

```csv
sample_id,workbook_path,treatment,sex
M01,data/m01.xlsx,vehicle,male
M02,data/m02.xlsx,vehicle,male
M03,data/m03.xlsx,aldosterone,male
M04,data/m04.xlsx,aldosterone,male
F01,data/f01.xlsx,vehicle,female
F02,data/f02.xlsx,vehicle,female
F03,data/f03.xlsx,aldosterone,female
F04,data/f04.xlsx,aldosterone,female
```

Key implementation rule:
- In replicate-aware mode, filenames are just file paths. The manifest defines the experimental design.
- In legacy exploratory mode, filenames may still encode the group label because there is no manifest.
- The config file is the one place users should edit paths, manifest locations, thresholds, and comparisons.

## Execution Flow

### Legacy exploratory mode

1. Load the configured group-level analysis entry.
2. Read exactly one workbook per configured group.
3. Average duplicate membrane spots within each workbook.
4. Join protocol-derived analyte labels from the extracted workbook.
5. Build one wide analyte-by-group dataset.
6. Optionally inspect the low-signal reference panel with `find_ref_thresh.R`.
7. Generate fold-change waterfall and barplot outputs.
8. Keep all outputs labeled as exploratory, not inferential.

### Replicate-aware inferential mode

1. Load the configured analysis entry and sample manifest.
2. Validate manifest columns, uniqueness of `sample_id`, existence of workbook paths, and presence of the configured subgroup column.
3. Read one LI-COR workbook per biological sample listed in the manifest.
4. Average duplicate membrane spots within each workbook.
5. Join protocol-derived analyte labels from the extracted workbook.
6. Build one long-format table with one row per `analyte x sample_id`, carrying treatment and subgroup metadata from the manifest.
7. For each configured subgroup level, run the configured treatment-versus-control contrast on log2-transformed signals.
8. Compute raw p-values, then adjust them within each comparison family using the configured `p.adjust` method, default `BH`.
9. Record low-signal status separately from inferential significance.
10. Write one canonical results file per comparison, a run-level index, and a combined workbook convenience export.
11. Keep exploratory threshold plots separate from inferential result tables.

### 3. Use per-analyte linear modeling on log2 signals for within-stratum comparisons

Within each configured stratum (for example `male` or `female`), the pipeline SHALL test the configured treatment-versus-control contrast on log2-transformed analyte signals using a per-analyte linear model. For the target `vehicle vs aldosterone` use case, this is equivalent to a two-group replicate-aware comparison while remaining extensible to future covariates.

Why this decision:
- The current ratio-only approach has no variance model.
- Linear models are transparent, easy to validate, and more extensible than hard-coding pairwise t-tests everywhere.
- Log2 scale aligns the effect-size calculation with existing fold-change interpretation.

Alternatives considered:
- Per-analyte Welch t-tests only.
  Viable for the immediate two-group case, but less extensible once subgrouping or extra covariates appear.
- Limma or empirical Bayes moderation.
  Potentially attractive later, but not necessary to satisfy the immediate requirement and would add conceptual and implementation overhead.

### 4. Apply multiplicity correction within each comparison family using `p.adjust`

For each stratum-specific treatment-versus-control comparison, the pipeline SHALL adjust analyte-level p-values across the tested analytes using `p.adjust`. The default method SHALL be `BH`, but the method MUST be configurable.

Why this decision:
- The user explicitly needs a significant-hit list after multiple-testing correction.
- `p.adjust` is native, well-understood, and sufficient for methods such as `BH`, `bonferroni`, `holm`, and `BY`.
- Correction by comparison family matches the actual question being asked, e.g. `male: vehicle vs aldosterone`.

Alternatives considered:
- Global adjustment across all strata and all contrasts.
  Rejected as the default because it mixes distinct inferential families and can be overly conservative for the stated within-sex questions.

### 5. Separate exploratory and inferential outputs

The output tree SHALL keep threshold-based exploratory plots separate from inferential result tables and any inferential plots. Any output labeled "significant" in the inferential branch must refer to adjusted statistical significance, not fold-change threshold passing.

Why this decision:
- The current terminology is misleading.
- Collaborators need to distinguish descriptive screening from tested findings quickly.

Alternatives considered:
- Preserve the current naming and document the caveat.
  Rejected because the current names encode the wrong semantics directly into the artifact tree.

### 6. Allow inference at `n = 2` per arm, but warn when `n < 3`

The inferential path SHALL allow analyte-level testing when each arm has at least two biological replicates. The pipeline SHALL mark comparisons with fewer than three replicates per arm as low-replication analyses and emit an explicit warning in the results metadata.

Why this decision:
- Two replicates per arm is the technical minimum for estimating within-group variance in a simple two-group comparison.
- Blocking `n = 2` would make the pipeline unusable for small pilot studies that still need a formal analysis.
- Carrying a warning preserves usability without pretending those results are stable.

Alternatives considered:
- Require three replicates per arm to run any inference.
  Rejected because it is too strict for the immediate collaborator use case and small exploratory animal studies.

### 7. Run inferential testing on all analytes, and flag low-signal analytes separately

The inferential path SHALL test all analytes that meet the replicate requirements. Low-signal status SHALL be recorded as a separate flag derived from the reference/background rule, rather than being used as a mandatory pre-test exclusion by default.

Why this decision:
- Hard pre-filtering changes the tested analyte family in a way that can be easy to overlook.
- Low signal is partly an assay-quality concern, not the same thing as inferential non-significance.
- A separate flag allows users to inspect adjusted significance together with assay-confidence information.

Alternatives considered:
- Drop low-signal analytes before testing by default.
  Rejected because it silently narrows the inferential family and can hide analytes that users may still want to review.

### 8. Write one results file per comparison, plus an index and a convenience workbook

The canonical inferential outputs SHALL be one machine-readable results file per stratum-specific comparison plus an index file summarizing the run. The pipeline SHOULD also emit a combined Excel workbook as a collaborator-friendly convenience export.

Why this decision:
- Per-comparison files are easier to inspect programmatically and scale cleanly as comparisons increase.
- An index file gives one place to see what was run and where each result file lives.
- A combined workbook is still useful for collaborators who prefer a single handoff file.

Alternatives considered:
- Only one combined workbook per run.
  Rejected as the canonical output because it is less robust for programmatic consumption and becomes awkward as analyses grow.

## Risks / Trade-offs

- [Manifest friction] → Mitigation: provide a simple template and validate required columns with direct error messages.
- [Model assumptions on small n] → Mitigation: use log2-transformed signals, report group sample sizes, and document the exact test family in the results.
- [Backward compatibility confusion] → Mitigation: keep the exploratory legacy path available for true one-workbook-per-group analyses, but fail fast when users attempt replicate-aware input without the required metadata.
- [Terminology drift between old and new outputs] → Mitigation: rename inferential outputs and update README/script comments to define each branch explicitly.
- [Future designs may need more than within-stratum two-group tests] → Mitigation: structure the analysis around formulas and manifest metadata rather than hard-coded ratio logic.

## Testing Strategy

The repo currently has little or no formal test coverage, so the replicate-aware refactor must be introduced behind a deliberate test harness rather than relying on ad hoc script runs.

### Test layers

1. Unit tests
   Cover path/config helpers, manifest validation, LI-COR import invariants, threshold helper behavior, multiplicity-method validation, and model-output semantics.
2. Integration tests
   Cover dataset assembly from fixture workbooks, mode dispatch, within-stratum comparisons, and inferential output generation.
3. Regression tests
   Lock down the current exploratory behavior so the refactor does not break one-workbook-per-group analyses.
4. Smoke tests
   Run the entrypoint scripts on tiny fixtures and assert output files and TSV payloads exist. Plot-image pixels should not be snapshot-tested.

### Test frameworks

- Use `testthat` for R code and R entrypoint behaviors.
- Use `pytest` for the Python protocol extractor.

### Package installation

- Use `pak` as the standard R package installer in `scripts/install_packages.R`.
- The install script should bootstrap `pak` if it is missing, then install any missing packages via `pak::pkg_install()`.

### Priority test categories

- P0: input labeling and import correctness
- P0: legacy exploratory workflow stability
- P0: replicate-aware inferential core
- P0: output semantics and routing
- P1: path/config resolution
- P1: robustness and perturbation checks such as row-order invariance
- P2: plot/render smoke tests

### Baseline tests required before the refactor

- Protocol extractor argument/path handling, including preset resolution, `--pages`, missing-PDF errors, and output-path creation.
- Path/config helper tests for repo-root detection, path resolution, output-root construction, and protocol-workbook errors.
- Legacy import tests for duplicate-spot averaging, non-`.xlsx` file ignoring, unknown group prefixes, duplicate group labels, `S001...` ordering validation, and `bad_analytes.xlsx` behavior.
- Threshold-helper tests for `ref_coords`, numeric threshold mode, combined mode, and empty-selection handling.
- Fold-change workflow tests for pair-specific filtering, `Reference`/`Negative` analyte exclusion, and threshold-hit marking.
- Output-tree tests for `main_analysis/`, `threshold_diagnostics/`, and `select_analytes/` generation.
- Legacy end-to-end smoke tests for `find_ref_thresh.R`, `main.R`, and `select-analytes-analysis.R`.

### Replicate-aware tests required before trusting inferential output

- Manifest validation tests for missing required columns, duplicate `sample_id`, missing workbook paths, missing subgroup columns, and unsupported `p_adjust_method`.
- Mode-dispatch tests to prove manifest-driven configs enter inferential mode and non-manifest one-workbook-per-group configs stay in exploratory mode.
- Sample-level dataset-assembly tests that verify one row per `analyte x sample_id` and prevent accidental collapse back to one value per group.
- Within-stratum comparison tests that run `vehicle vs aldosterone` separately in `male` and `female` and do not generate cross-sex contrasts unless requested.
- Minimum-replicate tests for `n = 2` warning behavior and `n = 1` untestable behavior.
- Multiplicity-correction tests confirming `BH` default behavior and correct comparison-family scoping.
- Low-signal semantics tests confirming low-signal analytes are flagged but still tested by default when inferentially eligible.
- Inferential output-structure tests confirming one result file per comparison, one run index, and one combined workbook.
- Permanent regression coverage proving the legacy VEGFRi/Dox-style exploratory path still runs and remains explicitly non-inferential.

### Synthetic fixtures

The test suite should use small synthetic fixtures rather than real collaborator data:

- `tests/fixtures/protocol/mock_protocol.xlsx`
- `tests/fixtures/licor/legacy_valid/`
- `tests/fixtures/licor/legacy_broken/`
- `tests/fixtures/licor/replicate_valid/`
- `tests/fixtures/manifests/`

The replicate-aware fixture should include at least:
- one male-only effect analyte
- one female-only effect analyte
- one null analyte
- one low-signal-but-changed analyte
- one analyte that is untestable because of insufficient data

### Silent-failure modes the test suite must catch

- analyte mislabeling caused by protocol extraction or LI-COR ordering drift
- accidental collapse of biological replicates into one group value
- leakage of pairwise filtering across unrelated treatment arms
- low-signal analytes being silently excluded from the inferential family
- multiplicity correction being applied across the wrong comparison family
- inferential outputs mislabeled as exploratory or exploratory outputs mislabeled as significant
- writes landing in the wrong user/analysis subtree
- row-order dependence in import or inferential results
- underpowered comparisons returning numeric p-values without explicit warning/status
- collaborator-facing plots rendering with clipped titles, truncated annotations, or literal placeholder text such as `NA`

## Final Multi-Lane Review And Release Gate

The automated test plan is necessary but not sufficient. Before this feature is accepted for collaborator use, the implementation revision must go through a final Research Partner review that inspects the actual code, fixtures, test evidence, and emitted outputs.

### Review workflow

1. Freeze the review target.
   Record the exact revision under review, the active config shape, the fixture set, and the scripts/spec artifacts being audited.
2. Run a Research Partner preflight.
   Reconstruct the code path, runtime/output roots, and artifact locations before specialist review begins.
3. Run implementation-accuracy review.
   Use an implementation-auditor lane to verify that code behavior matches the OpenSpec contract before judging whether the chosen method is statistically appropriate.
4. Run statistical-validity review.
   Use a stats-reviewer lane to verify replicate preservation, within-stratum contrasts, minimum-replicate handling, multiplicity correction, and low-signal semantics.
5. Run documentation-consistency review.
   Use a documentation-wizard lane to verify that `README.md`, methods text, config instructions, and output descriptions match the implemented interface and behavior.
6. Run robustness-and-testing review.
   Use a robustness-test-designer lane to verify that `testthat`, `pytest`, fixtures, smoke tests, and silent-failure coverage are sufficient and that the claimed release gate is evidence-backed.
7. Run optional scientific-interpretation review when needed.
   If manuscript-facing claims or biological interpretation are in scope, add a scientific-reviewer lane so claim boundaries are checked separately from implementation and statistics.
8. Synthesize findings.
   Use a review-synthesizer lane to combine the specialist outputs into one ranked review with blocking issues, non-blocking follow-ups, and a final acceptance decision.

### Required review artifacts

The final review must inspect, at minimum:

- `scripts/config/analysis_config.R`
- `scripts/main.R`, `scripts/find_ref_thresh.R`, `scripts/select-analytes-analysis.R`
- `scripts/helpers/project_paths.R`, `scripts/helpers/array_helper_scripts.R`, and any new replicate-aware statistical helpers
- `README.md`
- OpenSpec design, task, and requirement files for this change
- manifest templates and manifest-validation tests
- `testthat` inventory and latest results
- `pytest` inventory and latest results for the Python extractor
- synthetic fixture inputs for valid legacy runs, broken legacy runs, valid replicate-aware runs, and underpowered replicate-aware runs
- fixture-backed outputs from one legacy exploratory smoke run and one replicate-aware smoke run, including per-comparison result files, run index, combined workbook, and exploratory threshold artifacts
- a release review bundle containing one short lane report each for implementation, statistics, documentation, and robustness/testing, plus the synthesized final review

### Required pass/fail criteria

Release is blocked unless all of the following are true:

- no unresolved P0 or P1 findings remain in the synthesized final review
- `testthat` and `pytest` both pass on the same revision being reviewed
- fixture-backed tests cover all P0 categories and all explicitly named silent-failure modes in this design
- the replicate-aware path preserves one row per `analyte x sample_id`, keeps biological replicates separate, and uses manifest metadata rather than filenames as the source of truth
- within-stratum inferential outputs match the configured comparison set exactly, with no unintended subgroup contrasts
- comparisons with fewer than 2 replicates per arm are marked not testable, and `n = 2` per arm emits an explicit low-replication warning
- adjusted p-values are computed within the correct comparison family, default to `BH`, and inferential significance is derived from adjusted p-value plus alpha rather than fold-change thresholds
- low-signal analytes are flagged separately from inferential significance and are not silently excluded from testing by default
- exploratory and inferential outputs are clearly separated in naming and output routing, and inferential outputs include one canonical file per comparison plus a run-level index
- the legacy one-workbook-per-group exploratory path still runs on fixtures and remains explicitly non-inferential
- `README.md`, config docs, manifest docs, methods templates, and output descriptions match the implemented interface exactly
- collaborator-facing example plots, especially shortlist/waterfall outputs, render without clipped titles and without literal placeholder text such as `NA` in captions or annotations

The acceptance decision must be evidence-based. "The scripts ran once" or "the output looks plausible" is not enough for sign-off.

## Migration Plan

1. Add a manifest-driven input path alongside the legacy one-workbook-per-group exploratory path.
2. Refactor dataset assembly to produce a sample-level long-format table before any group collapsing.
3. Add a statistical analysis layer that emits per-comparison results tables with multiplicity correction.
4. Update output paths and labels so inferential and exploratory products are separated.
5. Update README and config examples to show both the original exploratory example and the new replicate-aware collaborator workflow.
6. Keep legacy exploratory scripts runnable during transition and after rollout, but prevent them from silently consuming replicate-aware datasets in an invalid way.

Rollback strategy:
- If the replicate-aware path is not ready, retain the legacy exploratory scripts without exposing inferential output claims.
