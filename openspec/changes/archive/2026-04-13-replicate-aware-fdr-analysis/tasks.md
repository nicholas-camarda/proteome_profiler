## 1. Test Harness And Fixtures

- [x] 1.1 Add `testthat` infrastructure for R helpers and script-level regression tests
- [x] 1.2 Add `pytest` infrastructure for `scripts/setup/extract_analyte_table.py`
- [x] 1.3 Add synthetic fixture files for protocol mapping, valid legacy LI-COR inputs, broken legacy LI-COR inputs, replicate-aware LI-COR inputs, and manifests
- [x] 1.4 Add regression coverage for the current one-workbook-per-group exploratory path so it remains runnable throughout the refactor
- [x] 1.5 Standardize `scripts/install_packages.R` on `pak`, including bootstrap when `pak` is not yet installed

## 2. Replicate-Aware Input Model

- [x] 2.1 Move user-editable analysis metadata into a dedicated config script such as `scripts/config/analysis_config.R`
- [x] 2.2 Keep `scripts/helpers/project_paths.R` focused on path resolution and validation rather than user-edited metadata
- [x] 2.3 Add a sample-manifest configuration path and validation rules to the user-editable config contract
- [x] 2.4 Add a documented manifest template with required columns `sample_id`, `workbook_path`, `treatment`, plus optional subgroup/covariate columns
- [x] 2.5 Refactor LI-COR import code to build a sample-level long-format dataset with `sample_id`, treatment, subgroup metadata, and analyte signal
- [x] 2.6 Keep technical duplicate spot averaging within each workbook while preserving biological replicates as separate observations
- [x] 2.7 Add fail-fast validation for missing manifest columns, duplicate sample IDs, missing workbook files, and ambiguous legacy multi-file group layouts
- [x] 2.8 Preserve the existing one-workbook-per-group exploratory import path for analyses without biological replicates

## 3. Differential-Analysis Engine

- [x] 3.1 Introduce a statistical-analysis helper that runs per-analyte within-stratum treatment-versus-control models on log2-transformed signals
- [x] 3.2 Add configuration for subgroup variable, contrast definitions, minimum replicate requirements, and alpha threshold
- [x] 3.3 Emit analyte-level effect estimates, sample counts, raw p-values, and test-status fields for each requested comparison
- [x] 3.4 Add explicit handling for untestable analytes or underpowered comparisons so the pipeline reports why inference was skipped

## 4. Multiplicity Correction And Significant-Hit Outputs

- [x] 4.1 Add configurable p-value adjustment with `BH` as the default and validation of supported methods
- [x] 4.2 Apply multiplicity correction separately within each comparison family and compute inferential significance from adjusted p-values
- [x] 4.3 Write machine-readable per-comparison results tables that include raw p-values, adjusted p-values, adjustment method, alpha, and significant-hit flags
- [x] 4.4 Rename or separate exploratory threshold-hit outputs so they cannot be mistaken for inferential significance

## 5. Script And Output Refactor

- [x] 5.1 Update `scripts/main.R` to dispatch explicitly between manifest-driven replicate-aware mode and legacy exploratory mode
- [x] 5.2 Keep `scripts/main.R` able to run the legacy exploratory path when only one workbook per group is provided
- [x] 5.3 Update `scripts/find_ref_thresh.R` so its role in the replicate-aware workflow is explicit and does not imply inferential significance
- [x] 5.4 Update `scripts/select-analytes-analysis.R` or replace it with a replicate-aware shortlist workflow that operates on inferential results where appropriate
- [x] 5.5 Reorganize output directories so exploratory plots, inferential tables, and shortlist outputs are separated and labeled by user, analysis, and comparison
- [x] 5.6 Fix collaborator-facing plot rendering so long titles/subtitles are not clipped and optional annotations do not print placeholder text such as `NA`

## 6. Verification And Documentation

- [x] 6.1 Add unit tests for protocol extraction, path/config helpers, and legacy LI-COR import invariants
- [x] 6.2 Add regression tests for `find_filter_thresh()`, `make_wf_data()`, and `make_graphs()` to catch silent semantic drift
- [x] 6.3 Add manifest-validation and mode-dispatch tests for the replicate-aware path
- [x] 6.4 Add integration tests for sample-level dataset assembly, within-stratum modeling, low-replication warnings, and multiplicity correction
- [x] 6.5 Add output-structure tests for per-user roots, exploratory vs inferential separation, per-comparison result files, run index, and combined workbook export
- [x] 6.6 Add smoke tests for `find_ref_thresh.R`, `main.R`, and `select-analytes-analysis.R` on tiny fixtures, checking file presence and TSV payloads instead of PNG pixels
- [x] 6.6.1 Add a shortlist-waterfall rendering check that verifies title/caption text is fully populated and does not contain literal `NA`
- [x] 6.7 Rewrite the README to document exact script order, manifest format, mode dispatch, replicate-aware statistics, and output semantics
- [x] 6.8 Add inline docstrings/comments for the new statistical and manifest-handling helpers so users can inspect the code quickly

## 7. Final Multi-Lane Review Gate

- [x] 7.1 Add a release-review bundle format that records the exact revision, fixture set, test results, and runtime artifacts under review
- [x] 7.2 Run a Research Partner preflight to reconstruct path truth, entrypoints, and artifact locations before final acceptance
- [x] 7.3 Run an implementation-accuracy review lane against the implementation revision and resolve or document all blocking findings
- [x] 7.4 Run a statistical-validity review lane against the implementation revision and resolve or document all blocking findings
- [x] 7.5 Run a documentation-consistency review lane against `README.md`, methods text, config docs, and output semantics
- [x] 7.6 Run a robustness-and-testing review lane against the final `testthat`/`pytest` suite, fixtures, smoke runs, and silent-failure coverage
- [x] 7.7 Run an optional scientific-interpretation review lane when manuscript-facing claims or biological interpretation are in scope
- [x] 7.8 Run a synthesized final review that merges lane outputs, ranks findings, and records blocking issues versus non-blocking follow-ups
- [x] 7.9 Treat unresolved P0/P1 findings, missing lane reports, or missing evidence artifacts as release blockers
