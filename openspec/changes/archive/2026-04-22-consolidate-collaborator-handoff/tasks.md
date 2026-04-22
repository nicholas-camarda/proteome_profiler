## 1. Baseline And Output Contract

- [x] 1.1 Capture the current replicate-aware output tree for a fixture run and identify any duplicate or stale summary artifacts.
- [x] 1.2 Add tests that define the expected top-level inferential artifacts: `comparison_workbook.xlsx`, method result workbooks, `methods_overview.md`, and `run_index.tsv`.
- [x] 1.3 Add tests that reject creation of an additional standalone run-summary workbook or report directory.

## 2. Consolidated Comparison Workbook

- [x] 2.1 Add workbook summary builders for `summary`, `input_qc_summary`, `method_summary`, and `significance_summary`.
- [x] 2.2 Write summary sheets before comparison result sheets in `inferential_results/comparison_workbook.xlsx`.
- [x] 2.3 Ensure summary values are derived from the executed result objects and corrected multi-method run index.
- [x] 2.4 Preserve method-specific result workbooks for detailed per-method results.
- [x] 2.5 Add workbook tests for sheet presence, sheet order, hit counts, tested-analyte counts, comparison slugs, and method names.

## 3. Selected-Analyte Summary Integration

- [x] 3.1 Keep selected-analyte inputs optional for setup validation and main inferential analysis.
- [x] 3.2 Add or refresh a selected-analyte summary sheet only when selected-analyte analysis is configured and executed.
- [x] 3.3 Include selected analyte names, comparison slugs, methods, plot paths or plot availability, QC status, and no-plot reasons where applicable.
- [x] 3.4 Preserve detailed selected-analyte outputs under their existing comparison-scoped folders.
- [x] 3.5 Add tests for runs with no selected analytes and runs with selected-analyte follow-up enabled.

## 4. Documentation And Handoff Flow

- [x] 4.1 Update `README.md` so it is the primary start-here guide with copyable clone, setup, validation, analysis, and output-review commands.
- [x] 4.2 Explain that `comparison_workbook.xlsx` is the first replicate-aware result to open and `run_index.tsv` is for exact artifact paths/provenance.
- [x] 4.3 Document selected analytes as an optional follow-up workflow, not a required primary-analysis setting.
- [x] 4.4 Update `docs/README.md` to show current method-specific 25-per-page barplots, method-specific waterfall plots, and selected-analyte outputs.
- [x] 4.5 Keep documentation present-tense and remove historical or obsolete-output framing from user-facing docs.

## 5. Replicate-Analysis Maintainability

- [x] 5.1 Inventory every function in `scripts/helpers/replicate_analysis.R`, identify its callers, and assign exactly one target owner module.
- [x] 5.2 Add baseline parity tests before moving functions, covering representative raw-log2 and normalized-t-test p-values, fold changes, FDR flags, run-index rows, workbook sheets, selected-analyte outputs, and plot path creation.
- [x] 5.3 Add concise roxygen-style docstrings and dependency notes for exported helpers and functions called directly by entry scripts or tests while functions still live in `replicate_analysis.R`.
- [x] 5.4 Create `scripts/helpers/replicate_method_specs.R` for method definitions, p-adjust validation, safe math helpers, and fold-change/SE transformation helpers.
- [x] 5.5 Create `scripts/helpers/replicate_inputs.R` for manifest parsing, workbook path resolution, sheet reading, coordinate cleanup, sample-level dataset construction, and input cleanliness summaries.
- [x] 5.6 Create `scripts/helpers/replicate_models.R` for comparison slugs, comparison validation, method-specific model fitting, multi-method execution, low-signal flags, and p-value/FDR result columns.
- [x] 5.7 Create `scripts/helpers/replicate_plotting.R` for waterfall plots, barplot data builders, pagination, selected-analyte plot renderers, common y-axis handling, star/bracket annotations, and SE whisker plotting.
- [x] 5.8 Create `scripts/helpers/replicate_outputs.R` for output cleanup, run-index construction, workbook writing, comparison-scoped tables, and methods overview writing.
- [x] 5.9 Create `scripts/helpers/replicate_selected_analytes.R` for selected-analyte config parsing, name validation and suggestions, comparison selection, selected-result extraction, selected QC summaries, and selected-output orchestration.
- [x] 5.10 Convert `scripts/helpers/replicate_analysis.R` into the public loader that sources the helper modules in dependency order and contains no analysis logic.
- [x] 5.11 After each module move or small batch of moves, run the parity tests and fix source-order or ownership issues before moving the next batch.
- [x] 5.12 Remove duplicated or obsolete helper paths rather than preserving fallback implementations.
- [x] 5.13 Record per-sample reference-spot normalization denominators and write reference-spot QC outputs under `input_qc/`.
- [x] 5.14 Add setup and main-run messages that summarize complete preferred, partial preferred, and protocol-table fallback reference-spot usage.

## 6. Verification

- [x] 6.1 Run the R test suite for workbook, inferential-output, selected-analyte, and script smoke tests.
- [x] 6.2 Run the Python protocol extraction tests if touched by documentation or setup checks.
- [x] 6.3 Run `Rscript scripts/check_setup.R` and at least one fixture or collaborator-safe demo analysis.
- [x] 6.4 Inspect `comparison_workbook.xlsx` for expected sheet names, sheet order, non-empty summary tabs, method/comparison labels, hit counts, QC totals, and selected-analyte summary content when applicable.
- [x] 6.5 Inspect representative generated plots, including full waterfall plots, significant-hit waterfall plots, 25-per-page barplots, significant-hit barplots, and selected-analyte bargraphs/waterfalls.
- [x] 6.6 Confirm inspected plots have correct titles, legends, control-first ordering, method labels, y-axis ranges, SE whiskers where applicable, bracket/star annotations where applicable, blank panel behavior, and no clipped labels or annotations.
- [x] 6.7 Record the inspected artifact paths and any issues found in the implementation notes or final handoff response.
- [x] 6.8 Run `openspec validate consolidate-collaborator-handoff` and fix any spec validation issues.
