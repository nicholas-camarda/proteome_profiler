## ADDED Requirements

### Requirement: Comparison workbook SHALL be the primary collaborator-facing replicate-aware result
Replicate-aware inferential analysis SHALL write `inferential_results/comparison_workbook.xlsx` as the primary human-readable workbook for reviewing configured comparisons and methods.

#### Scenario: Workbook includes run summary tabs
- **WHEN** replicate-aware inferential analysis completes
- **THEN** `comparison_workbook.xlsx` SHALL contain `summary`, `input_qc_summary`, `method_summary`, and `significance_summary` sheets before comparison result sheets

#### Scenario: Summary tabs are derived from executed results
- **WHEN** summary sheets are written
- **THEN** their counts, comparison slugs, method names, significance totals, and QC totals SHALL be derived from the same in-memory results and run index used to write method outputs

#### Scenario: Workbook remains the only collaborator-facing summary workbook
- **WHEN** replicate-aware inferential analysis writes its output tree
- **THEN** it SHALL NOT create an additional standalone run-summary workbook or report directory that duplicates the comparison workbook summary role

### Requirement: Output artifacts SHALL have distinct roles
The output tree SHALL preserve one clear purpose for each major artifact so users can navigate results without duplicate summary surfaces.

#### Scenario: Machine-readable provenance remains separate
- **WHEN** replicate-aware inferential analysis completes
- **THEN** `inferential_results/run_index.tsv` SHALL contain one row per comparison and method with generated output paths and SHALL be treated as machine-readable provenance rather than the primary collaborator-facing results table

#### Scenario: Method workbooks remain detailed result artifacts
- **WHEN** multiple primary methods are run
- **THEN** method-specific result workbooks SHALL remain available for full per-method tables and SHALL NOT be replaced by summary-only workbook sheets

#### Scenario: Methods overview remains prose only
- **WHEN** `inferential_results/methods_overview.md` is written
- **THEN** it SHALL explain estimands, fold-change definitions, raw p-value and BH interpretation, and major caveats without duplicating full result tables

#### Scenario: Comparison folders contain plots and comparison-scoped tables
- **WHEN** comparison-scoped outputs are written
- **THEN** `inferential_results/comparisons/<comparison_slug>/` SHALL contain method-organized plots and comparison-scoped tables, not a second copy of top-level run summaries

### Requirement: Generated workbooks and plots SHALL be inspected for content validity
Verification SHALL include inspection of representative generated results and plots, not only checks that expected files exist.

#### Scenario: Workbook content is inspected
- **WHEN** replicate-aware inferential analysis is regenerated for a fixture or collaborator-safe analysis
- **THEN** verification SHALL inspect `comparison_workbook.xlsx` sheet names, sheet order, non-empty summary tables, method/comparison labels, hit counts, QC totals, and selected-analyte summary content when selected analytes are configured

#### Scenario: Plot content is inspected
- **WHEN** representative waterfall, 25-per-page barplot, significant-hit barplot, and selected-analyte plot files are regenerated
- **THEN** verification SHALL inspect rendered plot dimensions, titles, legends, control-first ordering, explicit linear fold-change-ratio y-axis labeling, y-axis ranges, group-level SE whiskers on both control and treatment bars where applicable, bracket/star annotations where applicable, blank panel behavior, and absence of clipped labels or annotations

#### Scenario: Barplot SE construction is documented
- **WHEN** replicate-aware barplots or method notes are generated
- **THEN** documentation SHALL state that barplot whiskers are visual `+/- 1 SE` bars for each plotted group mean on the linear fold-change-ratio scale, not confidence intervals for the treatment-vs-control contrast
- **AND** documentation SHALL define `normalized_signal` as the per-sample, per-analyte averaged duplicate raw signal divided by that sample's reference-spot denominator, using preferred Reference Spots pairs `A1,2` and `J1,2` when present and otherwise available Reference Spots rows from that sample
- **AND** documentation SHALL state that `normalized_t_test` barplot whiskers use `normalized_signal` group SE divided by the control-arm `normalized_signal` mean
- **AND** documentation SHALL state that `raw_log2_lm` barplot whiskers use log2-scale group-mean SE bounds converted back to the linear ratio scale and therefore can appear asymmetric
- **AND** documentation SHALL distinguish barplot group-mean SE whiskers from waterfall treatment-effect SE whiskers
- **AND** documentation SHALL state that `low_signal_flag` is computed from the configured raw-signal low-signal threshold for every replicate-aware method and is separate from the `normalized_t_test` reference-spot denominator

#### Scenario: Reference-spot normalization is auditable
- **WHEN** sample-level replicate-aware input QC is generated
- **THEN** `input_qc/reference_spot_qc.tsv` SHALL contain one row per sample with the reference-spot denominator, reference row source, reference QC status, reference spot identifiers, protocol names, reference coordinates, and raw reference signals used to compute `normalized_signal`
- **AND** `input_qc/reference_spot_summary.tsv` SHALL summarize sample counts using complete preferred `A1,2`/`J1,2` reference spots, partial preferred reference spots, and protocol-table fallback reference rows
- **AND** setup checks SHALL parse the configured sample sheets and report the same reference-spot QC counts before the analysis run
- **AND** a missing, non-finite, or nonpositive reference denominator SHALL fail clearly before inferential results are written

#### Scenario: Inspection findings block completion
- **WHEN** inspection identifies stale tabs, wrong plot families, duplicated summaries, missing SE whiskers, clipped annotations, unreadable facet titles, or mismatched method labels
- **THEN** the change SHALL NOT be considered complete until the generated artifacts are corrected and regenerated

### Requirement: Selected-analyte summaries SHALL be conditional and comparison-scoped
Selected-analyte outputs SHALL remain optional follow-up artifacts and SHALL add summary content only when selected-analyte analysis is configured and executed.

#### Scenario: No selected analytes configured
- **WHEN** main inferential analysis runs without selected analytes configured
- **THEN** the analysis SHALL complete without selected-analyte summary sheets or selected-analyte output folders

#### Scenario: Selected-analyte follow-up writes summaries
- **WHEN** selected-analyte analysis is configured and executed
- **THEN** the workflow SHALL write or refresh a selected-analyte summary sheet in `comparison_workbook.xlsx` with selected analyte names, comparison slugs, methods, plot availability, and QC status

#### Scenario: Selected-analyte artifacts stay in selected-analyte tree
- **WHEN** selected-analyte plots and tables are generated
- **THEN** the detailed selected-analyte artifacts SHALL remain under `select_analytes/<comparison_slug>/` for exploratory mode or `select_analytes/<comparison_slug>/<method>/` for replicate-aware mode

### Requirement: Documentation SHALL describe the current handoff path
User-facing documentation SHALL explain the current clone-and-run workflow and output review path directly, without historical framing or references to obsolete plot types as the recommended output.

#### Scenario: README points to the primary output path
- **WHEN** a user reads the README after setup and analysis commands
- **THEN** the README SHALL direct them to open `inferential_results/comparison_workbook.xlsx` first for replicate-aware results and to use `run_index.tsv` only when they need exact artifact paths

#### Scenario: Example output docs show current plot families
- **WHEN** docs show example outputs for replicate-aware analyses
- **THEN** the examples SHALL use current method-specific 25-per-page barplots, method-specific waterfall plots, and selected-analyte plots rather than exploratory-only plots as the main examples

#### Scenario: No historical migration wording
- **WHEN** README or docs are updated for this change
- **THEN** they SHALL describe the supported behavior in present tense and SHALL NOT use migration framing such as "now supported" unless explicitly requested

### Requirement: Replicate-analysis helpers SHALL be auditable by responsibility
The replicate-aware implementation SHALL be organized and documented so reviewers can identify input/QC, modeling, plotting, output writing, and selected-analyte responsibilities without tracing unrelated code paths.

#### Scenario: Exported helpers have docstrings
- **WHEN** a helper function is intended to be called by an entry script or test
- **THEN** it SHALL have a concise roxygen-style docstring describing inputs, outputs, and statistical or plotting assumptions when relevant

#### Scenario: Responsibilities are separated
- **WHEN** replicate-aware helper code is refactored
- **THEN** method specs/math, input/QC, model fitting, plotting, workbook/path output, and selected-analyte logic SHALL be separated into clearly named helper files without duplicated alternate implementations

#### Scenario: Public source entrypoint remains stable
- **WHEN** entry scripts or tests need replicate-aware helpers
- **THEN** they SHALL source `scripts/helpers/replicate_analysis.R`, which SHALL load the owned helper modules in deterministic dependency order and SHALL NOT contain analysis logic itself

#### Scenario: Module ownership is explicit
- **WHEN** functions are moved out of `replicate_analysis.R`
- **THEN** each function SHALL have exactly one owning helper module and shared helpers SHALL be placed in the earliest module needed by the dependency order rather than copied

#### Scenario: Refactor preserves numerical results
- **WHEN** the helper organization changes without an intended statistical-method change
- **THEN** tests SHALL verify representative p-values, fold changes, FDR flags, run-index rows, workbook sheets, and selected-analyte outputs remain unchanged for fixture analyses

#### Scenario: Refactor is staged separately from behavior changes
- **WHEN** helper functions are moved between files
- **THEN** function bodies SHALL NOT change statistical calculations, filtering rules, plotting semantics, output paths, or workbook contents in the same step
