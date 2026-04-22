## Context

The repository has two audiences that need different interfaces to the same analysis: lab users who need a simple clone-and-run path, and reviewers/developers who need auditable method outputs and provenance. The current output tree already contains the major ingredients, but the collaborator-facing interpretation path is spread across method workbooks, comparison folders, plot folders, `run_index.tsv`, and prose notes. Adding a new report directory would make that worse.

The implementation also concentrates much of the replicate-aware workflow in `scripts/helpers/replicate_analysis.R`. That is workable for execution, but it is becoming hard to review because input validation, modeling, plotting, workbook writing, selected-analyte behavior, and path indexing are interleaved.

## Goals / Non-Goals

**Goals:**

- Make `inferential_results/comparison_workbook.xlsx` the first artifact a collaborator should open for replicate-aware results.
- Add concise summary tabs to that workbook so run status, input QC, method-level hit counts, and selected-analyte outputs are visible in one place.
- Keep method-specific full result files and plot folders available without duplicating their content into multiple summary artifacts.
- Preserve the organized output tree: workbook summaries at the inferential level, method-specific plots under comparison folders, selected-analyte follow-up under `select_analytes/`.
- Keep README as the start-here guide and avoid creating another user-facing setup guide unless the README becomes structurally unmanageable.
- Refactor replicate-analysis helper code through a staged, behavior-preserving module split with explicit ownership boundaries and parity checks.

**Non-Goals:**

- No GitHub Actions or CI setup.
- No new standalone report directory or second collaborator-facing workbook.
- No change to statistical estimands, p-value calculation, FDR calculation, low-signal filtering, plot thresholds, or selected-analyte inclusion rules.
- No backwards-compatibility branches or legacy rescue paths.

## Decisions

### Consolidate summaries into `comparison_workbook.xlsx`

`comparison_workbook.xlsx` should become the human-readable summary and comparison surface. It should include run-level summary tabs before comparison tables, then method/comparison result sheets. This avoids adding a new artifact while giving collaborators a clear opening point.

Alternative considered: add `run_report.xlsx` or a report folder. Rejected because it duplicates `comparison_workbook.xlsx`, makes artifact navigation harder, and violates the goal of reducing output sprawl.

### Keep `run_index.tsv` machine-readable

`run_index.tsv` should remain a complete path/provenance index with one row per comparison and method. It should not be treated as the collaborator-facing results table. Human-readable summaries should be derived from the same data and written into workbook tabs.

Alternative considered: remove `run_index.tsv` and rely on workbook sheets. Rejected because scripts and tests need a stable machine-readable output index.

### Keep prose method notes short

`methods_overview.md` should explain what the methods estimate, why raw-log2 and normalized t-test fold changes can differ, how BH is computed over usable tests, and what the significance flags mean. It should not become a full analysis report or duplicate workbook tables.

Alternative considered: expand the methods note into a narrative report. Rejected because the repo needs fewer user-facing artifacts, not another one.

### Make plot scales and SE bars explicit

Replicate-aware barplots should state that the y-axis is a linear fold-change ratio relative to the control group. Barplot whiskers should be described as visual `+/- 1 SE` bars for each plotted group mean on that method-specific linear ratio scale, not as confidence intervals for the treatment-vs-control contrast.

For `normalized_t_test`, the documentation should define `normalized_signal` before using it in formulas: one sample/analyte `normalized_signal` equals the averaged duplicate raw analyte signal divided by that same sample's reference-spot denominator. The denominator uses preferred Reference Spots pairs `A1,2` and `J1,2` when present and otherwise available Reference Spots rows from that sample. Control and treatment barplot whiskers are built from `normalized_signal` group standard errors divided by the control-arm `normalized_signal` mean. These whiskers should be symmetric on the plotted linear ratio scale unless a lower bound is clipped at zero. For `raw_log2_lm`, group SE bounds are computed on the log2 scale and converted back to the linear ratio scale, so the rendered whiskers can be asymmetric. Waterfall plots have a different role: they show the method-specific treatment effect estimate and use the treatment-effect SE.

Documentation should also distinguish the configured low-signal reference panel from membrane reference-spot normalization. `low_signal_flag` is computed from the raw-signal threshold and reported for every replicate-aware method. It is separate from the `normalized_t_test` denominator and does not remove analytes before p-value adjustment.

Alternative considered: label all whiskers generically as `+/- 1 SE`. Rejected because the same visual term refers to group-mean SE in barplots and treatment-effect SE in waterfalls, and the distinction matters for interpretation.

### Make reference-spot normalization auditable

The replicate-aware parser should record the exact reference rows and raw reference signals used to compute each sample's `normalized_signal` denominator. Setup checks should parse the configured sample sheets and report how many samples used complete preferred `A1,2`/`J1,2` raw reference spots, partial preferred reference spots, or protocol-table fallback reference rows. Main analysis should write `input_qc/reference_spot_qc.tsv` and `input_qc/reference_spot_summary.tsv` so the denominator source is reviewable without opening the raw workbook.

Alternative considered: only describe the reference-spot rule in documentation. Rejected because this would not detect workbook-layout, manifest-sheet, or protocol-join mistakes that silently change the denominator used by `normalized_t_test`.

### Treat selected-analyte summaries as conditional workbook content

Selected analytes are optional follow-up outputs. When configured and run, concise selected-analyte summary tabs should be included in the comparison workbook or refreshed by the selected-analyte workflow. When not configured, the main analysis should not create empty selected-analyte tabs.

Alternative considered: require selected analytes in `.env` for every run. Rejected because selected analytes are explicitly a post hoc plotting/review workflow, not required input for primary analysis.

### Split replicate-analysis helpers through a staged module boundary

`scripts/helpers/replicate_analysis.R` should remain the single public source entrypoint used by `scripts/main.R`, `scripts/find_ref_thresh.R`, `scripts/check_setup.R`, `scripts/select-analytes-analysis.R`, and tests. Its responsibility should shrink to loading the new helper modules in a deterministic order and defining no analysis logic itself.

The module split should be behavior-preserving and staged:

1. Add baseline parity tests and capture representative fixture outputs before moving functions.
2. Add docstrings and dependency notes to the existing functions while they still live in `replicate_analysis.R`.
3. Move functions into modules without changing function bodies except for comments, whitespace, and import/source ordering.
4. Run parity tests after each module move or small batch of module moves.
5. Inspect representative output workbooks and rendered plots after the refactor, not just file paths.
6. Only after parity and artifact inspection pass, make any strictly necessary call-site cleanup.

The target module boundaries are:

- `scripts/helpers/replicate_method_specs.R`: method definitions, p-adjust validation, safe math helpers, fold-change and SE transformation helpers.
- `scripts/helpers/replicate_inputs.R`: sample manifest parsing, workbook path resolution, raw sheet reading, coordinate cleanup, sample-level dataset construction, and input cleanliness summaries.
- `scripts/helpers/replicate_models.R`: comparison slug construction, comparison validation, method-specific model fitting, multi-method execution, low-signal flags, and p-value/FDR result columns.
- `scripts/helpers/replicate_plotting.R`: waterfall plots, barplot data builders, barplot pagination, selected-analyte plot renderers, common y-axis handling, star/bracket annotations, and SE whisker plotting.
- `scripts/helpers/replicate_outputs.R`: output-directory cleanup, run-index construction, workbook writing, comparison-scoped table writing, and methods overview writing.
- `scripts/helpers/replicate_selected_analytes.R`: selected-analyte config parsing, selected-name validation and suggestions, comparison selection, selected-result extraction, selected QC summaries, and selected-output orchestration.

The split must not create duplicate implementations or fallback paths. A function should have one owning module. If a helper is shared, it should move to the earliest appropriate module in the source order rather than being copied.

Alternative considered: leave `replicate_analysis.R` monolithic and only add tests. Rejected because the code is already difficult to inspect and future method changes will be riskier.

Alternative considered: move logic opportunistically while adding workbook features. Rejected because refactoring and behavior changes in the same step would make numerical discrepancies hard to diagnose.

## Risks / Trade-offs

- Workbook tabs can become too large if every method table is duplicated. Mitigation: summary tabs should aggregate only run-level and hit-count information; full method sheets remain the detailed tables.
- Selected-analyte workflow may run after the main inferential workbook is created. Mitigation: selected-analyte execution should refresh or append only the selected-analyte summary tabs without rewriting unrelated result content inconsistently.
- Refactoring `replicate_analysis.R` could accidentally change results. Mitigation: freeze representative fixture expectations first, move functions by module with no body-level behavior edits, and verify p-values, fold changes, FDR flags, run-index rows, workbook sheet names, selected-analyte outputs, and plot path creation after each stage.
- Module dependencies can become cyclic if plotting, output writing, and selected-analyte orchestration are split carelessly. Mitigation: use a one-way source order of method specs/math, inputs, models, plotting, outputs, selected-analyte orchestration; shared helpers move upward in that order rather than being duplicated.
- File-existence tests can miss bad handoff artifacts, including blank sheets, stale tabs, clipped labels, wrong plot families, or visually broken significance annotations. Mitigation: inspect representative generated workbooks and rendered PNGs/PDFs as part of verification, including workbook sheet content, plot dimensions, titles, legends, explicit linear fold-change-ratio y-axis labels, y-axis ranges, control and treatment SE whiskers, and selected-analyte outputs.
- README can become too long if it absorbs all setup detail. Mitigation: keep README as the linear quick start and use `docs/README.md` only for example output interpretation and visual examples.
