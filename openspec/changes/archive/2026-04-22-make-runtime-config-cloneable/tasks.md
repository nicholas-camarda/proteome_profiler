## 1. Environment Loading

- [x] 1.1 Add `dotenv` to the R dependency list and package installer.
- [x] 1.2 Call the `.env` loader before config sourcing in `scripts/main.R`, `scripts/find_ref_thresh.R`, and `scripts/select-analytes-analysis.R`.
- [x] 1.3 Add tests proving explicit process environment variables override `.env` values.

## 2. Single-File Run Configuration

- [x] 2.1 Add `.env.example` with commented sections for runtime paths, analysis identity, protocol inputs, manifest/data inputs, comparisons, thresholds, methods, and optional selected analytes.
- [x] 2.2 Implement `.env` parsing that converts delimiter-based values into the existing internal analysis-config shape.
- [x] 2.3 Define and validate simple `.env` conventions for lists, numeric vectors, comparisons, and optional selected analytes.
- [x] 2.4 Remove `scripts/config/analysis_config.R` as a supported run-configuration source.
- [x] 2.5 Remove `PROTEOME_PROFILER_CONFIG`-based run configuration from entry scripts and tests.
- [x] 2.6 Remove collaborator-specific active defaults from tracked runtime files.
- [x] 2.7 Add tests proving a fresh clone without private defaults does not select Nicole or any collaborator analysis implicitly.
- [x] 2.8 Add tests proving routine analysis configuration can be parsed from `.env` without editing or creating R scripts.
- [x] 2.9 Add tests proving omitted selected-analyte fields do not block setup validation or main analysis.
- [x] 2.10 Add tests proving selected-analyte follow-up stops clearly when selected analytes are omitted.
- [x] 2.11 Document and validate exact selected-analyte comparison slug syntax, including replicate-aware `<subgroup value>_<control label>_vs_<treatment label>` slugs and `|`-delimited multiple slugs.

## 3. Path Resolution

- [x] 3.1 Update path resolution to always search repo and runtime roots and search cloud roots only when configured.
- [x] 3.2 Ensure `PROTEOME_PROFILER_RUNTIME_ROOT` controls new output paths.
- [x] 3.3 Document and validate the recommended input layout: `manifests/`, `workbooks/`, and `protocols/`.
- [x] 3.4 Add tests for path resolution with no cloud root, with a cloud root, absolute workbook paths, relative workbook paths, and missing required paths.

## 4. Manifest Workbook Inputs

- [x] 4.1 Preserve support for one workbook per biological sample using `sample_id` and `workbook_path`.
- [x] 4.2 Preserve support for one workbook with multiple sample sheets using `sample_id`, `workbook_path`, and `sheet_name`.
- [x] 4.3 Ensure only manifest-listed sheets are read from multi-sheet workbooks.
- [x] 4.4 Add setup validation that fails clearly when a listed workbook or listed sheet is missing.
- [x] 4.5 Add README examples for both one-workbook-per-sample manifests and one-workbook-with-multiple-sheets manifests.

## 5. Setup Validation And Dependencies

- [x] 5.1 Add `scripts/check_setup.R` that loads `.env`, parses the selected run configuration, validates config shape, and prints resolved paths.
- [x] 5.2 Make `load_analysis_packages()` stop with actionable missing-package messages.
- [x] 5.3 Add Python dependency metadata or documentation for protocol extraction dependencies.
- [x] 5.4 Make setup validation report missing R packages, missing Python dependencies, missing Java/tabula requirements, missing protocol files, invalid manifests, missing workbooks, missing workbook sheets, and unwritable output roots.
- [x] 5.5 Add tests covering successful setup validation and major setup failure modes.
- [x] 5.6 Add Java setup guidance that gives macOS users a single-command install path when Java is unavailable for `tabula-py`.

## 6. Demo Or Smoke Workflow

- [x] 6.1 Provide a cloneable demo or generated fixture workflow that does not require private OneDrive data.
- [x] 6.2 Ensure the demo exercises env loading, config loading, path resolution, setup validation, and at least one analysis script.
- [x] 6.3 Add a smoke test for the documented demo command sequence.

## 7. Documentation

- [x] 7.1 Rewrite README quick start so the primary flow is copy `.env.example`, edit `.env`, run setup validation, then run the analysis scripts.
- [x] 7.2 Document every `.env` field with examples for legacy and replicate-aware analyses.
- [x] 7.3 State clearly that routine users do not edit, copy, or create R scripts.
- [x] 7.4 Update docs references so example outputs and setup commands match the runtime workflow.
- [x] 7.5 Write user-facing docs in present-state language only; do not reference removed config workflows, removed artifact layouts, or "now supported" language.
- [x] 7.6 Add beginner-level README setup instructions with copyable commands for opening Terminal, cloning `https://github.com/nicholas-camarda/proteome_profiler.git`, entering the repo, creating `.env`, installing dependencies, checking setup, and running the analysis sequence.
- [x] 7.7 Document replicate-aware waterfall and barplot SE whiskers, including ratio-scale barplot whiskers and annotation-safe y-axis behavior.

## 8. Replicate-Aware Plot Uncertainty

- [x] 8.1 Require replicate-aware waterfall plots to include finite method-specific `effect_se_log2` values and draw `+/- 1 SE` whiskers for all supported methods.
- [x] 8.2 Add replicate-aware barplot SE whiskers for treatment fold-change bars using ratio-scale bounds from `effect_estimate_log2 +/- effect_se_log2`.
- [x] 8.3 Ensure common y-axis limits include SE upper whiskers plus annotation headroom, and disable clipping for labels/brackets.
- [x] 8.4 Ensure replicate-aware selected-analyte bargraphs reuse the same SE-aware renderer and common y-axis logic.
- [x] 8.5 Add tests covering barplot SE bounds and y-axis headroom.
- [x] 8.6 Ensure multi-method `run_index.tsv` includes every configured method and all generated method-specific output paths.
- [x] 8.7 Ensure selected-analyte comparison resolution de-duplicates multi-method run-index rows before writing comparison folders.

## 9. Verification

- [x] 9.1 Run R unit tests for `.env` config/path/setup behavior.
- [x] 9.2 Run script smoke tests for legacy and replicate-aware workflows.
- [x] 9.3 Run Python protocol extraction tests.
- [x] 9.4 Run targeted replicate-aware plot and inferential output tests after adding SE whiskers.
- [x] 9.5 Rerun Nicole's configured workflow end-to-end and audit generated run indexes, result tables, plot paths, selected-analyte outputs, and significance summaries.
- [x] 9.6 Run the full R test suite and Python protocol extraction tests after the live audit fixes.
- [x] 9.7 Run `openspec status --change make-runtime-config-cloneable` and confirm the change is apply-ready.
