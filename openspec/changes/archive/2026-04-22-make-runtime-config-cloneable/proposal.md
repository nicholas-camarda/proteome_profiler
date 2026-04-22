## Why

The repository can run the current analyses, but a fresh clone still inherits machine-specific paths, collaborator-specific defaults, and split dependency assumptions that make it too fragile for future lab users. This change makes the runtime setup explicit, local-first, and fail-fast so a new user can configure a small `.env` file, run documented commands, and get clear setup errors before analysis begins.

## What Changes

- Add a first-class `.env` workflow that serves as the single collaborator-facing run configuration file, including machine paths, active analysis identity, input paths, comparisons, thresholds, analysis methods, and optional selected analytes.
- Use an R `.env` reader such as `dotenv` so environment loading is standard and readable rather than hand-rolled.
- Retire `scripts/config/analysis_config.R` as a supported run configuration source; the pipeline SHALL build run configuration from `.env`.
- Replace checked-in collaborator-specific defaults with clone-safe defaults and `.env.example` templates.
- Make cloud backup paths optional rather than required by path resolution.
- Add setup validation that checks selected analysis, config path, required R packages, Python protocol-extraction dependencies, protocol files, manifests, workbook paths, and output writability before long-running analysis starts.
- Add or document a cloneable demo/smoke path that can run without private OneDrive data.
- Make runtime package loading fail early with actionable missing-package messages.
- Add Python dependency documentation or package metadata for protocol workbook extraction.
- Add Java setup guidance for `tabula-py` so a missing Java runtime has a single-command macOS installation path.
- Document exact selected-analyte comparison slug syntax, including replicate-aware `subgroup_control_vs_treatment` slugs and `|`-delimited multiple slugs.
- Ensure replicate-aware waterfall and barplot outputs include method-specific standard-error whiskers, with barplot y-axis headroom that accounts for SE whiskers and annotations.
- Update README/docs so future users edit one `.env` file and then follow a command sequence, without creating a new R script for each analysis.

## Capabilities

### New Capabilities

- `cloneable-runtime-configuration`: Covers single-file `.env` run configuration, local-safe defaults, optional cloud roots, preflight setup validation, dependency checks, and a documented fresh-clone workflow.

### Modified Capabilities

- None.

## Impact

- Affected files include `.gitignore`, `.env.example`, README/docs, config parsing helpers, `scripts/helpers/project_paths.R`, `scripts/helpers/runtime_setup.R`, entry scripts, and new setup/check scripts.
- Runtime behavior changes to a single supported run-configuration source: one `.env` file defines the active run and internal R helpers parse/validate it.
- Scientific analysis logic, fold-change estimands, and inferential methods are not changed by this proposal.
- Replicate-aware plotting semantics include explicit `+/- 1 SE` whiskers for waterfall and barplot outputs when biological replicate methods are used.
