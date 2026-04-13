## Why

The current proteome profiler workflow is only valid for exploratory, one-workbook-per-group analyses. It averages duplicate spots within a membrane, collapses each treatment arm to a single value, and then labels fold-change threshold hits as "significant" without any biological replicate model or multiplicity correction.

This needs to change now because the next collaborator-facing use case is a factorial mouse-plasma study with biological replicates (`sex x treatment x replicate`). That study requires within-sex vehicle-versus-aldosterone comparisons, explicit replicate-aware statistics, and a defensible list of significant analytes after multiple-testing correction.

## What Changes

- Add support for sample-level biological replicate metadata instead of assuming one workbook per analysis group.
- Add replicate-aware analyte summaries that preserve per-sample values rather than collapsing each group to one number before testing.
- Add within-stratum differential testing so analyses can compare treatment arms inside a user-defined subgroup such as sex.
- Add multiplicity-corrected hit calling with a configurable adjustment method and `BH` as the default.
- Add explicit significant-hit tables as analysis outputs, separate from exploratory plots.
- Preserve the current one-workbook-per-group exploratory workflow for analyses that do not have biological replicates.
- Rename or relabel current threshold-based "significant" outputs so they are not conflated with statistical significance.
- Reject ambiguous multi-workbook layouts only when users attempt replicate-aware inference without the required sample metadata.

## Capabilities

### New Capabilities
- `replicate-aware-inputs`: ingest LI-COR workbooks as sample-level observations with explicit metadata for biological replicate, treatment, and optional stratification variables such as sex.
- `within-stratum-differential-analysis`: compute analyte-level treatment-versus-control comparisons within a configured subgroup, preserving biological replicate variance.
- `multiplicity-corrected-hit-reporting`: produce per-comparison results tables with raw p-values, adjusted p-values, configured adjustment method, effect estimates, and a clear significant-hit list.

### Modified Capabilities
- None.

## Impact

- Affected code: [scripts/helpers/array_helper_scripts.R](/Users/ncamarda/Projects/proteome_profiler/scripts/helpers/array_helper_scripts.R), [scripts/main.R](/Users/ncamarda/Projects/proteome_profiler/scripts/main.R), [scripts/find_ref_thresh.R](/Users/ncamarda/Projects/proteome_profiler/scripts/find_ref_thresh.R), [scripts/select-analytes-analysis.R](/Users/ncamarda/Projects/proteome_profiler/scripts/select-analytes-analysis.R), [scripts/helpers/project_paths.R](/Users/ncamarda/Projects/proteome_profiler/scripts/helpers/project_paths.R), [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R), [scripts/install_packages.R](/Users/ncamarda/Projects/proteome_profiler/scripts/install_packages.R), and [README.md](/Users/ncamarda/Projects/proteome_profiler/README.md).
- New outputs: replicate-aware analysis tables, multiplicity-corrected hit tables, and clarified output labels separating exploratory fold-change screens from inferential results.
- New configuration surface: sample metadata, subgroup definitions, statistical test selection, alpha threshold, and p-value adjustment method, while preserving the current exploratory config path for non-replicate analyses.
- Dependencies/systems: likely continued use of base R or tidyverse-compatible statistical tooling; no external service impact.
