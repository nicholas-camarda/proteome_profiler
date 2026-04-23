# Docs

This folder holds repo-local output examples that can be embedded directly in the documentation.

## Embedded Output Examples

These are local rendered example outputs from the current workspace. They are example outputs, not defaults for a new project.

## Replicate-Aware Examples

These examples show the biological-replicate outputs produced under `inferential_results/` and optional selected-analyte follow-up outputs produced under `select_analytes/`.

Open `inferential_results/comparison_workbook.xlsx` first for replicate-aware results. Use `inferential_results/run_index.tsv` when you need exact generated paths.

### Replicate-aware inferential waterfall

![Replicate-aware inferential waterfall example](output-examples/replicate-aware-inferential-waterfall-se.png)

Source output type: `inferential_results/comparisons/<comparison_slug>/waterfall_plots/<method>/<method>_waterfall*.png`

All replicate-aware inferential waterfall methods show method-specific effect estimates with `+/- 1 SE` whiskers for each plotted analyte.

### Replicate-aware inferential barplots

![Replicate-aware inferential barplot example](output-examples/replicate-aware-inferential-barplot-fdr-page-1.png)

All-tested source output type: `inferential_results/comparisons/<comparison_slug>/barplots/<method>/all_tested/<method>_barplot_all_tested_page_<n>.png`

Significant-hit source output type: `inferential_results/comparisons/<comparison_slug>/barplots/<method>/significant_hits/<threshold>/<method>_barplot_<threshold>_page_<n>.png`

This example shows fixed-size analyte title strips, a linear fold-change-ratio y-axis, `+/- 1 SE` whiskers on both control and treatment bars, bracketed `*` annotations, a page caption defining the significance threshold, and fixed y-axis limits within the barplot set. For `normalized_t_test`, the plotted values use per-sample `normalized_signal`, defined as averaged duplicate raw analyte signal divided by that sample's reference-spot denominator. The denominator and exact raw reference signals used for each sample are reported in `input_qc/reference_spot_qc.tsv`. The low-signal threshold creates `low_signal_flag` for every replicate-aware method; it is separate from reference-spot normalization and does not remove analytes before p-value adjustment. The treatment bar is the treatment-arm mean `normalized_signal` divided by the control-arm mean `normalized_signal`; whiskers are symmetric on the plotted linear ratio scale unless the lower bound is clipped at zero. For `raw_log2_lm`, whiskers may appear asymmetric because log2-scale group-mean SE bounds are converted back to the linear ratio scale.

### Replicate-aware selected-analyte outputs

Selected-analyte outputs are optional follow-up outputs. They are created only when `scripts/select-analytes-analysis.R` is run with selected-analyte coordinates configured.

Source output types: `select_analytes/<comparison_slug>/<method>/selected_results.tsv`, `select_analytes/<comparison_slug>/<method>/selected_analyte_qc.tsv`, `select_analytes/<comparison_slug>/<method>/selected_waterfall.png`, and `select_analytes/<comparison_slug>/<method>/selected_bargraphs/<Analyte>.png`.

For the broader explanation of what each output means and where it is generated, see the root [README.md](../README.md).

## Exploratory Examples

### Threshold diagnostics

![Threshold diagnostics example](output-examples/threshold-diagnostics-region-stats.png)

Source output type: `threshold_diagnostics/region_stats.png`

### Main analysis waterfall

![Main analysis waterfall example](output-examples/main-waterfall-sorafenib.png)

Source output type: `main_analysis/.../all_comparisons/waterfalls/...png`

### Main analysis combined barplots

![Main analysis combined barplot example](output-examples/main-barplot-combined-page-1.png)

Source output type: `main_analysis/.../fold_change_hits/.../barplots/...png`

### Selected-analytes waterfall

![Selected analytes waterfall example](output-examples/select-analytes-waterfall.png)

Source output type: `select_analytes/<comparison_slug>/selected_waterfall.png`

Selected-analyte bargraphs are written as one PNG per configured analyte under `select_analytes/<comparison_slug>/selected_bargraphs/`.

## Method Notes

- [Worked Example: `raw_log2_lm` vs `normalized_t_test`](method-comparison-worked-example.md)
