# Artifact Review Checklist

Use this checklist before sharing collaborator-facing results. It is a manual review guide; it does not create additional output files.

## Start With The Workbook

- Open `inferential_results/comparison_workbook.xlsx` first for replicate-aware analyses.
- Confirm the summary, input QC, method summary, significance summary, and comparison sheets match the intended run.
- Confirm `selected_analytes_summary` appears only when selected-analyte follow-up was run.
- Check that raw p-value significance and FDR thresholds are interpreted separately.
- Confirm `n_tested`, `n_p_adjust_hypotheses`, low-signal flags, and hit counts are plausible for each comparison and method.

## Check Inputs And QC

- Review `input_qc/reference_spot_qc.tsv` for reference-spot denominators, raw reference signals, source rows, and QC status.
- Confirm no sample has unexplained missing, nonpositive, or invalid reference denominators.
- Review `input_qc/summary.tsv` and any setup warnings for missing sheets, nonpositive signals, or skipped analytes.
- Confirm `threshold_diagnostics/` matches the configured low-signal reference panel and threshold.

## Check Methods And Paths

- Read `inferential_results/methods_overview.md` and confirm the method definitions match the analysis being shared.
- Use `inferential_results/run_index.tsv` to verify generated paths exist and point to the expected comparison folders.
- Confirm the comparison slugs match the manifest labels and `.env` comparison settings.

## Check Plots

- Inspect representative waterfall plots for every method being shared and confirm effect-direction ordering, labels, and `+/- 1 SE` whiskers are visible when biological replicates exist.
- Inspect representative 25-per-page barplots and confirm the y-axis is labeled as a linear fold-change ratio relative to the control group.
- Confirm control or vehicle appears first in each barplot.
- Confirm SE whiskers, bracket/star annotations, page captions, and fixed y-axis limits are readable and not clipped.
- Confirm selected-analyte plots are present only for configured selected analytes and the configured comparison/method folders.
- Review `selected_analyte_qc.tsv` for selected analytes that are low-signal flagged or not plotted.

## Handoff Decision

- Treat low-signal flags as QC context, not automatic proof that an analyte should be ignored.
- Do not describe raw p < alpha hits as FDR-significant.
- Do not describe fold-change magnitude as statistical evidence.
- If any QC flag, missing plot, or unexpected hit-count discrepancy remains unexplained, resolve it before sharing results.
