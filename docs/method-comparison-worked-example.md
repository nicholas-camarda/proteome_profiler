# Worked Example: `raw_log2_lm` vs `normalized_t_test`

This note explains what it means to say that `raw_log2_lm` models treatment effects on a scale where multiplicative changes are natural.

It uses one worked example comparison (`male: control vs treated`) and walks through the arithmetic for two analytes:

- `CXCL9/MIG`
- `IL-6`

## The Core Difference

`raw_log2_lm`:

- takes one averaged raw-signal value per biological sample
- keeps only finite positive signals
- applies `log2()` to those sample-level raw signals
- fits `lm(log2(signal) ~ treatment)`

So it asks:

> after log-transforming the raw measured signals, is the treated group higher or lower than control?

On the `log2` scale:

- `+1` means `2x`
- `-1` means `0.5x`
- `+0.585` means about `1.5x`

That is why multiplicative changes are natural there. Equal fold changes become equal distances on the `log2` scale.

`normalized_t_test`:

- recomputes one `normalized_signal` value for each sample and analyte from the raw sample sheets
- defines `normalized_signal` as averaged duplicate raw analyte signal divided by that sample's reference-spot denominator
- uses the mean raw signal of preferred Reference Spots pairs `A1,2` and `J1,2` as the denominator when those pairs are present; otherwise it uses available Reference Spots rows from that sample
- reports each sample's reference-spot denominator, reference rows, and raw reference signals in `input_qc/reference_spot_qc.tsv`
- runs a two-sample equal-variance t-test on those per-sample `normalized_signal` values
- reports effect size as `log2(mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values))`

So it asks:

> after reference-spot normalization, do the treated and control groups differ on the normalized scale?

These are different questions. The same analyte can look more separated after normalization, less separated after normalization, or roughly unchanged.

## Important Terminology: Two Different "Reference" Concepts

There are two different reference concepts in this pipeline.

1. Membrane reference spots used for normalization

- These are the vendor protocol reference spots on the membrane.
- `normalized_t_test` uses them directly.
- For each sample, the raw-data normalization denominator is the mean raw signal of preferred Reference Spots pairs `A1,2` and `J1,2` when they are available.
- If those preferred pairs are not available in a sample sheet, the analysis uses available Reference Spots rows from that sample.
- The resulting per-sample, per-analyte value is named `normalized_signal`.

2. User-chosen `thresholds$ref_coords` and `thresholds$ref_signal`

- These are not the normalization denominator.
- They are the low-signal reference panel chosen during `find_ref_thresh.R`.
- In replicate-aware analysis, they are used to create `low_signal_flag` for every configured method.
- They do not rescale the data before `raw_log2_lm` or `normalized_t_test`.
- They do not remove analytes before p-value adjustment.

So:

- `raw_log2_lm` does **not** normalize by the membrane reference spots
- `normalized_t_test` **does**
- both methods still carry the separate low-signal flag from `thresholds$ref_coords` / `ref_signal`

## Example 1: `CXCL9/MIG`

Male raw signals:

- control: `1450, 2785, 2730, 1765`
- treated: `2010, 3035, 5185, 2585`

Male normalized values:

- control: `0.0349, 0.0293, 0.0257, 0.0209`
- treated: `0.0386, 0.0422, 0.0456, 0.0417`

Here "normalized values" means the sample-level `normalized_signal` values for this analyte: each raw CXCL9/MIG signal divided by that same sample's reference-spot denominator.

### `raw_log2_lm`

Take `log2` of each raw signal:

- control:
  - `log2(1450) ≈ 10.50`
  - `log2(2785) ≈ 11.44`
  - `log2(2730) ≈ 11.41`
  - `log2(1765) ≈ 10.79`
- treated:
  - `log2(2010) ≈ 10.97`
  - `log2(3035) ≈ 11.57`
  - `log2(5185) ≈ 12.34`
  - `log2(2585) ≈ 11.34`

Group means on the `log2` scale:

- mean control log2 signal ≈ `11.05`
- mean treated log2 signal ≈ `11.57`

Difference:

- `11.57 - 11.05 ≈ 0.52`

Convert back to a fold change:

- `2^0.52 ≈ 1.43`

That matches the pipeline result:

- `effect_estimate_log2 ≈ 0.518`
- `fold_change_ratio ≈ 1.432`
- raw p-value ≈ `0.213`

### `normalized_t_test`

Group means on the normalized scale:

- mean control normalized ≈ `0.0277`
- mean treated normalized ≈ `0.0420`

Ratio:

- `0.0420 / 0.0277 ≈ 1.517`

Log2 ratio:

- `log2(1.517) ≈ 0.602`

That matches the pipeline result:

- `effect_estimate_log2 ≈ 0.602`
- `fold_change_ratio ≈ 1.517`
- raw p-value ≈ `0.00469`

### Why the p-values differ so much

After normalization, the groups separate much more cleanly:

- all normalized treated values are around `0.0386` to `0.0456`
- all normalized control values are around `0.0209` to `0.0349`

So the normalization reduces the apparent overlap for this analyte, and the t-test gets much stronger evidence for a group difference.

## Example 2: `IL-6`

Male raw signals:

- control: `2175, 4450, 5295, 6875`
- treated: `5810, 16400, 10950, 21550`

Male normalized values:

- control: `0.0523, 0.0469, 0.0498, 0.0814`
- treated: `0.1116, 0.2279, 0.0963, 0.3480`

### `raw_log2_lm`

The fitted result is:

- `effect_estimate_log2 ≈ 1.499`
- `fold_change_ratio ≈ 2.826`
- raw p-value ≈ `0.0329`

### `normalized_t_test`

The fitted result is:

- `effect_estimate_log2 ≈ 1.766`
- `fold_change_ratio ≈ 3.402`
- raw p-value ≈ `0.0579`

### Why this goes the other way

Here the normalized effect size is larger, but the normalized treated values are also much more spread out:

- approximately `0.096` to `0.348`

That larger within-group spread weakens the t-test, even though the mean difference is larger.

So normalization helps some analytes and hurts others. It is not automatically better or worse. It changes the question and can change the apparent separation between groups.

## What Is Statistically Different About the Two Methods

`raw_log2_lm`:

- works on raw measured intensity
- converts multiplicative raw-signal changes into additive differences on the `log2` scale
- effect size is treated minus control on the `log2(raw signal)` scale
- `2^(effect)` gives a raw-signal fold change

`normalized_t_test`:

- first divides each sample by its membrane reference-spot intensity
- then compares the per-sample `normalized_signal` values directly
- effect size is `log2(mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values))`

## Plot Scales And SE Bars

Replicate-aware barplots always use a linear fold-change-ratio y-axis relative to the control group.

For `normalized_t_test`:

- the control bar is `1`
- the treatment bar is `mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`
- the control whisker is `1 +/- SE(control-arm normalized_signal values) / mean(control-arm normalized_signal values)`
- the treatment whisker is `treatment fold change +/- SE(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`
- whiskers are symmetric on the plotted linear ratio scale unless the lower bound is clipped at zero

For `raw_log2_lm`:

- the control bar is `1`
- the treatment bar is `2^(mean log2 treatment - mean log2 control)`
- the control whisker is `2^(+/- SE(control log2 values))`
- the treatment whisker is `2^(effect_estimate_log2 +/- SE(treatment log2 values))`
- whiskers can look asymmetric because the bounds are computed on the log2 scale and converted back to the linear ratio scale

These barplot whiskers are visual group-mean SE bars on the plotted fold-change scale. They are not confidence intervals for the treatment-vs-control contrast, and they are not the same quantity as the standard error used for the treatment-effect test. Waterfall plots use `effect_estimate_log2 +/- effect_se_log2` because waterfalls visualize the method-specific treatment effect estimate.

## Practical Interpretation

Use `raw_log2_lm` when you want the question to be:

> are treated and control different on the raw measured signal scale, after log-transforming to interpret changes as fold-like differences?

Use `normalized_t_test` when you want the question to be:

> are treated and control different after each membrane has been normalized to its reference spots?

Neither method should silently replace the other. They are best treated as parallel first-class analyses that answer related but different questions.
