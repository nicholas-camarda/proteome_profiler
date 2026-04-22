# Worked Example: `raw_log2_lm` vs `normalized_t_test`

This note explains what it means to say that `raw_log2_lm` models treatment effects on a scale where multiplicative changes are natural.

It uses one real comparison from the Nicole `male: vehicle vs aldosterone` run and walks through the arithmetic for two analytes:

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

- recomputes workbook-style normalized values from the raw sample sheets
- uses the membrane reference spots as the denominator for each sample
- runs a two-sample equal-variance t-test on those normalized replicate values
- reports effect size as `log2(mean(normalized treatment) / mean(normalized control))`

So it asks:

> after reference-spot normalization, do the treated and control groups differ on the normalized scale?

These are different questions. The same analyte can look more separated after normalization, less separated after normalization, or roughly unchanged.

## Important Terminology: Two Different "Reference" Concepts

There are two different reference concepts in this pipeline.

1. Membrane reference spots used for normalization

- These are the vendor protocol reference spots on the membrane.
- `normalized_t_test` uses them directly.
- In the current code, the raw-data normalization denominator is computed from the preferred reference spot pairs `A1,2` and `J1,2` when they are available.

2. User-chosen `thresholds$ref_coords` and `thresholds$ref_signal`

- These are not the normalization denominator.
- They are the low-signal reference panel chosen during `find_ref_thresh.R`.
- In replicate-aware analysis, they are used to create `low_signal_flag`.
- They do not rescale the data before `raw_log2_lm` or `normalized_t_test`.

So:

- `raw_log2_lm` does **not** normalize by the membrane reference spots
- `normalized_t_test` **does**
- both methods still carry the separate low-signal flag from `thresholds$ref_coords` / `ref_signal`

## Example 1: `CXCL9/MIG`

Male raw signals:

- vehicle: `1450, 2785, 2730, 1765`
- aldosterone: `2010, 3035, 5185, 2585`

Male normalized values:

- vehicle: `0.0349, 0.0293, 0.0257, 0.0209`
- aldosterone: `0.0386, 0.0422, 0.0456, 0.0417`

### `raw_log2_lm`

Take `log2` of each raw signal:

- vehicle:
  - `log2(1450) ≈ 10.50`
  - `log2(2785) ≈ 11.44`
  - `log2(2730) ≈ 11.41`
  - `log2(1765) ≈ 10.79`
- aldosterone:
  - `log2(2010) ≈ 10.97`
  - `log2(3035) ≈ 11.57`
  - `log2(5185) ≈ 12.34`
  - `log2(2585) ≈ 11.34`

Group means on the `log2` scale:

- mean vehicle log2 signal ≈ `11.05`
- mean aldosterone log2 signal ≈ `11.57`

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

- mean vehicle normalized ≈ `0.0277`
- mean aldosterone normalized ≈ `0.0420`

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

- all normalized aldosterone values are around `0.0386` to `0.0456`
- all normalized vehicle values are around `0.0209` to `0.0349`

So the normalization reduces the apparent overlap for this analyte, and the t-test gets much stronger evidence for a group difference.

## Example 2: `IL-6`

Male raw signals:

- vehicle: `2175, 4450, 5295, 6875`
- aldosterone: `5810, 16400, 10950, 21550`

Male normalized values:

- vehicle: `0.0523, 0.0469, 0.0498, 0.0814`
- aldosterone: `0.1116, 0.2279, 0.0963, 0.3480`

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

Here the normalized effect size is larger, but the normalized aldosterone values are also much more spread out:

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
- then compares the normalized replicate values directly
- effect size is `log2(mean normalized treatment / mean normalized control)`

## Practical Interpretation

Use `raw_log2_lm` when you want the question to be:

> are treated and control different on the raw measured signal scale, after log-transforming to interpret changes as fold-like differences?

Use `normalized_t_test` when you want the question to be:

> are treated and control different after each membrane has been normalized to its reference spots?

Neither method should silently replace the other. They are best treated as parallel first-class analyses that answer related but different questions.
