# Statistical Review

## Resolved

- Minimum-replication gating still counts only finite positive signals toward `min_reps`.
- Underpowered analytes continue to surface as `not_testable` instead of being silently modeled.
- The stronger replicate-aware fixture now demonstrates multiple comparison families in one run and includes a real BH edge case where `raw_p_value < 0.05` but `adjusted_p_value > 0.05`.
- Uneven subgroup sizes are now directly exercised in the test suite, including run-level low-replication warnings.
- Partial missingness/nonpositive values are now covered in the stronger fixture and still produce the expected `not_testable` behavior for the affected analyte/comparison.

## Remaining Caveats

- Inference at `n = 2` per arm remains technically estimable but weak; this is still a design choice, not a bug.
- Low-signal status remains a QC flag rather than a hard inferential exclusion, which is consistent with the intended analysis path.

## Assessment

No current statistical blocker remains. The previous robustness concern about untested multi-comparison FDR and uneven/missing replicate structures is now materially addressed by the stronger fixture coverage.
