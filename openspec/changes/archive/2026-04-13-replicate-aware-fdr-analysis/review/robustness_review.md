# Robustness Review

## Previously Blocking Issue

- The earlier blocker was fixture realism for replicate-aware inference.

## Closed By Current Tree

- A larger replicate-aware fixture now exercises multiple treatment contrasts in one run.
- The stronger fixture includes a nontrivial BH/FDR case where a raw-significant analyte is no longer significant after multiplicity correction.
- Uneven subgroup sizes are now part of the tested fixture rather than assumed away.
- Partial missingness/nonpositive values are now covered in a way that drives a concrete `not_testable` outcome for one analyte/comparison.
- Manifest-driven exploratory mode is now smoke-tested separately from replicate-aware inference, reducing the chance of input-model drift.

## Remaining Follow-Up Backlog

- Add a workbook-level large replicate-aware import fixture if deeper end-to-end LI-COR import stress coverage becomes necessary.
- Add explicit tests for malformed LI-COR ordering in a larger mixed-quality fixture, not just the current targeted invariant checks.
- Expand Python extractor coverage beyond path resolution and API import compatibility.

## Assessment

The earlier release blocker is closed. The current fixture and smoke evidence is now strong enough to support release and spec archival for this change.
