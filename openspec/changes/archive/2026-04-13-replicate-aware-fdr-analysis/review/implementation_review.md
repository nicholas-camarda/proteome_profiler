# Implementation Review

## Resolved

- Input declaration is no longer split unnecessarily: exploratory and replicate-aware analyses can both use `input$manifest`, while `mode` controls analysis semantics.
- `find_ref_thresh.R` and `select-analytes-analysis.R` now follow that same input model instead of treating manifest-driven input as inherently inferential-only.
- Replicate-aware shortlist generation supports multiple comparison slugs in one run.
- The Python protocol extractor now runs against the real protocol PDF in this environment by handling the installed `tabula.io` import path.
- The shipped config surface is more compact and user-facing, while helper code normalizes it internally for backwards compatibility.

## Remaining Follow-Ups

- A dedicated image-content assertion for inferential waterfall title clipping would still be stronger than existence-only render checks.

## Assessment

No current implementation blocker remains. The code path now matches the intended architecture closely enough to support release and archival of this change.
