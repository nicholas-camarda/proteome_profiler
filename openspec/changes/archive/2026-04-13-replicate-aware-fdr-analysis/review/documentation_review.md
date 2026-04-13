# Documentation Review

## Resolved

- README now makes the manual interpretation step between `find_ref_thresh.R` and `main.R` explicit instead of implying a direct script-to-script handoff.
- The threshold diagnostic example is embedded where that manual decision actually happens.
- README now documents the unified input model: manifest-driven inputs are the documented interface for both exploratory and replicate-aware analyses, with `mode` controlling analysis semantics.
- README documents multi-comparison shortlist configuration through `shortlist$comparisons`.
- Root README now shows example output images directly instead of burying them only in `docs/README.md`.

## Remaining Follow-Ups

- Documentation still contains machine-local absolute file links and some maintainer-specific path examples. That is portability debt, not a methodological blocker.
- The output root still uses `output/plots/...` even though non-plot artifacts live there too. This is naming debt, not a correctness blocker.

## Assessment

No documentation blocker remains for this change. The major workflow confusions that previously affected correctness-of-use have been addressed.
