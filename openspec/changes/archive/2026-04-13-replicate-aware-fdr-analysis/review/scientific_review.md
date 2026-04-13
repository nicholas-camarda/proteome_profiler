# Scientific Review

## Key Risks

- The replicate-aware workflow is associational within subgroup, not causal and not interaction modeling.
- Shortlist outputs can fall back to top-ranked tested analytes when `selection_top_n` is configured and no significant hits exist.
- `n = 2` per arm is supported but should be treated as weak evidence even when adjusted p-values are reported.

## Guidance

- Report subgroup analyses as separate within-subgroup treatment comparisons unless an interaction model is added.
- Do not describe shortlist outputs as statistically significant without checking the underlying `results.tsv`.
- Add sample counts and low-replication warnings to any collaborator or manuscript-facing summary.

## Assessment

No causal overclaim is encoded in the current repo. Remaining risk is overinterpretation by downstream users rather than a direct scientific-logic bug in the implementation.
