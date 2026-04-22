## Why

The analysis now produces defensible replicate-aware results, but the collaborator-facing output surface is still too easy to misread: summaries, method results, plots, and provenance are split across many artifacts without a single clear review path. The next change should make the default handoff workflow cloneable, organized, and maintainable without adding another report tree or requiring users to understand the internals.

## What Changes

- Add a collaborator-facing output contract that makes `inferential_results/comparison_workbook.xlsx` the primary review workbook for replicate-aware analyses.
- Consolidate run-level summaries into workbook tabs rather than generating another standalone summary artifact.
- Keep `run_index.tsv` as the machine-readable provenance and path index, not as the primary human-readable result.
- Keep `methods_overview.md` as a short prose note explaining estimands, fold-change definitions, p-value adjustment, and major interpretation caveats.
- Preserve method-specific result workbooks and comparison-specific plot folders, but define exactly what each artifact is for.
- Add selected-analyte summary content to the comparison workbook only when selected analytes are configured and results exist.
- Harden README and docs so a new lab user can run the workflow from `.env`, manifest, workbook, and protocol inputs without editing R scripts.
- Refactor and document the replicate-analysis helper code by responsibility using a staged, behavior-preserving module split with parity checks before and after movement.
- Exclude GitHub Actions from this change.

## Capabilities

### New Capabilities

- `collaborator-facing-output-contract`: Defines the user-facing output roles, consolidated comparison workbook contents, documentation expectations, and maintainability requirements for replicate-aware analysis handoff.

### Modified Capabilities

- `cloneable-runtime-configuration`: Clarifies that the README remains the primary start-here guide for clone-and-run usage and that selected analytes are optional follow-up inputs.

## Impact

- Affected documentation: `README.md`, `docs/README.md`, and `.env.example` comments if needed for clarity.
- Affected output generation: `scripts/helpers/replicate_analysis.R` and any helper modules split from it.
- Affected tests: workbook-sheet tests, artifact-contract tests, selected-analyte summary tests, documentation consistency checks, and smoke tests for clean output organization.
- No GitHub Actions or CI setup is included in this change.
