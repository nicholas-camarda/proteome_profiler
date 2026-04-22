## Context

The current repository has a working analysis pipeline, but its default runtime state is not clone-safe. `scripts/config/analysis_config.R` contains live analysis entries, a collaborator-specific default analysis, and local OneDrive/runtime paths. Package installation is partially covered for R, while Python protocol extraction depends on undeclared Python/Java setup. Runtime package loading also does not stop immediately when packages are unavailable.

Future lab users need a simple path: clone the repo, create `.env`, prepare a manifest/data files, run a setup check, then run the documented commands. They should not edit, copy, or create R scripts for routine use. The implementation must not change the scientific analysis behavior.

## Goals / Non-Goals

**Goals:**

- Make `.env` the single user-facing run configuration file for routine use.
- Parse structured run settings from `.env` using clear conventions for lists and key/value fields, then normalize them into the existing internal analysis-config shape.
- Keep the manifest as the source of truth for sample-level workbook inputs, including both one-workbook-per-sample and one-workbook-with-multiple-sheets layouts.
- Remove collaborator-specific active defaults from tracked runtime config.
- Make cloud storage optional; repo and runtime roots must be sufficient for a normal local run.
- Add a fail-fast setup check that reports missing packages, paths, manifests, workbooks, and writable output roots before long-running analysis.
- Provide a cloneable demo or smoke workflow that does not require private OneDrive data.
- Update documentation so a non-programmer can follow a command sequence without editing, copying, or creating R scripts.
- Include beginner-operational README steps with copyable commands starting from opening Terminal and cloning the repository.
- Make selected-analyte comparison slug syntax explicit enough that collaborators can configure multiple selected comparisons without guessing.
- Ensure replicate-aware waterfall and barplot outputs show method-specific standard-error uncertainty with annotation-safe plot limits.

**Non-Goals:**

- Do not require collaborators to edit R scripts or create one R script per analysis.
- Do not change fold-change calculations, p-value adjustments, or inferential estimands.
- Do not add Shiny or another GUI in this change.
- Do not support silent fallback paths when setup is invalid; fail clearly instead.

## Decisions

1. **Use `dotenv` for standard `.env` loading.**
   - Decision: Add `dotenv` as an R dependency and load repo-root `.env` before resolving run configuration in all entry scripts.
   - Rationale: A standard package is more readable and less error-prone than custom parsing, and this dependency is small relative to the usability benefit.
   - Alternatives considered: Hand-roll `.env` parsing; rejected because correctness around quoting, blank values, and comments is not worth maintaining here.

2. **Make `.env` the collaborator-facing run sheet.**
   - Decision: `.env` SHALL contain both runtime paths and analysis run settings such as mode, protocol workbook/PDF, manifest path, treatment/subgroup columns, comparisons, thresholds, and methods. Selected analytes are optional follow-up settings.
   - Rationale: This satisfies the actual usability requirement: one editable file plus manifest/data files, no R script edits.
   - Alternatives considered: Per-analysis R files or YAML analysis files; rejected for this change because they add another artifact and another syntax users must learn.

3. **Remove R config as a supported run-configuration source.**
   - Decision: `scripts/config/analysis_config.R` should be removed or converted into internal implementation code that does not contain user-specific analyses. Entry scripts and setup checks should build the run configuration from `.env`.
   - Rationale: Maintaining two supported configuration paths creates ambiguity and makes the user model harder to explain. This change should leave `.env` as the only supported run-configuration path.
   - Alternatives considered: Keep R config loading as a second run path; rejected because no external users have adopted the repository and a second configuration source would add confusion without benefit.

4. **Preflight setup becomes a required documented command.**
   - Decision: Add `scripts/check_setup.R` that loads `.env`, resolves the selected config, validates dependencies and inputs, and prints resolved paths.
   - Rationale: Future users should see one actionable setup report before any analysis writes outputs.
   - Alternatives considered: Let each entry script fail independently; rejected because failures are later, less clear, and harder for non-programmers.

5. **Cloud path is optional.**
   - Decision: Path resolution should search repo root and runtime root always, and cloud root only when configured.
   - Rationale: A lab mate should not need your OneDrive layout or any cloud provider to run a local analysis.

6. **Manifest workbook inputs remain first-class.**
   - Decision: The manifest remains responsible for `sample_id`, `workbook_path`, optional `sheet_name`, treatment, subgroup columns, and any other sample-level metadata. The README should recommend `workbooks/` for raw sample workbooks and `protocols/` for protocol files.
   - Rationale: Sample metadata is tabular and belongs in CSV, not `.env`. This also cleanly supports collaborator workbooks where many raw sample sheets live in one Excel file.

7. **Selected-analyte comparison slugs are explicit config values.**
   - Decision: `PROTEOME_PROFILER_SHORTLIST_COMPARISONS` should use `|` between selected comparison slugs. Replicate-aware slugs use `<subgroup value>_<control label>_vs_<treatment label>`; non-stratified and legacy slugs use `<control label>_vs_<treatment label>`.
   - Rationale: This keeps the selected-analyte workflow deterministic when an analysis has multiple strata/comparisons, and it lets validation report the exact available slugs when a user mistypes one.

8. **Replicate-aware plots show method-specific uncertainty.**
   - Decision: Replicate-aware waterfall plots show `+/- 1 SE` whiskers on the log2 effect scale, and replicate-aware barplots show treatment-bar SE whiskers converted to the fold-change ratio scale from `effect_estimate_log2 +/- effect_se_log2`.
   - Rationale: Biological-replicate analyses estimate both an effect and uncertainty. Showing fold-change bars without uncertainty would overstate visual precision and make method-specific plots less defensible.
   - Implementation note: Ratio-scale whiskers are asymmetric when transformed back from symmetric log2 SE bounds.

9. **Provide a demo/smoke path.**
   - Decision: Include or generate small fixture data/config sufficient to run the workflow without private files.
   - Rationale: The fastest way to validate a clone is to run a known-good minimal analysis.

## Risks / Trade-offs

- **Risk: `.env` can become long.** → Mitigation: Use commented sections in `.env.example`, strict validation, and simple delimiter conventions; keep per-sample information in the manifest CSV.
- **Risk: Structured values in `.env` can be mistyped.** → Mitigation: `check_setup.R` validates parsed comparisons, thresholds, methods, optional analyte lists, and column names before analysis.
- **Risk: Adding setup checks duplicates validation logic.** → Mitigation: Reuse existing config/path/manifest validation helpers where possible; do not create a parallel validation universe.
- **Risk: Removing R config loading may require broad test updates.** → Mitigation: Update fixture/test helpers to generate `.env` run sheets and validate the same internal config objects from the new source.
