## ADDED Requirements

### Requirement: Runtime run settings SHALL be loadable from `.env`
The system SHALL load the active run configuration from a repo-root `.env` file before selecting an analysis mode or resolving project paths.

#### Scenario: Active analysis selected from `.env`
- **WHEN** `.env` defines `PROTEOME_PROFILER_ANALYSIS`
- **THEN** entry scripts SHALL use that value unless the process environment already provides an explicit override

#### Scenario: Analysis metadata selected from `.env`
- **WHEN** `.env` defines analysis mode, protocol paths, input manifest, comparison settings, thresholds, and methods
- **THEN** the pipeline SHALL parse those values into the internal analysis configuration used by the analysis scripts

#### Scenario: Optional selected analytes omitted
- **WHEN** `.env` does not define selected analytes
- **THEN** setup validation and the main analysis SHALL still run if all required analysis settings are valid

#### Scenario: Optional selected analytes provided
- **WHEN** `.env` defines selected analytes
- **THEN** the selected-analyte follow-up workflow SHALL use those analytes for user-requested focused outputs

#### Scenario: Selected-analyte comparison slugs configured
- **WHEN** `.env` defines `PROTEOME_PROFILER_SHORTLIST_COMPARISONS`
- **THEN** the pipeline SHALL parse multiple selected comparison slugs using `|` as the delimiter and SHALL document replicate-aware slugs as `<subgroup value>_<control label>_vs_<treatment label>`

#### Scenario: Non-stratified selected-analyte comparison slugs configured
- **WHEN** `.env` defines `PROTEOME_PROFILER_SHORTLIST_COMPARISONS` for a non-stratified or legacy analysis
- **THEN** the pipeline SHALL document and validate slugs as `<control label>_vs_<treatment label>`

#### Scenario: Invalid selected-analyte comparison slug
- **WHEN** a configured selected-analyte comparison slug does not match an available comparison
- **THEN** the workflow SHALL fail clearly with the expected slug format and the available comparison slugs

#### Scenario: Selected-analyte workflow without selected analytes
- **WHEN** the selected-analyte follow-up script is run without selected analytes configured
- **THEN** the script SHALL stop with a clear message explaining that selected analytes are optional for main analysis but required for selected-analyte outputs

#### Scenario: Runtime root selected from `.env`
- **WHEN** `.env` defines `PROTEOME_PROFILER_RUNTIME_ROOT`
- **THEN** path resolution SHALL use that value as the runtime output root

#### Scenario: Missing `.env` file
- **WHEN** no `.env` file exists
- **THEN** setup validation SHALL report that the file is missing and point the user to `.env.example`

### Requirement: `.env` SHALL be the supported run-configuration source
The pipeline SHALL build the active run configuration from `.env` and SHALL NOT require or support user-authored R config files as run configuration.

#### Scenario: Fresh clone without private config
- **WHEN** a user clones the repository without private data
- **THEN** the documented setup check SHALL fail only with actionable missing local setup messages, not by trying to run a collaborator-specific analysis by default

#### Scenario: R config override not supported
- **WHEN** a user attempts to configure a run by editing or pointing to an R config file
- **THEN** the supported workflow SHALL direct them to `.env` instead

#### Scenario: No collaborator-specific tracked defaults
- **WHEN** setup validation starts from a fresh clone
- **THEN** no tracked R config file SHALL select Nicole, Nick, or any collaborator-specific analysis implicitly

### Requirement: Routine use SHALL NOT require R script edits
The documented collaborator workflow SHALL NOT require users to edit, copy, or create R scripts for their own analyses.

#### Scenario: New collaborator configures an analysis
- **WHEN** a collaborator needs to configure a new analysis
- **THEN** the README SHALL direct them to edit `.env` and their manifest/input files, not an R script

#### Scenario: Entry scripts remain stable
- **WHEN** a collaborator switches from one analysis to another
- **THEN** they SHALL change `.env` values rather than modifying `scripts/main.R`, `scripts/find_ref_thresh.R`, `scripts/select-analytes-analysis.R`, or config helper code

### Requirement: Cloud storage roots SHALL be optional
The system SHALL resolve project-relative paths using the repository root and runtime root, and SHALL include a cloud root only when one is configured.

#### Scenario: Cloud root omitted
- **WHEN** `PROTEOME_PROFILER_CLOUD_PARENT` is empty or unset
- **THEN** path resolution SHALL not search or require a cloud storage location

#### Scenario: Cloud root configured
- **WHEN** `PROTEOME_PROFILER_CLOUD_PARENT` is configured
- **THEN** path resolution MAY include the configured cloud project root after repo and runtime candidates

### Requirement: Manifest workbook paths SHALL support documented input layouts
The manifest SHALL remain the source of truth for sample-level workbook inputs, and the documentation SHALL describe where users should place workbook files and how workbook paths are resolved.

#### Scenario: Recommended workbook folder layout
- **WHEN** a user follows the README for a new analysis
- **THEN** the README SHALL show a recommended layout with manifests under `manifests/`, raw sample workbooks under `workbooks/`, and protocol files under `protocols/`

#### Scenario: Relative workbook paths
- **WHEN** the manifest contains a relative `workbook_path` such as `workbooks/M_VEH_01.xlsx`
- **THEN** setup validation and analysis scripts SHALL resolve it against the repository root, runtime root, and configured cloud root candidates

#### Scenario: Absolute workbook paths
- **WHEN** the manifest contains an absolute `workbook_path`
- **THEN** setup validation and analysis scripts SHALL use that path directly and fail clearly if it does not exist

### Requirement: Multi-sheet sample workbooks SHALL be supported
The manifest SHALL support analyses where multiple raw sample sheets are stored in one Excel workbook, and the pipeline SHALL only read sheets explicitly listed in the manifest.

#### Scenario: One workbook per biological sample
- **WHEN** each biological sample has its own workbook and the manifest provides `sample_id`, `workbook_path`, treatment, and optional subgroup columns
- **THEN** the pipeline SHALL read each listed workbook as one sample without requiring `sheet_name`

#### Scenario: One workbook with multiple sample sheets
- **WHEN** multiple biological samples are stored as sheets in one workbook and the manifest provides `sample_id`, `workbook_path`, `sheet_name`, treatment, and optional subgroup columns
- **THEN** the pipeline SHALL read each listed sheet as one sample

#### Scenario: Analysis tabs excluded by omission
- **WHEN** a multi-sheet workbook also contains summary or analysis tabs such as `ALL data` or `Benjamini Hochberg`
- **THEN** the pipeline SHALL ignore those tabs unless they are explicitly listed in the manifest

#### Scenario: Missing listed sheet
- **WHEN** the manifest lists a `sheet_name` that does not exist in the referenced workbook
- **THEN** setup validation SHALL fail before analysis and name the missing workbook sheet

### Requirement: Setup validation SHALL fail fast with actionable diagnostics
The system SHALL provide a setup validation command that checks dependencies, selected analysis configuration, protocol files, manifests, workbook paths, and output writability before long-running analysis commands are run.

#### Scenario: Missing R package
- **WHEN** a required R package is not installed
- **THEN** setup validation SHALL fail with the missing package names and the documented installation command

#### Scenario: Missing Python protocol dependency
- **WHEN** protocol extraction dependencies are unavailable
- **THEN** setup validation SHALL report the missing dependency and the documented Python setup command

#### Scenario: Missing Java dependency for Tabula
- **WHEN** Java is unavailable and protocol extraction may require `tabula-py`
- **THEN** setup validation and documentation SHALL provide a single-command Java installation path for macOS users

#### Scenario: Missing analysis input
- **WHEN** the selected analysis references a missing manifest, protocol workbook, protocol PDF, data directory, or sample workbook
- **THEN** setup validation SHALL fail with the unresolved path and the candidate locations checked

#### Scenario: Valid setup
- **WHEN** dependencies and selected analysis inputs are valid
- **THEN** setup validation SHALL print the active analysis name, runtime output root, resolved input paths, parsed comparison settings, and next recommended command

### Requirement: Runtime package loading SHALL stop on missing packages
The analysis entry scripts SHALL stop before analysis execution when required R packages are missing.

#### Scenario: Entry script missing package
- **WHEN** a user runs `scripts/main.R`, `scripts/find_ref_thresh.R`, or `scripts/select-analytes-analysis.R` without required R packages installed
- **THEN** the script SHALL stop with a clear missing-package message instead of continuing to an obscure downstream error

### Requirement: A cloneable demo or smoke workflow SHALL be available
The repository SHALL provide a documented workflow that validates a clone without requiring private collaborator data or a user-specific OneDrive directory.

#### Scenario: New user runs demo check
- **WHEN** a new user follows the README demo/setup commands on a clean clone
- **THEN** the commands SHALL exercise config loading, dependency checks, path resolution, and at least one small analysis fixture or generated smoke dataset

### Requirement: Replicate-aware plots SHALL show standard-error uncertainty
Replicate-aware inferential and selected-analyte plotting SHALL include method-specific standard-error uncertainty when biological replicate methods are used.

#### Scenario: Replicate-aware waterfall standard errors
- **WHEN** replicate-aware inferential waterfall plots are generated for any supported method
- **THEN** each plotted effect estimate SHALL show `+/- 1 SE` whiskers derived from the method-specific `effect_se_log2`

#### Scenario: Replicate-aware barplot standard errors
- **WHEN** replicate-aware inferential barplots are generated for any supported method
- **THEN** treatment fold-change bars SHALL show `+/- 1 SE` whiskers converted from log2 effect scale to fold-change ratio scale

#### Scenario: Replicate-aware selected-analyte barplot standard errors
- **WHEN** replicate-aware selected-analyte bargraphs are generated
- **THEN** the selected-analyte bargraphs SHALL use the same fold-change barplot renderer and SHALL include the method-specific treatment SE whiskers when finite SE values are available

#### Scenario: Multi-method inferential run index
- **WHEN** replicate-aware inferential analysis runs multiple methods
- **THEN** `inferential_results/run_index.tsv` SHALL contain one row per comparison and method, with output paths for every generated method-specific table, waterfall plot, and barplot set

#### Scenario: Selected-analyte comparisons resolve independently of method rows
- **WHEN** selected-analyte follow-up reads a multi-method `run_index.tsv`
- **THEN** selected comparison slugs SHALL resolve to unique comparison folders rather than rerunning once per method row

#### Scenario: Ratio-scale SE whisker bounds
- **WHEN** a method reports `effect_estimate_log2` and `effect_se_log2`
- **THEN** treatment barplot whisker bounds SHALL be computed as `2^(effect_estimate_log2 - effect_se_log2)` and `2^(effect_estimate_log2 + effect_se_log2)`

#### Scenario: Annotation and number clipping avoidance
- **WHEN** replicate-aware barplots include bars, SE whiskers, value labels, and significance brackets
- **THEN** y-axis limits SHALL include the maximum finite bar or SE-whisker upper bound plus annotation headroom, and plot clipping SHALL not truncate bracket or label annotations

### Requirement: Documentation SHALL present `.env` as the primary setup surface
The README SHALL guide users through copying `.env.example`, editing one `.env` file, running setup validation, and then running the analysis scripts.

#### Scenario: Beginner clone instructions
- **WHEN** a user with no command-line experience reads the README
- **THEN** the README SHALL include explicit instructions to open Terminal, clone `https://github.com/nicholas-camarda/proteome_profiler.git`, enter the repository directory, copy `.env.example` to `.env`, and run the setup/analysis commands

#### Scenario: Copyable command blocks
- **WHEN** the README describes setup and execution steps
- **THEN** each step SHALL include copyable command blocks and brief plain-language context for what the command does

#### Scenario: New lab user follows README
- **WHEN** a lab user with no R config experience reads the quick start
- **THEN** the primary path SHALL be editing `.env` and following commands, with no R script edits required

#### Scenario: Documentation describes present behavior only
- **WHEN** user-facing documentation is updated for this runtime workflow
- **THEN** it SHALL describe the current supported workflow directly without historical comparisons or "now supported" language
