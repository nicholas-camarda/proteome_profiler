# replicate-aware-inputs Specification

## Purpose
TBD - created by archiving change replicate-aware-fdr-analysis. Update Purpose after archive.
## Requirements
### Requirement: Replicate-aware analyses SHALL use explicit sample metadata
The system SHALL support a manifest-driven input mode where each LI-COR workbook is linked to a unique biological sample identifier and the analysis metadata needed to define comparisons, including treatment and optional subgroup columns such as sex.

#### Scenario: Valid manifest-driven run
- **WHEN** the user provides a manifest with `sample_id`, workbook path, treatment, and subgroup metadata for all included workbooks
- **THEN** the pipeline SHALL build the analysis dataset from the manifest rather than inferring biology from filename prefixes

#### Scenario: Missing required metadata
- **WHEN** a manifest is missing a required column such as `sample_id`, workbook path, or treatment
- **THEN** the pipeline SHALL stop with a direct validation error naming the missing fields

### Requirement: Analysis mode SHALL be selected explicitly from configuration
The system SHALL choose between exploratory mode and replicate-aware inferential mode from the run configuration rather than from filename heuristics alone.

#### Scenario: Manifest-driven replicate-aware mode
- **WHEN** the run configuration includes a sample manifest
- **THEN** the pipeline SHALL enter replicate-aware mode and use the manifest as the experimental-design source of truth

#### Scenario: Exploratory mode
- **WHEN** the run configuration does not include a sample manifest and the input contains exactly one workbook per configured group
- **THEN** the pipeline SHALL enter exploratory mode

### Requirement: The manifest SHALL be the source of truth for experimental design
In replicate-aware mode, the manifest SHALL define sample identity and analysis-group metadata. Filenames MAY be human-readable, but the pipeline MUST NOT infer treatment, sex, or replicate identity from filenames when a manifest is present.

#### Scenario: Arbitrary workbook filenames
- **WHEN** the manifest points to workbooks with arbitrary filenames that do not encode treatment or subgroup labels
- **THEN** the pipeline SHALL use the manifest metadata rather than attempting to parse those filenames

#### Scenario: Configured subgroup column validation
- **WHEN** the analysis config requests a subgroup variable such as `sex`
- **THEN** the pipeline SHALL verify that the manifest includes that column before analysis begins

### Requirement: Replicate-aware manifests SHALL expose a minimal required schema
Replicate-aware input manifests SHALL contain at least `sample_id`, `workbook_path`, and `treatment` columns. Any subgrouped analysis SHALL also require the configured subgroup column, such as `sex`.

#### Scenario: Valid within-sex manifest
- **WHEN** the user requests within-sex treatment comparisons
- **THEN** the manifest SHALL include `sample_id`, `workbook_path`, `treatment`, and `sex`

#### Scenario: Unique sample identifiers
- **WHEN** the manifest contains duplicate `sample_id` values
- **THEN** the pipeline SHALL stop with a validation error rather than merging those rows implicitly

### Requirement: The analysis dataset SHALL preserve biological replicates
The system SHALL average the duplicate membrane spots within each workbook to obtain one analyte-level value per sample, but it MUST retain separate rows for different biological samples assigned to the same treatment arm or subgroup.

#### Scenario: Multiple mice in the same treatment arm
- **WHEN** four `male_vehicle` workbooks are present in the manifest
- **THEN** the canonical analysis dataset SHALL contain four distinct sample-level analyte measurements for that stratum and treatment arm rather than collapsing them to one group-level value

#### Scenario: Technical duplicate spot averaging
- **WHEN** a workbook contains the expected duplicate membrane spots for one analyte
- **THEN** the pipeline SHALL average those within-workbook duplicate spots before attaching the sample-level metadata

### Requirement: Ambiguous replicate layouts MUST be rejected
The system MUST fail fast when multiple workbooks would be merged into the same analysis group without explicit sample metadata.

#### Scenario: Multiple files share the same exploratory group label
- **WHEN** the exploratory import path sees more than one workbook for a group label such as `vehicle`
- **THEN** the pipeline SHALL stop with an error explaining that biological replicates require explicit replicate-aware metadata

### Requirement: Non-replicate exploratory analyses SHALL remain supported
The system SHALL continue to support the current one-workbook-per-group exploratory workflow for analyses that do not have biological replicates and therefore do not require a manifest-driven inferential path.

#### Scenario: Single workbook per group exploratory run
- **WHEN** the input contains exactly one workbook per configured analysis group
- **THEN** the pipeline SHALL allow the existing exploratory workflow to run without requiring sample-level biological replicate metadata

#### Scenario: Exploratory run output semantics
- **WHEN** a one-workbook-per-group exploratory analysis completes
- **THEN** the pipeline SHALL keep its outputs labeled as exploratory rather than inferentially significant
