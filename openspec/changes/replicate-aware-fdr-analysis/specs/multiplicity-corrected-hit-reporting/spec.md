## ADDED Requirements

### Requirement: P-value adjustment MUST be configurable and default to BH
The system MUST support a configurable p-value adjustment method for analyte-level differential testing, with `BH` as the default.

#### Scenario: Default adjustment method
- **WHEN** the user does not specify a p-value adjustment method
- **THEN** the pipeline SHALL apply the `BH` method to the analyte-level p-values for each comparison family

#### Scenario: Alternate adjustment method requested
- **WHEN** the user configures a supported adjustment method such as `holm`, `bonferroni`, or `BY`
- **THEN** the pipeline SHALL apply that method and record the selected method in the output metadata

### Requirement: Significant-hit calls SHALL be based on adjusted p-values
The system SHALL determine inferential significance from multiplicity-adjusted p-values and a configurable alpha threshold, rather than from fold-change thresholds alone.

#### Scenario: Significant analyte after correction
- **WHEN** an analyte has an adjusted p-value less than or equal to the configured alpha
- **THEN** the inferential results SHALL mark that analyte as significant for that comparison

#### Scenario: Exploratory fold-change hit without inferential significance
- **WHEN** an analyte passes an exploratory fold-change threshold but its adjusted p-value exceeds alpha
- **THEN** the pipeline SHALL not mark it as inferentially significant

### Requirement: The pipeline SHALL emit explicit multiplicity-corrected results tables
For each requested comparison family, the system SHALL write a machine-readable results table containing the analyte identifier, effect estimate, raw p-value, adjusted p-value, adjustment method, alpha threshold, and significant-hit flag.

#### Scenario: Per-comparison result export
- **WHEN** a within-stratum comparison completes successfully
- **THEN** the pipeline SHALL write a results table for that comparison with both raw and adjusted inferential results

#### Scenario: Output labeling clarity
- **WHEN** inferential results are written alongside exploratory outputs
- **THEN** the output tree and filenames SHALL distinguish adjusted-p-value significance from threshold-based exploratory hits

### Requirement: Inferential result exports SHALL include per-comparison files and run-level summaries
The system SHALL write one canonical machine-readable results file per stratum-specific comparison, plus a run-level index file that lists the available comparisons and their output locations. The system SHOULD also emit a combined workbook for collaborator-friendly review.

#### Scenario: Canonical per-comparison export
- **WHEN** the run includes `male: vehicle vs aldosterone` and `female: vehicle vs aldosterone`
- **THEN** the pipeline SHALL write separate inferential results files for each comparison and a run-level index file

#### Scenario: Convenience workbook export
- **WHEN** inferential results are generated successfully
- **THEN** the pipeline SHOULD write a combined workbook that collects the comparison tables into one collaborator-facing file
