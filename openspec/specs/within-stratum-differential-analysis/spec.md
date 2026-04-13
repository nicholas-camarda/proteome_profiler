# within-stratum-differential-analysis Specification

## Purpose
TBD - created by archiving change replicate-aware-fdr-analysis. Update Purpose after archive.
## Requirements
### Requirement: The pipeline SHALL support within-stratum treatment contrasts
The system SHALL allow the user to define treatment-versus-control contrasts within a configured subgroup variable, such as comparing `vehicle` versus `aldosterone` separately within `male` and within `female`.

#### Scenario: Within-sex treatment comparisons
- **WHEN** the analysis config defines `sex` as the stratification variable and `vehicle vs aldosterone` as the contrast
- **THEN** the pipeline SHALL run one analyte-level comparison for `male` and one analyte-level comparison for `female`

#### Scenario: No direct subgroup contrast requested
- **WHEN** the config only requests within-stratum treatment comparisons
- **THEN** the pipeline SHALL not generate direct `male vs female` analyte tests

### Requirement: Statistical testing SHALL use biological replicate variance
For each analyte and requested comparison, the system SHALL fit a replicate-aware statistical test on sample-level log2-transformed signals instead of using only a single collapsed group value.

#### Scenario: Replicate-aware analyte test
- **WHEN** a stratum contains multiple biological replicates in both treatment arms
- **THEN** the pipeline SHALL estimate an analyte-level treatment effect and p-value from the sample-level observations in that stratum

#### Scenario: Insufficient replicates for inference
- **WHEN** a requested comparison does not meet the minimum replicate requirement in one or both treatment arms
- **THEN** the pipeline SHALL mark that analyte comparison as not testable and report the reason instead of emitting a misleading p-value

#### Scenario: Low-replication comparison warning
- **WHEN** a requested comparison has exactly two biological replicates per arm
- **THEN** the pipeline SHALL run the analyte-level test and mark the comparison output as low-replication

### Requirement: Differential-analysis outputs SHALL include effect estimates and sample counts
Each analyte-level differential-analysis result SHALL report enough information to interpret the test, including the comparison label, subgroup label, per-arm sample counts, group summaries, and an effect estimate on the analysis scale.

#### Scenario: Results table export
- **WHEN** the pipeline writes a differential-analysis result table
- **THEN** each row SHALL include the analyte identifier, subgroup, comparison, per-arm sample counts, effect estimate, and raw p-value

### Requirement: Low-signal status SHALL be recorded separately from inferential significance
The system SHALL evaluate low-signal status for each analyte using the configured reference/background rule, but it MUST keep that status separate from the decision to run inferential testing by default.

#### Scenario: Low-signal analyte is still tested
- **WHEN** an analyte is flagged as low-signal but has sufficient biological replicates for a requested comparison
- **THEN** the pipeline SHALL still compute the inferential result and record the low-signal flag in the output

