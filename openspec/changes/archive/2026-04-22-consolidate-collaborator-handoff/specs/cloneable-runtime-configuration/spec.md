## MODIFIED Requirements

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

#### Scenario: Selected analytes documented as optional follow-up
- **WHEN** the README describes selected-analyte configuration
- **THEN** it SHALL state that selected analytes are optional for setup validation and main analysis, and are only required when the user runs the selected-analyte follow-up workflow

#### Scenario: Documentation describes present behavior only
- **WHEN** user-facing documentation is updated for this runtime workflow
- **THEN** it SHALL describe the current supported workflow directly without historical comparisons or "now supported" language
