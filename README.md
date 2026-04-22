# Proteome Profiler Array Analysis

This repository analyzes membrane-based Proteome Profiler arrays from LI-COR spot-intensity Excel exports.

It supports two analysis modes:

- `exploratory`: fold-change visualization when there are not enough biological replicates for inferential testing.
- `replicate`: biological-replicate-aware treatment-versus-control inference using a sample manifest.

Routine users edit one local `.env` file plus their manifest CSV. Users run the R scripts, but do not edit, copy, or create R scripts.

Rendered output examples are described in [docs/README.md](docs/README.md).

## Quick Start

### 1. Open Terminal

On macOS, open **Terminal** from Applications, Spotlight, or Launchpad.

### 2. Clone The Repository

Copy and paste:

```bash
git clone https://github.com/nicholas-camarda/proteome_profiler.git
cd proteome_profiler
```

### 3. Install Dependencies

Install the R packages:

```bash
Rscript scripts/install_packages.R
```

Install the Python packages used for protocol-table extraction:

```bash
python3 -m pip install -r requirements.txt
```

The protocol PDF extractor uses `tabula-py`, which requires Java. Confirm Java is available:

```bash
java -version
```

If that command says Java is missing, install a Java runtime on macOS with Homebrew:

```bash
brew install --cask temurin
```

Then run `java -version` again.

### 4. Create A Local Run Sheet

Copy the example `.env` file:

```bash
cp .env.example .env
```

Open `.env` in a text editor and edit it for the analysis you want to run.

Do not commit `.env`. It is ignored by git because it contains local paths and run settings.

### 5. Create The Demo Inputs

The example `.env` is configured for a small generated demo dataset. Create it with:

```bash
Rscript scripts/setup/create_demo_data.R
```

### 6. Check Setup

Run setup validation before running the analysis:

```bash
Rscript scripts/check_setup.R
```

This checks the `.env` file, required packages, protocol files, manifest columns, workbook paths, workbook sheets, and output folder writability.

### 7. Run The Analysis

Review low-signal threshold diagnostics:

```bash
Rscript scripts/find_ref_thresh.R
```

Run the main analysis:

```bash
Rscript scripts/main.R
```

For replicate-aware analyses, open this workbook first:

```text
PROTEOME_PROFILER_RUNTIME_ROOT/output/<user>/<slug>/inferential_results/comparison_workbook.xlsx
```

Use `inferential_results/run_index.tsv` when you need exact file paths for every generated table and plot.

Create selected-analyte follow-up plots if `PROTEOME_PROFILER_SHORTLIST_ANALYTES` is set in `.env`:

```bash
Rscript scripts/select-analytes-analysis.R
```

Selected analytes are optional. `find_ref_thresh.R` and `main.R` do not require selected analytes.

## Recommended File Layout

Use this layout for a project:

```text
proteome_profiler/
  .env
  manifests/
    my_samples.csv
  workbooks/
    M_VEH_01.xlsx
    M_VEH_02.xlsx
    M_ALDO_01.xlsx
    F_VEH_01.xlsx
    F_ALDO_01.xlsx
  protocols/
    cytokine_array_protocol.pdf
    cytokine_array_protocol.xlsx
```

Relative paths in `.env` and the manifest are resolved against:

```text
1. the repository folder
2. PROTEOME_PROFILER_RUNTIME_ROOT
3. PROTEOME_PROFILER_CLOUD_PARENT/proteome_profiler, if configured
```

Absolute paths are also allowed.

## The `.env` File

`.env` is the run sheet. It defines which analysis to run, where the inputs are, which comparisons to test, and where outputs go.

List values use `|`.

Comparison values use `control=treatment1|treatment2`.

Example replicate-aware run:

```bash
PROTEOME_PROFILER_ANALYSIS=nicole_mf_veh_vs_aldo
PROTEOME_PROFILER_MODE=replicate
PROTEOME_PROFILER_USER=nicole
PROTEOME_PROFILER_SLUG=nicole_mf_veh_vs_aldo

PROTEOME_PROFILER_RUNTIME_ROOT=~/ProjectsRuntime/proteome_profiler
PROTEOME_PROFILER_CLOUD_PARENT=

PROTEOME_PROFILER_PROTOCOL_PRESET=cytokine_xl
PROTEOME_PROFILER_PROTOCOL_WORKBOOK=protocols/cytokine_array_protocol.xlsx
PROTEOME_PROFILER_PROTOCOL_PDF=protocols/cytokine_array_protocol.pdf
PROTEOME_PROFILER_PROTOCOL_PAGES="17|18|19"

PROTEOME_PROFILER_INPUT_MANIFEST=manifests/my_samples.csv
PROTEOME_PROFILER_TREATMENT_COLUMN=treatment
PROTEOME_PROFILER_SUBGROUP_COLUMN=sex

PROTEOME_PROFILER_COMPARISONS="vehicle=aldosterone"
PROTEOME_PROFILER_REF_COORDS="E21,22|G13,14|E17,18"
PROTEOME_PROFILER_REF_SIGNAL=1500

PROTEOME_PROFILER_MIN_REPS_PER_ARM=2
PROTEOME_PROFILER_P_ADJUST_METHOD=BH
PROTEOME_PROFILER_ALPHA=0.05
PROTEOME_PROFILER_ANALYSIS_METHODS="raw_log2_lm|normalized_t_test"
```

Optional selected-analyte follow-up settings:

```bash
PROTEOME_PROFILER_SHORTLIST_COMPARISONS="male_vehicle_vs_aldosterone|female_vehicle_vs_aldosterone"
PROTEOME_PROFILER_SHORTLIST_METHODS="raw_log2_lm|normalized_t_test"
PROTEOME_PROFILER_SHORTLIST_ANALYTES="CCL3/CCL4/MIP-1α/β|IL-1α/IL-1F1|IL-10"
```

The selected-analyte settings are only required when running `Rscript scripts/select-analytes-analysis.R`.

### Writing `PROTEOME_PROFILER_SHORTLIST_COMPARISONS`

`PROTEOME_PROFILER_SHORTLIST_COMPARISONS` is only for selected-analyte follow-up plots. It tells `scripts/select-analytes-analysis.R` which finished comparison folders to use after `scripts/main.R` has already run.

For replicate-aware analyses with a subgroup column, write each selected comparison slug as:

```text
<subgroup value>_<control label>_vs_<treatment label>
```

Use the exact labels from the manifest and `.env`, converted to lowercase, with spaces, punctuation, and other non-letter/non-number characters replaced by `_`.

Example:

```text
PROTEOME_PROFILER_SUBGROUP_COLUMN=sex
PROTEOME_PROFILER_COMPARISONS="vehicle=aldosterone"
```

If the manifest has `sex` values `male` and `female`, the available comparison slugs are:

```text
male_vehicle_vs_aldosterone
female_vehicle_vs_aldosterone
```

To select both:

```bash
PROTEOME_PROFILER_SHORTLIST_COMPARISONS="male_vehicle_vs_aldosterone|female_vehicle_vs_aldosterone"
```

For a non-stratified or exploratory comparison, omit the subgroup prefix:

```bash
PROTEOME_PROFILER_SHORTLIST_COMPARISONS="vehicle_vs_aldosterone"
```

If you leave `PROTEOME_PROFILER_SHORTLIST_COMPARISONS` blank and the analysis has exactly one comparison, the selected-analyte script uses that comparison. If the analysis has multiple comparisons, set this field explicitly.

For `exploratory` mode without a manifest, use a data directory and group levels:

```bash
PROTEOME_PROFILER_MODE=exploratory
PROTEOME_PROFILER_INPUT_DATA_DIR=workbooks
PROTEOME_PROFILER_GROUP_LEVELS="vehicle|treated"
PROTEOME_PROFILER_COMPARISONS="vehicle=treated"
PROTEOME_PROFILER_FOLD_CHANGE=1.5
PROTEOME_PROFILER_GROUPS_PER_PAGE=25
```

## Manifest CSV

The manifest is the sample table. It tells the pipeline which Excel workbook belongs to each biological sample and which metadata labels define the analysis.

For replicate-aware analyses, include:

```text
sample_id
workbook_path
treatment
```

If the analysis is stratified, include the subgroup column named in `.env`, such as `sex`.

### One Workbook Per Sample

```csv
sample_id,workbook_path,treatment,sex
M_VEH_01,workbooks/M_VEH_01.xlsx,vehicle,male
M_VEH_02,workbooks/M_VEH_02.xlsx,vehicle,male
M_ALDO_01,workbooks/M_ALDO_01.xlsx,aldosterone,male
F_VEH_01,workbooks/F_VEH_01.xlsx,vehicle,female
F_ALDO_01,workbooks/F_ALDO_01.xlsx,aldosterone,female
```

### One Workbook With Multiple Sample Sheets

Use `sheet_name` when multiple samples are sheets in the same workbook:

```csv
sample_id,workbook_path,sheet_name,treatment,sex
M_VEH_01,workbooks/cytokine_array_041326.xlsx,MV1,vehicle,male
M_VEH_02,workbooks/cytokine_array_041326.xlsx,MV2,vehicle,male
M_ALDO_01,workbooks/cytokine_array_041326.xlsx,MA1,aldosterone,male
F_VEH_01,workbooks/cytokine_array_041326.xlsx,FV1,vehicle,female
F_ALDO_01,workbooks/cytokine_array_041326.xlsx,FA1,aldosterone,female
```

Only listed sheets are read. Summary tabs such as `ALL data` or `Benjamini Hochberg` are ignored unless they are listed in the manifest.

## Protocol Workbook

The R workflow needs a protocol workbook that maps array coordinates to analytes.

If you have the vendor protocol PDF, extract the workbook with:

```bash
python3 scripts/setup/extract_analyte_table.py --preset cytokine_xl
```

Then point `.env` to the workbook:

```bash
PROTEOME_PROFILER_PROTOCOL_WORKBOOK=output/cytoXL array kit - protocol.xlsx
```

## Outputs

Outputs are written under:

```text
PROTEOME_PROFILER_RUNTIME_ROOT/output/<user>/<slug>/
```

Common folders:

```text
threshold_diagnostics/
main_analysis/
inferential_results/
select_analytes/
```

`main_analysis/` is used for exploratory analyses without enough biological replicates for inferential testing. It visualizes raw-signal fold changes after the configured low-signal reference-panel threshold is applied.

`inferential_results/` is used for replicate-aware analyses. Open `comparison_workbook.xlsx` first. Its first sheets summarize the run, input/QC flags, method-level counts, and significance counts before the comparison-specific result sheets.

Detailed replicate-aware artifacts:

```text
inferential_results/
  comparison_workbook.xlsx
  raw_log2_lm_results.xlsx
  normalized_t_test_results.xlsx
  methods_overview.md
  run_index.tsv
  comparisons/
    <comparison_slug>/
      tables/<method>_results.tsv
      waterfall_plots/<method>/<method>_waterfall*.png
      barplots/<method>/all_tested/<method>_barplot_all_tested_page_<n>.png
      barplots/<method>/significant_hits/<threshold>/<method>_barplot_<threshold>_page_<n>.png
```

Artifact roles:

- `comparison_workbook.xlsx` is the first human-readable replicate-aware result.
- `input_qc/reference_spot_qc.tsv` lists, for each biological sample, the reference-spot denominator, reference rows, raw reference signals, and QC status used to compute `normalized_signal` for `normalized_t_test`.
- `<method>_results.xlsx` files contain full detailed result tables for one method.
- `methods_overview.md` explains the method estimands, fold-change definitions, and significance flags.
- `run_index.tsv` is the machine-readable path and provenance index.
- `comparisons/<comparison_slug>/` contains comparison-scoped method tables and plots.

`select_analytes/` is created by `scripts/select-analytes-analysis.R` when selected analytes are configured.

## Methods Text

### With Biological Replicates

Proteome Profiler array signals are processed in R by averaging duplicate membrane spots within each sample, mapping spot coordinates to analytes with a protocol-derived workbook, and performing per-analyte treatment-versus-control inference within each analysis stratum. The `raw_log2_lm` method fits `log2(signal) ~ treatment`, where `signal` is the averaged duplicate raw analyte signal for one biological sample. The treatment coefficient estimates the difference in mean log2 raw signal between arms and is reported as both a log2 effect and `2^(coefficient)` fold change.

The `normalized_t_test` method recomputes a per-sample, per-analyte `normalized_signal` from raw signals before testing. For one biological sample and one analyte, `normalized_signal = averaged duplicate raw analyte signal / sample reference-spot denominator`. The reference-spot denominator is the mean raw signal of that sample's preferred Reference Spots pairs `A1,2` and `J1,2` when present; if those preferred pairs are unavailable, the analysis uses available Reference Spots rows from that sample. The analysis records the exact reference rows and raw reference signals used for each sample in `input_qc/reference_spot_qc.tsv`, and `scripts/check_setup.R` reports the reference-spot QC counts before analysis. The method applies a two-sided equal-variance t-test to those per-sample `normalized_signal` values and reports ratio-scale fold change only when both group means are positive. Raw p-values are adjusted within each comparison family using the Benjamini-Hochberg method.

The low-signal threshold from `PROTEOME_PROFILER_REF_SIGNAL` is separate from reference-spot normalization. In replicate-aware results it creates `low_signal_flag` for every configured method, including both `raw_log2_lm` and `normalized_t_test`; it does not rescale values and does not remove analytes before p-value adjustment.

Replicate-aware barplots use a linear fold-change-ratio y-axis relative to the control group. For `normalized_t_test`, the control bar is `1`, the treatment bar is `mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`, the control whisker is `1 +/- SE(control-arm normalized_signal values) / mean(control-arm normalized_signal values)`, and the treatment whisker is `treatment fold change +/- SE(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`. For `raw_log2_lm`, bars are converted from log2-scale group means to the linear fold-change-ratio scale; these whiskers can look asymmetric after conversion back from log2. Barplot whiskers are visual group-mean SE bars, not confidence intervals for the treatment-vs-control contrast. Waterfall whiskers instead show `+/- 1 SE` of the method-specific treatment effect estimate.

### Without Biological Replicates

Proteome Profiler array signals are processed in R by averaging duplicate membrane spots within each membrane, mapping spot coordinates to analytes with a protocol-derived workbook, and visualizing treatment-versus-control fold changes. Because biological replicates are not available, no inferential statistical testing or multiplicity correction is performed.

## Troubleshooting

Run:

```bash
Rscript scripts/check_setup.R
```

This command reports missing packages, missing files, invalid manifest columns, missing workbook sheets, reference-spot denominator checks, and output-folder problems before the analysis runs.
