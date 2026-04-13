# Proteome Profiler Array Analysis

Exploratory analysis scripts for membrane-based proteome profiler arrays, plus an OpenSpec plan for a future replicate-aware inferential workflow.

The active workflow reads LI-COR spot-intensity exports, averages duplicate membrane spots into one analyte-level signal, normalizes each treatment arm to the control arm, and writes waterfall/barplot outputs for follow-up review.

If you want a docs entry point for the output examples, start with [docs/README.md](/Users/ncamarda/Projects/proteome_profiler/docs/README.md).

## Example Outputs

These are example PNGs copied from one real run into `docs/output-examples/` so they render directly in the repo. They illustrate the output structure and plotting style only; do not treat the specific analytes, groups, or thresholds shown here as defaults for a new project.

### Threshold Diagnostics

This is the diagnostic plot from `scripts/find_ref_thresh.R`. It is used to choose a practical raw-signal cutoff for low-intensity analytes.

![Threshold diagnostics example](docs/output-examples/threshold-diagnostics-region-stats.png)

To use this figure:

1. First review your own data across the currently configured groups and choose a provisional low-signal analyte panel. `scripts/find_ref_thresh.R` now writes `threshold_diagnostics/candidate_low_signal_analytes.tsv` as a programmatic starting point, but you still need to decide whether those analytes make biological and assay-level sense for your dataset.
2. Write the coordinates you chose into `ref_coords_to_make_filter` in [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R).
3. Re-run `Rscript scripts/find_ref_thresh.R` if you changed that panel.
4. Read the mean signal shown in the figure caption and the magenta vertical line in the lower-left histogram.
5. Round that mean to a practical integer cutoff and write it into `ref_thresh_to_filter` in [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R).
6. Re-run `Rscript scripts/main.R` so the main workflow filters each control-vs-treatment pair to analytes where both groups have raw signal above that cutoff.

If you want to compare multiple candidate cutoffs, set `ref_thresh_to_filter` to a vector such as `c(150, 200)`. `scripts/main.R` will generate one output subtree per threshold.

### Main Analysis Waterfall

This is one example `main_analysis/` waterfall from a control-versus-treatment comparison before restricting to fold-change hits.

![Main analysis waterfall example](docs/output-examples/main-waterfall-sorafenib.png)

### Main Analysis Combined Barplots

This is one faceted `main_analysis/` barplot page showing analytes that passed the configured fold-change threshold.

![Main analysis combined barplot example](docs/output-examples/main-barplot-combined-page-1.png)

### Select-Analytes Waterfall

This is an example focused shortlist waterfall from `scripts/select-analytes-analysis.R` for one configured comparison.

![Select analytes waterfall example](docs/output-examples/select-analytes-waterfall.png)

## Example Configuration

The user-editable config lives in [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R). Edit that file once for your machine and project layout. The active R scripts source that file explicitly at startup, and `scripts/helpers/project_paths.R` should usually not need to be edited.

The example entry in that file is illustrative only, not a required directory structure.

Example layout:

- Project repo: `/Users/ncamarda/Projects/proteome_profiler`
- Runtime/output root: `/Users/ncamarda/ProjectsRuntime/proteome_profiler`
- Cloud parent directory: `/Users/ncamarda/Library/CloudStorage/OneDrive-Personal/path/to/your/project-parent`
- Example dataset under that parent: `proteome_profiler/projects/example-array-project`

## Workflow Modes

There are two workflows to understand:

1. Implemented now: legacy exploratory mode
   This is the current one-workbook-per-group workflow already present in the repo.
2. Specified next: replicate-aware inferential mode
   This is the planned manifest-driven workflow described in `openspec/changes/replicate-aware-fdr-analysis/`. It is not implemented yet.

## What The Current Analysis Does

For a given config entry, the main pipeline currently computes the control-versus-treatment comparisons listed in that entry’s `comparisons` field.

This is an exploratory fold-change screen, not a replicated inferential analysis. The `significant` label in output paths means "passes the fold-change threshold", not "statistically significant".

## Quick Start

From the repo root:

```bash
python3 scripts/setup/extract_analyte_table.py --preset cytokine_xl
Rscript scripts/install_packages.R
Rscript scripts/find_ref_thresh.R
Rscript scripts/main.R
Rscript scripts/select-analytes-analysis.R
```

The Python extraction step comes first. The R scripts depend on the generated workbook at `output/cytoXL array kit - protocol.xlsx`.

Execution order and dependencies:

1. `python3 scripts/setup/extract_analyte_table.py --preset ...`
   This creates the analyte workbook that maps protocol positions to gene/protein names. The R scripts need that mapping to label the LI-COR spots correctly.
2. `Rscript scripts/install_packages.R`
   This is setup only and does not produce analysis output. It bootstraps `pak` if needed and then installs missing R packages with `pak::pkg_install()`.
3. `Rscript scripts/find_ref_thresh.R`
   Run this before `main.R` to inspect the manually chosen low-signal analyte panel and choose the raw-signal filter threshold for the active analysis config.
4. `Rscript scripts/main.R`
   This does not read the `find_ref_thresh.R` output file directly. It uses the threshold values already written into `scripts/config/analysis_config.R`, so `find_ref_thresh.R` is a decision aid, not a hard runtime dependency.
5. `Rscript scripts/select-analytes-analysis.R`
   Run this after `main.R` or after you have settled on the threshold settings you want, because it uses the same comparison and threshold configuration to create a focused shortlist view.

What each script does:

- [scripts/install_packages.R](/Users/ncamarda/Projects/proteome_profiler/scripts/install_packages.R): installs required R packages with `pak`.
- [scripts/find_ref_thresh.R](/Users/ncamarda/Projects/proteome_profiler/scripts/find_ref_thresh.R): reads the LI-COR workbooks plus the extracted protocol workbook, writes `threshold_diagnostics/candidate_low_signal_analytes.tsv` as a suggested starting panel, then summarizes the configured low-signal analyte panel in `threshold_diagnostics/region_stats.png`. Use that TSV to choose `ref_coords_to_make_filter`, then use the plot mean as the candidate cutoff to write into `ref_thresh_to_filter` in the config.
- [scripts/main.R](/Users/ncamarda/Projects/proteome_profiler/scripts/main.R): reads the same inputs, then runs the configured treatment-vs-control comparisons across the configured thresholds and writes the full waterfall/barplot output tree to `main_analysis/`.
- [scripts/select-analytes-analysis.R](/Users/ncamarda/Projects/proteome_profiler/scripts/select-analytes-analysis.R): reuses the same dataset and threshold configuration, but narrows the analysis to the configured shortlist comparison and writes those focused outputs to `select_analytes/`.

## Current Step-By-Step Workflow

This section describes the workflow that is implemented in the repo today.

1. Run `python3 scripts/setup/extract_analyte_table.py --preset ...` to create the protocol workbook that maps array coordinates to analyte names.
2. Edit [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R) and choose the analysis entry, group labels, comparisons, raw-signal thresholds, and fold-change thresholds.
3. Put exactly one LI-COR workbook per analysis group into the configured `data_dir`.
4. Run `Rscript scripts/find_ref_thresh.R`, review `threshold_diagnostics/candidate_low_signal_analytes.tsv`, and decide which analytes should define `ref_coords_to_make_filter` for your dataset.
5. Re-run `Rscript scripts/find_ref_thresh.R` if you changed `ref_coords_to_make_filter`, then update `ref_thresh_to_filter` in [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R) using the mean signal shown in `threshold_diagnostics/region_stats.png` as your candidate raw-signal cutoff.
6. Run `Rscript scripts/main.R` to build the analyte-by-group dataset, apply the configured raw-signal cutoff, compute treatment-vs-control fold changes, and write waterfall/barplot outputs.
7. Run `Rscript scripts/select-analytes-analysis.R` if you want a focused shortlist view for the configured selection comparison.

## How To Set `ref_thresh_to_filter`

`scripts/find_ref_thresh.R` does not edit the config for you. It writes a diagnostic figure so you can choose the threshold manually.

Use `threshold_diagnostics/region_stats.png` like this:

1. Run `Rscript scripts/find_ref_thresh.R` and inspect `threshold_diagnostics/candidate_low_signal_analytes.tsv`.
2. Review those suggested analytes against your data across the configured groups and decide which coordinates belong in `ref_coords_to_make_filter` in [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R).
3. Re-run `Rscript scripts/find_ref_thresh.R` after changing `ref_coords_to_make_filter`.
4. Open `threshold_diagnostics/region_stats.png`.
5. In the top panel, check whether the chosen analytes really cluster in the low-intensity range across groups. If they do not, revise `ref_coords_to_make_filter` and regenerate the diagnostic plot.
6. In the lower-left histogram, use the magenta vertical line and the reported mean signal as the candidate raw-signal cutoff.
7. Round that value to a practical integer and set `ref_thresh_to_filter = c(<value>)`.
8. If you want to compare sensitivity to stricter gating, provide multiple values such as `c(150, 200)`. `scripts/main.R` will create `main_analysis/ref_threshold_<value>/` outputs for each threshold.

In the current implementation, `main.R` keeps an analyte in a pairwise comparison only when both the control arm and the treatment arm have raw signal above the chosen `ref_thresh_to_filter` value.

## Planned Replicate-Aware Workflow

This section describes the manifest-driven workflow specified in [openspec/changes/replicate-aware-fdr-analysis/design.md](/Users/ncamarda/Projects/proteome_profiler/openspec/changes/replicate-aware-fdr-analysis/design.md). It is not implemented yet, but this is the intended design contract.

### Input Contract

In replicate-aware mode, the experimental design lives in a manifest file, not in the workbook filenames.

Required manifest columns:

- `sample_id`
- `workbook_path`
- `treatment`

Required additional columns depend on the analysis. For example, a subgrouped analysis may also need:

- `sex`

Example manifest:

```csv
sample_id,workbook_path,treatment,sex
M01,data/m01.xlsx,control,male
M02,data/m02.xlsx,control,male
M03,data/m03.xlsx,treated,male
M04,data/m04.xlsx,treated,male
F01,data/f01.xlsx,control,female
F02,data/f02.xlsx,control,female
F03,data/f03.xlsx,treated,female
F04,data/f04.xlsx,treated,female
```

In this mode:

- workbook filenames can be arbitrary
- the manifest is the source of truth for sample identity and experimental design
- the config chooses the subgroup column, treatment column, comparisons, alpha, and p-value adjustment method

Recommended replicate-aware config shape:

```r
list(
  user = "analyst_name",
  analysis_slug = "example_array_analysis",
  sample_manifest = "manifests/example_samples.csv",
  protocol_preset = "cytokine_xl",
  info_fn = "output/cytoXL array kit - protocol.xlsx",
  subgroup_var = "sex",
  treatment_var = "treatment",
  comparisons = list("control" = c("treated")),
  ref_thresh_to_filter = c(150),
  p_adjust_method = "BH",
  alpha = 0.05
)
```

### Planned Analysis Flow

1. Load the analysis config.
2. If `sample_manifest` is present, enter replicate-aware mode. If not, and there is exactly one workbook per group, stay in the legacy exploratory mode.
3. Validate the manifest columns, unique `sample_id` values, subgroup column presence, and workbook paths.
4. Read one LI-COR workbook per biological sample from the manifest.
5. Average duplicate membrane spots within each workbook to create one analyte-level value per sample.
6. Join the protocol workbook so array coordinates map to analyte names.
7. Build one long-format dataset with one row per `analyte x sample_id`.
8. For each subgroup level, run the configured treatment-vs-control comparison on log2-transformed signals.
9. Compute raw p-values, then adjust them within each comparison family using the configured method, default `BH`.
10. Record low-signal status as a separate flag rather than dropping analytes by default before testing.
11. Write one canonical results file per comparison, one run-level index, and one combined workbook for collaborator-friendly review.

## Methods Templates

Use these as short starting paragraphs for a paper methods section. Edit the protocol name, thresholds, and comparison labels to match the actual run.

### With Biological Replicates

Proteome Profiler array signals were processed in R by averaging duplicate membrane spots within each sample, mapping spot coordinates to analytes with a protocol-derived workbook, and modeling log2-transformed analyte signals for treatment-versus-control comparisons within each analysis stratum. Raw p-values were adjusted within each comparison family using the Benjamini-Hochberg method by default, and waterfall plots plus multiplicity-corrected results tables were generated. Low-signal analytes were flagged relative to the configured raw-intensity reference cutoff rather than excluded from testing by default.

### Without Biological Replicates

Proteome Profiler array signals were processed in R by averaging duplicate membrane spots within each membrane, mapping spot coordinates to analytes with a protocol-derived workbook, and visualizing treatment-versus-control fold changes with waterfall and barplot outputs. A dataset-specific raw-intensity cutoff was chosen from a manually selected low-signal analyte panel, and fold-change thresholds were used for exploratory hit nomination. Because biological replicates were not available, no inferential statistical testing or multiplicity correction was performed.

## Script Layout

- `scripts/main.R`, `scripts/find_ref_thresh.R`, `scripts/select-analytes-analysis.R`, `scripts/install_packages.R`: active entrypoints.
- `scripts/config/analysis_config.R`: the one config file users are expected to edit.
- `scripts/helpers/`: shared path-resolution, package loading, color setup, and R helper code.
- `scripts/setup/extract_analyte_table.py`: setup utility for turning protocol PDFs into analyte workbooks.
- `scripts/legacy/`: older project-specific and superseded scripts kept for reference.

## Protocol Extraction

The Python extractor is now a setup utility instead of a hard-coded one-off script.

Run it before the R scripts for a new machine or a new protocol workbook.

Examples:

```bash
python3 scripts/setup/extract_analyte_table.py --preset cytokine_xl
python3 scripts/setup/extract_analyte_table.py --preset angiogenesis
```

You can also pass explicit paths:

```bash
python3 scripts/setup/extract_analyte_table.py \
  --input "protocols/cytoXL array kit - protocol.pdf" \
  --output "output/cytoXL array kit - protocol.xlsx" \
  --pages 17,18,19
```

## Path Setup

Edit [scripts/config/analysis_config.R](/Users/ncamarda/Projects/proteome_profiler/scripts/config/analysis_config.R) and set:

- `runtime_root`
- `cloud_parent`
- the example or project-specific analysis entry you want to run

For each analysis entry, the important fields are:

- `user`: the analyst or collaborator who owns the output subtree, for example `analyst_name`
- `analysis_slug`: the folder name for that analysis under the user subtree
- `data_dir`: raw LI-COR export directory
- `info_fn`: generated analyte workbook from the Python extraction step
- thresholds, groups, and comparisons used by the R scripts
- `ref_coords_to_make_filter`: a manually chosen set of low-signal analyte coordinates used by `scripts/find_ref_thresh.R` to estimate a practical raw-signal floor. These are example low-signal analytes, not necessarily the protocol's true control spots.

That one script is the only place you need to change paths or analysis metadata for a new machine or project layout.

## Output Locations

Outputs are organized by user first, then by analysis:

```text
output/plots/<user>/<analysis_slug>/
```

For one example config, that root might be:

```text
output/plots/analyst_name/example_array_analysis/
```

Under the default path configuration, that resolves under `/Users/ncamarda/ProjectsRuntime/proteome_profiler/output/...`.

Inside that analysis root:

- `threshold_diagnostics/`: outputs from `scripts/find_ref_thresh.R`
- `main_analysis/`: outputs from `scripts/main.R`
- `select_analytes/`: outputs from `scripts/select-analytes-analysis.R`

Inside `main_analysis/`, folders are organized by threshold choice:

```text
main_analysis/
  ref_threshold_<value>/
    all_comparisons/
      waterfalls/
        <control comparison>/
    fold_change_hits/
      threshold_<value>/
        waterfalls/
          <control comparison>/
        barplots/
          <control comparison>/
            combined/
            <treatment comparison>/
              combined/
              individual/
```

`fold_change_hits` means "passes the configured fold-change cutoff." It does not mean statistical significance.

## Important Caveats

- This project is best described as exploratory membrane-array profiling used to nominate analytes for follow-up validation.
- The workflow compares each treatment arm to its configured control group; it does not currently estimate more complex contrast structures unless they are implemented explicitly.
- The pipeline currently assumes one LI-COR workbook per analysis group. It averages duplicate spots within a membrane, but it does not yet support multiple biological replicate workbooks per group.
- The code assumes LI-COR spot exports remain in the expected `S001...` order; the repo now validates that assumption and errors if the workbook order changes.
- The low-signal filter is comparison-specific in the current code, so one treatment arm no longer suppresses another arm's analytes during pairwise comparisons.
