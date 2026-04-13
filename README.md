# proteome_profiler

Exploratory analysis scripts for membrane-based proteome profiler arrays, with the current working example centered on the VEGFRi/Dox in-vivo mouse project.

The active workflow reads LI-COR spot-intensity exports, averages duplicate membrane spots into one analyte-level signal, normalizes each treatment arm to the control arm, and writes waterfall/barplot outputs for follow-up review.

## Current supported example

The tested handoff path in this repo is:

- Project repo: `/Users/ncamarda/Projects/proteome_profiler`
- Runtime/output root: `/Users/ncamarda/ProjectsRuntime/proteome_profiler`
- Cloud data root: `/Users/ncamarda/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects/proteome_profiler`
- Example dataset: `projects/Veh vs Sor Dox Lis - Cytokine XL`

## What The Analysis Does

For the VEGFRi/Dox cytokine and angiogenesis arrays, the main pipeline currently computes:

- `vehicle` vs `sorafenib`
- `vehicle` vs `sor + dox`
- `vehicle` vs `sor + lis`

This is an exploratory fold-change screen, not a replicated inferential analysis. The `significant` label in output paths means "passes the fold-change threshold", not "statistically significant".

## Quick Start

From the repo root:

```bash
Rscript scripts/install_packages.R
Rscript scripts/find_ref_thresh.R
Rscript scripts/main.R
Rscript scripts/select-analytes-analysis.R
```

What each script does:

- [scripts/install_packages.R](/Users/ncamarda/Projects/proteome_profiler/scripts/install_packages.R): installs required CRAN packages.
- [scripts/find_ref_thresh.R](/Users/ncamarda/Projects/proteome_profiler/scripts/find_ref_thresh.R): helps inspect low-signal regions for the VEGFRi/Dox cytokine example.
- [scripts/main.R](/Users/ncamarda/Projects/proteome_profiler/scripts/main.R): runs the main VEGFRi/Dox cytokine pipeline from the example `.env` file.
- [scripts/select-analytes-analysis.R](/Users/ncamarda/Projects/proteome_profiler/scripts/select-analytes-analysis.R): produces a sorafenib-focused shortlist view for the same dataset.

## Path Overrides

If a collaborator stores the runtime or cloud roots elsewhere, they can override the defaults with:

```bash
export PROTEOME_PROFILER_RUNTIME_ROOT="/path/to/ProjectsRuntime/proteome_profiler"
export PROTEOME_PROFILER_CLOUD_ROOT="/path/to/cloud/proteome_profiler"
export PROTEOME_PROFILER_ENV_FILE="/path/to/dox_cytokine_xl_array.env"
```

The main script defaults to the tested VEGFRi/Dox example env file when `PROTEOME_PROFILER_ENV_FILE` is not set.

## Output Locations

The current VEGFRi/Dox example writes to:

- `output/plots/nick/cytokine_xl_array`
- `output/plots/nick/cytokine_xl_array/select`

Under the default path configuration, those resolve to `/Users/ncamarda/ProjectsRuntime/proteome_profiler/output/...`.

## Important Caveats

- This project is best described as exploratory membrane-array profiling used to nominate analytes for follow-up validation.
- The workflow compares each treatment arm to `vehicle`; it does not currently estimate combo-versus-sorafenib rescue effects directly.
- The code assumes LI-COR spot exports remain in the expected `S001...` order; the repo now validates that assumption and errors if the workbook order changes.
- The low-signal filter is comparison-specific in the current code, so one treatment arm no longer suppresses another arm's analytes during pairwise comparisons.
