# proteome_profiler

Exploratory analysis scripts for membrane-based proteome profiler arrays.

The active workflow reads LI-COR spot-intensity exports, averages duplicate membrane spots into one analyte-level signal, normalizes each treatment arm to the control arm, and writes waterfall/barplot outputs for follow-up review.

## Current supported example

The shared config lives in [scripts/project_paths.R](/Users/ncamarda/Projects/proteome_profiler/scripts/project_paths.R). Edit that file once for your machine and project layout.

The VEGFRi/Dox entry in that file is an example configuration, not a required directory structure.

Example layout:

- Project repo: `/Users/ncamarda/Projects/proteome_profiler`
- Runtime/output root: `/Users/ncamarda/ProjectsRuntime/proteome_profiler`
- Cloud parent directory: `/Users/ncamarda/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects`
- Example dataset under that parent: `proteome_profiler/projects/Veh vs Sor Dox Lis - Cytokine XL`

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
- [scripts/main.R](/Users/ncamarda/Projects/proteome_profiler/scripts/main.R): runs the main pipeline using the example config in `scripts/project_paths.R`.
- [scripts/select-analytes-analysis.R](/Users/ncamarda/Projects/proteome_profiler/scripts/select-analytes-analysis.R): produces a sorafenib-focused shortlist view for the same dataset.

## Path Setup

Edit [scripts/project_paths.R](/Users/ncamarda/Projects/proteome_profiler/scripts/project_paths.R) and set:

- `runtime_root`
- `cloud_parent`
- the example or project-specific analysis entry you want to run

That one script is the only place you need to change paths for a new machine or project layout.

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
