# Preflight

## Scope

- Change: `replicate-aware-fdr-analysis`
- Review date: `2026-04-13`
- Revision under review: `6695924` plus uncommitted local changes in the working tree

## Path Truth

- Code root: `~/Projects/proteome_profiler`
- Runtime root: `proteome_profiler_config$runtime_root`
- Cloud fallback root: `proteome_profiler_config$cloud_parent/proteome_profiler`
- Path resolution order: `repo -> runtime_root -> cloud project root`
- New outputs are written under `runtime_root/output/plots/<user>/<analysis_slug>/`

## Active Entrypoints

- `scripts/main.R`: primary orchestrator; supports manifest-driven inputs for both exploratory and replicate-aware analyses
- `scripts/find_ref_thresh.R`: threshold diagnostic for both exploratory and replicate-aware analyses using the manifest-driven input model
- `scripts/select-analytes-analysis.R`: dual-mode shortlist writer with support for multiple inferential shortlist comparisons
- `scripts/install_packages.R`: `pak`-based dependency installer

## Primary Evidence Reviewed

- Source-of-truth config: `scripts/config/analysis_config.R`
- Path/config helpers: `scripts/helpers/project_paths.R`
- Replicate-aware implementation: `scripts/helpers/replicate_analysis.R`
- Legacy helper layer: `scripts/helpers/array_helper_scripts.R`
- Runtime docs: `README.md`, `docs/README.md`
- R tests: `tests/testthat/*`
- Python tests: `tests/python/test_extract_analyte_table.py`

## Verification Status

- `Rscript -e "testthat::test_dir('tests/testthat', reporter='summary')"`: passed
- `python3 -m pytest -q tests/python/test_extract_analyte_table.py`: passed

## Key Observations

- Input declaration is now unified around manifest-driven configuration; exploratory versus inferential behavior is decided by `mode`, not by whether a manifest exists.
- The new larger replicate-aware fixture directly exercises multiple comparison families, BH/FDR behavior that changes calls, uneven subgroup sizes, and partial missingness.
- The final review gate is still based on code and synthetic/fixture evidence rather than a live collaborator biological-replicate runtime, but the missing robustness cases from the earlier blocked review are now covered in tests.
