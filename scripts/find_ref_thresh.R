rm(list = ls())
source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = FALSE)

source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "replicate_analysis.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number1 protocol data
#' @note run scripts/setup/extract_analyte_table.py before this script if the protocol workbook has not been created yet.

initialize_runtime_config_from_env(required_env_file = TRUE)

# Load the shared analysis metadata, then write this script's outputs into the
# per-user analysis tree under `threshold_diagnostics/`.
# The configured coordinates are a manually chosen low-signal analyte panel
# used to estimate a practical raw-signal floor for this dataset.
analysis_name <- get_selected_analysis_name()
example_config <- get_analysis_config(analysis_name)
info_fn <- get_protocol_workbook_path(example_config)
analysis_output_root <- get_analysis_output_root(example_config)
output_dir <- file.path(analysis_output_root, "threshold_diagnostics")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ref_coords_to_make_filter <- example_config$ref_coords_to_make_filter



# read in analyte info
analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

## Build one analyte-level table for threshold diagnostics, regardless of
## whether the active analysis is exploratory one-workbook-per-group data or a
## manifest-driven replicate-aware run.
if (config_uses_sample_manifest(example_config)) {
    manifest_path <- resolve_project_path(example_config$sample_manifest, must_exist = TRUE)
    data_dir <- if (!is.null(example_config$data_dir)) {
        resolve_project_path(example_config$data_dir, must_exist = TRUE)
    } else {
        NULL
    }

    manifest_tbl <- read_sample_manifest(
        manifest_path = manifest_path,
        subgroup_var = example_config$subgroup_var,
        treatment_var = example_config$treatment_var
    ) %>%
        resolve_manifest_workbook_paths(data_dir = data_dir)

    sample_level_df <- build_sample_level_dataset(
        manifest = manifest_tbl,
        analyte_info = analyte_info_df,
        treatment_var = example_config$treatment_var,
        subgroup_var = example_config$subgroup_var,
        data_dir = data_dir
    )

    threshold_diag <- build_threshold_diagnostic_dataset(
        sample_data = sample_level_df,
        subgroup_var = example_config$subgroup_var
    )
    df <- threshold_diag$wide_df
    my_group_lvls <- threshold_diag$group_levels
} else {
    data_dir <- resolve_project_path(example_config$data_dir, must_exist = TRUE)
    my_group_lvls <- example_config$group_levels
    df <- make_plot_ready_dataset(data_dir, analyte_info = analyte_info_df, preview = TRUE, my_group_lvls = my_group_lvls)
}

# Suggest a reviewable starting panel of low-signal analytes so the user can
# choose `ref_coords_to_make_filter` with actual data in hand rather than by
# trial and error alone.
candidate_low_signal_df <- suggest_low_signal_panel(df)
write_tsv(candidate_low_signal_df, file.path(output_dir, "candidate_low_signal_analytes.tsv"))
message("Suggested low-signal analytes written to:")
message(file.path(output_dir, "candidate_low_signal_analytes.tsv"))

long_df <- df %>%
    pivot_longer(all_of(my_group_lvls), names_to = "group", values_to = "signal")

# CHeck this first if you don't want to look through the blots themselves and select regions
# ggplot(long_df, mapping = aes(x = signal)) +
#     geom_histogram(binwidth = 1000) +
#     facet_wrap(~group)


style_config <- build_group_style(my_group_lvls, scheme = "threshold")
my_group_lvls <- style_config$group_levels
my_colors <- style_config$fill


## This writes the diagnostic figure used to pick a reference threshold before
## running the full `main.R` workflow.
filt_lst <- find_filter_thresh(wide_df = df, ref_coords = ref_coords_to_make_filter, my_colors = my_colors)
# filt_lst <- find_filter_thresh(wide_df = df, filter_threshold = 100, my_colors = my_colors)

ggsave(
    filename = file.path(output_dir, "region_stats.png"), plot = filt_lst[[2]],
    width = 10, height = 10
)
