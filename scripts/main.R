# Main script for running the proteome profiler analysis and generating graphs.
# Nicholas Camarda
# 2026-06-17

source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = TRUE)
configure_progress_handlers()

# helper scripts
source(file.path("scripts", "config", "analysis_config.R"))
source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number0 This code requires that the analytes are measured in LICOR software starting in order of the protocol datasheet.
#' So the S001 corresponds to the very first item on the protocol coordinate system. Background should be measured last, after negative controls

#' @Number1 Parse the protocol data
#' Note: Run scripts/setup/extract_analyte_table.py before the R workflow if the protocol workbook has not been created yet.

# Load the named analysis example from scripts/config/analysis_config.R.
# This object defines:
# - who owns the output subtree (`user`)
# - the analysis folder name under that user (`analysis_slug`)
# - where the raw data and protocol workbook live
# - which groups, comparisons, and thresholds this script should run
example_config <- get_analysis_config("vegfri_dox_cytokine_xl")

info_fn <- get_protocol_workbook_path(example_config)
data_dir <- resolve_project_path(example_config$data_dir, must_exist = TRUE)
analysis_output_root <- get_analysis_output_root(example_config)
output_dir <- file.path(analysis_output_root, "main_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Pull the run knobs out of the config once so the body of the script reads as
# a straightforward orchestration layer rather than nested list indexing.
my_main_threshold <- example_config$main_threshold
my_ref_thresh_to_filter <- example_config$ref_thresh_to_filter
my_comparisons <- example_config$comparisons
my_groups_per_page <- example_config$groups_per_page
style_config <- build_group_style(example_config$group_levels, scheme = "main")
my_group_lvls <- style_config$group_levels

if (!is_perfect_square(n = my_groups_per_page)) {
    message(qq("Warning: @{my_groups_per_page} is not a perfect square. Facet plot results may be ugly.\n"))
}

my_colors <- style_config$fill
my_outline_colors <- style_config$outline

message("Using analysis config: vegfri_dox_cytokine_xl")
message(qq("Resolved data dir: @{data_dir}"))
message(qq("Resolved protocol table: @{info_fn}"))
message(qq("Resolved analysis output root: @{analysis_output_root}"))
message(qq("Resolved main-analysis output dir: @{output_dir}"))

analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

# Convert the raw LI-COR spot exports into one analyte-level dataframe that can
# be reused across all threshold combinations below.
my_initial_ready_df <- make_plot_ready_dataset(
    data_dir = data_dir,
    analyte_info = analyte_info_df,
    preview = FALSE,
    my_group_lvls = my_group_lvls
)

# Parallelism is applied across threshold combinations, not across the raw-data
# parsing step above, because the dataset only needs to be built once.
num_cores <- max(1, min(parallel::detectCores() - 1, length(my_main_threshold)))
message(qq("Switching to multisession mode. Using @{num_cores} cores..."))
plan(multisession, workers = num_cores)

total_iterations <- length(my_ref_thresh_to_filter) * length(my_main_threshold)

with_progress({
    message(qq("\nStarting processing of @{data_dir}"))
    message(qq("Output directory: @{output_dir}"))

    p_inner <- progressor(steps = total_iterations)
    for (i in seq_along(my_ref_thresh_to_filter)) {
        furrr::future_walk(.x = seq_along(my_main_threshold), .f = function(j, p_inner) {
            # Each worker writes one ref-threshold / fold-change-threshold slice
            # of the analysis tree under `main_analysis/`.
            p_inner(qq("Processing ref thresh @{my_ref_thresh_to_filter[i]} and main thresh @{my_main_threshold[j]}"))

            make_graphs(
                df = my_initial_ready_df,
                ref_thresh_to_filter = my_ref_thresh_to_filter[i],
                main_threshold = my_main_threshold[j],
                output_dir_full_path = output_dir,
                comparisons = my_comparisons,
                groups_per_page = my_groups_per_page
            )
        }, p_inner = p_inner, .options = furrr::furrr_options(seed = TRUE))
    }
})
plan(sequential)
