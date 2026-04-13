# Main script for running the proteome profiler analysis and generating graphs.
# Nicholas Camarda
# 2026-06-17

source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = TRUE)
configure_progress_handlers()

# helper scripts
config_path <- Sys.getenv("PROTEOME_PROFILER_CONFIG", unset = file.path("scripts", "config", "analysis_config.R"))
# Keep the user-editable config source explicit in the entry script while still
# allowing tests or collaborators to point at an alternate config file.
source(path.expand(config_path))
source(file.path("scripts", "helpers", "replicate_analysis.R"))
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
analysis_name <- get_selected_analysis_name(Sys.getenv("PROTEOME_PROFILER_ANALYSIS", unset = ""))
example_config <- get_analysis_config(analysis_name)
analysis_mode <- get_analysis_mode(example_config)

info_fn <- get_protocol_workbook_path(example_config)
analysis_output_root <- get_analysis_output_root(example_config)

# Pull the run knobs out of the config once so the body of the script reads as
# a straightforward orchestration layer rather than nested list indexing.
my_ref_thresh_to_filter <- example_config$ref_thresh_to_filter
my_comparisons <- example_config$comparisons
message(qq("Using analysis config: @{analysis_name}"))
message(qq("Resolved protocol table: @{info_fn}"))
message(qq("Resolved analysis output root: @{analysis_output_root}"))

analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

manifest_tbl <- NULL
sample_level_df <- NULL
manifest_path <- NULL
data_dir <- NULL

if (config_uses_sample_manifest(example_config)) {
    manifest_path <- resolve_project_path(example_config$sample_manifest, must_exist = TRUE)
    data_dir <- if (!is.null(example_config$data_dir) && !identical(example_config$data_dir, "")) {
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
}

if (analysis_mode == "replicate") {
    output_dir <- file.path(analysis_output_root, "inferential_results")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    message(qq("Resolved sample manifest: @{manifest_path}"))
    message(qq("Resolved inferential output dir: @{output_dir}"))

    low_signal_threshold <- if (length(my_ref_thresh_to_filter) > 0) my_ref_thresh_to_filter[[1]] else NA_real_
    if (length(my_ref_thresh_to_filter) > 1) {
        message("Replicate-aware inferential mode uses the first ref_thresh_to_filter value as the low-signal flag threshold.")
    }

    inferential_results <- run_within_stratum_differential_analysis(
        sample_data = sample_level_df,
        comparisons = my_comparisons,
        subgroup_var = example_config$subgroup_var,
        p_adjust_method = example_config$p_adjust_method,
        alpha = example_config$alpha,
        min_reps = ifelse(is.null(example_config$min_reps_per_arm), 2, example_config$min_reps_per_arm),
        low_signal_threshold = low_signal_threshold
    )

    write_inferential_outputs(
        inferential_results = inferential_results,
        output_dir = output_dir
    )
} else {
    if (!config_uses_sample_manifest(example_config)) {
        data_dir <- resolve_project_path(example_config$data_dir, must_exist = TRUE)
    }
    output_dir <- file.path(analysis_output_root, "main_analysis")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    my_main_threshold <- example_config$main_threshold
    my_groups_per_page <- example_config$groups_per_page
    exploratory_group_levels <- if (config_uses_sample_manifest(example_config)) {
        build_threshold_diagnostic_dataset(
            sample_data = sample_level_df,
            subgroup_var = example_config$subgroup_var
        )$group_levels
    } else {
        example_config$group_levels
    }
    style_config <- build_group_style(exploratory_group_levels, scheme = "main")
    my_group_lvls <- style_config$group_levels

    if (!is_perfect_square(n = my_groups_per_page)) {
        message(qq("Warning: @{my_groups_per_page} is not a perfect square. Facet plot results may be ugly.\n"))
    }

    my_colors <- style_config$fill
    my_outline_colors <- style_config$outline

    if (config_uses_sample_manifest(example_config)) {
        message(qq("Resolved exploratory sample manifest: @{manifest_path}"))
    } else {
        message(qq("Resolved data dir: @{data_dir}"))
    }
    message(qq("Resolved main-analysis output dir: @{output_dir}"))

    # Convert the raw LI-COR spot exports into one analyte-level dataframe that can
    # be reused across all threshold combinations below.
    my_initial_ready_df <- if (config_uses_sample_manifest(example_config)) {
        # Exploratory manifest-driven runs still collapse to one analyte-by-condition
        # table before plotting; the manifest only replaces filename-based input
        # declaration, not the legacy plotting semantics.
        build_threshold_diagnostic_dataset(
            sample_data = sample_level_df,
            subgroup_var = example_config$subgroup_var
        )$wide_df
    } else {
        make_plot_ready_dataset(
            data_dir = data_dir,
            analyte_info = analyte_info_df,
            preview = FALSE,
            my_group_lvls = my_group_lvls
        )
    }

    # Parallelism is applied across threshold combinations, not across the raw-data
    # parsing step above, because the dataset only needs to be built once.
    num_cores <- max(1, min(parallel::detectCores() - 1, length(my_main_threshold)))
    message(qq("Switching to multisession mode. Using @{num_cores} cores..."))
    plan(multisession, workers = num_cores)

    total_iterations <- length(my_ref_thresh_to_filter) * length(my_main_threshold)

    with_progress({
        if (config_uses_sample_manifest(example_config)) {
            message(qq("\nStarting exploratory manifest-driven processing of @{manifest_path}"))
        } else {
            message(qq("\nStarting processing of @{data_dir}"))
        }
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
}
