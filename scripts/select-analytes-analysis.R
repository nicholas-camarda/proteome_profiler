rm(list = ls())
source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = FALSE)

source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "replicate_analysis.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number1
#' This script writes explicit selected-analyte views under the per-user
#' analysis tree in `select_analytes/`. Both legacy and replicate-aware modes
#' use `shortlist$analytes` / `selection_analytes` as the source of truth and
#' write one comparison-scoped folder per selected comparison.
initialize_runtime_config_from_env(required_env_file = TRUE)

analysis_name <- get_selected_analysis_name()
example_config <- get_analysis_config(analysis_name)
info_fn <- get_protocol_workbook_path(example_config)
analysis_output_root <- get_analysis_output_root(example_config)
output_dir <- file.path(analysis_output_root, "select_analytes")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
selected_names <- get_selected_analyte_names(example_config)

if (is_replicate_aware_config(example_config)) {
    inferential_dir <- file.path(analysis_output_root, "inferential_results")
    run_index_path <- file.path(inferential_dir, "run_index.tsv")
    shortlist_methods <- resolve_inferential_shortlist_methods(example_config)
    stale_comparisons_dir <- file.path(output_dir, "comparisons")
    if (dir.exists(stale_comparisons_dir)) {
        unlink(stale_comparisons_dir, recursive = TRUE, force = TRUE)
    }

    if (!file.exists(run_index_path)) {
        stop(sprintf(
            paste(
                "Replicate-aware shortlist mode expects inferential outputs from scripts/main.R.",
                "Missing file: %s",
                "Run Rscript scripts/main.R first for analysis '%s'."
            ),
            run_index_path,
            analysis_name
        ))
    }

    run_index <- read_tsv(run_index_path, show_col_types = FALSE)
    selected_comparisons <- resolve_shortlist_comparisons(run_index, example_config)

    for (row_idx in seq_len(nrow(selected_comparisons))) {
        selected_comparison <- selected_comparisons[row_idx, , drop = FALSE]
        comparison_slug <- selected_comparison$comparison_slug[[1]]
        output_subdir <- file.path(output_dir, comparison_slug)
        dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
        clean_replicate_shortlist_comparison_root(output_subdir)

        method_index <- map_dfr(shortlist_methods, function(shortlist_method) {
            comparison_tbl <- read_inferential_comparison_results(
                inferential_dir = inferential_dir,
                comparison_slug = comparison_slug,
                analysis_method = shortlist_method
            )

            method_output_dir <- file.path(output_subdir, shortlist_method)
            dir.create(method_output_dir, recursive = TRUE, showWarnings = FALSE)
            clean_replicate_shortlist_comparison_dir(method_output_dir)

            validate_selected_analyte_names(selected_names, comparison_tbl$Name)
            selected_tbl <- comparison_tbl %>%
                filter(Name %in% selected_names) %>%
                arrange(match(Name, selected_names))

            selected_results_path <- file.path(method_output_dir, "selected_results.tsv")
            write_tsv(selected_tbl, selected_results_path)

            for (optional_column in c("fold_change_ratio", "effect_estimate_log2", "low_signal_flag", "raw_p_lt_alpha", "fdr_lt_0_20", "fdr_lt_0_25")) {
                if (!optional_column %in% names(selected_tbl)) {
                    selected_tbl[[optional_column]] <- NA
                }
            }

            qc_tbl <- selected_tbl %>%
                mutate(
                    plot_status = if_else(is.finite(fold_change_ratio), "plotted", "not_plotted"),
                    no_plot_reason = case_when(
                        !is.finite(fold_change_ratio) ~ "fold_change_ratio is not finite",
                        TRUE ~ NA_character_
                    )
                ) %>%
                select(
                    Name, Coordinate, analysis_method, analysis_method_label,
                    test_status, test_reason, low_signal_flag,
                    fold_change_ratio, effect_estimate_log2,
                    raw_p_value, adjusted_p_value, raw_p_lt_alpha, fdr_lt_0_20, fdr_lt_0_25,
                    plot_status, no_plot_reason
                )
            selected_qc_path <- file.path(method_output_dir, "selected_analyte_qc.tsv")
            write_tsv(qc_tbl, selected_qc_path)

            selected_waterfall_path <- file.path(method_output_dir, "selected_waterfall.png")
            waterfall_tbl <- selected_tbl %>%
                filter(test_status == "tested", is.finite(effect_estimate_log2))
            if (nrow(waterfall_tbl) > 0) {
                selected_waterfall <- plot_inferential_waterfall(
                    waterfall_tbl,
                    title = sprintf(
                        "Selected analytes\n%s\n%s",
                        unique(comparison_tbl$comparison_label),
                        unique(comparison_tbl$analysis_method_label)
                    )
                )
                ggsave(
                    filename = selected_waterfall_path,
                    plot = selected_waterfall,
                    width = max(8, 0.75 * nrow(waterfall_tbl)),
                    height = 8,
                    dpi = 150
                )
            } else if (file.exists(selected_waterfall_path)) {
                unlink(selected_waterfall_path, force = TRUE)
                selected_waterfall_path <- NA_character_
            } else {
                selected_waterfall_path <- NA_character_
            }

            raw_alpha_label <- sprintf("raw p < %.3g", unique(comparison_tbl$alpha)[[1]])
            significance_tbl <- selected_tbl %>%
                transmute(
                    Name,
                    selected_significant = !is.na(raw_p_lt_alpha) & raw_p_lt_alpha,
                    selected_significance_definition = raw_alpha_label
                )
            barplot_data <- build_inferential_fold_change_barplot_data(
                selected_tbl,
                mark_treatment_significant = FALSE
            )
            if (nrow(barplot_data) > 0) {
                barplot_data <- barplot_data %>%
                    left_join(significance_tbl, by = "Name") %>%
                    mutate(
                        significant = as.character(group) == unique(as.character(selected_tbl$treatment))[[1]] & selected_significant,
                        significance_label = if_else(significant, "*", NA_character_),
                        significance_definition = if_else(significant, selected_significance_definition, NA_character_)
                    ) %>%
                    select(-selected_significant, -selected_significance_definition)
            }

            bargraph_index <- write_selected_fold_change_bargraphs(
                barplot_data = barplot_data,
                output_dir = method_output_dir,
                title = sprintf(
                    "Selected analytes\n%s\n%s",
                    unique(comparison_tbl$comparison_label),
                    unique(comparison_tbl$analysis_method_label)
                )
            )
            bargraph_index_path <- file.path(method_output_dir, "selected_bargraph_index.tsv")
            write_tsv(bargraph_index, bargraph_index_path)

            tibble(
                analysis_method = unique(comparison_tbl$analysis_method),
                analysis_method_label = unique(comparison_tbl$analysis_method_label),
                method_output_dir = method_output_dir,
                selected_results_path = selected_results_path,
                selected_qc_path = selected_qc_path,
                selected_waterfall_path = selected_waterfall_path,
                selected_bargraph_index_path = bargraph_index_path
            )
        })

        write_tsv(method_index, file.path(output_subdir, "method_index.tsv"))

        message(qq(
            "Replicate-aware selected analytes written to @{output_subdir} ",
            "(methods=@{paste(shortlist_methods, collapse=', ')})"
        ))
    }
} else {
    data_dir <- if (config_uses_sample_manifest(example_config) ||
        is.null(example_config$data_dir) ||
        identical(example_config$data_dir, "")) {
        NULL
    } else {
        resolve_project_path(example_config$data_dir, must_exist = TRUE)
    }
    my_group_lvls <- example_config$group_levels

    # read in analyte info
    analyte_info_df <- read_excel(info_fn) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)
    validate_selected_analyte_names(selected_names, analyte_info_df$Name)

    df <- if (config_uses_sample_manifest(example_config)) {
        manifest_path <- resolve_project_path(example_config$sample_manifest, must_exist = TRUE)
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
        my_group_lvls <- threshold_diag$group_levels
        threshold_diag$wide_df
    } else {
        make_plot_ready_dataset(
            data_dir = data_dir,
            analyte_info = analyte_info_df,
            preview = TRUE,
            my_group_lvls = my_group_lvls
        )
    }

    style_config <- build_group_style(my_group_lvls, scheme = "main")
    my_group_lvls <- style_config$group_levels
    my_colors <- style_config$fill
    my_outline_colors <- style_config$outline

    my_ref_thresh <- example_config$ref_thresh_to_filter[[1]]
    sor_thresh <- example_config$selection_threshold

    stale_flat_outputs <- c(
        list.files(output_dir, pattern = "^waterfall_plot_.*\\.png$", full.names = TRUE),
        list.files(output_dir, pattern = "^barplot_select_.*\\.png$", full.names = TRUE),
        file.path(output_dir, "select_bargraphs")
    )
    invisible(walk(stale_flat_outputs, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))

    selected_pairs <- resolve_legacy_select_comparisons(example_config)
    for (pair_idx in seq_len(nrow(selected_pairs))) {
        selected_pair <- selected_pairs[pair_idx, , drop = FALSE]
        comparison_slug <- selected_pair$comparison_slug[[1]]
        control_label <- selected_pair$control[[1]]
        treatment_label <- selected_pair$treatment[[1]]
        output_subdir <- file.path(output_dir, comparison_slug)
        if (dir.exists(output_subdir)) {
            unlink(output_subdir, recursive = TRUE, force = TRUE)
        }
        dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)

        selected_lst_obj <- make_wf_data(
            df,
            my_main_threshold = sor_thresh,
            ref_thresh_to_filter_ = -Inf,
            comparisons = setNames(list(c(treatment_label)), control_label)
        )
        result_key <- qq("@{control_label} vs @{treatment_label}")
        selected_wf <- selected_lst_obj[[result_key]]$wf_dat[[treatment_label]] %>%
            filter(Name %in% selected_names) %>%
            arrange(match(Name, selected_names)) %>%
            mutate(order = factor(seq_along(order)))
        selected_dat <- selected_lst_obj[[result_key]]$dat_filtered[[treatment_label]] %>%
            filter(Name %in% selected_names) %>%
            arrange(match(Name, selected_names), group)

        selected_path <- file.path(output_subdir, "selected_analytes.tsv")
        write_tsv(selected_dat, selected_path)

        qc_tbl <- selected_dat %>%
            group_by(Name, Coordinate) %>%
            summarize(
                comparison_slug = comparison_slug,
                control = control_label,
                treatment = treatment_label,
                ref_signal_threshold = my_ref_thresh,
                min_signal = min(signal, na.rm = TRUE),
                passes_ref_threshold = all(signal > my_ref_thresh, na.rm = TRUE),
                low_signal_flag = !passes_ref_threshold,
                finite_relative_signal = all(is.finite(relative_signal)),
                plot_status = if_else(finite_relative_signal, "plotted", "not_plotted"),
                no_plot_reason = if_else(finite_relative_signal, NA_character_, "relative_signal is not finite"),
                .groups = "drop"
            ) %>%
            arrange(match(Name, selected_names))
        selected_qc_path <- file.path(output_subdir, "selected_analyte_qc.tsv")
        write_tsv(qc_tbl, selected_qc_path)

        selected_waterfall_path <- file.path(output_subdir, "selected_waterfall.png")
        if (nrow(selected_wf) > 0) {
            selected_waterfall <- plot_wf_graph(
                data = selected_wf,
                add_fc = FALSE,
                TITLE = qq("Selected analytes\n@{control_label} vs @{treatment_label}\nRef threshold > @{my_ref_thresh} reported in selected_analyte_qc.tsv"),
                main_cutoff = NULL,
                line_color = NA
            )
            ggsave(
                plot = selected_waterfall,
                filename = selected_waterfall_path,
                width = max(8, 0.75 * nrow(selected_wf)),
                height = 8,
                dpi = 150
            )
        }

        barplot_data <- build_legacy_selected_fold_change_barplot_data(
            selected_dat = selected_dat,
            mark_threshold_hits = FALSE,
            threshold = sor_thresh
        )
        bargraph_index <- write_selected_fold_change_bargraphs(
            barplot_data = barplot_data,
            output_dir = output_subdir,
            title = qq("Selected analytes\n@{control_label} vs @{treatment_label}\nRef threshold > @{my_ref_thresh} reported in selected_analyte_qc.tsv")
        )
        write_tsv(bargraph_index, file.path(output_subdir, "selected_bargraph_index.tsv"))

        message(qq("Legacy selected analytes written to @{output_subdir}"))
    }
}
