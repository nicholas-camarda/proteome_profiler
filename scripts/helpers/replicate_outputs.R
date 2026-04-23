# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

#' Write a cross-method comparison table for each inferential comparison
#'
#' @param method_results Named list returned under `$methods` from
#'   `run_replicate_analysis_methods()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written file paths.
build_workbook_sheet_names <- function(labels) {
    used_names <- character()

    vapply(seq_along(labels), function(i) {
        base_name <- labels[[i]] %>%
            str_replace_all("[\\\\/?*\\[\\]:]", "_") %>%
            str_replace_all("\\s+", "_") %>%
            str_sub(1, 31)

        if (base_name == "") {
            base_name <- sprintf("sheet_%02d", i)
        }

        candidate <- base_name
        suffix <- 1
        while (candidate %in% used_names) {
            suffix_text <- sprintf("_%02d", suffix)
            candidate <- str_sub(base_name, 1, 31 - nchar(suffix_text))
            candidate <- paste0(candidate, suffix_text)
            suffix <- suffix + 1
        }

        used_names <<- c(used_names, candidate)
        candidate
    }, character(1))
}

#' Write inferential result tables and plots to disk
#'
#' @param inferential_results Output of `run_within_stratum_differential_analysis()`.
#' @param output_dir Directory where inferential outputs should be written.
#'
#' @return Invisibly returns the output paths written.
build_inferential_workbook_filename <- function(run_index) {
    workbook_method <- run_index %>%
        distinct(analysis_method) %>%
        pull(analysis_method) %>%
        as.character()

    if (length(workbook_method) == 1 && !is.na(workbook_method[[1]]) && workbook_method[[1]] != "") {
        sprintf("%s_results.xlsx", workbook_method[[1]])
    } else {
        "results_workbook.xlsx"
    }
}

#' Return the table directory for one inferential comparison
#'
#' @param output_dir Root inferential output directory.
#' @param comparison_slug Comparison slug.
#'
#' @return Directory path.
inferential_comparison_tables_dir <- function(output_dir, comparison_slug) {
    file.path(output_dir, "comparisons", comparison_slug, "tables")
}

#' Return the waterfall directory for one comparison and method
#'
#' @param output_dir Root inferential output directory.
#' @param comparison_slug Comparison slug.
#' @param method_name Analysis method slug.
#'
#' @return Directory path.
inferential_comparison_waterfall_dir <- function(output_dir, comparison_slug, method_name) {
    file.path(output_dir, "comparisons", comparison_slug, "waterfall_plots", method_name)
}

#' Return the barplot directory for one comparison and method
#'
#' @param output_dir Root inferential output directory.
#' @param comparison_slug Comparison slug.
#' @param method_name Analysis method slug.
#'
#' @return Directory path.
inferential_comparison_barplot_dir <- function(output_dir, comparison_slug, method_name) {
    file.path(output_dir, "comparisons", comparison_slug, "barplots", method_name)
}

#' Collapse existing output paths for run-index cells
#'
#' @param paths Character vector of candidate paths.
#'
#' @return Semicolon-delimited existing paths, or `NA`.
collapse_existing_output_paths <- function(paths) {
    paths <- paths[file.exists(paths)]
    if (length(paths) == 0) {
        return(NA_character_)
    }
    paste(sort(paths), collapse = ";")
}

#' Build path columns for the multi-method inferential run index
#'
#' @param run_index Tibble with comparison and method rows.
#' @param output_dir Root inferential output directory.
#'
#' @return Tibble with output path columns.
build_inferential_output_paths_index <- function(run_index, output_dir) {
    if (nrow(run_index) == 0) {
        return(tibble())
    }

    run_index %>%
        transmute(
            comparison_slug,
            analysis_method
        ) %>%
        pmap_dfr(function(comparison_slug, analysis_method) {
            tables_dir <- inferential_comparison_tables_dir(output_dir, comparison_slug)
            waterfall_dir <- inferential_comparison_waterfall_dir(output_dir, comparison_slug, analysis_method)
            barplot_dir <- inferential_comparison_barplot_dir(output_dir, comparison_slug, analysis_method)

            tibble(
                comparison_slug = comparison_slug,
                analysis_method = analysis_method,
                result_path = collapse_existing_output_paths(file.path(tables_dir, sprintf("%s_results.tsv", analysis_method))),
                full_waterfall_path = collapse_existing_output_paths(file.path(waterfall_dir, sprintf("%s_waterfall.png", analysis_method))),
                raw_p_alpha_waterfall_path = collapse_existing_output_paths(file.path(waterfall_dir, sprintf("%s_waterfall_raw_p_lt_alpha.png", analysis_method))),
                fdr_0_20_waterfall_path = collapse_existing_output_paths(file.path(waterfall_dir, sprintf("%s_waterfall_fdr_lt_0_20.png", analysis_method))),
                fdr_0_25_waterfall_path = collapse_existing_output_paths(file.path(waterfall_dir, sprintf("%s_waterfall_fdr_lt_0_25.png", analysis_method))),
                all_tested_barplot_path = collapse_existing_output_paths(list.files(
                    file.path(barplot_dir, "all_tested"),
                    pattern = sprintf("^%s_barplot_all_tested_page_[0-9]+\\.png$", analysis_method),
                    full.names = TRUE
                )),
                raw_p_alpha_barplot_path = collapse_existing_output_paths(list.files(
                    file.path(barplot_dir, "significant_hits", "raw_p_lt_alpha"),
                    pattern = sprintf("^%s_barplot_raw_p_lt_alpha_page_[0-9]+\\.png$", analysis_method),
                    full.names = TRUE
                )),
                fdr_0_20_barplot_path = collapse_existing_output_paths(list.files(
                    file.path(barplot_dir, "significant_hits", "fdr_lt_0_20"),
                    pattern = sprintf("^%s_barplot_fdr_lt_0_20_page_[0-9]+\\.png$", analysis_method),
                    full.names = TRUE
                )),
                fdr_0_25_barplot_path = collapse_existing_output_paths(list.files(
                    file.path(barplot_dir, "significant_hits", "fdr_lt_0_25"),
                    pattern = sprintf("^%s_barplot_fdr_lt_0_25_page_[0-9]+\\.png$", analysis_method),
                    full.names = TRUE
                ))
            )
        })
}

#' Remove stale inferential comparison artifacts before writing current outputs
#'
#' @param comparison_dir Comparison output directory.
#' @param remove_organized_dirs Logical; remove current organized subfolders too.
#'
#' @return Invisibly returns `NULL`.
clean_inferential_comparison_dir <- function(comparison_dir, remove_organized_dirs = FALSE) {
    stale_flat_files <- list.files(
        comparison_dir,
        pattern = "^(results|raw_log2_lm_results|normalized_t_test_results|waterfall|raw_log2_lm_waterfall|normalized_t_test_waterfall|raw_log2_lm_barplot|normalized_t_test_barplot).*\\.(tsv|png)$",
        full.names = TRUE
    )
    stale_exploratory_dirs <- list.files(
        comparison_dir,
        pattern = "^(raw_log2_lm|normalized_t_test)_bargraphs_",
        full.names = TRUE
    )
    stale_paths <- c(stale_flat_files, stale_exploratory_dirs)
    if (remove_organized_dirs) {
        stale_paths <- c(stale_paths, file.path(comparison_dir, c("tables", "waterfall_plots", "barplots")))
    }
    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))
}

#' Write one method-specific inferential workbook
#'
#' @param inferential_results Output from one inferential method.
#' @param output_dir Root inferential output directory.
#' @param run_index Run-index rows for this method.
#'
#' @return Workbook path.
write_inferential_workbook <- function(inferential_results, output_dir, run_index = inferential_results$run_index) {
    workbook_sheets <- inferential_results$results
    workbook_sheets$run_index <- run_index

    workbook_path <- file.path(output_dir, build_inferential_workbook_filename(run_index))
    writexl::write_xlsx(workbook_sheets, workbook_path)
    workbook_path
}

#' Write all outputs for one inferential method
#'
#' @param inferential_results Output of `run_within_stratum_differential_analysis()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written output paths.
write_inferential_outputs <- function(inferential_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stale_workbook_paths <- file.path(output_dir, c(
        "combined_results.xlsx",
        "results_workbook.xlsx",
        "raw_log2_lm_results.xlsx",
        "normalized_t_test_results.xlsx"
    ))
    invisible(walk(stale_workbook_paths, function(path) {
        if (file.exists(path)) {
            unlink(path, force = TRUE)
            if (file.exists(path)) {
                stop(sprintf(
                    paste(
                        "Could not remove stale workbook '%s'.",
                        "Close the file in Excel or another viewer and rerun the analysis to avoid mixing old and new outputs."
                    ),
                    path
                ))
            }
        }
    }))

    output_rows <- map_dfr(names(inferential_results$results), function(comparison_slug) {
        comparison_dir <- file.path(output_dir, "comparisons", comparison_slug)
        dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
        clean_inferential_comparison_dir(comparison_dir, remove_organized_dirs = TRUE)

        result_tbl <- inferential_results$results[[comparison_slug]]
        method_name <- unique(result_tbl$analysis_method)
        if (length(method_name) != 1 || is.na(method_name) || identical(method_name, "")) {
            stop(sprintf(
                "Expected exactly one analysis_method in comparison '%s' when writing inferential outputs.",
                comparison_slug
            ))
        }

        tables_dir <- inferential_comparison_tables_dir(output_dir, comparison_slug)
        waterfall_dir <- inferential_comparison_waterfall_dir(output_dir, comparison_slug, method_name)
        barplot_dir <- inferential_comparison_barplot_dir(output_dir, comparison_slug, method_name)
        dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(waterfall_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)

        result_path <- file.path(tables_dir, sprintf("%s_results.tsv", method_name))
        write_tsv(result_tbl, result_path)

        full_waterfall_path <- file.path(waterfall_dir, sprintf("%s_waterfall.png", method_name))
        waterfall_plot <- plot_inferential_waterfall(
            result_tbl,
            title = sprintf("%s\n%s", unique(result_tbl$comparison_label), unique(result_tbl$analysis_method_label))
        )
        ggsave(
            filename = full_waterfall_path,
            plot = waterfall_plot,
            width = max(12, 0.4 * nrow(result_tbl)),
            height = 10
        )

        threshold_paths <- write_threshold_waterfall_set(
            result_tbl = result_tbl,
            output_dir = waterfall_dir,
            filename_prefix = method_name,
            title_suffix = unique(result_tbl$analysis_method_label)
        )
        bargraph_paths <- write_threshold_bargraph_set(
            result_tbl = result_tbl,
            output_dir = barplot_dir,
            filename_prefix = method_name
        )

        tibble(
            comparison_slug = comparison_slug,
            analysis_method = method_name,
            result_path = result_path,
            full_waterfall_path = full_waterfall_path,
            raw_p_alpha_waterfall_path = threshold_paths$raw_p_lt_alpha,
            fdr_0_20_waterfall_path = threshold_paths$fdr_lt_0_20,
            fdr_0_25_waterfall_path = threshold_paths$fdr_lt_0_25,
            all_tested_barplot_path = bargraph_paths$all_tested,
            raw_p_alpha_barplot_path = bargraph_paths$raw_p_lt_alpha,
            fdr_0_20_barplot_path = bargraph_paths$fdr_lt_0_20,
            fdr_0_25_barplot_path = bargraph_paths$fdr_lt_0_25
        )
    })

    run_index <- inferential_results$run_index %>%
        left_join(output_rows, by = c("comparison_slug", "analysis_method"))
    run_index_path <- file.path(output_dir, "run_index.tsv")
    write_tsv(run_index, run_index_path)

    workbook_path <- write_inferential_workbook(
        inferential_results = inferential_results,
        output_dir = output_dir,
        run_index = run_index
    )

    invisible(list(
        run_index = run_index_path,
        workbook = workbook_path,
        results = output_rows$result_path
    ))
}

#' Write non-primary method comparison tables and plots
#'
#' @param method_results Named list of method results.
#' @param output_dir Root inferential output directory.
#' @param primary_method Method already written by `write_inferential_outputs()`.
#'
#' @return Invisibly returns `NULL`.
write_method_specific_waterfalls <- function(method_results, output_dir, primary_method) {
    comparison_slugs <- method_results[[primary_method]]$run_index %>%
        pull(comparison_slug) %>%
        as.character()

    invisible(walk(comparison_slugs, function(comparison_slug) {
        comparison_dir <- file.path(output_dir, "comparisons", comparison_slug)
        dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
        clean_inferential_comparison_dir(comparison_dir)

        invisible(walk(names(method_results), function(method_name) {
            if (identical(method_name, primary_method)) {
                return(invisible(NULL))
            }

            result_tbl <- method_results[[method_name]]$results[[comparison_slug]]
            tables_dir <- inferential_comparison_tables_dir(output_dir, comparison_slug)
            waterfall_dir <- inferential_comparison_waterfall_dir(output_dir, comparison_slug, method_name)
            barplot_dir <- inferential_comparison_barplot_dir(output_dir, comparison_slug, method_name)
            dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(waterfall_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)

            result_path <- file.path(tables_dir, sprintf("%s_results.tsv", method_name))
            write_tsv(result_tbl, result_path)
            method_plot_path <- file.path(waterfall_dir, sprintf("%s_waterfall.png", method_name))
            method_plot <- plot_inferential_waterfall(
                result_tbl,
                title = sprintf("%s\n%s", unique(result_tbl$comparison_label), unique(result_tbl$analysis_method_label))
            )
            ggsave(
                filename = method_plot_path,
                plot = method_plot,
                width = max(12, 0.4 * nrow(result_tbl)),
                height = 10
            )

            write_threshold_waterfall_set(
                result_tbl = result_tbl,
                output_dir = waterfall_dir,
                filename_prefix = method_name,
                title_suffix = unique(result_tbl$analysis_method_label)
            )
            write_threshold_bargraph_set(
                result_tbl = result_tbl,
                output_dir = barplot_dir,
                filename_prefix = method_name
            )
        }))
    }))
}

#' Write the collaborator-facing cross-method comparison workbook
#'
#' @param method_results Named list returned under `$methods` from
#'   `run_replicate_analysis_methods()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written workbook path.
write_method_comparison_outputs <- function(method_results, output_dir) {
    combined_long <- map_dfr(names(method_results), function(method) {
        imap_dfr(method_results[[method]]$results, function(result_tbl, comparison_slug) {
            result_tbl %>%
                mutate(comparison_slug = comparison_slug)
        })
    })

    comparison_slugs <- combined_long %>%
        distinct(comparison_slug) %>%
        pull(comparison_slug) %>%
        as.character()

    comparison_tables <- map(comparison_slugs, function(comparison_slug) {
        comparison_long <- combined_long %>%
            filter(comparison_slug == .env$comparison_slug)

        comparison_long %>%
            select(
                Name,
                Coordinate,
                analysis_method,
                analysis_method_label,
                raw_p_value,
                raw_p_lt_alpha,
                adjusted_p_value,
                n_p_adjust_hypotheses,
                fdr_lt_0_20,
                fdr_lt_0_25,
                effect_estimate_log2,
                effect_se_log2,
                fold_change_ratio,
                fold_change_status,
                fold_change_note,
                test_status,
                control_n,
                treatment_n,
                low_signal_flag
            ) %>%
            pivot_wider(
                names_from = analysis_method,
                values_from = c(
                    analysis_method_label,
                    raw_p_value,
                    raw_p_lt_alpha,
                    adjusted_p_value,
                    n_p_adjust_hypotheses,
                    fdr_lt_0_20,
                    fdr_lt_0_25,
                    effect_estimate_log2,
                    effect_se_log2,
                    fold_change_ratio,
                    fold_change_status,
                    fold_change_note,
                    test_status,
                    control_n,
                    treatment_n,
                    low_signal_flag
                ),
                names_glue = "{analysis_method}_{.value}"
            ) %>%
            arrange(Name, Coordinate)
    })
    names(comparison_tables) <- comparison_slugs

    comparison_index <- map_dfr(comparison_slugs, function(comparison_slug) {
        comparison_long <- combined_long %>%
            filter(comparison_slug == .env$comparison_slug)

        comparison_long %>%
            group_by(comparison_slug, subgroup, control, treatment, analysis_method, analysis_method_label) %>%
            summarize(
                n_tested = sum(test_status == "tested", na.rm = TRUE),
                n_raw_p_below_0_05 = sum(raw_p_value < 0.05, na.rm = TRUE),
                n_raw_p_lt_alpha = sum(raw_p_lt_alpha, na.rm = TRUE),
                n_fdr_lt_0_20 = sum(fdr_lt_0_20, na.rm = TRUE),
                n_fdr_lt_0_25 = sum(fdr_lt_0_25, na.rm = TRUE),
                best_raw_p_value = suppressWarnings(min(raw_p_value, na.rm = TRUE)),
                .groups = "drop"
            ) %>%
            mutate(best_raw_p_value = if_else(is.infinite(best_raw_p_value), NA_real_, best_raw_p_value))
    })

    run_summary <- tibble(
        n_methods = n_distinct(combined_long$analysis_method),
        n_comparisons = n_distinct(combined_long$comparison_slug),
        n_method_comparison_rows = nrow(comparison_index),
        n_analyte_result_rows = nrow(combined_long),
        n_tested = sum(combined_long$test_status == "tested", na.rm = TRUE),
        n_not_testable = sum(combined_long$test_status != "tested", na.rm = TRUE),
        total_raw_p_lt_alpha = sum(combined_long$raw_p_lt_alpha, na.rm = TRUE),
        total_fdr_lt_0_20 = sum(combined_long$fdr_lt_0_20, na.rm = TRUE),
        total_fdr_lt_0_25 = sum(combined_long$fdr_lt_0_25, na.rm = TRUE)
    )

    input_qc_summary <- combined_long %>%
        group_by(comparison_slug, subgroup, control, treatment, analysis_method, analysis_method_label) %>%
        summarize(
            n_analyte_rows = n(),
            n_tested = sum(test_status == "tested", na.rm = TRUE),
            n_not_testable = sum(test_status != "tested", na.rm = TRUE),
            n_low_signal_flagged = sum(low_signal_flag, na.rm = TRUE),
            n_low_replication_warning = sum(low_replication_warning, na.rm = TRUE),
            n_missing_fold_change = sum(!is.finite(fold_change_ratio), na.rm = TRUE),
            n_tested_missing_effect_se = sum(test_status == "tested" & !is.finite(effect_se_log2), na.rm = TRUE),
            min_control_n = suppressWarnings(min(control_n, na.rm = TRUE)),
            min_treatment_n = suppressWarnings(min(treatment_n, na.rm = TRUE)),
            .groups = "drop"
        ) %>%
        mutate(
            min_control_n = if_else(is.infinite(min_control_n), NA_real_, min_control_n),
            min_treatment_n = if_else(is.infinite(min_treatment_n), NA_real_, min_treatment_n)
        ) %>%
        arrange(comparison_slug, analysis_method)

    method_summary <- combined_long %>%
        group_by(analysis_method, analysis_method_label) %>%
        summarize(
            n_comparisons = n_distinct(comparison_slug),
            n_analyte_rows = n(),
            n_tested = sum(test_status == "tested", na.rm = TRUE),
            n_not_testable = sum(test_status != "tested", na.rm = TRUE),
            total_raw_p_lt_alpha = sum(raw_p_lt_alpha, na.rm = TRUE),
            total_fdr_lt_0_20 = sum(fdr_lt_0_20, na.rm = TRUE),
            total_fdr_lt_0_25 = sum(fdr_lt_0_25, na.rm = TRUE),
            best_raw_p_value = suppressWarnings(min(raw_p_value, na.rm = TRUE)),
            best_adjusted_p_value = suppressWarnings(min(adjusted_p_value, na.rm = TRUE)),
            .groups = "drop"
        ) %>%
        mutate(
            best_raw_p_value = if_else(is.infinite(best_raw_p_value), NA_real_, best_raw_p_value),
            best_adjusted_p_value = if_else(is.infinite(best_adjusted_p_value), NA_real_, best_adjusted_p_value)
        ) %>%
        arrange(analysis_method)

    comparison_workbook_path <- file.path(output_dir, "comparison_workbook.xlsx")
    workbook_sheets <- c(
        list(
            summary = run_summary,
            input_qc_summary = input_qc_summary,
            method_summary = method_summary,
            significance_summary = comparison_index
        ),
        stats::setNames(comparison_tables, build_workbook_sheet_names(comparison_slugs))
    )
    writexl::write_xlsx(workbook_sheets, comparison_workbook_path)

    invisible(list(
        comparison_workbook = comparison_workbook_path
    ))
}

#' Write or refresh selected-analyte summary content in the comparison workbook
#'
#' @param inferential_dir Root `inferential_results/` directory containing
#'   `comparison_workbook.xlsx`.
#' @param selected_summary_tbl Tibble summarizing selected-analyte outputs.
#'
#' @return Path to the refreshed comparison workbook.
write_selected_analyte_workbook_summary <- function(inferential_dir, selected_summary_tbl) {
    comparison_workbook_path <- file.path(inferential_dir, "comparison_workbook.xlsx")
    if (!file.exists(comparison_workbook_path)) {
        stop(sprintf(
            paste(
                "Cannot write selected-analyte workbook summary because comparison_workbook.xlsx is missing:",
                "%s",
                "Run scripts/main.R before scripts/select-analytes-analysis.R."
            ),
            comparison_workbook_path
        ), call. = FALSE)
    }

    existing_sheet_names <- readxl::excel_sheets(comparison_workbook_path)
    existing_sheets <- setNames(
        map(existing_sheet_names, function(sheet_name) {
            readxl::read_excel(comparison_workbook_path, sheet = sheet_name)
        }),
        existing_sheet_names
    )
    existing_sheets <- existing_sheets[names(existing_sheets) != "selected_analytes_summary"]

    selected_summary_tbl <- selected_summary_tbl %>%
        arrange(comparison_slug, analysis_method, Name)

    insert_after <- match("significance_summary", names(existing_sheets))
    if (is.na(insert_after)) {
        insert_after <- match("summary", names(existing_sheets))
    }
    if (is.na(insert_after)) {
        insert_after <- 0
    }

    refreshed_sheets <- append(
        existing_sheets,
        list(selected_analytes_summary = selected_summary_tbl),
        after = insert_after
    )
    writexl::write_xlsx(refreshed_sheets, comparison_workbook_path)
    comparison_workbook_path
}

#' Write a collaborator-facing method overview for dual-method runs
#'
#' @param analysis_methods Character vector of method slugs.
#' @param output_dir Root inferential output directory.
#'
#' @return Absolute path to the written markdown file.
write_inferential_method_overview <- function(analysis_methods, output_dir) {
    overview_path <- file.path(output_dir, "methods_overview.md")
    workbook_names <- sprintf("`%s_results.xlsx`", analysis_methods)
    workbook_text <- paste(workbook_names, collapse = ", ")

    lines <- c(
        "# Inferential Methods Overview",
        "",
        "This file is a short guide to what each configured inferential method means and how to read the replicate-aware results.",
        "",
        "## Start Here",
        "",
        "- Open `comparison_workbook.xlsx` first. It is the main collaborator-facing summary across all configured methods.",
        sprintf("- The method-specific workbooks are %s. Open them only when you want one method's full tables in workbook form.", workbook_text),
        "- Use `run_index.tsv` when you need exact file locations for comparison-specific tables, waterfall plots, or barplots.",
        "- Use `input_qc/reference_spot_qc.tsv` when you want to inspect the per-sample reference-spot denominators used by `normalized_t_test`.",
        "",
        "## How To Read The Plots",
        "",
        "- Waterfall plots rank analytes by the method-specific effect estimate. Waterfall whiskers are +/- 1 standard error for that plotted effect estimate.",
        "- `all_tested` barplots show every analyte with a finite, ratio-scale-interpretable fold change for that method.",
        "- `raw_p_lt_alpha`, `fdr_lt_0_20`, and `fdr_lt_0_25` barplots and waterfall plots are filtered views of the same results table, not separate analyses.",
        "- Barplots use a linear fold-change-ratio y-axis relative to the control group.",
        "- Barplot whiskers are visual group-mean SE bars on the plotted fold-change scale, not confidence intervals for the treatment-vs-control contrast. The p-values come from the method's treatment-versus-control test, not from whether the whiskers overlap.",
        "- Barplot page captions define the star threshold, for example `* = FDR < 0.20`, next to the page `N`.",
        "- `low_signal_flag` is computed from the configured raw-signal low-signal threshold and is reported for every method. It is a QC flag, not an exclusion before p-value adjustment.",
        "",
        "## Compare Methods Carefully",
        "",
        "- Fold-change values are method-specific and should not be compared across methods as if they were the same quantity.",
        "- `raw_log2_lm` and `normalized_t_test` answer different questions, so different p-values and different fold changes are expected.",
        "- If ratio-scale fold change is not interpretable, the results table marks that row with `fold_change_status` and `fold_change_note`.",
        "- BH correction uses only the analytes that produced a valid p-value within that comparison family.",
        "",
        "## Shared Definitions",
        "",
        "- `normalized_signal` means one sample/analyte value computed as averaged duplicate raw analyte signal divided by that sample's reference-spot denominator.",
        "- For `normalized_t_test`, the denominator uses preferred Reference Spots pairs `A1,2` and `J1,2` when present, otherwise available Reference Spots rows from that sample.",
        "- For `normalized_t_test`, the control bar is `1`; the treatment bar is `mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`.",
        "- For `raw_log2_lm`, bar heights come from log2-scale group means converted back to the linear ratio scale by taking `2^x` of the plotted log2 values, which turns additive log2 differences back into ordinary fold-change ratios (so `1` means equal to control, `2` means twofold above control, and `0.5` means half of control). Because of that back-transformation, whiskers can appear asymmetric on the linear y-axis.",
        ""
    )

    for (method in analysis_methods) {
        method_spec <- get_inferential_method_spec(method)
        lines <- c(
            lines,
            sprintf("## %s (`%s`)", method_spec$label, method),
            "",
            "Statistical methods:",
            sprintf("- %s", method_spec$statistical_methods),
            "",
            "So what this means:",
            sprintf("- %s", method_spec$so_what),
            "",
            "Statistical strengths:",
            sprintf("- %s", method_spec$pros),
            "",
            "Statistical limitations:",
            sprintf("- %s", method_spec$cons),
            ""
        )
    }

    writeLines(lines, overview_path)
    overview_path
}

#' Remove stale multi-method inferential artifacts
#'
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns `NULL`.
clean_multi_method_output_dir <- function(output_dir) {
    stale_paths <- c(
        file.path(output_dir, "methods"),
        file.path(output_dir, "method_comparison"),
        file.path(output_dir, "method_index.tsv"),
        file.path(output_dir, "comparison_index.tsv"),
        file.path(output_dir, "all_methods_long.tsv")
    )

    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))
}

#' Write one or more inferential methods to disk
#'
#' @param inferential_method_results Output of `run_replicate_analysis_methods()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written output paths.
write_multi_method_inferential_outputs <- function(inferential_method_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    clean_multi_method_output_dir(output_dir)

    primary_method <- inferential_method_results$primary_method
    method_results <- inferential_method_results$methods
    analysis_methods <- names(method_results)

    primary_paths <- write_inferential_outputs(
        inferential_results = method_results[[primary_method]],
        output_dir = output_dir
    )
    write_method_specific_waterfalls(
        method_results = method_results,
        output_dir = output_dir,
        primary_method = primary_method
    )
    combined_run_index <- map_dfr(method_results, "run_index") %>%
        left_join(
            build_inferential_output_paths_index(., output_dir),
            by = c("comparison_slug", "analysis_method")
        )
    write_tsv(combined_run_index, file.path(output_dir, "run_index.tsv"))

    method_workbooks <- setNames(
        map_chr(analysis_methods, function(method_name) {
            method_run_index <- combined_run_index %>%
                filter(.data$analysis_method == method_name)

            write_inferential_workbook(
                inferential_results = method_results[[method_name]],
                output_dir = output_dir,
                run_index = method_run_index
            )
        }),
        analysis_methods
    )

    comparison_paths <- write_method_comparison_outputs(
        method_results = method_results,
        output_dir = output_dir
    )
    overview_path <- write_inferential_method_overview(
        analysis_methods = analysis_methods,
        output_dir = output_dir
    )

    invisible(list(
        primary = primary_paths,
        workbooks = method_workbooks,
        comparison = comparison_paths,
        overview = overview_path
    ))
}
