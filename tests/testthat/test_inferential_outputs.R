test_that("inferential outputs include per-comparison files, index, and workbook", {
    tmp_dir <- tempfile("inferential-output-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")
    output_dir <- file.path(tmp_dir, "inferential_results")

    write_protocol_fixture(protocol_path)
    create_replicate_fixture_dir(data_dir)
    create_replicate_manifest(manifest_path, data_dir)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    manifest_tbl <- read_sample_manifest(manifest_path, subgroup_var = "sex", treatment_var = "treatment") %>%
        resolve_manifest_workbook_paths(data_dir = data_dir)

    sample_df <- build_sample_level_dataset(
        manifest = manifest_tbl,
        analyte_info = analyte_info,
        treatment_var = "treatment",
        subgroup_var = "sex",
        data_dir = data_dir
    )

    results <- run_within_stratum_differential_analysis(
        sample_data = sample_df,
        comparisons = list("control" = c("treated")),
        subgroup_var = "sex",
        p_adjust_method = "BH",
        alpha = 0.05,
        min_reps = 2,
        low_signal_threshold = 70
    )

    out_paths <- write_inferential_outputs(results, output_dir)

    expect_true(file.exists(out_paths$run_index))
    expect_true(file.exists(out_paths$workbook))
    expect_true(all(file.exists(out_paths$results)))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_treated", "tables", "raw_log2_lm_results.tsv")))
    expect_false(file.exists(file.path(output_dir, "comparisons", "male_control_vs_treated", "results.tsv")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_treated", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_treated", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall.png")))

    male_results <- readr::read_tsv(
        file.path(output_dir, "comparisons", "male_control_vs_treated", "tables", "raw_log2_lm_results.tsv"),
        show_col_types = FALSE
    )
    expect_true(all(c("fdr_lt_0_20", "fdr_lt_0_25") %in% names(male_results)))
})

test_that("inferential outputs include explicit raw-p and FDR waterfall tiers when hits exist", {
    output_dir <- file.path(tempfile("inferential-sig-"), "inferential_results")
    dir.create(output_dir, recursive = TRUE)

    sample_df <- create_complex_replicate_sample_data()
    results <- run_within_stratum_differential_analysis(
        sample_data = sample_df,
        comparisons = list("control" = c("drug_a", "drug_b")),
        subgroup_var = "sex",
        p_adjust_method = "BH",
        alpha = 0.05,
        min_reps = 2,
        low_signal_threshold = 70
    )

    male_hits <- results$results[["male_control_vs_drug_a"]] %>%
        filter(raw_p_lt_alpha)
    barplot_data <- build_inferential_fold_change_barplot_data(male_hits)
    expect_equal(levels(barplot_data$group), c("control", "drug_a"))
    expect_equal(levels(barplot_data$short_group), c("control", "drug_a"))

    write_inferential_outputs(results, output_dir)

    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_raw_p_lt_alpha.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_fdr_lt_0_25.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "all_tested", "raw_log2_lm_barplot_all_tested_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "significant_hits", "raw_p_lt_alpha", "raw_log2_lm_barplot_raw_p_lt_alpha_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "significant_hits", "fdr_lt_0_20", "raw_log2_lm_barplot_fdr_lt_0_20_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "significant_hits", "fdr_lt_0_25", "raw_log2_lm_barplot_fdr_lt_0_25_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_raw_p_lt_alpha.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_fdr_lt_0_25.png")))
})

test_that("multi-method inferential outputs write method-specific workbooks plus comparison outputs", {
    output_dir <- file.path(tempfile("inferential-multi-"), "inferential_results")
    dir.create(output_dir, recursive = TRUE)

    sample_df <- create_complex_replicate_sample_data()
    results <- run_replicate_analysis_methods(
        sample_data = sample_df,
        comparisons = list("control" = c("drug_a")),
        subgroup_var = "sex",
        analysis_methods = c("raw_log2_lm", "normalized_t_test"),
        p_adjust_method = "BH",
        alpha = 0.05,
        min_reps = 2,
        low_signal_threshold = 70
    )

    walk(results$methods, function(method_result) {
        walk(method_result$results, function(comparison_result) {
            plotted_rows <- comparison_result %>%
                filter(test_status == "tested", is.finite(effect_estimate_log2))
            expect_true("effect_se_log2" %in% names(plotted_rows))
            expect_true(all(is.finite(plotted_rows$effect_se_log2)))
        })
    })

    out_paths <- write_multi_method_inferential_outputs(results, output_dir)

    expect_true(file.exists(file.path(output_dir, "run_index.tsv")))
    expect_true(file.exists(file.path(output_dir, "raw_log2_lm_results.xlsx")))
    expect_true(file.exists(file.path(output_dir, "normalized_t_test_results.xlsx")))
    expect_true(file.exists(file.path(output_dir, "comparison_workbook.xlsx")))
    expect_true(file.exists(file.path(output_dir, "methods_overview.md")))
    expect_false(file.exists(file.path(output_dir, "run_summary.xlsx")))
    expect_false(dir.exists(file.path(output_dir, "run_report")))
    expect_false(dir.exists(file.path(output_dir, "summary_report")))

    run_index <- readr::read_tsv(file.path(output_dir, "run_index.tsv"), show_col_types = FALSE)
    expect_equal(sort(unique(run_index$analysis_method)), c("normalized_t_test", "raw_log2_lm"))
    expect_equal(nrow(run_index), 4)
    expect_true(all(!is.na(run_index$result_path)))
    expect_true(all(file.exists(as.character(run_index$result_path))))
    expect_true(all(!is.na(run_index$full_waterfall_path)))
    expect_true(all(file.exists(as.character(run_index$full_waterfall_path))))

    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "tables", "raw_log2_lm_results.tsv")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "tables", "normalized_t_test_results.tsv")))
    expect_false(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "results.tsv")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "normalized_t_test", "normalized_t_test_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "raw_log2_lm", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "waterfall_plots", "normalized_t_test", "normalized_t_test_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "all_tested", "raw_log2_lm_barplot_all_tested_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "normalized_t_test", "all_tested", "normalized_t_test_barplot_all_tested_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "raw_log2_lm", "significant_hits", "fdr_lt_0_20", "raw_log2_lm_barplot_fdr_lt_0_20_page_1.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "barplots", "normalized_t_test", "significant_hits", "fdr_lt_0_20", "normalized_t_test_barplot_fdr_lt_0_20_page_1.png")))
    expect_true(file.exists(out_paths$primary$run_index))
    expect_true(all(file.exists(unname(out_paths$workbooks))))
    expect_true(file.exists(out_paths$comparison$comparison_workbook))

    normalized_sheet_names <- openxlsx::getSheetNames(file.path(output_dir, "normalized_t_test_results.xlsx"))
    expect_true("run_index" %in% normalized_sheet_names)
    expect_true(any(normalized_sheet_names != "run_index"))
    normalized_workbook_index <- readxl::read_excel(file.path(output_dir, "normalized_t_test_results.xlsx"), sheet = "run_index")
    expect_true(all(normalized_workbook_index$analysis_method == "normalized_t_test"))
    expect_true(all(!is.na(normalized_workbook_index$result_path)))

    comparison_workbook_path <- file.path(output_dir, "comparison_workbook.xlsx")
    comparison_sheet_names <- openxlsx::getSheetNames(comparison_workbook_path)
    expect_equal(
        comparison_sheet_names[1:4],
        c("summary", "input_qc_summary", "method_summary", "significance_summary")
    )
    expect_false("selected_analytes_summary" %in% comparison_sheet_names)

    comparison_summary <- readxl::read_excel(comparison_workbook_path, sheet = "summary")
    input_qc_summary <- readxl::read_excel(comparison_workbook_path, sheet = "input_qc_summary")
    method_summary <- readxl::read_excel(comparison_workbook_path, sheet = "method_summary")
    significance_summary <- readxl::read_excel(comparison_workbook_path, sheet = "significance_summary")

    expect_equal(comparison_summary$n_methods[[1]], 2)
    expect_equal(comparison_summary$n_comparisons[[1]], 2)
    expect_equal(comparison_summary$n_method_comparison_rows[[1]], 4)
    expect_equal(nrow(input_qc_summary), 4)
    expect_true(all(c("n_low_signal_flagged", "n_not_testable", "n_low_replication_warning") %in% names(input_qc_summary)))
    expect_equal(sort(method_summary$analysis_method), c("normalized_t_test", "raw_log2_lm"))
    expect_equal(nrow(significance_summary), 4)
    expect_true(all(c("n_tested", "n_raw_p_lt_alpha", "n_fdr_lt_0_20", "n_fdr_lt_0_25") %in% names(significance_summary)))

    raw_written <- readr::read_tsv(
        file.path(output_dir, "comparisons", "male_control_vs_drug_a", "tables", "raw_log2_lm_results.tsv"),
        show_col_types = FALSE
    )
    raw_memory <- results$methods$raw_log2_lm$results$male_control_vs_drug_a
    raw_written_row <- raw_written %>% filter(Name == "Analyte Strong A")
    raw_memory_row <- raw_memory %>% filter(Name == "Analyte Strong A")
    expect_equal(raw_written_row$raw_p_value, raw_memory_row$raw_p_value, tolerance = 1e-12)
    expect_equal(raw_written_row$adjusted_p_value, raw_memory_row$adjusted_p_value, tolerance = 1e-12)
    expect_equal(raw_written_row$fold_change_ratio, raw_memory_row$fold_change_ratio, tolerance = 1e-12)
    expect_equal(raw_written_row$fdr_lt_0_20, raw_memory_row$fdr_lt_0_20)

    normalized_written <- readr::read_tsv(
        file.path(output_dir, "comparisons", "male_control_vs_drug_a", "tables", "normalized_t_test_results.tsv"),
        show_col_types = FALSE
    )
    normalized_memory <- results$methods$normalized_t_test$results$male_control_vs_drug_a
    normalized_written_row <- normalized_written %>% filter(Name == "Analyte Strong A")
    normalized_memory_row <- normalized_memory %>% filter(Name == "Analyte Strong A")
    expect_equal(normalized_written_row$raw_p_value, normalized_memory_row$raw_p_value, tolerance = 1e-12)
    expect_equal(normalized_written_row$effect_se_log2, normalized_memory_row$effect_se_log2, tolerance = 1e-12)
    expect_equal(normalized_written_row$fold_change_ratio, normalized_memory_row$fold_change_ratio, tolerance = 1e-12)

    expect_false(dir.exists(file.path(output_dir, "methods")))
    expect_false(dir.exists(file.path(output_dir, "method_comparison")))
    expect_false(file.exists(file.path(output_dir, "method_index.tsv")))
})

test_that("sample-level cleanliness outputs summarize nonpositive rows", {
    output_dir <- file.path(tempfile("input-qc-"), "input_qc")
    dir.create(output_dir, recursive = TRUE)

    sample_df <- create_complex_replicate_sample_data()
    qc_paths <- write_sample_level_cleanliness_outputs(sample_df, output_dir)

    expect_true(file.exists(qc_paths$summary))
    expect_true(file.exists(qc_paths$sample_summary))
    expect_true(file.exists(qc_paths$analyte_summary))
    expect_true(file.exists(qc_paths$issue_rows))
    expect_true(file.exists(qc_paths$reference_spot_summary))
    expect_true(file.exists(qc_paths$reference_spot_qc))

    summary_tbl <- readr::read_tsv(qc_paths$summary, show_col_types = FALSE)
    issue_tbl <- readr::read_tsv(qc_paths$issue_rows, show_col_types = FALSE)
    reference_summary_tbl <- readr::read_tsv(qc_paths$reference_spot_summary, show_col_types = FALSE)
    reference_qc_tbl <- readr::read_tsv(qc_paths$reference_spot_qc, show_col_types = FALSE)

    expect_gt(summary_tbl$n_nonpositive_signal[[1]], 0)
    expect_true(all(issue_tbl$signal_issue %in% c("negative", "zero", "missing", "non_finite")))
    expect_equal(reference_summary_tbl$n_samples[[1]], dplyr::n_distinct(sample_df$sample_id))
    expect_true(all(c(
        "normalization_denominator",
        "normalization_reference_status",
        "reference_qc_issue"
    ) %in% names(reference_qc_tbl)))
    expect_true(all(reference_qc_tbl$reference_qc_issue == "ok"))
})
