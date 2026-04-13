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
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_treated", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_treated", "raw_log2_lm_waterfall.png")))

    male_results <- readr::read_tsv(
        file.path(output_dir, "comparisons", "male_control_vs_treated", "results.tsv"),
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

    write_inferential_outputs(results, output_dir)

    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall_raw_p_lt_alpha.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall_fdr_lt_0_25.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "raw_log2_lm_waterfall_raw_p_lt_alpha.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_drug_a", "raw_log2_lm_waterfall_fdr_lt_0_25.png")))
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

    out_paths <- write_multi_method_inferential_outputs(results, output_dir)

    expect_true(file.exists(file.path(output_dir, "run_index.tsv")))
    expect_true(file.exists(file.path(output_dir, "raw_log2_lm_results.xlsx")))
    expect_true(file.exists(file.path(output_dir, "normalized_t_test_results.xlsx")))
    expect_true(file.exists(file.path(output_dir, "comparison_workbook.xlsx")))
    expect_true(file.exists(file.path(output_dir, "methods_overview.md")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "normalized_t_test_waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "raw_log2_lm_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_drug_a", "normalized_t_test_waterfall_fdr_lt_0_20.png")))
    expect_true(file.exists(out_paths$primary$run_index))
    expect_true(all(file.exists(unname(out_paths$workbooks))))
    expect_true(file.exists(out_paths$comparison$comparison_workbook))

    normalized_sheet_names <- openxlsx::getSheetNames(file.path(output_dir, "normalized_t_test_results.xlsx"))
    expect_true("run_index" %in% normalized_sheet_names)
    expect_true(any(normalized_sheet_names != "run_index"))

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

    summary_tbl <- readr::read_tsv(qc_paths$summary, show_col_types = FALSE)
    issue_tbl <- readr::read_tsv(qc_paths$issue_rows, show_col_types = FALSE)

    expect_gt(summary_tbl$n_nonpositive_signal[[1]], 0)
    expect_true(all(issue_tbl$signal_issue %in% c("negative", "zero", "missing", "non_finite")))
})
