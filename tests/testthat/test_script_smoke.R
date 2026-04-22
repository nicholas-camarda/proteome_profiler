expect_script_success <- function(script_result) {
    expect_equal(
        script_result$status,
        0,
        info = paste(script_result$output, collapse = "\n")
    )
}

test_that("exploratory entry scripts run end-to-end on fixture data", {
    tmp_dir <- tempfile("exploratory-smoke-")
    dir.create(tmp_dir)

    fixture_paths <- write_env_fixture(
        path = file.path(tmp_dir, ".env"),
        runtime_root = file.path(tmp_dir, "runtime"),
        analysis_name = "legacy_smoke"
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(shortlist_result)

    analysis_root <- file.path(fixture_paths$runtime_root, "output", "tester", "legacy_smoke")
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "candidate_low_signal_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "region_stats.png")))
    expect_true(file.exists(file.path(analysis_root, "main_analysis", "ref_threshold_80", "all_comparisons", "waterfalls", "vehicle vs treated", "vehicle vs treated--treated-main_waterfall-all.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "vehicle_vs_treated", "selected_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "vehicle_vs_treated", "selected_analyte_qc.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "vehicle_vs_treated", "selected_waterfall.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "vehicle_vs_treated", "selected_bargraphs", "Analyte A.png")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "select_bargraphs")))
})

test_that("manifest-driven exploratory entry scripts run end-to-end on fixture data", {
    tmp_dir <- tempfile("exploratory-manifest-smoke-")
    dir.create(tmp_dir)

    fixture_paths <- write_env_fixture(
        path = file.path(tmp_dir, ".env"),
        runtime_root = file.path(tmp_dir, "runtime"),
        analysis_name = "legacy_manifest_smoke"
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(shortlist_result)

    analysis_root <- file.path(fixture_paths$runtime_root, "output", "tester", "legacy_manifest_smoke")
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "candidate_low_signal_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "region_stats.png")))
    expect_true(file.exists(file.path(analysis_root, "main_analysis", "ref_threshold_70", "all_comparisons", "waterfalls", "control vs treated", "control vs treated--treated-main_waterfall-all.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "control_vs_treated", "selected_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "control_vs_treated", "selected_analyte_qc.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "control_vs_treated", "selected_waterfall.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "control_vs_treated", "selected_bargraphs", "Analyte A.png")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "select_bargraphs")))
})

test_that("replicate-aware entry scripts run end-to-end on fixture data", {
    tmp_dir <- tempfile("replicate-smoke-")
    dir.create(tmp_dir)

    fixture_paths <- write_env_fixture(
        path = file.path(tmp_dir, ".env"),
        runtime_root = file.path(tmp_dir, "runtime"),
        analysis_name = "replicate_smoke"
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        env_path = fixture_paths$env_path
    )
    expect_script_success(shortlist_result)

    analysis_root <- file.path(fixture_paths$runtime_root, "output", "tester", "replicate_smoke")
    threshold_candidates <- readr::read_tsv(
        file.path(analysis_root, "threshold_diagnostics", "candidate_low_signal_analytes.tsv"),
        show_col_types = FALSE
    )
    run_index <- readr::read_tsv(
        file.path(analysis_root, "inferential_results", "run_index.tsv"),
        show_col_types = FALSE
    )
    method_index <- readr::read_tsv(
        file.path(analysis_root, "select_analytes", "male_control_vs_treated", "method_index.tsv"),
        show_col_types = FALSE
    )
    raw_selected_qc <- readr::read_tsv(
        file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "selected_analyte_qc.tsv"),
        show_col_types = FALSE
    )
    normalized_selected_qc <- readr::read_tsv(
        file.path(analysis_root, "select_analytes", "male_control_vs_treated", "normalized_t_test", "selected_analyte_qc.tsv"),
        show_col_types = FALSE
    )
    comparison_workbook_path <- file.path(analysis_root, "inferential_results", "comparison_workbook.xlsx")

    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "candidate_low_signal_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "region_stats.png")))
    expect_true(nrow(threshold_candidates) > 0)
    expect_true(all(c("Name", "Coordinate", "suggestion_score") %in% names(threshold_candidates)))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "run_index.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "male_control_vs_treated", "tables", "raw_log2_lm_results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "female_control_vs_treated", "tables", "raw_log2_lm_results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "male_control_vs_treated", "tables", "normalized_t_test_results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "female_control_vs_treated", "tables", "normalized_t_test_results.tsv")))
    expect_false(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "male_control_vs_treated", "results.tsv")))
    expect_true(file.exists(comparison_workbook_path))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "method_index.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "selected_results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "selected_analyte_qc.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "selected_waterfall.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "selected_bargraphs", "Analyte A.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "normalized_t_test", "selected_results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "normalized_t_test", "selected_analyte_qc.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "normalized_t_test", "selected_bargraphs", "Analyte A.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "female_control_vs_treated", "normalized_t_test", "selected_results.tsv")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "comparisons")))
    expect_false(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "shortlist.tsv")))
    expect_false(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "shortlist_summary.tsv")))
    expect_false(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "shortlist_waterfall.png")))
    expect_false(file.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "bargraph_index.tsv")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "bargraphs")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_p_lt_alpha")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "raw_p_lt_alpha")))
    expect_false(dir.exists(file.path(analysis_root, "select_analytes", "male_control_vs_treated", "raw_log2_lm", "fdr_lt_0_20")))
    expect_true("result_path" %in% names(run_index))
    expect_equal(sort(unique(run_index$analysis_method)), c("normalized_t_test", "raw_log2_lm"))
    expect_true(all(!is.na(run_index$result_path)))
    expect_true(all(file.exists(as.character(run_index$result_path))))
    expect_true("any_low_replication_warning" %in% names(run_index))
    expect_true(all(c("analysis_method", "analysis_method_label", "method_output_dir", "selected_results_path", "selected_qc_path", "selected_bargraph_index_path") %in% names(method_index)))
    expect_equal(sort(method_index$analysis_method), c("normalized_t_test", "raw_log2_lm"))
    expect_true(all(file.exists(as.character(method_index$selected_results_path))))
    expect_true(all(file.exists(as.character(method_index$selected_qc_path))))
    expect_true(all(file.exists(as.character(method_index$selected_bargraph_index_path))))
    expect_equal(raw_selected_qc$analysis_method[[1]], "raw_log2_lm")
    expect_equal(normalized_selected_qc$analysis_method[[1]], "normalized_t_test")
    expect_true(all(raw_selected_qc$Name %in% c("Analyte A", "Analyte B")))
    expect_true(any(raw_selected_qc$low_signal_flag))
    expect_true(all(c("plot_status", "no_plot_reason") %in% names(normalized_selected_qc)))

    comparison_workbook_sheets <- openxlsx::getSheetNames(comparison_workbook_path)
    expect_equal(
        comparison_workbook_sheets[1:5],
        c("summary", "input_qc_summary", "method_summary", "significance_summary", "selected_analytes_summary")
    )
    selected_summary <- readxl::read_excel(comparison_workbook_path, sheet = "selected_analytes_summary")
    expect_equal(nrow(selected_summary), 8)
    expect_equal(sort(unique(selected_summary$comparison_slug)), c("female_control_vs_treated", "male_control_vs_treated"))
    expect_equal(sort(unique(selected_summary$analysis_method)), c("normalized_t_test", "raw_log2_lm"))
    expect_true(all(c("plot_status", "no_plot_reason", "bargraph_path", "selected_qc_path") %in% names(selected_summary)))
    expect_true(all(file.exists(as.character(stats::na.omit(selected_summary$selected_qc_path)))))
})

test_that("selected-analyte follow-up stops clearly when analytes are omitted", {
    tmp_dir <- tempfile("select-no-analytes-")
    dir.create(tmp_dir)

    fixture_paths <- write_env_fixture(
        path = file.path(tmp_dir, ".env"),
        runtime_root = file.path(tmp_dir, "runtime"),
        analysis_name = "replicate_smoke",
        include_selected_analytes = FALSE
    )

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        env_path = fixture_paths$env_path
    )

    expect_false(identical(shortlist_result$status, 0))
    expect_true(any(grepl("optional for setup validation and main analysis", shortlist_result$output)))
})
