expect_script_success <- function(script_result) {
    expect_equal(
        script_result$status,
        0,
        info = paste(script_result$output, collapse = "\n")
    )
}

test_that("legacy entry scripts run end-to-end on fixture data", {
    tmp_dir <- tempfile("legacy-smoke-")
    dir.create(tmp_dir)

    fixture_paths <- write_analysis_config_fixture(
        path = file.path(tmp_dir, "analysis_config.R"),
        runtime_root = file.path(tmp_dir, "runtime")
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        analysis_name = "legacy_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        analysis_name = "legacy_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        analysis_name = "legacy_smoke",
        config_path = fixture_paths$config_path
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
    tmp_dir <- tempfile("legacy-manifest-smoke-")
    dir.create(tmp_dir)

    fixture_paths <- write_analysis_config_fixture(
        path = file.path(tmp_dir, "analysis_config.R"),
        runtime_root = file.path(tmp_dir, "runtime")
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        analysis_name = "legacy_manifest_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        analysis_name = "legacy_manifest_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        analysis_name = "legacy_manifest_smoke",
        config_path = fixture_paths$config_path
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

    fixture_paths <- write_analysis_config_fixture(
        path = file.path(tmp_dir, "analysis_config.R"),
        runtime_root = file.path(tmp_dir, "runtime")
    )

    find_thresh_result <- run_repo_script(
        script_path = "scripts/find_ref_thresh.R",
        analysis_name = "replicate_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(find_thresh_result)

    main_result <- run_repo_script(
        script_path = "scripts/main.R",
        analysis_name = "replicate_smoke",
        config_path = fixture_paths$config_path
    )
    expect_script_success(main_result)

    shortlist_result <- run_repo_script(
        script_path = "scripts/select-analytes-analysis.R",
        analysis_name = "replicate_smoke",
        config_path = fixture_paths$config_path
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
})
