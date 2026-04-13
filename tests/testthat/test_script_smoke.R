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
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "waterfall_plot_REF-80_SOR-1_5.png")))
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
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "waterfall_plot_REF-70_SOR-1_5.png")))
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
    shortlist_summary <- readr::read_tsv(
        file.path(analysis_root, "select_analytes", "comparisons", "male_control_vs_treated", "shortlist_summary.tsv"),
        show_col_types = FALSE
    )
    female_shortlist_summary <- readr::read_tsv(
        file.path(analysis_root, "select_analytes", "comparisons", "female_control_vs_treated", "shortlist_summary.tsv"),
        show_col_types = FALSE
    )

    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "candidate_low_signal_analytes.tsv")))
    expect_true(file.exists(file.path(analysis_root, "threshold_diagnostics", "region_stats.png")))
    expect_true(nrow(threshold_candidates) > 0)
    expect_true(all(c("Name", "Coordinate", "suggestion_score") %in% names(threshold_candidates)))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "run_index.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "male_control_vs_treated", "results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "inferential_results", "comparisons", "female_control_vs_treated", "results.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "comparisons", "male_control_vs_treated", "shortlist.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "comparisons", "male_control_vs_treated", "shortlist_waterfall.png")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "comparisons", "female_control_vs_treated", "shortlist.tsv")))
    expect_true(file.exists(file.path(analysis_root, "select_analytes", "comparisons", "female_control_vs_treated", "shortlist_waterfall.png")))
    expect_true("result_path" %in% names(run_index))
    expect_true(all(!is.na(run_index$result_path)))
    expect_true(all(file.exists(as.character(run_index$result_path))))
    expect_true("any_low_replication_warning" %in% names(run_index))
    expect_equal(shortlist_summary$shortlist_basis[[1]], "fdr_lt_0_20")
    expect_equal(female_shortlist_summary$shortlist_basis[[1]], "top_n_fallback")
})
