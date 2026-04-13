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
    expect_true(file.exists(file.path(output_dir, "comparisons", "male_control_vs_treated", "waterfall.png")))
    expect_true(file.exists(file.path(output_dir, "comparisons", "female_control_vs_treated", "waterfall.png")))
})
