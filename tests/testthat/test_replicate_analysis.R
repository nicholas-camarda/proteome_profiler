test_that("sample manifest validation catches missing required columns", {
    tmp_dir <- tempfile("manifest-bad-")
    dir.create(tmp_dir)
    manifest_path <- file.path(tmp_dir, "bad_manifest.csv")
    readr::write_csv(tibble(sample_id = "S1", workbook_path = "x.xlsx"), manifest_path)

    expect_error(
        read_sample_manifest(manifest_path, subgroup_var = "sex", treatment_var = "treatment"),
        "missing required columns"
    )
})

test_that("sample manifest validation catches blank required values", {
    tmp_dir <- tempfile("manifest-blank-")
    dir.create(tmp_dir)
    manifest_path <- file.path(tmp_dir, "blank_manifest.csv")
    readr::write_csv(
        tibble(
            sample_id = c("S1", ""),
            workbook_path = c("one.xlsx", "two.xlsx"),
            treatment = c("control", "treated"),
            sex = c("male", "female")
        ),
        manifest_path
    )

    expect_error(
        read_sample_manifest(manifest_path, subgroup_var = "sex", treatment_var = "treatment"),
        "blank or NA values"
    )
})

test_that("multi-sheet workbooks require a manifest sheet_name", {
    tmp_dir <- tempfile("manifest-multisheet-missing-")
    dir.create(tmp_dir)
    workbook_path <- file.path(tmp_dir, "collaborator.xlsx")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

    create_multisheet_replicate_workbook(workbook_path)
    readr::write_csv(
        tibble(
            sample_id = c("M01", "M02"),
            workbook_path = rep(workbook_path, 2),
            treatment = c("control", "treated"),
            sex = c("male", "male")
        ),
        manifest_path
    )

    manifest_tbl <- read_sample_manifest(manifest_path, subgroup_var = "sex", treatment_var = "treatment")

    expect_error(
        resolve_manifest_workbook_paths(manifest_tbl),
        "contains multiple sheets"
    )
})

test_that("replicate-aware dataset can read one workbook with many sample sheets", {
    tmp_dir <- tempfile("replicate-multisheet-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    workbook_path <- file.path(tmp_dir, "collaborator.xlsx")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

    write_protocol_fixture(protocol_path)
    create_multisheet_replicate_workbook(workbook_path)
    create_multisheet_replicate_manifest(manifest_path, workbook_path)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    manifest_tbl <- read_sample_manifest(manifest_path, subgroup_var = "sex", treatment_var = "treatment") %>%
        resolve_manifest_workbook_paths()

    sample_df <- build_sample_level_dataset(
        manifest = manifest_tbl,
        analyte_info = analyte_info,
        treatment_var = "treatment",
        subgroup_var = "sex"
    )

    expect_equal(length(unique(sample_df$sample_id)), 4)
    expect_equal(nrow(sample_df %>% filter(Name == "Analyte A")), 4)
    expect_true(all(sample_df$sheet_name %in% c("M01_sheet", "M02_sheet", "M03_sheet", "M04_sheet")))
    expect_true(all(sample_df$workbook_path == normalizePath(workbook_path, winslash = "/", mustWork = TRUE)))
    expect_true("normalized_signal" %in% names(sample_df))
    expect_true(any(is.finite(sample_df$normalized_signal)))
})

test_that("replicate-aware dataset preserves one row per analyte by sample", {
    tmp_dir <- tempfile("replicate-valid-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

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

    expect_equal(length(unique(sample_df$sample_id)), 8)
    expect_equal(nrow(sample_df %>% filter(Name == "Analyte A")), 8)
})

test_that("within-stratum inferential analysis stays within subgroup and adjusts p-values", {
    tmp_dir <- tempfile("replicate-inferential-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

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

    expect_equal(nrow(results$run_index), 2)
    expect_true(all(c("male_control_vs_treated", "female_control_vs_treated") %in% results$run_index$comparison_slug))

    male_results <- results$results[["male_control_vs_treated"]]
    female_results <- results$results[["female_control_vs_treated"]]

    expect_true(all(male_results$subgroup == "male"))
    expect_true(all(female_results$subgroup == "female"))
    expect_true("adjusted_p_value" %in% names(male_results))
    expect_true("effect_se_log2" %in% names(male_results))
    expect_true(all(is.finite(male_results$effect_se_log2[male_results$test_status == "tested"])))
    expect_true("low_signal_flag" %in% names(male_results))
    expect_true(male_results$low_replication_warning[1])
    expect_equal(
        male_results$adjusted_p_value[male_results$test_status == "tested"],
        p.adjust(male_results$raw_p_value[male_results$test_status == "tested"], method = "BH")
    )
})

test_that("inferential analysis reports underpowered analytes as not testable", {
    tmp_dir <- tempfile("replicate-underpowered-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

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
        min_reps = 3,
        low_signal_threshold = 70
    )

    male_results <- results$results[["male_control_vs_treated"]]
    expect_true(all(male_results$test_status == "not_testable"))
    expect_true(all(is.na(male_results$raw_p_value)))
    expect_true(all(str_detect(male_results$test_reason, "requires >= 3 finite positive signals per arm")))
})

test_that("configured comparison labels must exist within each subgroup", {
    tmp_dir <- tempfile("replicate-bad-contrast-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

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

    expect_error(
        run_within_stratum_differential_analysis(
            sample_data = sample_df,
            comparisons = list("control" = c("treated_typo")),
            subgroup_var = "sex",
            p_adjust_method = "BH",
            alpha = 0.05,
            min_reps = 2,
            low_signal_threshold = 70
        ),
        "Configured treatment label"
    )
})

test_that("inferential gating counts only finite positive signals toward minimum replicates", {
    tmp_dir <- tempfile("replicate-nonpositive-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    manifest_path <- file.path(tmp_dir, "manifest.csv")

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
    ) %>%
        mutate(signal = if_else(sample_id == "M01" & Name == "Analyte A", 0, signal))

    results <- run_within_stratum_differential_analysis(
        sample_data = sample_df,
        comparisons = list("control" = c("treated")),
        subgroup_var = "sex",
        p_adjust_method = "BH",
        alpha = 0.05,
        min_reps = 2,
        low_signal_threshold = 70
    )

    analyte_a_row <- results$results[["male_control_vs_treated"]] %>%
        filter(Name == "Analyte A")

    expect_equal(analyte_a_row$control_rows_total[[1]], 2)
    expect_equal(analyte_a_row$control_n[[1]], 1)
    expect_equal(analyte_a_row$test_status[[1]], "not_testable")
    expect_true(str_detect(analyte_a_row$test_reason[[1]], "finite positive signals per arm"))
})

test_that("strong replicate-aware fixture exercises BH, uneven subgroups, and partial missingness", {
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

    expect_equal(
        sort(results$run_index$comparison_slug),
        sort(c(
            "male_control_vs_drug_a",
            "male_control_vs_drug_b",
            "female_control_vs_drug_a",
            "female_control_vs_drug_b"
        ))
    )

    male_drug_a <- results$results[["male_control_vs_drug_a"]]
    female_drug_b <- results$results[["female_control_vs_drug_b"]]

    expect_equal(nrow(male_drug_a), 12)

    moderate_a <- male_drug_a %>% filter(Name == "Analyte Moderate A")
    expect_lt(moderate_a$raw_p_value[[1]], 0.05)
    expect_gt(moderate_a$adjusted_p_value[[1]], 0.05)
    expect_true(moderate_a$raw_p_lt_alpha[[1]])
    expect_false(is.na(moderate_a$fdr_lt_0_20[[1]]))

    strong_a <- male_drug_a %>% filter(Name == "Analyte Strong A")
    expect_true(strong_a$fdr_lt_0_20[[1]])
    expect_true(strong_a$fdr_lt_0_25[[1]])

    female_missing <- female_drug_b %>% filter(Name == "Analyte Missing")
    expect_equal(female_missing$test_status[[1]], "not_testable")
    expect_equal(female_missing$treatment_rows_total[[1]], 2)
    expect_equal(female_missing$treatment_n[[1]], 0)

    low_signal_row <- female_drug_b %>% filter(Name == "Analyte Low")
    expect_true(low_signal_row$low_signal_flag[[1]])

    female_run_index <- results$run_index %>%
        filter(subgroup == "female")
    expect_true(all(female_run_index$any_low_replication_warning))
    expect_true(all(female_run_index$min_control_n == 2))
})

test_that("normalized t-test uses normalized sheet values", {
    sample_df <- create_complex_replicate_sample_data()

    results <- run_within_stratum_differential_analysis(
        sample_data = sample_df,
        comparisons = list("control" = c("drug_a")),
        subgroup_var = "sex",
        p_adjust_method = "BH",
        alpha = 0.05,
        min_reps = 2,
        low_signal_threshold = 70,
        analysis_method = "normalized_t_test"
    )

    male_drug_a <- results$results[["male_control_vs_drug_a"]]
    strong_a <- male_drug_a %>% filter(Name == "Analyte Strong A")
    expected_p <- with(
        sample_df %>% filter(sex == "male", treatment %in% c("control", "drug_a"), Name == "Analyte Strong A"),
        stats::t.test(normalized_signal ~ treatment, var.equal = TRUE)$p.value
    )

    expect_equal(strong_a$analysis_method[[1]], "normalized_t_test")
    expect_equal(strong_a$analysis_method_label[[1]], get_inferential_method_spec("normalized_t_test")$label)
    expect_equal(strong_a$raw_p_value[[1]], expected_p)
    expect_true(is.finite(strong_a$fold_change_ratio[[1]]))
    expected_se <- with(
        sample_df %>% filter(sex == "male", treatment %in% c("control", "drug_a"), Name == "Analyte Strong A"),
        pooled_log2_ratio_se(
            control_values = normalized_signal[treatment == "control"],
            treatment_values = normalized_signal[treatment == "drug_a"],
            control_mean = mean(normalized_signal[treatment == "control"]),
            treatment_mean = mean(normalized_signal[treatment == "drug_a"])
        )
    )
    expect_equal(strong_a$effect_se_log2[[1]], expected_se)
    expect_true(is.finite(strong_a$control_mean_normalized[[1]]))
    expect_true(is.finite(strong_a$treatment_mean_normalized[[1]]))
})
