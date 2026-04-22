test_that("config validation distinguishes legacy and replicate-aware modes", {
    legacy_config <- list(
        user = "nick",
        analysis_slug = "legacy_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        data_dir = "data/legacy",
        comparisons = list("control" = c("treated")),
        group_levels = c("control", "treated")
    )

    inferential_config <- list(
        mode = "replicate",
        user = "lauren",
        analysis_slug = "replicate_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        sample_manifest = "manifests/example_samples.csv",
        subgroup_var = "sex",
        treatment_var = "treatment",
        comparisons = list("control" = c("treated")),
        p_adjust_method = "BH",
        alpha = 0.05
    )

    expect_false(is_replicate_aware_config(legacy_config))
    expect_true(is_replicate_aware_config(inferential_config))
    expect_no_error(validate_analysis_config(legacy_config))
    expect_no_error(validate_analysis_config(inferential_config))
})

test_that("manifest-driven exploratory configs validate without inferential fields", {
    exploratory_manifest_config <- list(
        mode = "legacy",
        user = "nick",
        analysis_slug = "legacy_manifest_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        sample_manifest = "manifests/example_samples.csv",
        treatment_var = "treatment",
        comparisons = list("control" = c("treated"))
    )

    expect_false(is_replicate_aware_config(exploratory_manifest_config))
    expect_no_error(validate_analysis_config(exploratory_manifest_config))
})

test_that("replicate-aware configs default to BH when p_adjust_method is omitted", {
    config_env <- environment(get_analysis_config)
    config_env$proteome_profiler_config$analyses$replicate_missing_method <- list(
        user = "lauren",
        analysis_slug = "replicate_default_bh",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        sample_manifest = "manifests/example_samples.csv",
        subgroup_var = "sex",
        treatment_var = "treatment",
        comparisons = list("control" = c("treated")),
        alpha = 0.05
    )

    inferred_config <- get_analysis_config("replicate_missing_method")
    expect_equal(inferred_config$p_adjust_method, "BH")
    expect_no_error(validate_analysis_config(inferred_config))
})

test_that("grouped user-facing config shape is normalized into internal fields", {
    config_env <- environment(get_analysis_config)
    config_env$proteome_profiler_config$analyses$grouped_example <- list(
        user = "lauren",
        slug = "grouped_slug",
        protocol = list(
            preset = "cytokine_xl",
            workbook = "output/protocol.xlsx"
        ),
        input = list(
            manifest = "manifests/example_samples.csv",
            subgroup = "sex",
            treatment = "treatment"
        ),
        comparisons = list("control" = c("treated")),
        thresholds = list(
            ref_coords = c("A3,4"),
            ref_signal = c(150)
        ),
        stats = list(
            alpha = 0.05
        ),
        shortlist = list(
            comparison = "male_control_vs_treated",
            analytes = c("Analyte A", "Analyte B")
        )
    )

    normalized_config <- get_analysis_config("grouped_example")

    expect_equal(normalized_config$analysis_slug, "grouped_slug")
    expect_equal(normalized_config$protocol_preset, "cytokine_xl")
    expect_equal(normalized_config$info_fn, "output/protocol.xlsx")
    expect_equal(normalized_config$sample_manifest, "manifests/example_samples.csv")
    expect_equal(normalized_config$subgroup_var, "sex")
    expect_equal(normalized_config$treatment_var, "treatment")
    expect_equal(normalized_config$ref_coords_to_make_filter, c("A3,4"))
    expect_equal(normalized_config$ref_thresh_to_filter, c(150))
    expect_equal(
        normalized_config$selection_comparison_slugs,
        "male_control_vs_treated"
    )
    expect_equal(normalized_config$selection_comparison_slug, "male_control_vs_treated")
    expect_equal(normalized_config$selection_analytes, c("Analyte A", "Analyte B"))
    expect_equal(normalized_config$p_adjust_method, "BH")
})

test_that("grouped shortlist comparisons accept multiple comparison slugs", {
    config_env <- environment(get_analysis_config)
    config_env$proteome_profiler_config$analyses$multi_shortlist_example <- list(
        user = "lauren",
        slug = "multi_shortlist_slug",
        protocol = list(
            preset = "cytokine_xl",
            workbook = "output/protocol.xlsx"
        ),
        input = list(
            manifest = "manifests/example_samples.csv",
            subgroup = "sex",
            treatment = "treatment"
        ),
        comparisons = list("control" = c("treated")),
        thresholds = list(
            ref_coords = c("A3,4"),
            ref_signal = c(150)
        ),
        stats = list(
            alpha = 0.05
        ),
        shortlist = list(
            comparisons = c("male_control_vs_treated", "female_control_vs_treated")
        )
    )

    normalized_config <- get_analysis_config("multi_shortlist_example")

    expect_equal(
        normalized_config$selection_comparison_slugs,
        c("male_control_vs_treated", "female_control_vs_treated")
    )
})

test_that("select-analytes requires explicit analyte names", {
    expect_error(
        get_selected_analyte_names(list()),
        "shortlist\\$analytes"
    )
})

test_that("legacy select comparisons resolve to comparison-scoped slugs", {
    config <- list(
        comparisons = list(vehicle = c("sorafenib", "sor + dox")),
        selection_comparison_slugs = c("vehicle_vs_sorafenib", "vehicle_vs_sor_dox")
    )

    selected_pairs <- resolve_legacy_select_comparisons(config)

    expect_equal(
        selected_pairs$comparison_slug,
        c("vehicle_vs_sorafenib", "vehicle_vs_sor_dox")
    )
    expect_equal(selected_pairs$control, c("vehicle", "vehicle"))
    expect_equal(selected_pairs$treatment, c("sorafenib", "sor + dox"))
})

test_that("replicate-aware grouped configs reject unsupported shortlist knobs", {
    config_env <- environment(get_analysis_config)
    config_env$proteome_profiler_config$analyses$replicate_shortlist_basis <- list(
        user = "lauren",
        slug = "replicate_shortlist_basis",
        protocol = list(
            preset = "cytokine_xl",
            workbook = "output/protocol.xlsx"
        ),
        input = list(
            manifest = "manifests/example_samples.csv",
            subgroup = "sex",
            treatment = "treatment"
        ),
        comparisons = list("control" = c("treated")),
        stats = list(
            alpha = 0.05
        ),
        shortlist = list(
            comparisons = c("male_control_vs_treated"),
            basis = "fdr_lt_0_20",
            top_n = 10,
            write_bargraphs = TRUE
        )
    )

    expect_error(
        get_analysis_config("replicate_shortlist_basis"),
        "uses explicit `shortlist\\$analytes`"
    )
})

test_that("analysis output roots are grouped by user and analysis slug", {
    config <- list(
        user = "lauren",
        analysis_slug = "replicate_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        comparisons = list("control" = c("treated")),
        data_dir = "data/legacy",
        group_levels = c("control", "treated")
    )

    output_root <- get_analysis_output_root(config)
    expect_true(grepl("output/lauren/replicate_example$", output_root))
})

test_that("selected analysis name comes from env input or optional config default", {
    config_env <- environment(get_selected_analysis_name)
    config_env$proteome_profiler_config$default_analysis <- "legacy_example"

    expect_equal(get_selected_analysis_name("replicate_example"), "replicate_example")
    expect_equal(get_selected_analysis_name(""), "legacy_example")
})

test_that("selected analysis name errors when no explicit or configured default exists", {
    config_env <- environment(get_selected_analysis_name)
    old_default <- config_env$proteome_profiler_config$default_analysis
    old_analyses <- config_env$proteome_profiler_config$analyses
    config_env$proteome_profiler_config$default_analysis <- NULL
    config_env$proteome_profiler_config$analyses <- list(
        legacy_example = list(),
        replicate_example = list()
    )
    on.exit({
        config_env$proteome_profiler_config$default_analysis <- old_default
        config_env$proteome_profiler_config$analyses <- old_analyses
    }, add = TRUE)

    expect_error(
        get_selected_analysis_name(""),
        "No analysis selected"
    )
})
