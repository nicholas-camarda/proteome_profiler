proteome_env_names <- c(
    "PROTEOME_PROFILER_ENV_FILE",
    "PROTEOME_PROFILER_ANALYSIS",
    "PROTEOME_PROFILER_MODE",
    "PROTEOME_PROFILER_USER",
    "PROTEOME_PROFILER_SLUG",
    "PROTEOME_PROFILER_RUNTIME_ROOT",
    "PROTEOME_PROFILER_CLOUD_PARENT",
    "PROTEOME_PROFILER_PROTOCOL_PRESET",
    "PROTEOME_PROFILER_PROTOCOL_WORKBOOK",
    "PROTEOME_PROFILER_PROTOCOL_PDF",
    "PROTEOME_PROFILER_INPUT_MANIFEST",
    "PROTEOME_PROFILER_TREATMENT_COLUMN",
    "PROTEOME_PROFILER_SUBGROUP_COLUMN",
    "PROTEOME_PROFILER_COMPARISONS",
    "PROTEOME_PROFILER_REF_COORDS",
    "PROTEOME_PROFILER_REF_SIGNAL",
    "PROTEOME_PROFILER_ANALYSIS_METHODS",
    "PROTEOME_PROFILER_SHORTLIST_COORDS",
    "PROTEOME_PROFILER_SHORTLIST_ANALYTES"
)

restore_proteome_env_later <- function() {
    old_values <- Sys.getenv(proteome_env_names, unset = NA_character_)
    Sys.unsetenv(proteome_env_names)
    function() {
        Sys.unsetenv(proteome_env_names)
        for (idx in seq_along(proteome_env_names)) {
            if (!is.na(old_values[[idx]])) {
                do.call(Sys.setenv, as.list(stats::setNames(old_values[[idx]], proteome_env_names[[idx]])))
            }
        }
    }
}

test_that("config validation distinguishes exploratory and replicate-aware modes", {
    exploratory_config <- list(
        user = "nick",
        analysis_slug = "exploratory_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        data_dir = "data/exploratory",
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

    expect_false(is_replicate_aware_config(exploratory_config))
    expect_true(is_replicate_aware_config(inferential_config))
    expect_no_error(validate_analysis_config(exploratory_config))
    expect_no_error(validate_analysis_config(inferential_config))
})

test_that("manifest-driven exploratory configs validate without inferential fields", {
    exploratory_manifest_config <- list(
        mode = "exploratory",
        user = "nick",
        analysis_slug = "exploratory_manifest_example",
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
            coords = c("A1,2", "A3,4")
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
    expect_equal(normalized_config$selection_coords, c("A1,2", "A3,4"))
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

test_that(".env run sheet parses into the internal analysis config", {
    restore_env <- restore_proteome_env_later()
    on.exit(restore_env(), add = TRUE)

    env_path <- tempfile(".env-")
    writeLines(c(
        "PROTEOME_PROFILER_ANALYSIS=env_example",
        "PROTEOME_PROFILER_MODE=replicate",
        "PROTEOME_PROFILER_USER=lauren",
        "PROTEOME_PROFILER_SLUG=env_slug",
        "PROTEOME_PROFILER_RUNTIME_ROOT=/tmp/proteome-runtime",
        "PROTEOME_PROFILER_PROTOCOL_PRESET=cytokine_xl",
        "PROTEOME_PROFILER_PROTOCOL_WORKBOOK=output/protocol.xlsx",
        "PROTEOME_PROFILER_INPUT_MANIFEST=manifests/example_samples.csv",
        "PROTEOME_PROFILER_TREATMENT_COLUMN=treatment",
        "PROTEOME_PROFILER_SUBGROUP_COLUMN=sex",
        "PROTEOME_PROFILER_COMPARISONS=control=treated",
        "PROTEOME_PROFILER_REF_COORDS=A3,4",
        "PROTEOME_PROFILER_REF_SIGNAL=150",
        "PROTEOME_PROFILER_ANALYSIS_METHODS=raw_log2_lm|normalized_t_test",
        "PROTEOME_PROFILER_SHORTLIST_COORDS=A1,2|A3,4"
    ), env_path)
    Sys.setenv(PROTEOME_PROFILER_ENV_FILE = env_path)

    initialize_runtime_config_from_env(required_env_file = TRUE)
    normalized_config <- get_analysis_config("env_example")

    expect_equal(get_selected_analysis_name(), "env_example")
    expect_equal(normalized_config$analysis_slug, "env_slug")
    expect_equal(normalized_config$sample_manifest, "manifests/example_samples.csv")
    expect_equal(normalized_config$subgroup_var, "sex")
    expect_equal(normalized_config$comparisons, list(control = "treated"))
    expect_equal(normalized_config$ref_thresh_to_filter, 150)
    expect_equal(normalized_config$analysis_methods, c("raw_log2_lm", "normalized_t_test"))
    expect_equal(normalized_config$selection_coords, c("A1,2", "A3,4"))
})

test_that("explicit process environment values override .env file values", {
    restore_env <- restore_proteome_env_later()
    on.exit(restore_env(), add = TRUE)

    env_path <- tempfile(".env-")
    writeLines(c(
        "PROTEOME_PROFILER_ANALYSIS=env_override",
        "PROTEOME_PROFILER_MODE=exploratory",
        "PROTEOME_PROFILER_USER=tester",
        "PROTEOME_PROFILER_SLUG=file_slug",
        "PROTEOME_PROFILER_RUNTIME_ROOT=/tmp/proteome-runtime",
        "PROTEOME_PROFILER_PROTOCOL_WORKBOOK=output/protocol.xlsx",
        "PROTEOME_PROFILER_INPUT_DATA_DIR=data/exploratory",
        "PROTEOME_PROFILER_GROUP_LEVELS=vehicle|treated",
        "PROTEOME_PROFILER_COMPARISONS=vehicle=treated"
    ), env_path)
    Sys.setenv(
        PROTEOME_PROFILER_ENV_FILE = env_path,
        PROTEOME_PROFILER_SLUG = "process_slug"
    )

    initialize_runtime_config_from_env(required_env_file = TRUE)
    normalized_config <- get_analysis_config("env_override")

    expect_equal(normalized_config$analysis_slug, "process_slug")
})

test_that("selected analytes can be omitted from .env for main analysis config", {
    restore_env <- restore_proteome_env_later()
    on.exit(restore_env(), add = TRUE)

    env_path <- tempfile(".env-")
    writeLines(c(
        "PROTEOME_PROFILER_ANALYSIS=no_selected",
        "PROTEOME_PROFILER_MODE=replicate",
        "PROTEOME_PROFILER_USER=tester",
        "PROTEOME_PROFILER_SLUG=no_selected",
        "PROTEOME_PROFILER_RUNTIME_ROOT=/tmp/proteome-runtime",
        "PROTEOME_PROFILER_PROTOCOL_WORKBOOK=output/protocol.xlsx",
        "PROTEOME_PROFILER_INPUT_MANIFEST=manifests/example_samples.csv",
        "PROTEOME_PROFILER_TREATMENT_COLUMN=treatment",
        "PROTEOME_PROFILER_SUBGROUP_COLUMN=sex",
        "PROTEOME_PROFILER_COMPARISONS=control=treated",
        "PROTEOME_PROFILER_ALPHA=0.05"
    ), env_path)
    Sys.setenv(PROTEOME_PROFILER_ENV_FILE = env_path)

    initialize_runtime_config_from_env(required_env_file = TRUE)
    normalized_config <- get_analysis_config("no_selected")

    expect_no_error(validate_analysis_config(normalized_config))
    expect_null(normalized_config$selection_coords)
})

test_that("selected analyte coordinates are optional for main analysis but required for select-analytes", {
    expect_error(
        get_selected_analyte_coordinates(list()),
        "optional for setup validation and main analysis"
    )
})

test_that("selected analyte coordinates normalize common user-entered forms", {
    available_tbl <- tibble::tibble(
        Name = c("Analyte A", "Analyte B"),
        Coordinate = c("A1,2", "A3,4")
    )

    expect_equal(
        compact_coordinate_text(c("A1, A2", "A1,2", "A1,A2", "a1, 2")),
        rep("A1,2", 4)
    )

    selected_tbl <- filter_selected_analyte_coordinates(
        available_tbl,
        selected_coords = c("A3,A4", "A1, 2")
    )

    expect_equal(selected_tbl$Name, c("Analyte B", "Analyte A"))
    expect_equal(selected_tbl$Coordinate, c("A3,4", "A1,2"))
})

test_that("missing selected analyte coordinates fail with coordinate suggestions", {
    available_tbl <- tibble::tibble(
        Name = c("Analyte A", "Analyte B"),
        Coordinate = c("A1,2", "A3,4")
    )

    expect_error(
        validate_selected_analyte_coordinates(c("A9,10"), available_tbl),
        "A9,10 \\(closest available: A1,2, A3,4\\)"
    )
})

test_that("deprecated selected-analyte name config fails clearly", {
    restore_env <- restore_proteome_env_later()
    on.exit(restore_env(), add = TRUE)

    env_path <- tempfile(".env-")
    writeLines(c(
        "PROTEOME_PROFILER_ANALYSIS=deprecated_selected_names",
        "PROTEOME_PROFILER_MODE=replicate",
        "PROTEOME_PROFILER_USER=tester",
        "PROTEOME_PROFILER_SLUG=deprecated_selected_names",
        "PROTEOME_PROFILER_RUNTIME_ROOT=/tmp/proteome-runtime",
        "PROTEOME_PROFILER_PROTOCOL_WORKBOOK=output/protocol.xlsx",
        "PROTEOME_PROFILER_INPUT_MANIFEST=manifests/example_samples.csv",
        "PROTEOME_PROFILER_TREATMENT_COLUMN=treatment",
        "PROTEOME_PROFILER_SUBGROUP_COLUMN=sex",
        "PROTEOME_PROFILER_COMPARISONS=control=treated",
        "PROTEOME_PROFILER_ALPHA=0.05",
        "PROTEOME_PROFILER_SHORTLIST_ANALYTES=Analyte A|Analyte B"
    ), env_path)
    Sys.setenv(PROTEOME_PROFILER_ENV_FILE = env_path)

    expect_error(
        initialize_runtime_config_from_env(required_env_file = TRUE),
        "PROTEOME_PROFILER_SHORTLIST_ANALYTES is not supported"
    )
})

test_that("exploratory select comparisons resolve to comparison-scoped slugs", {
    config <- list(
        comparisons = list(vehicle = c("sorafenib", "sor + dox")),
        selection_comparison_slugs = c("vehicle_vs_sorafenib", "vehicle_vs_sor_dox")
    )

    selected_pairs <- resolve_exploratory_select_comparisons(config)

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
        "uses explicit `shortlist\\$coords`"
    )
})

test_that("analysis output roots are grouped by user and analysis slug", {
    config <- list(
        user = "lauren",
        analysis_slug = "replicate_example",
        protocol_preset = "cytokine_xl",
        info_fn = "output/protocol.xlsx",
        comparisons = list("control" = c("treated")),
        data_dir = "data/exploratory",
        group_levels = c("control", "treated")
    )

    output_root <- get_analysis_output_root(config)
    expect_true(grepl("output/lauren/replicate_example$", output_root))
})

test_that("selected analysis name comes from active .env-derived config", {
    config_env <- environment(get_selected_analysis_name)
    config_env$proteome_profiler_config$default_analysis <- "exploratory_example"

    expect_equal(get_selected_analysis_name("replicate_example"), "replicate_example")
    expect_equal(get_selected_analysis_name(""), "exploratory_example")
})

test_that("selected analysis name errors when active config has no analysis name", {
    config_env <- environment(get_selected_analysis_name)
    old_default <- config_env$proteome_profiler_config$default_analysis
    old_analyses <- config_env$proteome_profiler_config$analyses
    config_env$proteome_profiler_config$default_analysis <- NULL
    config_env$proteome_profiler_config$analyses <- list(
        exploratory_example = list(),
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
