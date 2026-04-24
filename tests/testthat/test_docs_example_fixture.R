test_that("synthetic docs example fixture regenerates public-safe docs outputs", {
    old_wd <- setwd(repo_root)
    on.exit(setwd(old_wd), add = TRUE)

    script_output <- system2(
        file.path(R.home("bin"), "Rscript"),
        "tests/generate_docs_example_fixture.R",
        stdout = TRUE,
        stderr = TRUE
    )
    script_status <- attr(script_output, "status")
    expect_true(
        is.null(script_status) || identical(script_status, 0L),
        info = paste(script_output, collapse = "\n")
    )

    docs_output_dir <- file.path(repo_root, "docs", "output-examples")
    expected_docs_files <- file.path(docs_output_dir, c(
        "replicate-aware-inferential-waterfall-example.png",
        "replicate-aware-inferential-barplot-example-page-1.png",
        "threshold-diagnostics-example.png",
        "main-waterfall-example.png",
        "main-barplot-example-page-1.png",
        "selected-analytes-waterfall-example.png"
    ))
    expect_true(all(file.exists(expected_docs_files)))

    fixture_root <- file.path(repo_root, "tests", "test_output", "docs_example_fixture")
    exploratory_root <- file.path(fixture_root, "output", "docs", "docs_exploratory_examples")
    replicate_root <- file.path(fixture_root, "output", "docs", "docs_replicate_examples")
    expect_true(dir.exists(exploratory_root))
    expect_true(dir.exists(replicate_root))

    run_index_path <- file.path(replicate_root, "inferential_results", "run_index.tsv")
    expect_true(file.exists(run_index_path))

    run_index <- readr::read_tsv(run_index_path, show_col_types = FALSE)
    expect_equal(
        sort(unique(run_index$comparison_slug)),
        c("stratum_1_group_a_vs_group_b", "stratum_2_group_a_vs_group_b")
    )
    expect_false(any(grepl("vehicle|sorafenib|aldosterone|nicole|nick", run_index$comparison_slug, ignore.case = TRUE)))

    expect_true(file.exists(file.path(
        replicate_root,
        "select_analytes",
        "stratum_1_group_a_vs_group_b",
        "raw_log2_lm",
        "selected_waterfall.png"
    )))
    expect_true(file.exists(file.path(
        replicate_root,
        "select_analytes",
        "stratum_2_group_a_vs_group_b",
        "normalized_t_test",
        "selected_bargraph_index.tsv"
    )))
    expect_true(file.exists(file.path(
        exploratory_root,
        "select_analytes",
        "group_a_vs_group_b",
        "selected_waterfall.png"
    )))

    docs_readme <- readLines(file.path(repo_root, "docs", "README.md"), warn = FALSE)
    docs_readme <- paste(docs_readme, collapse = "\n")
    expect_match(docs_readme, "replicate-aware-inferential-waterfall-example\\.png")
    expect_match(docs_readme, "main-waterfall-example\\.png")
    expect_match(docs_readme, "synthetic public-safe fixture")
    expect_false(grepl("main-waterfall-sorafenib\\.png", docs_readme))
    expect_false(grepl("replicate-aware-inferential-waterfall-se\\.png", docs_readme))
})
