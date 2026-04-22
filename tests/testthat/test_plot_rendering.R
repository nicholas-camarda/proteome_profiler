test_that("waterfall plot caption does not print literal NA when title has two lines", {
    plot_data <- tibble(
        Name = c("A", "B"),
        relative_signal = c(2, 0.5),
        log_relative_signal = log2(c(2, 0.5)),
        main_group = "treated",
        control = "control"
    )

    plt <- plot_wf_graph(
        data = plot_data,
        TITLE = "Short title\nSubtitle only",
        main_cutoff = NULL
    )

    expect_false(grepl("NA", plt$labels$caption, fixed = TRUE))
})

test_that("inferential waterfall captions standard-error whiskers when available", {
    results_tbl <- tibble(
        Name = c("A", "B"),
        effect_estimate_log2 = c(1, -0.5),
        effect_se_log2 = c(0.2, 0.1),
        test_status = "tested",
        waterfall_subtitle = "Effect estimate is treatment minus control"
    )

    plt <- plot_inferential_waterfall(results_tbl, title = "control vs treated")

    expect_true(grepl("whiskers = +/- 1 SE", plt$labels$caption, fixed = TRUE))
})

test_that("inferential waterfall requires standard errors for plotted effects", {
    results_tbl <- tibble(
        Name = c("A", "B"),
        effect_estimate_log2 = c(1, -0.5),
        effect_se_log2 = c(0.2, NA_real_),
        test_status = "tested",
        waterfall_subtitle = "Effect estimate is treatment minus control"
    )

    expect_error(
        plot_inferential_waterfall(results_tbl, title = "control vs treated"),
        "missing finite standard errors"
    )
})

test_that("faceted inferential barplots define star threshold in page caption", {
    result_tbl <- tibble(
        Name = c("A", "B"),
        Coordinate = c("A1", "A2"),
        control = "control",
        treatment = "drug",
        fold_change_ratio = c(1.8, 1.3),
        effect_estimate_log2 = log2(c(1.8, 1.3)),
        effect_se_log2 = c(0.2, 0.1),
        control_fold_change_ymin = c(0.9, 0.95),
        control_fold_change_ymax = c(1.1, 1.05),
        treatment_fold_change_ymin = c(1.6, 1.2),
        treatment_fold_change_ymax = c(2.0, 1.4)
    )

    barplot_data <- build_inferential_fold_change_barplot_data(
        result_tbl,
        mark_treatment_significant = TRUE,
        significance_label = "*",
        significance_definition = "FDR < 0.20"
    )
    pages <- plot_inferential_fold_change_barplot_pages(
        barplot_data,
        title = "FDR < 0.20 Hits",
        groups_per_page = 25
    )

    expect_true(all(na.omit(barplot_data$significance_label) == "*"))
    expect_equal(pages[[1]]$patches$annotation$caption, "* = FDR < 0.20; Analytes shown = 2 / 2")
})

test_that("faceted inferential barplots render fixed-size panel title strips", {
    result_tbl <- tibble(
        Name = c("LDL R", "C-Reactive Protein"),
        Coordinate = c("G23,24", "C15,16"),
        control = "vehicle",
        treatment = "aldosterone",
        fold_change_ratio = c(1.4, 1.5),
        effect_estimate_log2 = log2(c(1.4, 1.5)),
        effect_se_log2 = c(0.1, 0.1),
        control_fold_change_ymin = c(0.9, 0.9),
        control_fold_change_ymax = c(1.1, 1.1),
        treatment_fold_change_ymin = c(1.3, 1.4),
        treatment_fold_change_ymax = c(1.5, 1.6)
    )

    barplot_data <- build_inferential_fold_change_barplot_data(
        result_tbl,
        mark_treatment_significant = TRUE,
        significance_label = "*",
        significance_definition = "FDR < 0.25"
    )
    pages <- plot_inferential_fold_change_barplot_pages(
        barplot_data,
        title = "FDR < 0.25 Hits",
        groups_per_page = 25
    )
    out_path <- tempfile(fileext = ".png")

    ggsave(out_path, plot = pages[[1]], width = 10, height = 10, dpi = 72)

    expect_true(file.exists(out_path))
    expect_gt(file.info(out_path)$size, 0)
})

test_that("replicate-aware barplot data includes control and treatment SE whiskers on fold-change scale", {
    result_tbl <- tibble(
        Name = "A",
        Coordinate = "A1",
        control = "control",
        treatment = "drug",
        fold_change_ratio = 2,
        effect_estimate_log2 = 1,
        effect_se_log2 = 0.25,
        control_fold_change_ymin = 0.8,
        control_fold_change_ymax = 1.2,
        treatment_fold_change_ymin = 1.7,
        treatment_fold_change_ymax = 2.3
    )

    barplot_data <- build_inferential_fold_change_barplot_data(result_tbl)
    treatment_row <- barplot_data %>% filter(as.character(group) == "drug")
    control_row <- barplot_data %>% filter(as.character(group) == "control")

    expect_equal(control_row$relative_signal_ymin[[1]], 0.8)
    expect_equal(control_row$relative_signal_ymax[[1]], 1.2)
    expect_equal(treatment_row$relative_signal_ymin[[1]], 1.7)
    expect_equal(treatment_row$relative_signal_ymax[[1]], 2.3)
})

test_that("replicate-aware barplot y-axis includes SE whisker headroom", {
    result_tbl <- tibble(
        Name = "A",
        Coordinate = "A1",
        control = "control",
        treatment = "drug",
        fold_change_ratio = 2,
        effect_estimate_log2 = 1,
        effect_se_log2 = 1,
        control_fold_change_ymin = 0.5,
        control_fold_change_ymax = 1.5,
        treatment_fold_change_ymin = 1.1,
        treatment_fold_change_ymax = 4
    )

    barplot_data <- build_inferential_fold_change_barplot_data(result_tbl)
    pages <- plot_inferential_fold_change_barplot_pages(
        barplot_data,
        title = "All Tested",
        groups_per_page = 1
    )

    expect_true(length(pages) == 1)
    expect_gt(max(barplot_data$relative_signal_ymax, na.rm = TRUE), max(barplot_data$relative_signal, na.rm = TRUE))
})
