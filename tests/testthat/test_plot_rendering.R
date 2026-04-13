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
