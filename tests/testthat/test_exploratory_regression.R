test_that("find_filter_thresh uses selected low-signal coordinates as the cutoff reference", {
    tmp_dir <- tempfile("exploratory-thresh-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")

    write_protocol_fixture(protocol_path)
    create_exploratory_fixture_dir(data_dir)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    wide_df <- make_plot_ready_dataset(
        data_dir = data_dir,
        analyte_info = analyte_info,
        preview = FALSE,
        my_group_lvls = c("vehicle", "treated")
    )

    threshold_obj <- find_filter_thresh(
        wide_df = wide_df,
        ref_coords = "A3,4",
        my_colors = c(vehicle = "#4C78A8", treated = "#F58518")
    )

    expect_equal(threshold_obj[[1]]$mean[[1]], 67.5)
    expect_s3_class(threshold_obj[[2]], "patchwork")
})

test_that("make_wf_data keeps pairwise filtering local to the requested comparison", {
    tmp_dir <- tempfile("exploratory-wf-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")

    write_protocol_fixture(protocol_path)
    create_exploratory_fixture_dir(data_dir)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    wide_df <- make_plot_ready_dataset(
        data_dir = data_dir,
        analyte_info = analyte_info,
        preview = FALSE,
        my_group_lvls = c("vehicle", "treated")
    )

    assign("my_group_lvls", c("vehicle", "treated"), envir = environment(make_wf_data))

    wf_obj <- make_wf_data(
        data = wide_df,
        ref_thresh_to_filter_ = 80,
        my_main_threshold = 1.5,
        comparisons = list("vehicle" = c("treated"))
    )

    treated_wf <- wf_obj[["vehicle vs treated"]]$wf_dat[["treated"]]
    expect_equal(nrow(treated_wf), 1)
    expect_equal(as.character(treated_wf$Name[[1]]), "Analyte A")
    expect_true(treated_wf$significant[[1]])
    expect_equal(treated_wf$relative_signal[[1]], 2)
})

test_that("make_graphs writes waterfall and threshold-hit artifacts to the expected tree", {
    tmp_dir <- tempfile("exploratory-graphs-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")
    output_dir <- file.path(tmp_dir, "outputs")

    write_protocol_fixture(protocol_path)
    create_exploratory_fixture_dir(data_dir)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    wide_df <- make_plot_ready_dataset(
        data_dir = data_dir,
        analyte_info = analyte_info,
        preview = FALSE,
        my_group_lvls = c("vehicle", "treated")
    )

    assign("my_group_lvls", c("vehicle", "treated"), envir = environment(make_wf_data))
    assign("my_colors", c(vehicle = "#4C78A8", treated = "#F58518"), envir = environment(make_graphs))
    assign("my_outline_colors", c(vehicle = "#2F4A60", treated = "#8C4D11"), envir = environment(make_graphs))

    make_graphs(
        df = wide_df,
        ref_thresh_to_filter = 80,
        main_threshold = 1.5,
        output_dir_full_path = output_dir,
        groups_per_page = 4,
        comparisons = list("vehicle" = c("treated"))
    )

    expect_true(file.exists(file.path(output_dir, "ref_threshold_80", "all_comparisons", "waterfalls", "vehicle vs treated", "vehicle vs treated--treated-main_waterfall-all.tsv")))
    expect_true(file.exists(file.path(output_dir, "ref_threshold_80", "fold_change_hits", "threshold_1.5", "waterfalls", "vehicle vs treated", "vehicle vs treated--treated-main_waterfall-significant.tsv")))
    expect_true(file.exists(file.path(output_dir, "ref_threshold_80", "fold_change_hits", "threshold_1.5", "barplots", "vehicle vs treated", "combined", "vehicle vs treated-barplot_page_1__significant.png")))
})
