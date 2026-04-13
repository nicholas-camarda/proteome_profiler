test_that("legacy import averages duplicate spots into one wide analyte row per group", {
    tmp_dir <- tempfile("legacy-valid-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")

    write_protocol_fixture(protocol_path)
    create_legacy_fixture_dir(data_dir)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    out <- make_plot_ready_dataset(
        data_dir = data_dir,
        analyte_info = analyte_info,
        preview = FALSE,
        my_group_lvls = c("vehicle", "treated")
    )

    expect_true(all(c("vehicle", "treated") %in% names(out)))
    expect_equal(nrow(out), 4)
    expect_equal(out$vehicle[out$Name == "Analyte A"], 100)
    expect_equal(out$treated[out$Name == "Analyte A"], 200)
})

test_that("legacy import rejects duplicate group files", {
    tmp_dir <- tempfile("legacy-dup-")
    dir.create(tmp_dir)
    protocol_path <- file.path(tmp_dir, "protocol.xlsx")
    data_dir <- file.path(tmp_dir, "data")

    write_protocol_fixture(protocol_path)
    create_legacy_fixture_dir(data_dir, duplicate_group = TRUE)

    analyte_info <- readxl::read_excel(protocol_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    expect_error(
        make_plot_ready_dataset(
            data_dir = data_dir,
            analyte_info = analyte_info,
            preview = FALSE,
            my_group_lvls = c("vehicle", "treated")
        ),
        "multiple LI-COR workbooks"
    )
})
