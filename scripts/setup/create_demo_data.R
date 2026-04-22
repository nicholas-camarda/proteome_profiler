source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = FALSE)

demo_root <- file.path("demo")
protocol_dir <- file.path(demo_root, "protocols")
manifest_dir <- file.path(demo_root, "manifests")
workbook_dir <- file.path(demo_root, "workbooks")

dir.create(protocol_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(workbook_dir, recursive = TRUE, showWarnings = FALSE)

protocol_tbl <- tibble::tibble(
    `Analyte/Control` = c(
        "Analyte A",
        "Analyte B",
        "Reference Spots",
        "Negative Control"
    ),
    Coordinate = c("A1, A2", "A3, A4", "J1, J2", "J23, J24")
)
writexl::write_xlsx(protocol_tbl, file.path(protocol_dir, "demo_protocol.xlsx"))
writeLines("Demo placeholder protocol PDF for setup validation.", file.path(protocol_dir, "demo_protocol.pdf"))

write_demo_workbook <- function(path, analyte_signals) {
    signal_rows <- unlist(lapply(analyte_signals, function(value) c(value, value)), use.names = FALSE)
    workbook_df <- tibble::tibble(
        Name = c(sprintf("S%03d", seq_along(signal_rows)), "B001"),
        Signal = c(signal_rows, 1)
    )

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Sheet1")
    openxlsx::writeData(wb, "Sheet1", workbook_df, startRow = 4)
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

sample_tbl <- tibble::tribble(
    ~sample_id, ~treatment, ~sex, ~signals,
    "M01", "control", "male", c(100, 60, 500, 20),
    "M02", "control", "male", c(102, 61, 500, 21),
    "M03", "treated", "male", c(250, 62, 500, 19),
    "M04", "treated", "male", c(260, 63, 500, 18),
    "F01", "control", "female", c(100, 80, 500, 20),
    "F02", "control", "female", c(99, 81, 500, 19),
    "F03", "treated", "female", c(101, 79, 500, 18),
    "F04", "treated", "female", c(100, 82, 500, 17)
)

manifest_tbl <- sample_tbl %>%
    mutate(
        workbook_path = file.path(workbook_dir, paste0(sample_id, ".xlsx"))
    )

purrr::walk2(manifest_tbl$workbook_path, manifest_tbl$signals, write_demo_workbook)

manifest_tbl %>%
    transmute(sample_id, workbook_path, treatment, sex) %>%
    readr::write_csv(file.path(manifest_dir, "demo_replicate_samples.csv"))

message("Demo data written under demo/.")
message("Copy .env.example to .env, then run Rscript scripts/check_setup.R.")
