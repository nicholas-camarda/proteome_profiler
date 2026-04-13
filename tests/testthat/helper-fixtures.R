library(openxlsx)
library(readr)
library(tibble)
library(dplyr)

repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "scripts", "helpers", "runtime_setup.R"), local = TRUE)
load_analysis_packages(include_parallel = TRUE)

proteome_profiler_config <- list(
    runtime_root = tempdir(),
    cloud_parent = tempdir(),
    analyses = list()
)

source(file.path(repo_root, "scripts", "helpers", "replicate_analysis.R"), local = TRUE)
source(file.path(repo_root, "scripts", "helpers", "project_paths.R"), local = TRUE)
source(file.path(repo_root, "scripts", "helpers", "array_helper_scripts.R"), local = TRUE)

fixture_protocol_table <- function() {
    tibble(
        `Analyte/Control` = c(
            "Analyte A",
            "Analyte B",
            "Reference Spots",
            "Negative Control"
        ),
        Coordinate = c("A1, A2", "A3, A4", "J1, J2", "J23, J24")
    )
}

write_protocol_fixture <- function(path) {
    writexl::write_xlsx(fixture_protocol_table(), path)
}

write_licor_workbook <- function(path, analyte_signals) {
    stopifnot(length(analyte_signals) == 4)

    signal_rows <- unlist(map(analyte_signals, ~ c(.x, .x)), use.names = FALSE)
    workbook_df <- tibble(
        Name = c(sprintf("S%03d", seq_along(signal_rows)), "B001"),
        Signal = c(signal_rows, 1)
    )

    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    writeData(wb, "Sheet1", workbook_df, startRow = 4)
    saveWorkbook(wb, path, overwrite = TRUE)
}

create_legacy_fixture_dir <- function(root_dir, duplicate_group = FALSE) {
    dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
    write_licor_workbook(file.path(root_dir, "vehicle - mouse1.xlsx"), c(100, 90, 500, 10))
    write_licor_workbook(file.path(root_dir, "treated - mouse1.xlsx"), c(200, 45, 500, 10))

    if (duplicate_group) {
        write_licor_workbook(file.path(root_dir, "vehicle - mouse2.xlsx"), c(110, 92, 500, 10))
    }

    root_dir
}

create_replicate_fixture_dir <- function(root_dir) {
    dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

    write_licor_workbook(file.path(root_dir, "alpha.xlsx"), c(100, 60, 500, 20))
    write_licor_workbook(file.path(root_dir, "beta.xlsx"), c(102, 61, 500, 21))
    write_licor_workbook(file.path(root_dir, "gamma.xlsx"), c(250, 62, 500, 19))
    write_licor_workbook(file.path(root_dir, "delta.xlsx"), c(260, 63, 500, 18))
    write_licor_workbook(file.path(root_dir, "epsilon.xlsx"), c(100, 80, 500, 20))
    write_licor_workbook(file.path(root_dir, "zeta.xlsx"), c(99, 81, 500, 19))
    write_licor_workbook(file.path(root_dir, "eta.xlsx"), c(101, 79, 500, 18))
    write_licor_workbook(file.path(root_dir, "theta.xlsx"), c(100, 82, 500, 17))

    root_dir
}

create_replicate_manifest <- function(path, data_dir) {
    manifest <- tibble(
        sample_id = c("M01", "M02", "M03", "M04", "F01", "F02", "F03", "F04"),
        workbook_path = file.path(data_dir, c("alpha.xlsx", "beta.xlsx", "gamma.xlsx", "delta.xlsx", "epsilon.xlsx", "zeta.xlsx", "eta.xlsx", "theta.xlsx")),
        treatment = c("control", "control", "treated", "treated", "control", "control", "treated", "treated"),
        sex = c("male", "male", "male", "male", "female", "female", "female", "female")
    )
    write_csv(manifest, path)
    path
}

#' Create sample metadata for a larger replicate-aware stress fixture
#'
#' This fixture is intentionally more realistic than the smoke fixture: it uses
#' two subgroup levels, two treatment contrasts per subgroup, and uneven sample
#' counts across subgroup-treatment cells.
#'
#' @return Tibble with `sample_id`, `sex`, and `treatment`.
create_complex_replicate_sample_metadata <- function() {
    tibble(
        sample_id = c(
            "MC1", "MC2", "MC3",
            "MA1", "MA2", "MA3",
            "MB1", "MB2",
            "FC1", "FC2",
            "FA1", "FA2", "FA3",
            "FB1", "FB2"
        ),
        sex = c(rep("male", 8), rep("female", 7)),
        treatment = c(
            "control", "control", "control",
            "drug_a", "drug_a", "drug_a",
            "drug_b", "drug_b",
            "control", "control",
            "drug_a", "drug_a", "drug_a",
            "drug_b", "drug_b"
        )
    )
}

#' Create a larger synthetic sample-level replicate-aware dataset
#'
#' The signals are chosen so this fixture exercises multiple comparison
#' families, a raw-significant analyte that fails BH correction, uneven subgroup
#' sizes, low-signal flagging, and one analyte that becomes not testable due to
#' partial missingness in one subgroup/contrast.
#'
#' @return Tibble in the same shape expected by
#'   `run_within_stratum_differential_analysis()`.
create_complex_replicate_sample_data <- function() {
    sample_meta <- create_complex_replicate_sample_metadata()

    signals <- tribble(
        ~Name, ~Coordinate, ~MC1, ~MC2, ~MC3, ~MA1, ~MA2, ~MA3, ~MB1, ~MB2, ~FC1, ~FC2, ~FA1, ~FA2, ~FA3, ~FB1, ~FB2,
        "Analyte Strong A", "A1, A2", 100, 110, 90, 220, 230, 210, 100, 105, 100, 110, 210, 220, 205, 95, 100,
        "Analyte Moderate A", "A3, A4", 100, 110, 90, 135, 145, 120, 100, 105, 100, 103, 115, 140, 118, 100, 104,
        "Analyte Strong B", "A5, A6", 100, 105, 95, 100, 102, 98, 210, 220, 100, 100, 100, 102, 98, 150, 155,
        "Analyte Null", "A7, A8", 100, 102, 98, 101, 103, 99, 99, 100, 100, 99, 101, 100, 102, 99, 101,
        "Analyte Low", "A9, A10", 40, 42, 38, 55, 58, 52, 35, 36, 39, 41, 54, 56, 53, 34, 35,
        "Analyte Missing", "A11, A12", 100, 100, 100, 100, 100, 100, 150, 160, 100, 100, 100, 100, 100, 0, NA,
        "Analyte Down", "A13, A14", 100, 102, 98, 70, 72, 68, 100, 99, 100, 101, 80, 81, 79, 100, 102,
        "Null 2", "B1, B2", 100, 102, 98, 99, 101, 100, 101, 99, 100, 99, 101, 100, 102, 98, 101,
        "Null 3", "B3, B4", 100, 98, 102, 101, 99, 100, 100, 99, 100, 101, 99, 100, 98, 99, 100,
        "Null 4", "B5, B6", 95, 100, 105, 97, 101, 103, 95, 97, 99, 100, 97, 98, 102, 100, 99,
        "Null 5", "B7, B8", 110, 108, 112, 109, 110, 111, 111, 109, 112, 111, 110, 112, 111, 110, 109,
        "Null 6", "B9, B10", 85, 90, 88, 86, 89, 87, 86, 87, 84, 85, 83, 86, 84, 85, 86
    )

    signals %>%
        pivot_longer(-c(Name, Coordinate), names_to = "sample_id", values_to = "signal") %>%
        left_join(sample_meta, by = "sample_id") %>%
        mutate(
            Sname = Coordinate,
            workbook_path = "synthetic"
        )
}

write_analysis_config_fixture <- function(path, runtime_root, cloud_parent = tempfile("cloud-root-")) {
    dir.create(runtime_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(cloud_parent, recursive = TRUE, showWarnings = FALSE)

    protocol_path <- file.path(runtime_root, "output", "protocol.xlsx")
    legacy_data_dir <- file.path(runtime_root, "data", "legacy")
    replicate_data_dir <- file.path(runtime_root, "data", "replicate")
    replicate_manifest_path <- file.path(runtime_root, "manifests", "replicate.csv")

    dir.create(dirname(protocol_path), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(replicate_manifest_path), recursive = TRUE, showWarnings = FALSE)

    write_protocol_fixture(protocol_path)
    create_legacy_fixture_dir(legacy_data_dir)
    create_replicate_fixture_dir(replicate_data_dir)
    create_replicate_manifest(replicate_manifest_path, replicate_data_dir)

    config_obj <- list(
        runtime_root = runtime_root,
        cloud_parent = cloud_parent,
        analyses = list(
            legacy_smoke = list(
                mode = "legacy",
                user = "tester",
                analysis_slug = "legacy_smoke",
                protocol_preset = "cytokine_xl",
                info_fn = "output/protocol.xlsx",
                protocol_pdf = "protocols/mock_protocol.pdf",
                data_dir = "data/legacy",
                group_levels = c("vehicle", "treated"),
                comparisons = list("vehicle" = c("treated")),
                ref_thresh_to_filter = c(80),
                main_threshold = c(1.5),
                groups_per_page = 4,
                ref_coords_to_make_filter = c("A3,4"),
                selection_control = "vehicle",
                selection_group = "treated",
                selection_threshold = 1.5
            ),
            legacy_manifest_smoke = list(
                mode = "legacy",
                user = "tester",
                analysis_slug = "legacy_manifest_smoke",
                protocol_preset = "cytokine_xl",
                info_fn = "output/protocol.xlsx",
                protocol_pdf = "protocols/mock_protocol.pdf",
                sample_manifest = "manifests/replicate.csv",
                treatment_var = "treatment",
                comparisons = list("control" = c("treated")),
                ref_coords_to_make_filter = c("A3,4"),
                ref_thresh_to_filter = c(70),
                main_threshold = c(1.5),
                groups_per_page = 4,
                selection_control = "control",
                selection_group = "treated",
                selection_threshold = 1.5
            ),
            replicate_smoke = list(
                mode = "replicate",
                user = "tester",
                analysis_slug = "replicate_smoke",
                protocol_preset = "cytokine_xl",
                info_fn = "output/protocol.xlsx",
                protocol_pdf = "protocols/mock_protocol.pdf",
                sample_manifest = "manifests/replicate.csv",
                subgroup_var = "sex",
                treatment_var = "treatment",
                comparisons = list("control" = c("treated")),
                ref_coords_to_make_filter = c("A3,4"),
                ref_thresh_to_filter = c(70),
                min_reps_per_arm = 2,
                p_adjust_method = "BH",
                alpha = 0.05,
                selection_comparison_slugs = c("male_control_vs_treated", "female_control_vs_treated"),
                selection_top_n = 5
            )
        )
    )

    config_lines <- c(
        "proteome_profiler_config <-",
        capture.output(dput(config_obj))
    )
    writeLines(config_lines, path)

    invisible(list(
        config_path = path,
        runtime_root = runtime_root,
        legacy_data_dir = legacy_data_dir,
        replicate_data_dir = replicate_data_dir,
        replicate_manifest_path = replicate_manifest_path
    ))
}

run_repo_script <- function(script_path, analysis_name, config_path) {
    old_wd <- setwd(repo_root)
    on.exit(setwd(old_wd), add = TRUE)

    output <- system2(
        "Rscript",
        script_path,
        stdout = TRUE,
        stderr = TRUE,
        env = c(
            paste0("PROTEOME_PROFILER_ANALYSIS=", analysis_name),
            paste0("PROTEOME_PROFILER_CONFIG=", config_path)
        )
    )

    status <- attr(output, "status")
    list(
        output = output,
        status = if (is.null(status)) 0 else status
    )
}
