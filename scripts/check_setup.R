rm(list = ls())

source(file.path("scripts", "helpers", "runtime_setup.R"))

missing_packages <- missing_analysis_packages(include_parallel = TRUE)
if (length(missing_packages) > 0) {
    stop(sprintf(
        paste(
            "Missing required R package(s): %s",
            "Install them with: Rscript scripts/install_packages.R"
        ),
        paste(missing_packages, collapse = ", ")
    ), call. = FALSE)
}

load_analysis_packages(include_parallel = TRUE)

source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "replicate_analysis.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

check_command <- function(command, args = character(), missing_message) {
    result <- tryCatch(
        system2(command, args, stdout = TRUE, stderr = TRUE),
        error = function(error) structure(conditionMessage(error), status = 127)
    )
    status <- attr(result, "status")
    if (is.null(status)) {
        return(invisible(TRUE))
    }
    stop(missing_message, call. = FALSE)
}

message("Checking Proteome Profiler setup...")

check_command(
    "python3",
    c("-c", shQuote("import pandas, tabula, openpyxl")),
    paste(
        "Missing Python protocol-extraction dependency.",
        "Install Python requirements with: python3 -m pip install -r requirements.txt"
    )
)
check_command(
    "java",
    "-version",
    paste(
        "Java is required by tabula-py for protocol PDF extraction.",
        "On macOS with Homebrew, install Java with: brew install --cask temurin",
        "Then rerun: java -version"
    )
)

initialize_runtime_config_from_env(required_env_file = TRUE)
analysis_name <- get_selected_analysis_name()
example_config <- get_analysis_config(analysis_name)
analysis_mode <- get_analysis_mode(example_config)

validate_analysis_config(example_config)

protocol_workbook_path <- get_protocol_workbook_path(example_config)
protocol_pdf_path <- if (!is.null(example_config$protocol_pdf) && !identical(example_config$protocol_pdf, "")) {
    resolve_project_path(example_config$protocol_pdf, must_exist = TRUE)
} else {
    NA_character_
}

analysis_output_root <- get_analysis_output_root(example_config)
if (!dir.exists(analysis_output_root)) {
    output_created <- suppressWarnings(dir.create(analysis_output_root, recursive = TRUE, showWarnings = FALSE))
    if (!isTRUE(output_created) && !dir.exists(analysis_output_root)) {
        stop(sprintf(
            paste(
                "Could not create analysis output root:",
                "%s",
                "Check PROTEOME_PROFILER_RUNTIME_ROOT in .env and confirm the folder is writable."
            ),
            analysis_output_root
        ), call. = FALSE)
    }
}
write_test_path <- tempfile("write-test-", tmpdir = analysis_output_root)
tryCatch(
    writeLines("ok", write_test_path),
    error = function(error) {
        stop(sprintf(
            paste(
                "Analysis output root is not writable:",
                "%s",
                "Reason: %s"
            ),
            analysis_output_root,
            conditionMessage(error)
        ), call. = FALSE)
    }
)
unlink(write_test_path, force = TRUE)

manifest_path <- NA_character_
data_dir <- NA_character_
resolved_manifest <- NULL
reference_spot_summary <- NULL

if (config_uses_sample_manifest(example_config)) {
    manifest_path <- resolve_project_path(example_config$sample_manifest, must_exist = TRUE)
    data_dir <- if (!is.null(example_config$data_dir) && !identical(example_config$data_dir, "")) {
        resolve_project_path(example_config$data_dir, must_exist = TRUE)
    } else {
        NA_character_
    }

    resolved_manifest <- read_sample_manifest(
        manifest_path = manifest_path,
        subgroup_var = example_config$subgroup_var,
        treatment_var = example_config$treatment_var
    ) %>%
        resolve_manifest_workbook_paths(data_dir = if (is.na(data_dir)) NULL else data_dir)

    analyte_info_df <- read_excel(protocol_workbook_path) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    sample_level_df <- build_sample_level_dataset(
        manifest = resolved_manifest,
        analyte_info = analyte_info_df,
        treatment_var = example_config$treatment_var,
        subgroup_var = example_config$subgroup_var,
        data_dir = if (is.na(data_dir)) NULL else data_dir
    )

    reference_spot_summary <- summarize_reference_spot_qc(sample_level_df) %>%
        summarize_reference_spot_qc_counts()
} else {
    data_dir <- resolve_project_path(example_config$data_dir, must_exist = TRUE)
}

message("Setup check passed.")
message(sprintf("Active analysis: %s", analysis_name))
message(sprintf("Analysis mode: %s", analysis_mode))
message(sprintf("Runtime root: %s", path.expand(proteome_profiler_config$runtime_root)))
message(sprintf("Output root: %s", analysis_output_root))
message(sprintf("Protocol workbook: %s", protocol_workbook_path))
if (!is.na(protocol_pdf_path)) {
    message(sprintf("Protocol PDF: %s", protocol_pdf_path))
}
if (!is.na(manifest_path)) {
    message(sprintf("Manifest: %s", manifest_path))
    message(sprintf("Samples: %s", nrow(resolved_manifest)))
    if ("resolved_sheet_name" %in% names(resolved_manifest)) {
        message(sprintf("Workbook sheets checked: %s", nrow(distinct(resolved_manifest, resolved_workbook_path, resolved_sheet_name))))
    }
    message(sprintf(
        "Reference spot QC: %s/%s samples used complete preferred A1,2/J1,2 raw reference spots; %s partial preferred; %s protocol-table fallback.",
        reference_spot_summary$n_preferred_complete[[1]],
        reference_spot_summary$n_samples[[1]],
        reference_spot_summary$n_preferred_partial[[1]],
        reference_spot_summary$n_protocol_fallback[[1]]
    ))
    if (reference_spot_summary$n_reference_qc_issues[[1]] > 0) {
        message(sprintf(
            "Reference spot QC samples requiring review: %s",
            reference_spot_summary$n_reference_qc_issues[[1]]
        ))
    }
} else {
    message(sprintf("Data directory: %s", data_dir))
}
message(sprintf("Comparisons: %s", paste(names(example_config$comparisons), "=", vapply(example_config$comparisons, paste, collapse = "|", FUN.VALUE = character(1)), collapse = "; ")))
message("Next command: Rscript scripts/find_ref_thresh.R")
