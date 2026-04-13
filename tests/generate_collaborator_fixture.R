library(openxlsx)
library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(stringr)

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "scripts", "helpers", "runtime_setup.R"), local = TRUE)
load_analysis_packages(include_parallel = FALSE)

#' Return the canonical root used for generated collaborator fixture outputs
#'
#' Everything for this synthetic run lives under one ignored tree so the input
#' workbooks, manifest, config, and generated outputs can be inspected together.
#'
#' @return Absolute path to `tests/test_output`.
get_fixture_root <- function() {
    normalizePath(file.path(repo_root, "tests", "test_output"), winslash = "/", mustWork = FALSE)
}

#' Build enough coordinate pairs for a synthetic Proteome Profiler protocol table
#'
#' The pipeline expects one protocol row per duplicate membrane-spot pair. This
#' helper emits display coordinates like `"A1, A2"` in a stable order.
#'
#' @param n_pairs Integer number of coordinate pairs required.
#'
#' @return Character vector of coordinate-pair labels.
generate_coordinate_pairs <- function(n_pairs) {
    row_labels <- c(LETTERS, paste0("A", LETTERS))
    coordinate_pairs <- character()

    for (row_label in row_labels) {
        for (col_start in seq(1, 24, by = 2)) {
            coordinate_pairs <- c(coordinate_pairs, sprintf("%s%s, %s%s", row_label, col_start, row_label, col_start + 1))
            if (length(coordinate_pairs) >= n_pairs) {
                return(coordinate_pairs)
            }
        }
    }

    stop(sprintf("Could not generate %s coordinate pairs.", n_pairs))
}

#' Build the synthetic protocol workbook table for the collaborator mock design
#'
#' The collaborator design uses 111 analytes plus protocol-style reference and
#' negative-control rows so the active pipeline behaves like a real array run.
#'
#' @return Tibble ready to be written as the protocol workbook.
build_protocol_table <- function() {
    analyte_names <- sprintf("Gene %03d", seq_len(111))
    control_names <- c("Reference Spots", "Negative Control")
    coordinate_pairs <- generate_coordinate_pairs(length(analyte_names) + length(control_names))

    tibble(
        `Analyte/Control` = c(analyte_names, control_names),
        Coordinate = coordinate_pairs
    )
}

#' Return the sample manifest for the collaborator experimental design
#'
#' The design is two sexes by two treatment arms with four biological replicates
#' per arm. Each biological replicate is represented by one LI-COR workbook.
#'
#' @return Tibble with `sample_id`, `sex`, `treatment`, and `workbook_path`.
build_sample_manifest <- function() {
    crossing(
        sex = c("male", "female"),
        treatment = c("vehicle", "aldosterone"),
        replicate = seq_len(4)
    ) %>%
        mutate(
            sample_id = sprintf(
                "%s_%s_%02d",
                if_else(sex == "male", "M", "F"),
                if_else(treatment == "vehicle", "VEH", "ALDO"),
                replicate
            ),
            workbook_path = file.path("tests", "test_output", "workbooks", sprintf("%s.xlsx", sample_id))
        ) %>%
        select(sample_id, workbook_path, treatment, sex)
}

#' Assign each synthetic analyte to a signal-pattern category
#'
#' The categories create a collaborator-facing mock run with strong male-only and
#' female-only treatment effects, shared effects, low-signal analytes, and nulls.
#'
#' @param gene_names Character vector of analyte names.
#'
#' @return Named character vector of category labels.
assign_gene_categories <- function(gene_names) {
    stopifnot(length(gene_names) == 111)

    categories <- c(
        rep("low_signal", 8),
        rep("male_up_strong", 10),
        rep("male_down_strong", 8),
        rep("female_up_strong", 10),
        rep("female_down_strong", 8),
        rep("both_up", 8),
        rep("both_down", 8),
        rep("male_up_moderate", 8),
        rep("female_up_moderate", 8),
        rep("null", 35)
    )

    setNames(categories, gene_names)
}

#' Simulate sample-level analyte intensities for the collaborator mock design
#'
#' The generated signals are randomized but reproducible. They are designed to
#' yield a mixture of significant hits, null analytes, and low-signal analytes
#' under within-sex BH-adjusted inference.
#'
#' @param manifest Tibble returned by `build_sample_manifest()`.
#' @param protocol_tbl Tibble returned by `build_protocol_table()`.
#' @param seed Integer random seed.
#'
#' @return Named list mapping `sample_id` to per-analyte numeric signals.
simulate_sample_signals <- function(manifest, protocol_tbl, seed = 20260413) {
    set.seed(seed)

    gene_tbl <- protocol_tbl %>%
        filter(!str_starts(`Analyte/Control`, "Reference")) %>%
        filter(!str_starts(`Analyte/Control`, "Negative")) %>%
        transmute(
            gene = `Analyte/Control`,
            category = assign_gene_categories(`Analyte/Control`),
            baseline = case_when(
                category == "low_signal" ~ runif(n(), min = 28, max = 55),
                TRUE ~ runif(n(), min = 95, max = 220)
            )
        )

    effect_multiplier <- function(category, sex, treatment) {
        if (treatment == "vehicle") {
            return(1)
        }

        switch(
            category,
            low_signal = 1,
            male_up_strong = if (sex == "male") 2.3 else 1.02,
            male_down_strong = if (sex == "male") 0.48 else 1.02,
            female_up_strong = if (sex == "female") 2.2 else 1.03,
            female_down_strong = if (sex == "female") 0.52 else 1.03,
            both_up = 1.7,
            both_down = 0.62,
            male_up_moderate = if (sex == "male") 1.35 else 1.02,
            female_up_moderate = if (sex == "female") 1.32 else 1.02,
            null = 1.0,
            1.0
        )
    }

    sample_ids <- manifest$sample_id
    sample_signal_list <- vector("list", length(sample_ids))
    names(sample_signal_list) <- sample_ids

    for (row_idx in seq_len(nrow(manifest))) {
        sex_value <- manifest$sex[[row_idx]]
        treatment_value <- manifest$treatment[[row_idx]]
        sample_id_value <- manifest$sample_id[[row_idx]]

        noise_sd <- if (treatment_value == "vehicle") 0.11 else 0.13
        signals <- map_dbl(seq_len(nrow(gene_tbl)), function(gene_idx) {
            gene_row <- gene_tbl[gene_idx, , drop = FALSE]
            mean_signal <- gene_row$baseline[[1]] * effect_multiplier(
                category = gene_row$category[[1]],
                sex = sex_value,
                treatment = treatment_value
            )
            max(0.5, mean_signal * exp(rnorm(1, mean = 0, sd = noise_sd)))
        })

        names(signals) <- gene_tbl$gene
        sample_signal_list[[sample_id_value]] <- signals
    }

    sample_signal_list
}

#' Write one LI-COR-style workbook for a biological sample
#'
#' The active parser expects sequential `S001...` spot labels, duplicate rows per
#' analyte, and a background row after the analyte/control spots.
#'
#' @param path Output workbook path.
#' @param signal_vector Numeric vector with one entry per protocol row.
#' @param seed Integer seed used only for small duplicate-spot jitter.
#'
#' @return Invisibly returns `path`.
write_licor_workbook <- function(path, signal_vector, seed) {
    set.seed(seed)

    duplicated_rows <- unlist(map(signal_vector, function(signal_value) {
        jitter <- rnorm(2, mean = 0, sd = max(0.5, 0.02 * signal_value))
        pmax(0.5, signal_value + jitter)
    }), use.names = FALSE)

    workbook_df <- tibble(
        Name = c(sprintf("S%03d", seq_along(duplicated_rows)), "B001"),
        Signal = c(duplicated_rows, 1)
    )

    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    writeData(wb, "Sheet1", workbook_df, startRow = 4)
    saveWorkbook(wb, path, overwrite = TRUE)

    invisible(path)
}

quote_string <- function(text_value) {
    sprintf("\"%s\"", text_value)
}

#' Write the synthetic collaborator analysis config file
#'
#' @param fixture_root Absolute fixture-root path.
#' @param low_signal_coords Character vector of coordinate pairs to use as the
#'   low-signal reference panel.
#' @param user_name Character user name for the output subtree.
#'
#' @return Invisibly returns the config path.
write_collaborator_config <- function(fixture_root, low_signal_coords, user_name = "nicole") {
    config_path <- file.path(fixture_root, "analysis_config.R")

    config_lines <- c(
        "proteome_profiler_config <- list(",
        sprintf("    runtime_root = %s,", quote_string(fixture_root)),
        sprintf("    cloud_parent = %s,", quote_string(fixture_root)),
        "    default_analysis = \"collaborator_aldo_mock\",",
        "    analyses = list(",
        "        collaborator_aldo_mock = list(",
        "            mode = \"replicate\",",
        sprintf("            user = %s,", quote_string(user_name)),
        "            slug = \"aldo_plasma_mock_111_gene\",",
        "            protocol = list(",
        "                preset = \"cytokine_xl\",",
        "                workbook = \"tests/test_output/output/mock_protocol.xlsx\"",
        "            ),",
        "            input = list(",
        "                manifest = \"tests/test_output/manifests/collaborator_aldo_samples.csv\",",
        "                subgroup = \"sex\",",
        "                treatment = \"treatment\"",
        "            ),",
        "            comparisons = list(",
        "                vehicle = c(\"aldosterone\")",
        "            ),",
        "            thresholds = list(",
        sprintf("                ref_coords = c(%s),", paste(vapply(low_signal_coords, quote_string, character(1)), collapse = ", ")),
        "                ref_signal = c(70)",
        "            ),",
        "            stats = list(",
        "                min_reps_per_arm = 2,",
        "                p_adjust_method = \"BH\",",
        "                alpha = 0.05",
        "            ),",
        "            shortlist = list(",
        "                comparisons = c(\"male_vehicle_vs_aldosterone\", \"female_vehicle_vs_aldosterone\"),",
        "                top_n = 10",
        "            )",
        "        )",
        "    )",
        ")"
    )
    writeLines(config_lines, config_path)

    invisible(config_path)
}

#' Materialize all synthetic collaborator input files under `tests/test_output`
#'
#' @param fixture_root Absolute fixture-root path.
#'
#' @return Named list of important generated input paths and metadata.
write_collaborator_fixture_inputs <- function(fixture_root) {
    if (dir.exists(fixture_root)) {
        unlink(fixture_root, recursive = TRUE, force = TRUE)
    }

    dir.create(fixture_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "workbooks"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "manifests"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "output"), recursive = TRUE, showWarnings = FALSE)

    protocol_tbl <- build_protocol_table()
    manifest_tbl <- build_sample_manifest()
    signal_list <- simulate_sample_signals(manifest_tbl, protocol_tbl)

    protocol_path <- file.path(fixture_root, "output", "mock_protocol.xlsx")
    write_xlsx(protocol_tbl, protocol_path)

    for (sample_idx in seq_len(nrow(manifest_tbl))) {
        sample_id_value <- manifest_tbl$sample_id[[sample_idx]]
        workbook_abs_path <- file.path(repo_root, manifest_tbl$workbook_path[[sample_idx]])

        signal_vector <- c(
            signal_list[[sample_id_value]],
            "Reference Spots" = runif(1, min = 450, max = 560),
            "Negative Control" = runif(1, min = 8, max = 18)
        )

        write_licor_workbook(
            path = workbook_abs_path,
            signal_vector = unname(signal_vector),
            seed = 1000 + sample_idx
        )
    }

    manifest_path <- file.path(fixture_root, "manifests", "collaborator_aldo_samples.csv")
    write_csv(manifest_tbl, manifest_path)

    initial_low_signal_coords <- protocol_tbl %>%
        slice(1:8) %>%
        pull(Coordinate) %>%
        str_replace_all(", ", ",")

    config_path <- write_collaborator_config(
        fixture_root = fixture_root,
        low_signal_coords = initial_low_signal_coords,
        user_name = "nicole"
    )

    list(
        fixture_root = fixture_root,
        protocol_path = protocol_path,
        manifest_path = manifest_path,
        config_path = config_path,
        manifest_tbl = manifest_tbl,
        protocol_tbl = protocol_tbl,
        initial_low_signal_coords = initial_low_signal_coords
    )
}

#' Run one repo entry script against the generated collaborator fixture config
#'
#' @param script_path Relative repo path to the entry script.
#' @param config_path Absolute path to the generated config file.
#'
#' @return Character vector of combined stdout/stderr.
run_repo_script <- function(script_path, config_path) {
    output <- system2(
        "Rscript",
        script_path,
        stdout = TRUE,
        stderr = TRUE,
        env = c(
            "PROTEOME_PROFILER_ANALYSIS=collaborator_aldo_mock",
            paste0("PROTEOME_PROFILER_CONFIG=", config_path)
        )
    )

    status <- attr(output, "status")
    if (!is.null(status) && status != 0) {
        stop(sprintf(
            "Script %s failed with status %s\n%s",
            script_path,
            status,
            paste(output, collapse = "\n")
        ))
    }

    output
}

#' Write a collaborator-facing summary of the generated synthetic run
#'
#' @param fixture_root Absolute fixture-root path.
#'
#' @return Invisibly returns the summary path.
write_fixture_summary <- function(fixture_root, selected_low_signal_coords) {
    analysis_root <- file.path(fixture_root, "output", "nicole", "aldo_plasma_mock_111_gene")
    run_index_path <- file.path(analysis_root, "inferential_results", "run_index.tsv")
    run_index <- read_tsv(run_index_path, show_col_types = FALSE)

    male_shortlist_summary_path <- file.path(
        analysis_root,
        "select_analytes", "comparisons", "male_vehicle_vs_aldosterone", "shortlist_summary.tsv"
    )
    female_shortlist_summary_path <- file.path(
        analysis_root,
        "select_analytes", "comparisons", "female_vehicle_vs_aldosterone", "shortlist_summary.tsv"
    )

    male_shortlist_summary <- read_tsv(male_shortlist_summary_path, show_col_types = FALSE)
    female_shortlist_summary <- read_tsv(female_shortlist_summary_path, show_col_types = FALSE)

    summary_lines <- c(
        "# Collaborator Mock Output",
        "",
        "Synthetic design:",
        "- 2 sexes (`male`, `female`)",
        "- 2 treatments (`vehicle`, `aldosterone`)",
        "- 4 biological replicates per sex-treatment cell",
        "- 111 analytes plus protocol-style reference/control rows",
        sprintf("- selected low-signal coordinates: %s", paste(selected_low_signal_coords, collapse = ", ")),
        "",
        sprintf("Config: `%s`", file.path(fixture_root, "analysis_config.R")),
        sprintf("Manifest: `%s`", file.path(fixture_root, "manifests", "collaborator_aldo_samples.csv")),
        sprintf("Protocol workbook: `%s`", file.path(fixture_root, "output", "mock_protocol.xlsx")),
        sprintf("Analysis root: `%s`", analysis_root),
        "",
        "Key outputs:",
        sprintf("- threshold diagnostics: `%s`", file.path(analysis_root, "threshold_diagnostics", "region_stats.png")),
        sprintf("- run index: `%s`", run_index_path),
        sprintf("- combined inferential workbook: `%s`", file.path(analysis_root, "inferential_results", "combined_results.xlsx")),
        sprintf("- male comparison waterfall: `%s`", file.path(analysis_root, "inferential_results", "comparisons", "male_vehicle_vs_aldosterone", "waterfall.png")),
        sprintf("- female comparison waterfall: `%s`", file.path(analysis_root, "inferential_results", "comparisons", "female_vehicle_vs_aldosterone", "waterfall.png")),
        sprintf("- male shortlist waterfall: `%s`", file.path(analysis_root, "select_analytes", "comparisons", "male_vehicle_vs_aldosterone", "shortlist_waterfall.png")),
        sprintf("- female shortlist waterfall: `%s`", file.path(analysis_root, "select_analytes", "comparisons", "female_vehicle_vs_aldosterone", "shortlist_waterfall.png")),
        "",
        "Run index preview:",
        paste(capture.output(print(run_index)), collapse = "\n"),
        "",
        "Male shortlist summary:",
        paste(capture.output(print(male_shortlist_summary)), collapse = "\n"),
        "",
        "Female shortlist summary:",
        paste(capture.output(print(female_shortlist_summary)), collapse = "\n")
    )

    summary_path <- file.path(fixture_root, "README.md")
    writeLines(summary_lines, summary_path)

    invisible(summary_path)
}

fixture_root <- get_fixture_root()
fixture_inputs <- write_collaborator_fixture_inputs(fixture_root)

old_wd <- setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

run_repo_script("scripts/find_ref_thresh.R", fixture_inputs$config_path)
candidate_path <- file.path(
    fixture_root,
    "output", "nicole", "aldo_plasma_mock_111_gene",
    "threshold_diagnostics", "candidate_low_signal_analytes.tsv"
)
suggested_low_signal_coords <- read_tsv(candidate_path, show_col_types = FALSE) %>%
    slice_head(n = 8) %>%
    pull(Coordinate)
write_collaborator_config(
    fixture_root = fixture_root,
    low_signal_coords = suggested_low_signal_coords,
    user_name = "nicole"
)
run_repo_script("scripts/find_ref_thresh.R", fixture_inputs$config_path)
run_repo_script("scripts/main.R", fixture_inputs$config_path)
run_repo_script("scripts/select-analytes-analysis.R", fixture_inputs$config_path)
write_fixture_summary(fixture_root, suggested_low_signal_coords)

cat("Generated collaborator mock fixture and outputs under:\n")
cat(fixture_root, "\n")
