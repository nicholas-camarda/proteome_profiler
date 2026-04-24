library(openxlsx)
library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(stringr)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "scripts", "helpers", "runtime_setup.R"), local = TRUE)
load_analysis_packages(include_parallel = FALSE)

docs_fixture_root <- function() {
    normalizePath(
        file.path(repo_root, "tests", "test_output", "docs_example_fixture"),
        winslash = "/",
        mustWork = FALSE
    )
}

compact_coordinate <- function(value) {
    value %>%
        str_replace_all("\\s+", "") %>%
        str_replace_all("([A-Za-z])(\\d+),([A-Za-z])(\\d+)", "\\1\\2,\\4") %>%
        str_to_upper()
}

generate_coordinate_pairs <- function(n_pairs) {
    row_labels <- c(LETTERS, paste0("A", LETTERS))
    coordinate_pairs <- character()

    for (row_label in row_labels) {
        for (col_start in seq(1, 24, by = 2)) {
            coordinate_pairs <- c(
                coordinate_pairs,
                sprintf("%s%s, %s%s", row_label, col_start, row_label, col_start + 1)
            )
            if (length(coordinate_pairs) >= n_pairs) {
                return(coordinate_pairs)
            }
        }
    }

    stop(sprintf("Could not generate %s coordinate pairs.", n_pairs))
}

build_protocol_table <- function() {
    analyte_names <- sprintf("Marker %03d", seq_len(112))
    control_names <- c("Reference Spots", "Negative Control")
    coordinate_pairs <- generate_coordinate_pairs(length(analyte_names) + length(control_names))

    tibble(
        `Analyte/Control` = c(analyte_names, control_names),
        Coordinate = coordinate_pairs
    )
}

assign_exploratory_categories <- function(gene_names) {
    categories <- c(
        rep("low_signal", 10),
        rep("group_b_up_strong", 18),
        rep("group_b_down_strong", 14),
        rep("group_b_up_moderate", 12),
        rep("group_b_down_moderate", 10),
        rep("group_b_up_noisy", 8),
        rep("null", length(gene_names) - 72)
    )

    setNames(categories, gene_names)
}

assign_replicate_categories <- function(gene_names) {
    categories <- c(
        rep("low_signal", 8),
        rep("stratum_1_up_strong", 12),
        rep("stratum_1_down_strong", 10),
        rep("stratum_2_up_strong", 12),
        rep("stratum_2_down_strong", 10),
        rep("both_up", 10),
        rep("both_down", 8),
        rep("stratum_1_up_moderate", 8),
        rep("stratum_2_up_moderate", 8),
        rep("null", length(gene_names) - 86)
    )

    setNames(categories, gene_names)
}

build_exploratory_manifest <- function(fixture_root) {
    crossing(
        treatment = c("Group A", "Group B"),
        replicate = seq_len(4)
    ) %>%
        mutate(
            sample_id = sprintf(
                "%s_%02d",
                if_else(treatment == "Group A", "GA", "GB"),
                replicate
            ),
            workbook_path = file.path(
                "tests",
                "test_output",
                "docs_example_fixture",
                "workbooks",
                "exploratory",
                sprintf("%s.xlsx", sample_id)
            )
        ) %>%
        select(sample_id, workbook_path, treatment)
}

build_replicate_manifest <- function(fixture_root) {
    crossing(
        stratum = c("Stratum 1", "Stratum 2"),
        treatment = c("Group A", "Group B"),
        replicate = seq_len(4)
    ) %>%
        mutate(
            sample_id = sprintf(
                "%s_%s_%02d",
                if_else(stratum == "Stratum 1", "S1", "S2"),
                if_else(treatment == "Group A", "GA", "GB"),
                replicate
            ),
            workbook_path = file.path(
                "tests",
                "test_output",
                "docs_example_fixture",
                "workbooks",
                "replicate",
                sprintf("%s.xlsx", sample_id)
            )
        ) %>%
        select(sample_id, workbook_path, treatment, stratum)
}

simulate_exploratory_signals <- function(manifest, protocol_tbl, seed = 20260424) {
    set.seed(seed)

    gene_tbl <- protocol_tbl %>%
        filter(!`Analyte/Control` %in% c("Reference Spots", "Negative Control")) %>%
        transmute(
            gene = `Analyte/Control`,
            category = assign_exploratory_categories(`Analyte/Control`),
            baseline = case_when(
                category == "low_signal" ~ runif(n(), min = 55, max = 165),
                TRUE ~ runif(n(), min = 260, max = 980)
            )
        )

    effect_multiplier <- function(category, treatment_label) {
        if (identical(treatment_label, "Group A")) {
            return(1)
        }

        switch(
            category,
            low_signal = 1,
            group_b_up_strong = 2.35,
            group_b_down_strong = 0.46,
            group_b_up_moderate = 1.55,
            group_b_down_moderate = 0.7,
            group_b_up_noisy = 1.75,
            null = 1,
            1
        )
    }

    sample_signal_list <- vector("list", nrow(manifest))
    names(sample_signal_list) <- manifest$sample_id

    for (row_idx in seq_len(nrow(manifest))) {
        treatment_value <- manifest$treatment[[row_idx]]
        sample_id_value <- manifest$sample_id[[row_idx]]

        signals <- map_dbl(seq_len(nrow(gene_tbl)), function(gene_idx) {
            gene_row <- gene_tbl[gene_idx, , drop = FALSE]
            mean_signal <- gene_row$baseline[[1]] * effect_multiplier(gene_row$category[[1]], treatment_value)
            noise_sd <- if (identical(gene_row$category[[1]], "group_b_up_noisy")) 0.22 else 0.12
            max(0.5, mean_signal * exp(rnorm(1, mean = 0, sd = noise_sd)))
        })

        names(signals) <- gene_tbl$gene
        sample_signal_list[[sample_id_value]] <- signals
    }

    sample_signal_list
}

simulate_replicate_signals <- function(manifest, protocol_tbl, seed = 20260425) {
    set.seed(seed)

    gene_tbl <- protocol_tbl %>%
        filter(!`Analyte/Control` %in% c("Reference Spots", "Negative Control")) %>%
        transmute(
            gene = `Analyte/Control`,
            category = assign_replicate_categories(`Analyte/Control`),
            baseline = case_when(
                category == "low_signal" ~ runif(n(), min = 24, max = 52),
                TRUE ~ runif(n(), min = 100, max = 245)
            )
        )

    effect_multiplier <- function(category, stratum_label, treatment_label) {
        if (identical(treatment_label, "Group A")) {
            return(1)
        }

        switch(
            category,
            low_signal = 1,
            stratum_1_up_strong = if (identical(stratum_label, "Stratum 1")) 2.4 else 1.03,
            stratum_1_down_strong = if (identical(stratum_label, "Stratum 1")) 0.44 else 1.03,
            stratum_2_up_strong = if (identical(stratum_label, "Stratum 2")) 2.25 else 1.04,
            stratum_2_down_strong = if (identical(stratum_label, "Stratum 2")) 0.48 else 1.04,
            both_up = 1.72,
            both_down = 0.63,
            stratum_1_up_moderate = if (identical(stratum_label, "Stratum 1")) 1.45 else 1.02,
            stratum_2_up_moderate = if (identical(stratum_label, "Stratum 2")) 1.42 else 1.02,
            null = 1,
            1
        )
    }

    sample_signal_list <- vector("list", nrow(manifest))
    names(sample_signal_list) <- manifest$sample_id

    for (row_idx in seq_len(nrow(manifest))) {
        stratum_value <- manifest$stratum[[row_idx]]
        treatment_value <- manifest$treatment[[row_idx]]
        sample_id_value <- manifest$sample_id[[row_idx]]

        signals <- map_dbl(seq_len(nrow(gene_tbl)), function(gene_idx) {
            gene_row <- gene_tbl[gene_idx, , drop = FALSE]
            mean_signal <- gene_row$baseline[[1]] * effect_multiplier(
                gene_row$category[[1]],
                stratum_value,
                treatment_value
            )
            max(0.5, mean_signal * exp(rnorm(1, mean = 0, sd = 0.13)))
        })

        names(signals) <- gene_tbl$gene
        sample_signal_list[[sample_id_value]] <- signals
    }

    sample_signal_list
}

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

quoted_env_value <- function(value) {
    if (is.null(value) || length(value) == 0 || all(is.na(value)) || all(value == "")) {
        return(NULL)
    }

    text_value <- paste(as.character(value), collapse = "")
    if (grepl("[[:space:]|;=#]", text_value)) {
        return(sprintf("\"%s\"", gsub("\"", "\\\\\"", text_value)))
    }

    text_value
}

write_env_file <- function(path, values) {
    env_lines <- imap_chr(values, function(value, name) {
        rendered_value <- quoted_env_value(value)
        if (is.null(rendered_value)) {
            return(NA_character_)
        }
        sprintf("%s=%s", name, rendered_value)
    })

    writeLines(stats::na.omit(env_lines), path)
    invisible(path)
}

build_exploratory_env_values <- function(fixture_root, low_signal_coords, selected_coords) {
    list(
        PROTEOME_PROFILER_ANALYSIS = "docs_exploratory_examples",
        PROTEOME_PROFILER_MODE = "exploratory",
        PROTEOME_PROFILER_USER = "docs",
        PROTEOME_PROFILER_SLUG = "docs_exploratory_examples",
        PROTEOME_PROFILER_RUNTIME_ROOT = fixture_root,
        PROTEOME_PROFILER_CLOUD_PARENT = fixture_root,
        PROTEOME_PROFILER_PROTOCOL_PRESET = "cytokine_xl",
        PROTEOME_PROFILER_PROTOCOL_WORKBOOK = "tests/test_output/docs_example_fixture/input/mock_protocol.xlsx",
        PROTEOME_PROFILER_INPUT_MANIFEST = "tests/test_output/docs_example_fixture/manifests/exploratory_samples.csv",
        PROTEOME_PROFILER_TREATMENT_COLUMN = "treatment",
        PROTEOME_PROFILER_COMPARISONS = "Group A=Group B",
        PROTEOME_PROFILER_REF_COORDS = paste(low_signal_coords, collapse = "|"),
        PROTEOME_PROFILER_REF_SIGNAL = "70",
        PROTEOME_PROFILER_FOLD_CHANGE = "1.5",
        PROTEOME_PROFILER_GROUPS_PER_PAGE = "25",
        PROTEOME_PROFILER_SHORTLIST_COMPARISONS = "group_a_vs_group_b",
        PROTEOME_PROFILER_SHORTLIST_FOLD_CHANGE = "1.5",
        PROTEOME_PROFILER_SHORTLIST_COORDS = paste(selected_coords, collapse = "|")
    )
}

build_replicate_env_values <- function(fixture_root, low_signal_coords, selected_coords) {
    list(
        PROTEOME_PROFILER_ANALYSIS = "docs_replicate_examples",
        PROTEOME_PROFILER_MODE = "replicate",
        PROTEOME_PROFILER_USER = "docs",
        PROTEOME_PROFILER_SLUG = "docs_replicate_examples",
        PROTEOME_PROFILER_RUNTIME_ROOT = fixture_root,
        PROTEOME_PROFILER_CLOUD_PARENT = fixture_root,
        PROTEOME_PROFILER_PROTOCOL_PRESET = "cytokine_xl",
        PROTEOME_PROFILER_PROTOCOL_WORKBOOK = "tests/test_output/docs_example_fixture/input/mock_protocol.xlsx",
        PROTEOME_PROFILER_INPUT_MANIFEST = "tests/test_output/docs_example_fixture/manifests/replicate_samples.csv",
        PROTEOME_PROFILER_TREATMENT_COLUMN = "treatment",
        PROTEOME_PROFILER_SUBGROUP_COLUMN = "stratum",
        PROTEOME_PROFILER_COMPARISONS = "Group A=Group B",
        PROTEOME_PROFILER_REF_COORDS = paste(low_signal_coords, collapse = "|"),
        PROTEOME_PROFILER_REF_SIGNAL = "70",
        PROTEOME_PROFILER_MIN_REPS_PER_ARM = "2",
        PROTEOME_PROFILER_P_ADJUST_METHOD = "BH",
        PROTEOME_PROFILER_ALPHA = "0.05",
        PROTEOME_PROFILER_ANALYSIS_METHODS = "raw_log2_lm|normalized_t_test",
        PROTEOME_PROFILER_SHORTLIST_COMPARISONS = "stratum_1_group_a_vs_group_b|stratum_2_group_a_vs_group_b",
        PROTEOME_PROFILER_SHORTLIST_METHODS = "raw_log2_lm|normalized_t_test",
        PROTEOME_PROFILER_SHORTLIST_COORDS = paste(selected_coords, collapse = "|")
    )
}

pick_category_coordinates <- function(protocol_tbl, category_lookup, category_names, n = 1) {
    selected_tbl <- protocol_tbl %>%
        filter(`Analyte/Control` %in% names(category_lookup)) %>%
        mutate(category = unname(category_lookup[`Analyte/Control`]))

    coords <- map(category_names, function(category_name) {
        selected_tbl %>%
            filter(category == category_name) %>%
            slice_head(n = n) %>%
            pull(Coordinate)
    }) %>%
        unlist(use.names = FALSE)

    compact_coordinate(coords)
}

write_fixture_inputs <- function(fixture_root) {
    if (dir.exists(fixture_root)) {
        unlink(fixture_root, recursive = TRUE, force = TRUE)
    }

    dir.create(file.path(fixture_root, "input"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "manifests"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "workbooks", "exploratory"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fixture_root, "workbooks", "replicate"), recursive = TRUE, showWarnings = FALSE)

    protocol_tbl <- build_protocol_table()
    exploratory_manifest <- build_exploratory_manifest(fixture_root)
    replicate_manifest <- build_replicate_manifest(fixture_root)
    exploratory_signals <- simulate_exploratory_signals(exploratory_manifest, protocol_tbl)
    replicate_signals <- simulate_replicate_signals(replicate_manifest, protocol_tbl)

    protocol_path <- file.path(fixture_root, "input", "mock_protocol.xlsx")
    writexl::write_xlsx(protocol_tbl, protocol_path)

    invisible(walk2(
        exploratory_manifest$sample_id,
        seq_along(exploratory_manifest$sample_id),
        function(sample_id_value, idx) {
            write_licor_workbook(
                path = file.path(repo_root, exploratory_manifest$workbook_path[[idx]]),
                signal_vector = unname(c(
                    exploratory_signals[[sample_id_value]],
                    "Reference Spots" = runif(1, min = 460, max = 560),
                    "Negative Control" = runif(1, min = 8, max = 18)
                )),
                seed = 5000 + idx
            )
        }
    ))

    invisible(walk2(
        replicate_manifest$sample_id,
        seq_along(replicate_manifest$sample_id),
        function(sample_id_value, idx) {
            write_licor_workbook(
                path = file.path(repo_root, replicate_manifest$workbook_path[[idx]]),
                signal_vector = unname(c(
                    replicate_signals[[sample_id_value]],
                    "Reference Spots" = runif(1, min = 455, max = 565),
                    "Negative Control" = runif(1, min = 8, max = 18)
                )),
                seed = 7000 + idx
            )
        }
    ))

    exploratory_manifest_path <- file.path(fixture_root, "manifests", "exploratory_samples.csv")
    replicate_manifest_path <- file.path(fixture_root, "manifests", "replicate_samples.csv")
    write_csv(exploratory_manifest, exploratory_manifest_path)
    write_csv(replicate_manifest, replicate_manifest_path)

    exploratory_categories <- assign_exploratory_categories(sprintf("Marker %03d", seq_len(112)))
    replicate_categories <- assign_replicate_categories(sprintf("Marker %03d", seq_len(112)))

    list(
        protocol_tbl = protocol_tbl,
        exploratory_low_signal_coords = pick_category_coordinates(
            protocol_tbl,
            exploratory_categories,
            category_names = "low_signal",
            n = 8
        ),
        exploratory_selected_coords = pick_category_coordinates(
            protocol_tbl,
            exploratory_categories,
            category_names = c("group_b_up_strong", "group_b_down_strong", "group_b_up_moderate"),
            n = 1
        ),
        replicate_low_signal_coords = pick_category_coordinates(
            protocol_tbl,
            replicate_categories,
            category_names = "low_signal",
            n = 8
        ),
        replicate_selected_coords = pick_category_coordinates(
            protocol_tbl,
            replicate_categories,
            category_names = c("stratum_1_up_strong", "stratum_1_down_strong", "stratum_2_up_strong"),
            n = 1
        )
    )
}

run_repo_script <- function(script_path, env_path) {
    old_wd <- setwd(repo_root)
    on.exit(setwd(old_wd), add = TRUE)

    output <- system2(
        file.path(R.home("bin"), "Rscript"),
        script_path,
        stdout = TRUE,
        stderr = TRUE,
        env = c(sprintf("PROTEOME_PROFILER_ENV_FILE=%s", env_path))
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

    invisible(output)
}

analysis_root <- function(fixture_root, slug) {
    file.path(fixture_root, "output", "docs", slug)
}

remove_stale_docs_examples <- function(docs_output_dir) {
    stale_files <- file.path(docs_output_dir, c(
        "replicate-aware-inferential-waterfall-example.png",
        "replicate-aware-inferential-barplot-example-page-1.png",
        "threshold-diagnostics-example.png",
        "main-waterfall-example.png",
        "main-barplot-example-page-1.png",
        "selected-analytes-waterfall-example.png",
        "main-barplot-combined-page-1.png",
        "main-waterfall-sorafenib.png",
        "replicate-aware-inferential-barplot-fdr-page-1.png",
        "replicate-aware-inferential-waterfall-se.png",
        "select-analytes-waterfall.png",
        "threshold-diagnostics-region-stats.png"
    ))

    invisible(walk(stale_files, function(path) {
        if (file.exists(path)) {
            unlink(path, force = TRUE)
        }
    }))
}

copy_required_file <- function(from, to) {
    if (!file.exists(from)) {
        stop(sprintf("Expected fixture output is missing: %s", from))
    }

    copied <- file.copy(from, to, overwrite = TRUE)
    if (!isTRUE(copied)) {
        stop(sprintf("Could not copy fixture output to docs path: %s", to))
    }
}

copy_docs_examples <- function(fixture_root) {
    docs_output_dir <- file.path(repo_root, "docs", "output-examples")
    dir.create(docs_output_dir, recursive = TRUE, showWarnings = FALSE)
    remove_stale_docs_examples(docs_output_dir)

    exploratory_root <- analysis_root(fixture_root, "docs_exploratory_examples")
    replicate_root <- analysis_root(fixture_root, "docs_replicate_examples")
    replicate_run_index <- read_tsv(
        file.path(replicate_root, "inferential_results", "run_index.tsv"),
        show_col_types = FALSE
    )

    raw_log2_row <- replicate_run_index %>%
        filter(
            comparison_slug == "stratum_1_group_a_vs_group_b",
            analysis_method == "raw_log2_lm"
        ) %>%
        slice_head(n = 1)
    normalized_row <- replicate_run_index %>%
        filter(
            comparison_slug == "stratum_1_group_a_vs_group_b",
            analysis_method == "normalized_t_test"
        ) %>%
        slice_head(n = 1)

    if (nrow(raw_log2_row) != 1 || nrow(normalized_row) != 1) {
        stop("Could not identify expected replicate comparison rows in run_index.tsv.")
    }

    normalized_fdr_barplots <- str_split(normalized_row$fdr_0_25_barplot_path[[1]], pattern = ";", simplify = TRUE)
    normalized_fdr_barplot_path <- normalized_fdr_barplots[[1]]
    if (is.na(normalized_fdr_barplot_path) || identical(normalized_fdr_barplot_path, "")) {
        stop("Expected normalized_t_test FDR < 0.25 barplot path is missing from run_index.tsv.")
    }

    copy_required_file(
        file.path(exploratory_root, "threshold_diagnostics", "region_stats.png"),
        file.path(docs_output_dir, "threshold-diagnostics-example.png")
    )
    copy_required_file(
        file.path(
            exploratory_root,
            "main_analysis",
            "ref_threshold_70",
            "all_comparisons",
            "waterfalls",
            "Group A vs Group B",
            "Group A vs Group B--Group B-main_waterfall-all.png"
        ),
        file.path(docs_output_dir, "main-waterfall-example.png")
    )
    copy_required_file(
        file.path(
            exploratory_root,
            "main_analysis",
            "ref_threshold_70",
            "fold_change_hits",
            "threshold_1.5",
            "barplots",
            "Group A vs Group B",
            "combined",
            "Group A vs Group B-barplot_page_1__significant.png"
        ),
        file.path(docs_output_dir, "main-barplot-example-page-1.png")
    )
    copy_required_file(
        file.path(exploratory_root, "select_analytes", "group_a_vs_group_b", "selected_waterfall.png"),
        file.path(docs_output_dir, "selected-analytes-waterfall-example.png")
    )
    copy_required_file(
        raw_log2_row$full_waterfall_path[[1]],
        file.path(docs_output_dir, "replicate-aware-inferential-waterfall-example.png")
    )
    copy_required_file(
        normalized_fdr_barplot_path,
        file.path(docs_output_dir, "replicate-aware-inferential-barplot-example-page-1.png")
    )
}

main <- function() {
    fixture_root <- docs_fixture_root()
    fixture_inputs <- write_fixture_inputs(fixture_root)

    exploratory_env_path <- file.path(fixture_root, ".env.docs_exploratory")
    replicate_env_path <- file.path(fixture_root, ".env.docs_replicate")

    write_env_file(
        exploratory_env_path,
        build_exploratory_env_values(
            fixture_root = fixture_root,
            low_signal_coords = fixture_inputs$exploratory_low_signal_coords,
            selected_coords = fixture_inputs$exploratory_selected_coords
        )
    )
    write_env_file(
        replicate_env_path,
        build_replicate_env_values(
            fixture_root = fixture_root,
            low_signal_coords = fixture_inputs$replicate_low_signal_coords,
            selected_coords = fixture_inputs$replicate_selected_coords
        )
    )

    run_repo_script("scripts/find_ref_thresh.R", exploratory_env_path)
    run_repo_script("scripts/main.R", exploratory_env_path)
    run_repo_script("scripts/select-analytes-analysis.R", exploratory_env_path)

    run_repo_script("scripts/find_ref_thresh.R", replicate_env_path)
    run_repo_script("scripts/main.R", replicate_env_path)
    run_repo_script("scripts/select-analytes-analysis.R", replicate_env_path)

    copy_docs_examples(fixture_root)

    cat("Generated synthetic docs example fixture under:\n")
    cat(fixture_root, "\n")
}

main()
