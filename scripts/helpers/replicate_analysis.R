#' Identify whether an analysis config requests replicate-aware inference
#'
#' Input declaration and analysis mode are now separate. A manifest may be used
#' for either exploratory or inferential runs, while replicate-aware inference
#' is enabled only when `mode = "replicate"`.
#'
#' @param example_config Named list describing one analysis run.
#'
#' @return Logical scalar.
is_replicate_aware_config <- function(example_config) {
    get_analysis_mode(example_config) == "replicate"
}

#' Validate the p-value adjustment method for inferential analysis
#'
#' @param method Character scalar naming a supported `p.adjust()` method.
#'
#' @return The validated method, invisibly.
validate_p_adjust_method <- function(method) {
    supported_methods <- p.adjust.methods
    if (is.null(method) || !method %in% supported_methods) {
        stop(sprintf(
            "Unsupported p-value adjustment method '%s'. Supported methods: %s",
            method,
            paste(supported_methods, collapse = ", ")
        ))
    }
    invisible(method)
}

#' Return metadata for one supported replicate-aware inferential method
#'
#' @param method Character scalar naming the method.
#'
#' @return Named list describing the method.
get_inferential_method_spec <- function(method) {
    switch(method,
        raw_log2_lm = list(
            slug = "raw_log2_lm",
            label = "Raw signal log2 linear model",
            value_column = "signal",
            waterfall_subtitle = "Effect estimate is treatment minus control on the log2 raw-signal scale",
            statistical_methods = c(
                "Fits one per-analyte linear model `log2(signal) ~ treatment` within each subgroup.",
                "Tests the treatment coefficient against zero, so the null hypothesis is equal mean log2 raw signal between treatment and control.",
                "Reports effect size on the log2 raw-signal scale and derives fold change as `2^(treatment coefficient)`."
            ),
            so_what = c(
                "Use this when you want to ask whether treatment changed the analyte on its raw signal scale after putting the signals on a log2 scale.",
                "In practice, this is usually the better choice if you want effect sizes that behave like fold changes and a method tied directly to the raw measured intensities."
            ),
            pros = c(
                "Targets a multiplicative treatment effect on the raw-signal scale, with the coefficient equal to a log2 ratio of geometric means.",
                "The log2 transform makes up- and down-regulation symmetric and can reduce scale-dependent variance when variability increases with signal intensity.",
                "The linear-model formulation gives one treatment effect estimate on the modeled scale even when replicate counts are unbalanced across arms."
            ),
            cons = c(
                "Requires finite positive raw signals, so zero or negative background-corrected values are excluded and can make some analytes untestable.",
                "Relies on linear-model assumptions on the log2 scale, especially approximate homoscedasticity and roughly normal residual behavior, which are hard to assess with small n.",
                "Answers a different estimand than the normalized-value workflow, so p-values and effect sizes need not match workbook-style mean comparisons."
            )
        ),
        normalized_t_test = list(
            slug = "normalized_t_test",
            label = "Normalized t-test",
            value_column = "normalized_signal",
            waterfall_subtitle = paste(
                "Bars show log2(mean normalized treatment / mean normalized control);",
                "p-values come from two-sided equal-variance t-tests on raw-data normalized replicate values."
            ),
            statistical_methods = c(
                "Recomputes workbook-style normalized replicate values from the raw signals for each analyte.",
                "Runs a two-sided equal-variance two-sample t-test on those normalized replicate values within each subgroup.",
                "Tests equality of mean normalized signal between treatment and control; when both group means are positive, reports fold change as `mean(treatment normalized) / mean(control normalized)`."
            ),
            so_what = c(
                "Use this when you want each analyte signal to be divided by the sample's reference-spot intensity first, and then you want to compare treatment and control on that normalized scale.",
                "In practice, this is the better choice if you want the question to be 'do the groups differ after reference-spot normalization?' rather than 'do the raw signals differ on a log2 scale?'"
            ),
            pros = c(
                "Targets differences in mean normalized signal on the same scale used by the normalized replicate workflow.",
                "The test is simple and transparent: per analyte, it compares the two group means of the normalized replicate values directly.",
                "Can still produce a p-value when normalized values are finite even if a ratio-scale fold change is not interpretable because one mean is nonpositive."
            ),
            cons = c(
                "Inference depends on the normalization step, so any artifacts introduced by normalization propagate into the means, standard errors, and p-values.",
                "Uses an equal-variance two-sample t-test; if group variances differ, p-values can be miscalibrated, especially with small n.",
                "Targets arithmetic means of normalized values rather than a multiplicative effect on the raw-signal scale, so results can diverge from the raw log2 model."
            )
        ),
        stop(sprintf(
            "Unsupported inferential analysis method '%s'. Supported methods: %s",
            method,
            paste(c("raw_log2_lm", "normalized_t_test"), collapse = ", ")
        ))
    )
}

#' Validate and normalize configured inferential methods
#'
#' @param methods Character vector of requested method slugs.
#'
#' @return Unique validated method vector.
validate_analysis_methods <- function(methods) {
    if (is.null(methods) || length(methods) == 0) {
        return("raw_log2_lm")
    }

    methods <- unique(unname(as.character(methods)))
    walk(methods, get_inferential_method_spec)
    methods
}

#' Mean helper that returns `NA` instead of `NaN` when no finite values exist
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar.
safe_mean_or_na <- function(x) {
    finite_x <- as.numeric(x)[is.finite(as.numeric(x))]
    if (length(finite_x) == 0) {
        return(NA_real_)
    }

    mean(finite_x)
}

safe_var_or_na <- function(x) {
    finite_x <- as.numeric(x)[is.finite(as.numeric(x))]
    if (length(finite_x) < 2) {
        return(NA_real_)
    }

    stats::var(finite_x)
}

pooled_log2_ratio_se <- function(control_values, treatment_values, control_mean, treatment_mean) {
    control_values <- as.numeric(control_values)[is.finite(as.numeric(control_values))]
    treatment_values <- as.numeric(treatment_values)[is.finite(as.numeric(treatment_values))]
    control_n <- length(control_values)
    treatment_n <- length(treatment_values)

    if (
        control_n < 2 ||
        treatment_n < 2 ||
        !is.finite(control_mean) ||
        !is.finite(treatment_mean) ||
        control_mean <= 0 ||
        treatment_mean <= 0
    ) {
        return(NA_real_)
    }

    control_var <- safe_var_or_na(control_values)
    treatment_var <- safe_var_or_na(treatment_values)
    if (!is.finite(control_var) || !is.finite(treatment_var)) {
        return(NA_real_)
    }

    pooled_df <- control_n + treatment_n - 2
    if (pooled_df <= 0) {
        return(NA_real_)
    }

    pooled_var <- ((control_n - 1) * control_var + (treatment_n - 1) * treatment_var) / pooled_df
    sqrt(pooled_var * ((1 / (treatment_n * treatment_mean^2)) + (1 / (control_n * control_mean^2)))) / log(2)
}

collapse_workbook_normalized_signal <- function(x) {
    numeric_x <- suppressWarnings(as.numeric(x))
    finite_nonzero <- numeric_x[is.finite(numeric_x) & numeric_x != 0]
    if (length(finite_nonzero) > 0) {
        return(finite_nonzero[[1]])
    }

    safe_mean_or_na(numeric_x)
}

compute_raw_sheet_normalized_signal <- function(sample_pair_tbl) {
    pair_tbl <- sample_pair_tbl %>%
        mutate(
            compact_coordinate = compact_coordinate_text(Coordinate),
            compact_sname = compact_coordinate_text(str_replace_all(Sname, "([A-Z]{1,2})0([0-9])", "\\1\\2"))
        )

    preferred_reference_rows <- pair_tbl %>%
        filter(compact_sname %in% c("A1,2", "J1,2"))

    reference_rows <- if (nrow(preferred_reference_rows) > 0) {
        preferred_reference_rows
    } else {
        fallback_reference_rows <- pair_tbl %>%
            filter(Name == "Reference Spots") %>%
            filter(str_detect(compact_coordinate, "^[A-Z]{1,2}1,2$"))

        if (nrow(fallback_reference_rows) > 0) {
            fallback_reference_rows
        } else {
            pair_tbl %>%
                filter(Name == "Reference Spots")
        }
    }

    if (nrow(reference_rows) == 0) {
        stop("Workbook-style normalized t-test requires at least one 'Reference Spots' pair in the protocol table.")
    }

    normalization_denominator <- safe_mean_or_na(reference_rows$signal)
    if (!is.finite(normalization_denominator) || normalization_denominator == 0) {
        stop("Workbook-style normalized t-test could not compute a finite nonzero reference normalization denominator from raw signals.")
    }

    pair_tbl %>%
        mutate(
            normalized_signal = signal / normalization_denominator
        ) %>%
        select(-compact_coordinate, -compact_sname)
}

normalize_pair_coordinate_key <- function(x) {
    compact_coordinate_text(str_replace_all(x, "([A-Z]{1,2})0([0-9])", "\\1\\2"))
}

#' Build a stable comparison slug for inferential outputs
#'
#' @param subgroup_label Character scalar naming the subgroup level.
#' @param control_label Character scalar naming the control arm.
#' @param treatment_label Character scalar naming the treatment arm.
#'
#' @return Character scalar slug safe for output paths.
make_comparison_slug <- function(subgroup_label, control_label, treatment_label) {
    str_replace_all(
        str_to_lower(sprintf("%s__%s_vs_%s", subgroup_label, control_label, treatment_label)),
        "[^a-z0-9]+",
        "_"
    ) %>%
        str_replace_all("^_|_$", "")
}

#' Build a stable comparison slug for a pairwise exploratory comparison
#'
#' Legacy exploratory analyses do not have a subgroup label, but their
#' select-analytes outputs still need comparison-scoped folders. This helper
#' uses the same normalization rules as inferential comparison slugs while
#' omitting the subgroup component.
#'
#' @param control_label Character scalar naming the control arm.
#' @param treatment_label Character scalar naming the treatment arm.
#'
#' @return Character scalar slug safe for output paths.
make_pairwise_comparison_slug <- function(control_label, treatment_label) {
    str_replace_all(
        str_to_lower(sprintf("%s_vs_%s", control_label, treatment_label)),
        "[^a-z0-9]+",
        "_"
    ) %>%
        str_replace_all("^_|_$", "")
}

#' Return configured explicit analyte names for the select-analytes workflow
#'
#' The select-analytes workflow is intentionally explicit: users name the
#' analytes they want to follow up, and the script writes focused plots for
#' those names rather than deriving a shortlist from significance or
#' fold-change thresholds.
#'
#' @param example_config Named analysis config returned by `get_analysis_config()`.
#'
#' @return Character vector of configured analyte names.
get_selected_analyte_names <- function(example_config) {
    selected_names <- example_config$selection_analytes
    if (is.null(selected_names) || length(selected_names) == 0 || all(is.na(selected_names)) || all(selected_names == "")) {
        stop("Select-analytes requires `shortlist$analytes` / `selection_analytes` to list explicit analyte names.")
    }

    unique(as.character(selected_names))
}

#' Suggest close protocol/result names for missing selected analytes
#'
#' @param missing_name Character scalar requested by the user.
#' @param available_names Character vector of available analyte names.
#' @param n Integer number of suggestions to return.
#'
#' @return Character vector of nearby available names.
suggest_selected_analyte_names <- function(missing_name, available_names, n = 5) {
    available_names <- unique(as.character(available_names))
    if (length(available_names) == 0) {
        return(character())
    }

    distances <- utils::adist(tolower(missing_name), tolower(available_names))
    available_names[order(as.numeric(distances), available_names)][seq_len(min(n, length(available_names)))]
}

#' Validate selected analyte names against available analysis names
#'
#' @param selected_names Character vector requested in the config.
#' @param available_names Character vector present in the data or results.
#'
#' @return Invisibly returns `selected_names`; errors if any are missing.
validate_selected_analyte_names <- function(selected_names, available_names) {
    missing_names <- setdiff(selected_names, unique(as.character(available_names)))
    if (length(missing_names) == 0) {
        return(invisible(selected_names))
    }

    suggestion_text <- map_chr(missing_names, function(missing_name) {
        suggestions <- suggest_selected_analyte_names(missing_name, available_names)
        sprintf(
            "%s (closest available: %s)",
            missing_name,
            paste(suggestions, collapse = ", ")
        )
    })

    stop(sprintf(
        "Configured selected analyte name(s) were not found: %s",
        paste(suggestion_text, collapse = "; ")
    ))
}

#' Resolve legacy select-analytes comparisons into control/treatment pairs
#'
#' Legacy configs can define many exploratory comparisons, but select-analytes
#' writes one focused output folder per selected pair. If explicit comparison
#' slugs are configured, they are matched against all configured pairwise
#' comparison slugs; otherwise the legacy `selection_control` and
#' `selection_group` fields identify the selected pair(s).
#'
#' @param example_config Named analysis config returned by `get_analysis_config()`.
#'
#' @return Tibble with `comparison_slug`, `control`, and `treatment` columns.
resolve_legacy_select_comparisons <- function(example_config) {
    configured_comparisons <- example_config$comparisons
    if (is.null(configured_comparisons) || length(configured_comparisons) == 0) {
        stop("Legacy select-analytes requires at least one configured comparison.")
    }

    all_pairs <- imap_dfr(configured_comparisons, function(treatments, control_label) {
        tibble(
            control = as.character(control_label),
            treatment = as.character(treatments)
        )
    }) %>%
        mutate(comparison_slug = make_pairwise_comparison_slug(control, treatment)) %>%
        select(comparison_slug, control, treatment)

    configured_slugs <- example_config$selection_comparison_slugs
    if (!is.null(configured_slugs) && length(configured_slugs) > 0) {
        selected_pairs <- all_pairs %>%
            filter(comparison_slug %in% configured_slugs)
        missing_slugs <- setdiff(configured_slugs, selected_pairs$comparison_slug)
        if (length(missing_slugs) > 0) {
            stop(sprintf(
                "Configured legacy select comparison slug(s) were not found: %s. Available slugs: %s",
                paste(missing_slugs, collapse = ", "),
                paste(all_pairs$comparison_slug, collapse = ", ")
            ))
        }
        return(selected_pairs)
    }

    if (!is.null(example_config$selection_control) && !is.null(example_config$selection_group)) {
        selected_pairs <- all_pairs %>%
            filter(
                control == example_config$selection_control,
                treatment %in% as.character(example_config$selection_group)
            )
        if (nrow(selected_pairs) == 0) {
            stop(sprintf(
                "Legacy select comparison '%s' vs '%s' was not found. Available slugs: %s",
                example_config$selection_control,
                paste(example_config$selection_group, collapse = ", "),
                paste(all_pairs$comparison_slug, collapse = ", ")
            ))
        }
        return(selected_pairs)
    }

    all_pairs
}

#' Read and validate a replicate-aware sample manifest
#'
#' @param manifest_path Absolute or project-relative path to the manifest CSV/TSV.
#' @param subgroup_var Optional subgroup column required for within-stratum runs.
#' @param treatment_var Name of the treatment column.
#'
#' @return Tibble containing validated sample metadata.
read_sample_manifest <- function(manifest_path, subgroup_var = NULL, treatment_var = "treatment") {
    manifest_ext <- tools::file_ext(manifest_path)
    manifest <- if (manifest_ext %in% c("tsv", "txt")) {
        read_tsv(manifest_path, show_col_types = FALSE)
    } else {
        read_csv(manifest_path, show_col_types = FALSE)
    }

    required_columns <- c("sample_id", "workbook_path", treatment_var)
    if (!is.null(subgroup_var)) {
        required_columns <- c(required_columns, subgroup_var)
    }

    missing_columns <- setdiff(required_columns, names(manifest))
    if (length(missing_columns) > 0) {
        stop(sprintf(
            "Sample manifest is missing required columns: %s",
            paste(missing_columns, collapse = ", ")
        ))
    }

    blank_value_columns <- required_columns[vapply(required_columns, function(column_name) {
        any(is.na(manifest[[column_name]]) | str_trim(as.character(manifest[[column_name]])) == "")
    }, logical(1))]
    if (length(blank_value_columns) > 0) {
        stop(sprintf(
            "Sample manifest contains blank or NA values in required columns: %s",
            paste(blank_value_columns, collapse = ", ")
        ))
    }

    duplicate_ids <- manifest %>%
        count(sample_id, name = "n_rows") %>%
        filter(n_rows > 1)
    if (nrow(duplicate_ids) > 0) {
        stop(sprintf(
            "Sample manifest contains duplicate sample_id values: %s",
            paste(duplicate_ids$sample_id, collapse = ", ")
        ))
    }

    if ("sheet_name" %in% names(manifest)) {
        manifest <- manifest %>%
            mutate(sheet_name = na_if(str_trim(as.character(sheet_name)), ""))
    }

    manifest
}

#' Resolve workbook paths declared in a sample manifest
#'
#' @param manifest Tibble returned by `read_sample_manifest()`.
#' @param data_dir Optional absolute data-directory root to resolve relative
#'   workbook paths against before falling back to project path resolution.
#'
#' @return Tibble with an added `resolved_workbook_path` column.
resolve_manifest_workbook_paths <- function(manifest, data_dir = NULL) {
    resolved_paths <- map_chr(manifest$workbook_path, function(path_value) {
        if (grepl("^~|^/", path_value)) {
            expanded <- path.expand(path_value)
            if (!file.exists(expanded)) {
                stop(sprintf("Workbook path does not exist: %s", path_value))
            }
            return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
        }

        data_dir_candidate <- if (!is.null(data_dir)) file.path(data_dir, path_value) else NULL
        if (!is.null(data_dir_candidate) && file.exists(data_dir_candidate)) {
            return(normalizePath(data_dir_candidate, winslash = "/", mustWork = TRUE))
        }

        resolve_project_path(path_value, must_exist = TRUE)
    })

    manifest_with_paths <- manifest %>%
        mutate(
            workbook_path = as.character(workbook_path),
            resolved_workbook_path = resolved_paths
        )

    workbook_sheet_map <- manifest_with_paths %>%
        distinct(resolved_workbook_path) %>%
        mutate(available_sheets = map(resolved_workbook_path, readxl::excel_sheets))

    if (!"sheet_name" %in% names(manifest_with_paths)) {
        manifest_with_paths$sheet_name <- NA_character_
    }

    manifest_with_paths %>%
        left_join(workbook_sheet_map, by = "resolved_workbook_path") %>%
        mutate(
            resolved_sheet_name = pmap_chr(
                list(workbook_path, resolved_workbook_path, sheet_name, available_sheets),
                function(workbook_path, resolved_workbook_path, sheet_name, available_sheets) {
                    if (is.na(sheet_name) || identical(sheet_name, "")) {
                        if (length(available_sheets) == 1) {
                            return(available_sheets[[1]])
                        }

                        stop(sprintf(
                            paste(
                                "Workbook '%s' contains multiple sheets.",
                                "Add a non-blank sheet_name column to the manifest for each sample.",
                                "Available sheets: %s"
                            ),
                            workbook_path,
                            paste(available_sheets, collapse = ", ")
                        ))
                    }

                    if (!sheet_name %in% available_sheets) {
                        stop(sprintf(
                            "Workbook '%s' does not contain manifest sheet_name '%s'. Available sheets: %s",
                            workbook_path,
                            sheet_name,
                            paste(available_sheets, collapse = ", ")
                        ))
                    }

                    sheet_name
                }
            )
        ) %>%
        dplyr::select(-available_sheets)
}

#' Read one workbook sheet and normalize it to spot names plus raw signal
#'
#' This accepts both the original LI-COR export layout where the header begins
#' on row 4 with columns `Name` and `Signal`, and collaborator workbooks where
#' each sample sheet starts on row 1 with columns such as `Spot Name` and
#' `Signal`.
#'
#' @param workbook_path Absolute path to an Excel workbook.
#' @param sheet_name Character scalar naming the worksheet to read.
#' @param sample_id Optional sample identifier. When a collaborator workbook
#'   contains a matching normalized column, it is carried through as
#'   `NormalizedSignal`.
#' @param normalized_column_candidates Optional character vector of fallback
#'   column names to search for collaborator-style normalized values.
#'
#' @return Tibble with `Sname`, `Signal`, optional `NormalizedSignal`, and
#'   `sname_grouping`.
read_licor_signal_sheet <- function(workbook_path, sheet_name, sample_id = NULL, normalized_column_candidates = NULL) {
    format_workbook_ref <- function(path, sheet) {
        sprintf("%s [%s]", path, sheet)
    }

    build_sequential_layout <- function(signal_tbl) {
        signal_tbl_no_background <- signal_tbl %>%
            filter(!startsWith(Sname, "B"))

        expected_snames <- sprintf("S%03d", seq_len(nrow(signal_tbl_no_background)))
        if (
            nrow(signal_tbl_no_background) > 0 &&
            nrow(signal_tbl_no_background) %% 2 == 0 &&
            identical(signal_tbl_no_background$Sname, expected_snames)
        ) {
            return(signal_tbl_no_background %>%
                mutate(
                    sname_grouping = rep(seq_len(nrow(.) / 2), each = 2),
                    layout_type = "sequential"
                ))
        }

        NULL
    }

    build_coordinate_layout <- function(signal_tbl) {
        if (nrow(signal_tbl) == 0 || nrow(signal_tbl) %% 2 != 0) {
            return(NULL)
        }

        row_labels <- str_extract(signal_tbl$Sname, "^[A-Z]{1,2}")
        col_labels <- suppressWarnings(as.integer(str_extract(signal_tbl$Sname, "[0-9]{2}$")))

        if (
            any(is.na(row_labels)) ||
            any(is.na(col_labels)) ||
            !all(col_labels[c(TRUE, FALSE)] %% 2 == 1) ||
            !all(col_labels[c(FALSE, TRUE)] == col_labels[c(TRUE, FALSE)] + 1) ||
            !all(row_labels[c(TRUE, FALSE)] == row_labels[c(FALSE, TRUE)])
        ) {
            return(NULL)
        }

        signal_tbl %>%
            mutate(
                sname_grouping = rep(seq_len(nrow(.) / 2), each = 2),
                layout_type = "coordinate"
            )
    }

    extract_normalized_signal <- function(raw_tbl) {
        candidate_names <- c(sample_id, normalized_column_candidates)
        candidate_names <- candidate_names[!is.na(candidate_names) & candidate_names != ""]
        candidate_names <- unique(candidate_names)

        if (length(candidate_names) == 0) {
            return(rep(NA_real_, nrow(raw_tbl)))
        }

        for (candidate_name in candidate_names) {
            matching_cols <- which(names(raw_tbl) == candidate_name)
            if (length(matching_cols) == 1) {
                return(suppressWarnings(as.numeric(raw_tbl[[matching_cols[[1]]]])))
            }
        }

        rep(NA_real_, nrow(raw_tbl))
    }

    normalize_signal_table <- function(raw_tbl) {
        if ("Name" %in% names(raw_tbl) && "Signal" %in% names(raw_tbl)) {
            signal_tbl <- tibble(
                Sname = as.character(raw_tbl[["Name"]]),
                Signal = suppressWarnings(as.numeric(raw_tbl[["Signal"]])),
                NormalizedSignal = rep(NA_real_, nrow(raw_tbl))
            )
        } else if ("Spot Name" %in% names(raw_tbl) && "Signal" %in% names(raw_tbl)) {
            signal_tbl <- tibble(
                Sname = as.character(raw_tbl[["Spot Name"]]),
                Signal = suppressWarnings(as.numeric(raw_tbl[["Signal"]])),
                NormalizedSignal = extract_normalized_signal(raw_tbl)
            )
        } else {
            return(NULL)
        }

        signal_tbl <- signal_tbl %>%
            filter(!is.na(Sname), str_trim(Sname) != "") %>%
            mutate(Sname = str_trim(Sname))

        sequential_layout <- build_sequential_layout(signal_tbl)
        if (!is.null(sequential_layout)) {
            return(sequential_layout)
        }

        coordinate_layout <- build_coordinate_layout(signal_tbl)
        if (!is.null(coordinate_layout)) {
            return(coordinate_layout)
        }

        NULL
    }

    for (skip_rows in c(3, 0)) {
        raw_tbl <- suppressMessages(suppressWarnings(
            read_excel(
                workbook_path,
                sheet = sheet_name,
                skip = skip_rows,
                .name_repair = "minimal"
            )
        ))
        normalized_tbl <- normalize_signal_table(raw_tbl)
        if (!is.null(normalized_tbl)) {
            return(normalized_tbl)
        }
    }

    stop(sprintf(
        paste(
            "Could not parse LI-COR signal data from %s.",
            "Expected either row-4 columns `Name`/`Signal` or row-1 columns `Spot Name`/`Signal`.",
            "Multi-sheet collaborator workbooks must point each sample to the correct sheet_name."
        ),
        format_workbook_ref(workbook_path, sheet_name)
    ))
}

#' Convert protocol coordinates into the compact display format used downstream
#'
#' @param text_vector Character vector like `"A1, A2"`.
#'
#' @return Character vector like `"A1,2"`.
compact_coordinate_text <- function(text_vector) {
    text_list <- str_split(text_vector, ", ")

    output <- sapply(text_list, function(words) {
        letters <- gsub("\\d+", "", words)
        numbers <- gsub("\\D+", "", words)

        unique_letters <- unique(letters)
        pieces <- character(length(unique_letters))
        for (i in seq_along(unique_letters)) {
            matching_numbers <- numbers[letters == unique_letters[i]]
            pieces[i] <- paste(unique_letters[i], paste(matching_numbers, collapse = ","), sep = "")
        }

        paste(pieces, collapse = ", ")
    })

    unname(output)
}

#' Read one LI-COR workbook and average duplicate membrane spots within sample
#'
#' @param workbook_path Absolute path to one LI-COR workbook.
#' @param analyte_info Protocol-derived analyte lookup table.
#' @param sample_metadata One-row tibble containing sample metadata.
#' @param remove_dat Optional tibble of technical artifacts to drop. Uses columns
#'   `Sname`, `sample_id`, and `group` when present.
#' @param treatment_var Name of the treatment column in `sample_metadata`.
#'
#' @return Sample-level analyte tibble with one row per analyte.
read_licor_sample_workbook <- function(workbook_path, analyte_info, sample_metadata, remove_dat = NULL, treatment_var = "treatment", sheet_name = NULL) {
    sample_id_value <- sample_metadata$sample_id[[1]]
    treatment_value <- sample_metadata[[treatment_var]][[1]]
    sheet_name_value <- if ("resolved_sheet_name" %in% names(sample_metadata)) sample_metadata$resolved_sheet_name[[1]] else sheet_name
    if (is.null(sheet_name_value) || is.na(sheet_name_value) || identical(sheet_name_value, "")) {
        stop(sprintf("Missing resolved sheet_name for sample_id '%s'.", sample_id_value))
    }

    res_temp_pre <- read_licor_signal_sheet(
        workbook_path = workbook_path,
        sheet_name = sheet_name_value,
        sample_id = sample_id_value,
        normalized_column_candidates = c(sheet_name_value)
    )

    if (is.null(remove_dat)) {
        remove_dat <- tibble(Sname = character(), sample_id = character(), group = character())
    }

    sample_remove_dat <- remove_dat %>%
        filter((is.na(sample_id) | .data$sample_id == .env$sample_id_value)) %>%
        filter((is.na(group) | .data$group == .env$treatment_value)) %>%
        distinct(Sname)

    res_temp <- res_temp_pre %>%
        anti_join(sample_remove_dat, by = "Sname")

    sample_metadata_clean <- sample_metadata %>%
        dplyr::select(-any_of(c("resolved_workbook_path", "resolved_sheet_name")))
    sheet_layout_type <- unique(res_temp$layout_type)
    if (length(sheet_layout_type) != 1) {
        stop(sprintf(
            "Expected exactly one sheet layout type for sample_id '%s'; found: %s",
            sample_id_value,
            paste(sheet_layout_type, collapse = ", ")
        ))
    }

    analyte_info_with_keys <- analyte_info %>%
        mutate(coordinate_key = compact_coordinate_text(Coordinate))

    result_tbl <- res_temp %>%
        group_by(sname_grouping) %>%
        reframe(
            signal = safe_mean_or_na(Signal),
            cached_normalized_signal = collapse_workbook_normalized_signal(NormalizedSignal),
            Sname = str_c(Sname, collapse = ", ")
        ) %>%
        mutate(raw_coordinate_key = normalize_pair_coordinate_key(Sname))

    metadata_columns <- setdiff(names(analyte_info), "sname_grouping")
    coordinate_join_tbl <- analyte_info_with_keys %>%
        dplyr::select(-sname_grouping) %>%
        rename_with(~ paste0(.x, "_coord"), -coordinate_key)
    grouping_join_tbl <- analyte_info %>%
        rename_with(~ paste0(.x, "_group"), all_of(metadata_columns))

    result_tbl <- result_tbl %>%
        left_join(coordinate_join_tbl, by = c("raw_coordinate_key" = "coordinate_key"))

    if (identical(sheet_layout_type, "sequential")) {
        result_tbl <- result_tbl %>%
            left_join(grouping_join_tbl, by = "sname_grouping")

        for (column_name in metadata_columns) {
            result_tbl[[column_name]] <- dplyr::coalesce(
                result_tbl[[paste0(column_name, "_coord")]],
                result_tbl[[paste0(column_name, "_group")]]
            )
        }
    } else {
        for (column_name in metadata_columns) {
            result_tbl[[column_name]] <- result_tbl[[paste0(column_name, "_coord")]]
        }
    }

    result_tbl <- result_tbl %>%
        dplyr::select(
            sname_grouping,
            signal,
            cached_normalized_signal,
            Sname,
            raw_coordinate_key,
            all_of(metadata_columns)
        )

    result_tbl <- result_tbl %>%
        dplyr::select(-raw_coordinate_key) %>%
        compute_raw_sheet_normalized_signal()

    bind_cols(
        result_tbl,
        sample_metadata_clean[rep(1, nrow(result_tbl)), , drop = FALSE]
    ) %>%
        mutate(
            sample_id = .env$sample_id_value,
            treatment = .env$treatment_value,
            workbook_path = .env$workbook_path,
            sheet_name = .env$sheet_name_value
        )
}

#' Read optional technical-artifact exclusions for legacy or sample-level inputs
#'
#' @param remove_fn Path to an optional `bad_analytes.xlsx` workbook.
#'
#' @return Tibble with any available exclusion columns.
read_bad_analytes_file <- function(remove_fn) {
    if (!file.exists(remove_fn)) {
        return(tibble(Sname = character(), group = character(), sample_id = character()))
    }

    remove_dat <- read_excel(remove_fn)
    names(remove_dat) <- c("Sname", names(remove_dat)[-1])
    remove_dat <- remove_dat %>%
        mutate(Sname = as.character(Sname))

    if (!"group" %in% names(remove_dat)) {
        remove_dat$group <- NA_character_
    }
    if (!"sample_id" %in% names(remove_dat)) {
        remove_dat$sample_id <- NA_character_
    }

    remove_dat %>%
        mutate(
            group = as.character(group),
            sample_id = as.character(sample_id)
        )
}

#' Build the canonical sample-level dataset for replicate-aware analyses
#'
#' @param manifest Validated sample manifest with resolved workbook paths.
#' @param analyte_info Protocol-derived analyte lookup table.
#' @param treatment_var Name of the treatment column.
#' @param subgroup_var Optional subgroup column to carry into the output.
#' @param data_dir Optional directory containing raw LI-COR workbooks.
#'
#' @return Tibble with one row per `analyte x sample_id`.
build_sample_level_dataset <- function(manifest, analyte_info, treatment_var = "treatment", subgroup_var = NULL, data_dir = NULL) {
    remove_fn <- if (!is.null(data_dir)) file.path(data_dir, "..", "bad_analytes.xlsx") else ""
    remove_dat <- read_bad_analytes_file(remove_fn)

    manifest %>%
        mutate(sample_tbl = pmap(
            list(resolved_workbook_path, resolved_sheet_name, sample_id, workbook_path),
            function(resolved_workbook_path, resolved_sheet_name, sample_id_value, workbook_path) {
                sample_metadata <- manifest %>%
                    filter(.data$sample_id == .env$sample_id_value) %>%
                    slice_head(n = 1)
                read_licor_sample_workbook(
                    workbook_path = resolved_workbook_path,
                    analyte_info = analyte_info,
                    sample_metadata = sample_metadata,
                    remove_dat = remove_dat,
                    treatment_var = treatment_var,
                    sheet_name = resolved_sheet_name
                )
            }
        )) %>%
        dplyr::select(sample_tbl) %>%
        unnest(sample_tbl) %>%
        mutate(Coordinate = compact_coordinate_text(Coordinate)) %>%
        {
            if (!is.null(subgroup_var)) {
                mutate(., subgroup = .data[[subgroup_var]])
            } else {
                .
            }
        } %>%
        mutate(
            treatment = as.character(.data[[treatment_var]])
        )
}

#' Summarize sample-level input cleanliness before downstream analysis
#'
#' This preflight report is descriptive only. It surfaces rows that will not be
#' usable on the log2 scale and highlights whether those issues cluster within
#' particular samples or analytes.
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#'
#' @return Named list of tibbles for run-level, sample-level, analyte-level, and
#'   row-level issue summaries.
summarize_sample_level_cleanliness <- function(sample_data) {
    if (!"workbook_path" %in% names(sample_data)) {
        sample_data$workbook_path <- NA_character_
    }
    if (!"sheet_name" %in% names(sample_data)) {
        sample_data$sheet_name <- NA_character_
    }
    if (!"sex" %in% names(sample_data)) {
        sample_data$sex <- NA_character_
    }

    issue_rows <- sample_data %>%
        mutate(
            signal_issue = case_when(
                is.na(signal) ~ "missing",
                !is.finite(signal) ~ "non_finite",
                signal < 0 ~ "negative",
                signal == 0 ~ "zero",
                TRUE ~ "ok"
            )
        )

    summary_tbl <- tibble(
        n_rows = nrow(issue_rows),
        n_samples = n_distinct(issue_rows$sample_id),
        n_analytes = n_distinct(issue_rows$Name),
        n_missing_signal = sum(issue_rows$signal_issue == "missing"),
        n_non_finite_signal = sum(issue_rows$signal_issue == "non_finite"),
        n_zero_signal = sum(issue_rows$signal_issue == "zero"),
        n_negative_signal = sum(issue_rows$signal_issue == "negative"),
        n_nonpositive_signal = sum(issue_rows$signal_issue %in% c("zero", "negative"))
    )

    sample_summary_tbl <- issue_rows %>%
        group_by(sample_id, sex, treatment, workbook_path, sheet_name) %>%
        summarize(
            n_rows = n(),
            n_missing_signal = sum(signal_issue == "missing"),
            n_non_finite_signal = sum(signal_issue == "non_finite"),
            n_zero_signal = sum(signal_issue == "zero"),
            n_negative_signal = sum(signal_issue == "negative"),
            n_nonpositive_signal = sum(signal_issue %in% c("zero", "negative")),
            min_signal = min(signal, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(min_signal = if_else(is.infinite(min_signal), NA_real_, min_signal)) %>%
        arrange(desc(n_nonpositive_signal), sample_id)

    analyte_summary_tbl <- issue_rows %>%
        group_by(Name, Coordinate) %>%
        summarize(
            n_rows = n(),
            n_missing_signal = sum(signal_issue == "missing"),
            n_non_finite_signal = sum(signal_issue == "non_finite"),
            n_zero_signal = sum(signal_issue == "zero"),
            n_negative_signal = sum(signal_issue == "negative"),
            n_nonpositive_signal = sum(signal_issue %in% c("zero", "negative")),
            min_signal = min(signal, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(min_signal = if_else(is.infinite(min_signal), NA_real_, min_signal)) %>%
        arrange(desc(n_nonpositive_signal), min_signal, Name, Coordinate)

    issue_rows_tbl <- issue_rows %>%
        filter(signal_issue != "ok") %>%
        select(sample_id, sex, treatment, workbook_path, sheet_name, Name, Coordinate, Sname, signal, signal_issue) %>%
        arrange(signal_issue, signal, sample_id, Name)

    list(
        summary = summary_tbl,
        sample_summary = sample_summary_tbl,
        analyte_summary = analyte_summary_tbl,
        issue_rows = issue_rows_tbl
    )
}

#' Write sample-level cleanliness summaries to disk
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#' @param output_dir Directory where cleanliness outputs should be written.
#'
#' @return Invisibly returns paths to the written summary files.
write_sample_level_cleanliness_outputs <- function(sample_data, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    cleanliness <- summarize_sample_level_cleanliness(sample_data)

    summary_path <- file.path(output_dir, "summary.tsv")
    sample_summary_path <- file.path(output_dir, "sample_summary.tsv")
    analyte_summary_path <- file.path(output_dir, "analyte_summary.tsv")
    issue_rows_path <- file.path(output_dir, "issue_rows.tsv")

    write_tsv(cleanliness$summary, summary_path)
    write_tsv(cleanliness$sample_summary, sample_summary_path)
    write_tsv(cleanliness$analyte_summary, analyte_summary_path)
    write_tsv(cleanliness$issue_rows, issue_rows_path)

    invisible(list(
        summary = summary_path,
        sample_summary = sample_summary_path,
        analyte_summary = analyte_summary_path,
        issue_rows = issue_rows_path
    ))
}

#' Collapse sample-level data into a wide analyte table for threshold diagnostics
#'
#' `find_ref_thresh.R` needs one column per displayed experimental condition so
#' users can inspect low-signal analytes across the configured replicate
#' conditions. In replicate-aware mode we derive those condition columns from the
#' sample-level dataset instead of assuming one workbook per group.
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#' @param subgroup_var Optional subgroup column such as `sex`.
#'
#' @return Named list with `wide_df` and `group_levels`.
build_threshold_diagnostic_dataset <- function(sample_data, subgroup_var = NULL) {
    if (!is.null(subgroup_var) && subgroup_var %in% names(sample_data)) {
        sample_data <- sample_data %>%
            mutate(
                threshold_group = str_c(.data[[subgroup_var]], treatment, sep = " | ")
            )
    } else {
        sample_data <- sample_data %>%
            mutate(threshold_group = treatment)
    }

    group_levels <- sample_data %>%
        distinct(threshold_group) %>%
        pull(threshold_group) %>%
        as.character()

    wide_df <- sample_data %>%
        group_by(Name, Sname, Coordinate, threshold_group) %>%
        summarize(signal = mean(signal, na.rm = TRUE), .groups = "drop") %>%
        mutate(threshold_group = factor(threshold_group, levels = group_levels)) %>%
        pivot_wider(names_from = threshold_group, values_from = signal) %>%
        arrange(Name, Coordinate)

    list(
        wide_df = wide_df,
        group_levels = group_levels
    )
}

#' Summarize low-signal status for inferential comparisons
#'
#' @param comparison_df Sample-level analyte data for one subgroup/comparison.
#' @param threshold Numeric raw-signal threshold.
#'
#' @return Tibble with one row per analyte and a `low_signal_flag`.
summarize_low_signal_flags <- function(comparison_df, threshold) {
    comparison_df %>%
        group_by(Name, Coordinate) %>%
        summarize(
            # Flag the analyte if any sample in the current comparison falls at
            # or below the raw-signal cutoff. Inference still runs; this is a
            # QC marker rather than a hard exclusion rule.
            low_signal_flag = any(signal <= threshold, na.rm = TRUE),
            .groups = "drop"
        )
}

#' Fit one per-analyte within-stratum treatment-versus-control model
#'
#' @param comparison_df Sample-level analyte data for a single subgroup/contrast.
#' @param control_label Character scalar naming the control arm.
#' @param treatment_label Character scalar naming the treatment arm.
#' @param subgroup_label Character scalar naming the subgroup level.
#' @param subgroup_var Character scalar naming the subgroup column.
#' @param alpha Numeric alpha threshold.
#' @param min_reps Integer minimum replicates per arm to run inference.
#' @param low_signal_threshold Numeric threshold used for low-signal flagging.
#'
#' @return Tibble of analyte-level inferential results with raw p-values.
fit_one_comparison_raw_log2_lm <- function(comparison_df, control_label, treatment_label, subgroup_label, subgroup_var, alpha, min_reps = 2, low_signal_threshold = NULL) {
    low_signal_flags <- if (!is.null(low_signal_threshold) && !is.na(low_signal_threshold)) {
        summarize_low_signal_flags(comparison_df, low_signal_threshold)
    } else {
        comparison_df %>%
            distinct(Name, Coordinate) %>%
            mutate(low_signal_flag = FALSE)
    }

    comparison_df %>%
        group_by(Name, Coordinate) %>%
        group_modify(function(dat, keys) {
            dat <- dat %>%
                mutate(
                    treatment = factor(treatment, levels = c(control_label, treatment_label)),
                    log2_signal = NA_real_
                )

            valid_signal_rows <- dat$signal > 0 & is.finite(dat$signal)
            dat$log2_signal[valid_signal_rows] <- log2(dat$signal[valid_signal_rows])

            control_rows_total <- sum(dat$treatment == control_label)
            treatment_rows_total <- sum(dat$treatment == treatment_label)
            control_n <- sum(dat$treatment == control_label & is.finite(dat$log2_signal))
            treatment_n <- sum(dat$treatment == treatment_label & is.finite(dat$log2_signal))
            control_mean_signal <- mean(dat$signal[dat$treatment == control_label], na.rm = TRUE)
            treatment_mean_signal <- mean(dat$signal[dat$treatment == treatment_label], na.rm = TRUE)

            if (control_n < min_reps || treatment_n < min_reps) {
                return(tibble(
                    control = control_label,
                    treatment = treatment_label,
                    subgroup = subgroup_label,
                    subgroup_var = subgroup_var,
                    control_rows_total = control_rows_total,
                    treatment_rows_total = treatment_rows_total,
                    control_n = control_n,
                    treatment_n = treatment_n,
                    control_mean_signal = control_mean_signal,
                    treatment_mean_signal = treatment_mean_signal,
                    control_mean_normalized = safe_mean_or_na(dat$normalized_signal[dat$treatment == control_label]),
                    treatment_mean_normalized = safe_mean_or_na(dat$normalized_signal[dat$treatment == treatment_label]),
                    effect_estimate_log2 = NA_real_,
                    effect_se_log2 = NA_real_,
                    fold_change_ratio = NA_real_,
                    fold_change_status = "not_available",
                    fold_change_note = "Fold change unavailable because the analyte was not testable in the raw log2 model.",
                    raw_p_value = NA_real_,
                    test_status = "not_testable",
                    test_reason = sprintf("requires >= %s finite positive signals per arm", min_reps),
                    low_replication_warning = control_n < 3 || treatment_n < 3,
                    alpha = alpha
                ))
            }

            model_dat <- dat %>%
                filter(is.finite(log2_signal))

            model_fit <- lm(log2_signal ~ treatment, data = model_dat)
            fit_summary <- summary(model_fit)$coefficients
            coef_name <- sprintf("treatment%s", treatment_label)
            raw_p_value <- if (coef_name %in% rownames(fit_summary)) fit_summary[coef_name, "Pr(>|t|)"] else NA_real_
            effect_estimate_log2 <- if (coef_name %in% rownames(fit_summary)) unname(coef(model_fit)[coef_name]) else NA_real_
            effect_se_log2 <- if (coef_name %in% rownames(fit_summary)) fit_summary[coef_name, "Std. Error"] else NA_real_

            tibble(
                control = control_label,
                treatment = treatment_label,
                subgroup = subgroup_label,
                subgroup_var = subgroup_var,
                control_rows_total = control_rows_total,
                treatment_rows_total = treatment_rows_total,
                control_n = control_n,
                treatment_n = treatment_n,
                control_mean_signal = control_mean_signal,
                treatment_mean_signal = treatment_mean_signal,
                control_mean_normalized = safe_mean_or_na(dat$normalized_signal[dat$treatment == control_label]),
                treatment_mean_normalized = safe_mean_or_na(dat$normalized_signal[dat$treatment == treatment_label]),
                effect_estimate_log2 = effect_estimate_log2,
                effect_se_log2 = effect_se_log2,
                fold_change_ratio = if_else(is.finite(effect_estimate_log2), 2^effect_estimate_log2, NA_real_),
                fold_change_status = if_else(is.finite(effect_estimate_log2), "reported", "not_available"),
                fold_change_note = if_else(
                    is.finite(effect_estimate_log2),
                    NA_character_,
                    "Fold change unavailable because the raw log2 model did not yield a finite effect estimate."
                ),
                raw_p_value = raw_p_value,
                test_status = "tested",
                test_reason = NA_character_,
                low_replication_warning = control_n < 3 || treatment_n < 3,
                alpha = alpha
            )
        }) %>%
        ungroup() %>%
        left_join(low_signal_flags, by = c("Name", "Coordinate"))
}

#' Fit one per-analyte workbook-style normalized t-test within subgroup
#'
#' @param comparison_df Sample-level analyte data for a single subgroup/contrast.
#' @param control_label Character scalar naming the control arm.
#' @param treatment_label Character scalar naming the treatment arm.
#' @param subgroup_label Character scalar naming the subgroup level.
#' @param subgroup_var Character scalar naming the subgroup column.
#' @param alpha Numeric alpha threshold.
#' @param min_reps Integer minimum replicates per arm to run inference.
#' @param low_signal_threshold Numeric threshold used for low-signal flagging.
#'
#' @return Tibble of analyte-level inferential results with raw p-values.
fit_one_comparison_normalized_t_test <- function(comparison_df, control_label, treatment_label, subgroup_label, subgroup_var, alpha, min_reps = 2, low_signal_threshold = NULL) {
    low_signal_flags <- if (!is.null(low_signal_threshold) && !is.na(low_signal_threshold)) {
        summarize_low_signal_flags(comparison_df, low_signal_threshold)
    } else {
        comparison_df %>%
            distinct(Name, Coordinate) %>%
            mutate(low_signal_flag = FALSE)
    }

    comparison_df %>%
        group_by(Name, Coordinate) %>%
        group_modify(function(dat, keys) {
            dat <- dat %>%
                mutate(
                    treatment = factor(treatment, levels = c(control_label, treatment_label))
                )

            valid_rows <- is.finite(dat$normalized_signal)
            control_rows_total <- sum(dat$treatment == control_label)
            treatment_rows_total <- sum(dat$treatment == treatment_label)
            control_n <- sum(dat$treatment == control_label & valid_rows)
            treatment_n <- sum(dat$treatment == treatment_label & valid_rows)
            control_mean_signal <- safe_mean_or_na(dat$signal[dat$treatment == control_label])
            treatment_mean_signal <- safe_mean_or_na(dat$signal[dat$treatment == treatment_label])
            control_mean_normalized <- safe_mean_or_na(dat$normalized_signal[dat$treatment == control_label & valid_rows])
            treatment_mean_normalized <- safe_mean_or_na(dat$normalized_signal[dat$treatment == treatment_label & valid_rows])

            if (control_n < min_reps || treatment_n < min_reps) {
                return(tibble(
                    control = control_label,
                    treatment = treatment_label,
                    subgroup = subgroup_label,
                    subgroup_var = subgroup_var,
                    control_rows_total = control_rows_total,
                    treatment_rows_total = treatment_rows_total,
                    control_n = control_n,
                    treatment_n = treatment_n,
                    control_mean_signal = control_mean_signal,
                    treatment_mean_signal = treatment_mean_signal,
                    control_mean_normalized = control_mean_normalized,
                    treatment_mean_normalized = treatment_mean_normalized,
                    effect_estimate_log2 = NA_real_,
                    effect_se_log2 = NA_real_,
                    fold_change_ratio = NA_real_,
                    fold_change_status = "not_available",
                    fold_change_note = "Fold change unavailable because the analyte was not testable in the normalized t-test.",
                    raw_p_value = NA_real_,
                    test_status = "not_testable",
                    test_reason = sprintf("requires >= %s finite normalized values per arm", min_reps),
                    low_replication_warning = control_n < 3 || treatment_n < 3,
                    alpha = alpha
                ))
            }

            model_dat <- dat %>%
                filter(valid_rows)

            t_test <- tryCatch(
                stats::t.test(normalized_signal ~ treatment, data = model_dat, var.equal = TRUE, alternative = "two.sided"),
                error = function(e) e
            )

            if (inherits(t_test, "error")) {
                return(tibble(
                    control = control_label,
                    treatment = treatment_label,
                    subgroup = subgroup_label,
                    subgroup_var = subgroup_var,
                    control_rows_total = control_rows_total,
                    treatment_rows_total = treatment_rows_total,
                    control_n = control_n,
                    treatment_n = treatment_n,
                    control_mean_signal = control_mean_signal,
                    treatment_mean_signal = treatment_mean_signal,
                    control_mean_normalized = control_mean_normalized,
                    treatment_mean_normalized = treatment_mean_normalized,
                    effect_estimate_log2 = NA_real_,
                    effect_se_log2 = NA_real_,
                    fold_change_ratio = NA_real_,
                    fold_change_status = "not_available",
                    fold_change_note = "Fold change unavailable because the normalized t-test could not be fit for this analyte.",
                    raw_p_value = NA_real_,
                    test_status = "not_testable",
                    test_reason = conditionMessage(t_test),
                    low_replication_warning = control_n < 3 || treatment_n < 3,
                    alpha = alpha
                ))
            }

            fold_change_ratio <- if (
                is.finite(control_mean_normalized) &&
                is.finite(treatment_mean_normalized) &&
                control_mean_normalized > 0 &&
                treatment_mean_normalized > 0
            ) {
                treatment_mean_normalized / control_mean_normalized
            } else {
                NA_real_
            }
            effect_estimate_log2 <- if (!is.na(fold_change_ratio) && is.finite(fold_change_ratio) && fold_change_ratio > 0) {
                log2(fold_change_ratio)
            } else {
                NA_real_
            }
            fold_change_status <- if (is.finite(fold_change_ratio) && fold_change_ratio > 0) {
                "reported"
            } else {
                "not_ratio_scale"
            }
            fold_change_note <- if (fold_change_status == "reported") {
                NA_character_
            } else {
                paste(
                    "Fold change omitted because one or both group means were nonpositive after normalization,",
                    "so a ratio-scale fold change is not interpretable."
                )
            }
            effect_se_log2 <- pooled_log2_ratio_se(
                control_values = model_dat$normalized_signal[model_dat$treatment == control_label],
                treatment_values = model_dat$normalized_signal[model_dat$treatment == treatment_label],
                control_mean = control_mean_normalized,
                treatment_mean = treatment_mean_normalized
            )

            tibble(
                control = control_label,
                treatment = treatment_label,
                subgroup = subgroup_label,
                subgroup_var = subgroup_var,
                control_rows_total = control_rows_total,
                treatment_rows_total = treatment_rows_total,
                control_n = control_n,
                treatment_n = treatment_n,
                control_mean_signal = control_mean_signal,
                treatment_mean_signal = treatment_mean_signal,
                control_mean_normalized = control_mean_normalized,
                treatment_mean_normalized = treatment_mean_normalized,
                effect_estimate_log2 = effect_estimate_log2,
                effect_se_log2 = effect_se_log2,
                fold_change_ratio = fold_change_ratio,
                fold_change_status = fold_change_status,
                fold_change_note = fold_change_note,
                raw_p_value = unname(t_test$p.value),
                test_status = "tested",
                test_reason = NA_character_,
                low_replication_warning = control_n < 3 || treatment_n < 3,
                alpha = alpha
            )
        }) %>%
        ungroup() %>%
        left_join(low_signal_flags, by = c("Name", "Coordinate"))
}

#' Validate that configured treatment contrasts exist within each subgroup
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#' @param comparisons Named list of control-to-treatment contrasts.
#' @param subgroup_var Character scalar naming the subgroup column.
#'
#' @return Invisibly returns `NULL`; called for side effects.
validate_requested_comparisons <- function(sample_data, comparisons, subgroup_var) {
    subgroup_levels <- sample_data %>%
        distinct(.data[[subgroup_var]]) %>%
        pull(1) %>%
        as.character()

    for (subgroup_level in subgroup_levels) {
        subgroup_treatments <- sample_data %>%
            filter(.data[[subgroup_var]] == subgroup_level) %>%
            distinct(treatment) %>%
            pull(treatment) %>%
            as.character()

        for (control_label in names(comparisons)) {
            if (!control_label %in% subgroup_treatments) {
                stop(sprintf(
                    "Configured control label '%s' was not found in subgroup '%s'. Available treatment labels: %s",
                    control_label,
                    subgroup_level,
                    paste(subgroup_treatments, collapse = ", ")
                ))
            }

            for (treatment_label in comparisons[[control_label]]) {
                if (!treatment_label %in% subgroup_treatments) {
                    stop(sprintf(
                        "Configured treatment label '%s' was not found in subgroup '%s'. Available treatment labels: %s",
                        treatment_label,
                        subgroup_level,
                        paste(subgroup_treatments, collapse = ", ")
                    ))
                }
            }
        }
    }

    invisible(NULL)
}

#' Run replicate-aware inferential analysis within subgroup strata
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#' @param comparisons Named list of control-to-treatment contrasts.
#' @param subgroup_var Character scalar naming the subgroup column.
#' @param p_adjust_method Character scalar naming the `p.adjust()` method.
#' @param alpha Numeric alpha threshold.
#' @param min_reps Integer minimum replicates per arm to run inference.
#' @param low_signal_threshold Optional raw-signal threshold for flagging.
#'
#' @return Named list of per-comparison result tibbles plus a run index.
run_within_stratum_differential_analysis <- function(sample_data, comparisons, subgroup_var, p_adjust_method = "BH", alpha = 0.05, min_reps = 2, low_signal_threshold = NULL, analysis_method = "raw_log2_lm") {
    validate_p_adjust_method(p_adjust_method)
    validate_requested_comparisons(sample_data, comparisons, subgroup_var)
    method_spec <- get_inferential_method_spec(analysis_method)

    if (!method_spec$value_column %in% names(sample_data)) {
        stop(sprintf(
            "Inferential method '%s' requires sample_data column '%s'.",
            analysis_method,
            method_spec$value_column
        ))
    }
    if (analysis_method == "normalized_t_test" &&
        all(!is.finite(sample_data[[method_spec$value_column]]))) {
        stop(paste(
            "Normalized t-test requires finite normalized sheet values.",
            "No finite normalized_signal values were found in the imported sample data."
        ))
    }

    subgroup_levels <- sample_data %>%
        distinct(.data[[subgroup_var]]) %>%
        pull(1) %>%
        as.character()

    result_tables <- list()
    run_index <- list()

    for (subgroup_level in subgroup_levels) {
        subgroup_df <- sample_data %>%
            filter(.data[[subgroup_var]] == subgroup_level)

        for (control_label in names(comparisons)) {
            for (treatment_label in comparisons[[control_label]]) {
                comparison_df <- subgroup_df %>%
                    filter(treatment %in% c(control_label, treatment_label)) %>%
                    mutate(treatment = factor(treatment, levels = c(control_label, treatment_label))) %>%
                    filter(!startsWith(Name, "Reference")) %>%
                    filter(!startsWith(Name, "Negative"))

                fit_fun <- switch(
                    analysis_method,
                    raw_log2_lm = fit_one_comparison_raw_log2_lm,
                    normalized_t_test = fit_one_comparison_normalized_t_test
                )

                result_tbl <- fit_fun(
                    comparison_df = comparison_df,
                    control_label = control_label,
                    treatment_label = treatment_label,
                    subgroup_label = subgroup_level,
                    subgroup_var = subgroup_var,
                    alpha = alpha,
                    min_reps = min_reps,
                    low_signal_threshold = low_signal_threshold
                )

                tested_mask <- result_tbl$test_status == "tested" & is.finite(result_tbl$raw_p_value)
                adjusted_p_values <- rep(NA_real_, nrow(result_tbl))
                adjusted_p_values[tested_mask] <- p.adjust(result_tbl$raw_p_value[tested_mask], method = p_adjust_method)
                n_p_adjust_hypotheses <- sum(tested_mask)

                result_tbl <- result_tbl %>%
                    mutate(
                        adjusted_p_value = adjusted_p_values,
                        p_adjust_method = p_adjust_method,
                        n_p_adjust_hypotheses = n_p_adjust_hypotheses,
                        raw_p_lt_alpha = if_else(!is.na(raw_p_value), raw_p_value < alpha, FALSE),
                        fdr_lt_0_20 = if_else(!is.na(adjusted_p_value), adjusted_p_value < 0.20, FALSE),
                        fdr_lt_0_25 = if_else(!is.na(adjusted_p_value), adjusted_p_value < 0.25, FALSE),
                        comparison_label = sprintf("%s: %s vs %s", subgroup_level, control_label, treatment_label),
                        analysis_method = analysis_method,
                        analysis_method_label = method_spec$label,
                        value_column_used = method_spec$value_column,
                        waterfall_subtitle = method_spec$waterfall_subtitle
                    ) %>%
                    arrange(raw_p_value, desc(abs(effect_estimate_log2)))

                comparison_slug <- make_comparison_slug(
                    subgroup_label = subgroup_level,
                    control_label = control_label,
                    treatment_label = treatment_label
                )

                result_tables[[comparison_slug]] <- result_tbl
                run_index[[comparison_slug]] <- tibble(
                    comparison_slug = comparison_slug,
                    subgroup = subgroup_level,
                    control = control_label,
                    treatment = treatment_label,
                    analysis_method = analysis_method,
                    analysis_method_label = method_spec$label,
                    n_analytes = nrow(result_tbl),
                    n_raw_p_lt_alpha = sum(result_tbl$raw_p_lt_alpha, na.rm = TRUE),
                    n_fdr_lt_0_20 = sum(result_tbl$fdr_lt_0_20, na.rm = TRUE),
                    n_fdr_lt_0_25 = sum(result_tbl$fdr_lt_0_25, na.rm = TRUE),
                    any_low_replication_warning = any(result_tbl$low_replication_warning, na.rm = TRUE),
                    min_control_n = min(result_tbl$control_n, na.rm = TRUE),
                    min_treatment_n = min(result_tbl$treatment_n, na.rm = TRUE)
                )
            }
        }
    }

    list(
        results = result_tables,
        run_index = bind_rows(run_index)
    )
}

#' Run one or more replicate-aware inferential methods
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#' @param comparisons Named list of control-to-treatment contrasts.
#' @param subgroup_var Character scalar naming the subgroup column.
#' @param analysis_methods Character vector of supported method slugs.
#' @param p_adjust_method Character scalar naming the `p.adjust()` method.
#' @param alpha Numeric alpha threshold.
#' @param min_reps Integer minimum replicates per arm to run inference.
#' @param low_signal_threshold Optional raw-signal threshold for flagging.
#'
#' @return Named list containing per-method inferential outputs.
run_replicate_analysis_methods <- function(sample_data, comparisons, subgroup_var, analysis_methods = "raw_log2_lm", p_adjust_method = "BH", alpha = 0.05, min_reps = 2, low_signal_threshold = NULL) {
    analysis_methods <- validate_analysis_methods(analysis_methods)

    method_results <- setNames(
        map(analysis_methods, function(method) {
            run_within_stratum_differential_analysis(
                sample_data = sample_data,
                comparisons = comparisons,
                subgroup_var = subgroup_var,
                p_adjust_method = p_adjust_method,
                alpha = alpha,
                min_reps = min_reps,
                low_signal_threshold = low_signal_threshold,
                analysis_method = method
            )
        }),
        analysis_methods
    )

    method_index <- map_dfr(analysis_methods, function(method) {
        method_spec <- get_inferential_method_spec(method)
        run_index <- method_results[[method]]$run_index

        tibble(
            analysis_method = method,
            analysis_method_label = method_spec$label,
            n_comparisons = nrow(run_index),
            total_raw_p_lt_alpha = sum(run_index$n_raw_p_lt_alpha, na.rm = TRUE),
            total_fdr_lt_0_20 = sum(run_index$n_fdr_lt_0_20, na.rm = TRUE),
            total_fdr_lt_0_25 = sum(run_index$n_fdr_lt_0_25, na.rm = TRUE)
        )
    })

    list(
        primary_method = analysis_methods[[1]],
        methods = method_results,
        method_index = method_index
    )
}

#' Resolve which inferential comparison(s) should feed the shortlist workflow
#'
#' @param run_index Tibble describing available inferential comparisons.
#' @param example_config Analysis config list.
#'
#' @return Tibble with one row per selected comparison.
resolve_shortlist_comparisons <- function(run_index, example_config) {
    configured_slugs <- example_config$selection_comparison_slugs
    if (is.null(configured_slugs) && !is.null(example_config$selection_comparison_slug)) {
        configured_slugs <- example_config$selection_comparison_slug
    }
    if (!is.null(configured_slugs)) {
        configured_slugs <- unique(unname(as.character(configured_slugs)))
        selected <- run_index %>%
            filter(comparison_slug %in% configured_slugs)

        missing_slugs <- setdiff(configured_slugs, selected$comparison_slug)
        if (length(missing_slugs) > 0) {
            stop(sprintf(
                paste(
                    "Configured shortlist comparison slug(s) were not found: %s",
                    "Available slugs: %s"
                ),
                paste(missing_slugs, collapse = ", "),
                paste(run_index$comparison_slug, collapse = ", ")
            ))
        }

        return(
            tibble(comparison_slug = configured_slugs) %>%
                left_join(selected, by = "comparison_slug")
        )
    }

    if (!is.null(example_config$selection_subgroup) &&
        !is.null(example_config$selection_control) &&
        !is.null(example_config$selection_group)) {
        derived_slug <- make_comparison_slug(
            subgroup_label = example_config$selection_subgroup,
            control_label = example_config$selection_control,
            treatment_label = example_config$selection_group
        )
        selected <- run_index %>%
            filter(comparison_slug == derived_slug)
        if (nrow(selected) == 0) {
            stop(sprintf(
                "Derived shortlist comparison '%s' was not found. Available slugs: %s",
                derived_slug,
                paste(run_index$comparison_slug, collapse = ", ")
            ))
        }
        return(slice_head(selected, n = 1))
    }

    if (nrow(run_index) == 1) {
        return(slice_head(run_index, n = 1))
    }

    stop(sprintf(
        paste(
            "Shortlist selection is ambiguous because multiple inferential comparisons are available.",
            "Set shortlist$comparisons or selection_subgroup/selection_control/selection_group in the analysis config.",
            "Available slugs: %s"
        ),
        paste(run_index$comparison_slug, collapse = ", ")
    ))
}

resolve_inferential_shortlist_methods <- function(example_config) {
    configured_method <- example_config$selection_method
    if (is.null(configured_method) || length(configured_method) == 0 || all(is.na(configured_method)) || all(configured_method == "")) {
        configured_method <- example_config$analysis_methods
    }
    validate_analysis_methods(configured_method)
}

read_inferential_comparison_results <- function(inferential_dir, comparison_slug, analysis_method) {
    workbook_path <- file.path(inferential_dir, sprintf("%s_results.xlsx", analysis_method))
    if (!file.exists(workbook_path)) {
        stop(sprintf(
            paste(
                "Missing workbook for shortlist method '%s': %s",
                "Run scripts/main.R first or choose a method that has been written."
            ),
            analysis_method,
            workbook_path
        ))
    }

    available_sheets <- readxl::excel_sheets(workbook_path)
    if (!comparison_slug %in% available_sheets) {
        stop(sprintf(
            "Workbook '%s' does not contain comparison sheet '%s'. Available sheets: %s",
            workbook_path,
            comparison_slug,
            paste(available_sheets, collapse = ", ")
        ))
    }

    readxl::read_excel(workbook_path, sheet = comparison_slug)
}

inferential_plot_tiers <- function(alpha) {
    list(
        list(
            flag_column = "raw_p_lt_alpha",
            file_stub = "raw_p_lt_alpha",
            title_prefix = sprintf("Raw p < %.3g Hits", alpha),
            barplot_label = sprintf("raw p < %.3g", alpha)
        ),
        list(
            flag_column = "fdr_lt_0_20",
            file_stub = "fdr_lt_0_20",
            title_prefix = "FDR < 0.20 Hits",
            barplot_label = "FDR < 0.20"
        ),
        list(
            flag_column = "fdr_lt_0_25",
            file_stub = "fdr_lt_0_25",
            title_prefix = "FDR < 0.25 Hits",
            barplot_label = "FDR < 0.25"
        )
    )
}

write_threshold_waterfall_set <- function(result_tbl, output_dir, filename_prefix = NULL, title_suffix = NULL) {
    plot_specs <- inferential_plot_tiers(unique(result_tbl$alpha))
    output_paths <- list()
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    invisible(walk(plot_specs, function(spec) {
        filtered_tbl <- result_tbl %>%
            filter(test_status == "tested", .data[[spec$flag_column]])

        if (is.null(filename_prefix)) {
            plot_filename <- sprintf("waterfall_%s.png", spec$file_stub)
        } else {
            plot_filename <- sprintf("%s_waterfall_%s.png", filename_prefix, spec$file_stub)
        }

        plot_path <- file.path(output_dir, plot_filename)

        if (nrow(filtered_tbl) > 0) {
            plot_title <- c(spec$title_prefix, unique(result_tbl$comparison_label), title_suffix) %>%
                discard(~ is.null(.x) || identical(.x, "")) %>%
                paste(collapse = "\n")

            plot_obj <- plot_inferential_waterfall(
                filtered_tbl,
                title = plot_title
            )
            ggsave(
                filename = plot_path,
                plot = plot_obj,
                width = max(12, 0.6 * nrow(filtered_tbl)),
                height = 10
            )
            output_paths[[spec$file_stub]] <<- plot_path
        } else {
            if (file.exists(plot_path)) {
                unlink(plot_path, force = TRUE)
            }
            output_paths[[spec$file_stub]] <<- NA_character_
        }
    }))

    output_paths
}

build_inferential_fold_change_barplot_data <- function(hit_tbl, mark_treatment_significant = TRUE, significance_label = "*", significance_definition = "threshold hit") {
    hit_tbl <- hit_tbl %>%
        filter(is.finite(fold_change_ratio))

    if (nrow(hit_tbl) == 0) {
        return(tibble())
    }

    control_label <- unique(as.character(hit_tbl$control))
    treatment_label <- unique(as.character(hit_tbl$treatment))
    if (length(control_label) != 1 || length(treatment_label) != 1) {
        stop("Expected one control label and one treatment label when building inferential fold-change barplots.")
    }
    group_levels <- c(control_label, treatment_label)

    bind_rows(
        hit_tbl %>%
            transmute(
                Name,
                Coordinate,
                group = as.character(control),
                relative_signal = 1,
                significant = FALSE,
                significance_label = NA_character_,
                significance_definition = NA_character_
            ),
        hit_tbl %>%
            transmute(
                Name,
                Coordinate,
                group = as.character(treatment),
                relative_signal = fold_change_ratio,
                significant = mark_treatment_significant,
                significance_label = if_else(mark_treatment_significant, significance_label, NA_character_),
                significance_definition = if_else(mark_treatment_significant, significance_definition, NA_character_)
            )
    ) %>%
        mutate(
            group = factor(group, levels = group_levels),
            short_group = factor(as.character(group), levels = levels(group)),
            display_name = str_split(as.character(Name), pattern = "/", simplify = TRUE)[, 1],
            display_name = ifelse(display_name == "Complement Component C5", "C5a", display_name),
            Name_Coordinate = str_c(display_name, "\n(", Coordinate, ")")
        )
}

plot_inferential_fold_change_barplot_pages <- function(barplot_data, title, groups_per_page = 25, ncol = NULL, nrow = NULL, common_y_max = NULL) {
    if (nrow(barplot_data) == 0) {
        return(list())
    }

    if (is.null(ncol) || is.null(nrow)) {
        ncol <- ceiling(sqrt(groups_per_page))
        nrow <- ceiling(groups_per_page / ncol)
    }
    groups_per_page <- ncol * nrow

    treatment_group <- setdiff(levels(barplot_data$group), levels(barplot_data$group)[[1]])[[1]]
    ordered_names <- barplot_data %>%
        filter(group == treatment_group) %>%
        arrange(desc(relative_signal)) %>%
        pull(Name_Coordinate)
    if (length(ordered_names) == 0) {
        ordered_names <- unique(barplot_data$Name_Coordinate)
    }

    barplot_data <- barplot_data %>%
        mutate(Name_Coordinate = factor(Name_Coordinate, levels = ordered_names))

    total_pages <- ceiling(length(ordered_names) / groups_per_page)
    fill_values <- setNames(c("#FFFFFF", "#737373")[seq_along(levels(barplot_data$group))], levels(barplot_data$group))
    color_values <- setNames(c("#000000", "#000000")[seq_along(levels(barplot_data$group))], levels(barplot_data$group))
    control_label <- levels(barplot_data$group)[[1]]

    y_axis_title <- sprintf("Fold Change to %s", control_label)
    significance_definitions <- barplot_data %>%
        filter(significant, !is.na(significance_label), !is.na(significance_definition)) %>%
        distinct(significance_label, significance_definition) %>%
        arrange(significance_label, significance_definition) %>%
        mutate(caption_text = sprintf("%s = %s", significance_label, significance_definition)) %>%
        pull(caption_text)
    if (is.null(common_y_max)) {
        finite_set_signals <- barplot_data$relative_signal[is.finite(barplot_data$relative_signal)]
        set_signal_max <- max(finite_set_signals, 1, na.rm = TRUE)
        if (!is.finite(set_signal_max) || set_signal_max <= 0) {
            set_signal_max <- 1
        }
        common_y_max <- max(pretty(c(0, set_signal_max * 1.22), n = 5), na.rm = TRUE)
    }

    map(seq_len(total_pages), function(page_i) {
        start_idx <- (page_i - 1) * groups_per_page + 1
        end_idx <- min(page_i * groups_per_page, length(ordered_names))
        page_names <- ordered_names[start_idx:end_idx]
        spacer_names <- character()
        if (length(page_names) < groups_per_page) {
            spacer_names <- sprintf("IGNORE_SPACER_%02d", seq_len(groups_per_page - length(page_names)))
        }

        page_data <- barplot_data %>%
            filter(Name_Coordinate %in% page_names) %>%
            mutate(
                Name_Coordinate = factor(Name_Coordinate, levels = page_names)
            )

        facet_labeler <- function(values) {
            wrapped_values <- str_wrap(as.character(values), width = 10)
            if_else(str_detect(as.character(values), "^IGNORE_SPACER_"), "", wrapped_values)
        }

        make_panel_plot <- function(panel_name) {
            if (str_detect(panel_name, "^IGNORE_SPACER_")) {
                return(patchwork::plot_spacer())
            }

            panel_title <- facet_labeler(panel_name)
            panel_data <- page_data %>%
                filter(as.character(Name_Coordinate) == panel_name, is.finite(relative_signal))
            labels_df <- panel_data %>%
                mutate(
                    label = as.character(round(relative_signal, 1)),
                    label_y = pmax(relative_signal, 0) / 6
                )
            finite_relative_signals <- panel_data$relative_signal[is.finite(panel_data$relative_signal)]
            panel_max <- max(finite_relative_signals, 1, na.rm = TRUE)
            significance_df <- panel_data %>%
                filter(significant, !is.na(significance_label)) %>%
                summarize(
                    group1 = as.character(levels(short_group)[[1]]),
                    group2 = as.character(levels(short_group)[[2]]),
                    y.position = panel_max * 1.12,
                    label = first(significance_label),
                    .groups = "drop"
                ) %>%
                filter(!is.na(label))

            title_plot <- ggplot() +
                annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = panel_title,
                    size = 4,
                    fontface = "bold",
                    lineheight = 0.9
                ) +
                xlim(0, 1) +
                ylim(0, 1) +
                theme_void() +
                theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

            panel_plot <- ggplot(panel_data, aes(x = short_group, y = relative_signal, fill = group, color = group)) +
                geom_col(linewidth = rel(1.1)) +
                geom_text(
                    data = labels_df,
                    aes(x = short_group, y = label_y, label = label),
                    size = 6,
                    hjust = -0.5,
                    angle = 90,
                    color = "black",
                    inherit.aes = FALSE
                ) +
                scale_y_continuous(
                    limits = c(0, common_y_max),
                    expand = expansion(mult = c(0, 0.02))
                ) +
                theme_prism(base_size = 16) +
                scale_fill_manual(values = fill_values) +
                scale_color_manual(values = color_values) +
                labs(
                    x = NULL,
                    y = NULL
                ) +
                theme(
                    plot.title = element_blank(),
                    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5),
                    axis.text = element_text(size = rel(1.75)),
                    panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
                    panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
                    panel.spacing = unit(1.25, "lines")
                ) +
                scale_x_discrete(labels = NULL)

            if (nrow(significance_df) > 0) {
                panel_plot <- panel_plot +
                    ggpubr::stat_pvalue_manual(
                        significance_df,
                        label = "label",
                        xmin = "group1",
                        xmax = "group2",
                        y.position = "y.position",
                        tip.length = 0.01,
                        bracket.size = 0.7,
                        size = 12,
                        inherit.aes = FALSE
                    )
            }

            title_plot / panel_plot +
                patchwork::plot_layout(heights = c(0.26, 1))
        }

        panel_plots <- map(c(as.character(page_names), spacer_names), make_panel_plot)
        plot_grid <- patchwork::wrap_plots(
            panel_plots,
            ncol = ncol,
            nrow = nrow,
            guides = "collect"
        ) &
            theme(legend.position = "right")

        plot_grid +
            patchwork::plot_annotation(
                title = title,
                subtitle = y_axis_title,
                caption = paste(c(significance_definitions, sprintf("N = %s / %s", end_idx, length(ordered_names))), collapse = "; ")
            ) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.35)),
                plot.subtitle = element_text(hjust = 0.5, face = "bold", size = rel(1.05)),
                plot.caption = element_text(hjust = 1)
            )
    })
}

#' Return a filesystem-safe filename stem for one analyte
#'
#' @param analyte_name Character scalar naming one analyte.
#'
#' @return Character scalar safe for use as a PNG filename stem.
sanitize_selected_analyte_filename <- function(analyte_name) {
    analyte_name %>%
        str_split(pattern = "/", simplify = TRUE) %>%
        .[, 1] %>%
        str_replace_all("[\\\\/?*\\[\\]:]", "_") %>%
        str_replace_all("\\s+", " ") %>%
        str_trim()
}

#' Build barplot-ready fold-change rows from legacy selected-analyte data
#'
#' @param selected_dat Tibble with legacy pairwise rows from `make_wf_data()`.
#' @param mark_threshold_hits Logical indicating whether treatment bars should
#'   be marked when they exceed the configured fold-change threshold.
#' @param threshold Numeric fold-change threshold used only when
#'   `mark_threshold_hits = TRUE`.
#'
#' @return Tibble in the shape expected by
#'   `plot_inferential_fold_change_barplot_pages()`.
build_legacy_selected_fold_change_barplot_data <- function(selected_dat, mark_threshold_hits = FALSE, threshold = NULL) {
    if (nrow(selected_dat) == 0) {
        return(tibble())
    }

    control_label <- unique(as.character(selected_dat$control))
    if (length(control_label) != 1) {
        stop("Expected one control label when building legacy selected-analyte barplots.")
    }

    selected_dat %>%
        mutate(
            group = factor(as.character(group), levels = unique(as.character(group))),
            short_group = factor(as.character(group), levels = levels(group)),
            display_name = str_split(as.character(Name), pattern = "/", simplify = TRUE)[, 1],
            display_name = ifelse(display_name == "Complement Component C5", "C5a", display_name),
            Name_Coordinate = str_c(display_name, "\n(", Coordinate, ")"),
            significant = if (isTRUE(mark_threshold_hits) && !is.null(threshold)) {
                as.character(group) != control_label & is.finite(relative_signal) &
                    (relative_signal > threshold | relative_signal < 1 / threshold)
            } else {
                FALSE
            },
            significance_label = NA_character_,
            significance_definition = NA_character_
        ) %>%
        select(Name, Coordinate, group, short_group, relative_signal, significant, significance_label, significance_definition, display_name, Name_Coordinate)
}

#' Write one selected-analyte fold-change bargraph per analyte
#'
#' This uses the same renderer as the inferential 25-per-page barplots, but
#' renders each selected analyte as a one-panel page. The y-axis maximum is
#' fixed across the whole selected set so the individual PNGs remain
#' comparable within a comparison/method.
#'
#' @param barplot_data Tibble returned by a selected-analyte barplot-data
#'   builder.
#' @param output_dir Directory where `selected_bargraphs/` should be written.
#' @param title Character title to use on each single-analyte plot.
#' @param height Numeric plot height in inches.
#' @param width Numeric plot width in inches.
#'
#' @return Tibble with one row per written PNG.
write_selected_fold_change_bargraphs <- function(barplot_data, output_dir, title, height = 6.5, width = 6.5) {
    bargraph_dir <- file.path(output_dir, "selected_bargraphs")
    if (dir.exists(bargraph_dir)) {
        unlink(bargraph_dir, recursive = TRUE, force = TRUE)
    }
    dir.create(bargraph_dir, recursive = TRUE, showWarnings = FALSE)

    plot_data <- barplot_data %>%
        filter(is.finite(relative_signal))
    if (nrow(plot_data) == 0) {
        return(tibble(Name = character(), bargraph_path = character()))
    }

    finite_set_signals <- plot_data$relative_signal[is.finite(plot_data$relative_signal)]
    set_signal_max <- max(finite_set_signals, 1, na.rm = TRUE)
    if (!is.finite(set_signal_max) || set_signal_max <= 0) {
        set_signal_max <- 1
    }
    selected_common_y_max <- max(pretty(c(0, set_signal_max * 1.22), n = 5), na.rm = TRUE)

    unique_names <- plot_data %>%
        distinct(Name, display_name) %>%
        arrange(match(Name, unique(barplot_data$Name)))

    map_dfr(seq_len(nrow(unique_names)), function(row_idx) {
        selected_name <- unique_names$Name[[row_idx]]
        filename_stem <- sanitize_selected_analyte_filename(unique_names$display_name[[row_idx]])
        analyte_plot_data <- plot_data %>%
            filter(Name == selected_name)
        plot_page <- plot_inferential_fold_change_barplot_pages(
            barplot_data = analyte_plot_data,
            title = title,
            groups_per_page = 1,
            ncol = 1,
            nrow = 1,
            common_y_max = selected_common_y_max
        )[[1]]
        plot_path <- file.path(bargraph_dir, sprintf("%s.png", filename_stem))
        ggsave(
            filename = plot_path,
            plot = plot_page,
            height = height,
            width = width,
            dpi = 150
        )

        tibble(Name = selected_name, bargraph_path = plot_path)
    })
}

write_threshold_bargraph_set <- function(result_tbl, output_dir, filename_prefix, groups_per_page = 25) {
    plot_specs <- c(
        list(list(
            flag_column = NULL,
            file_stub = "all_tested",
            title_prefix = "All Tested Analytes",
            mark_treatment_significant = FALSE,
            barplot_label = NA_character_
        )),
        map(inferential_plot_tiers(unique(result_tbl$alpha)), function(spec) {
            spec$mark_treatment_significant <- TRUE
            spec
        })
    )
    output_paths <- list()
    ncol <- ceiling(sqrt(groups_per_page))
    nrow <- ceiling(groups_per_page / ncol)
    groups_per_page <- ncol * nrow
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stale_flat_paths <- list.files(
        output_dir,
        pattern = sprintf("^%s_barplot_.*_page_[0-9]+\\.png$", filename_prefix),
        full.names = TRUE
    )
    stale_legacy_dirs <- list.files(
        output_dir,
        pattern = sprintf("^%s_bargraphs_", filename_prefix),
        full.names = TRUE
    )
    invisible(walk(c(stale_flat_paths, stale_legacy_dirs), function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))

    invisible(walk(plot_specs, function(spec) {
        plot_output_dir <- if (identical(spec$file_stub, "all_tested")) {
            file.path(output_dir, "all_tested")
        } else {
            file.path(output_dir, "significant_hits", spec$file_stub)
        }
        if (dir.exists(plot_output_dir)) {
            unlink(plot_output_dir, recursive = TRUE, force = TRUE)
        }

        hit_tbl <- result_tbl %>%
            filter(test_status == "tested")
        if (!is.null(spec$flag_column)) {
            hit_tbl <- hit_tbl %>%
                filter(.data[[spec$flag_column]])
        }

        barplot_data <- build_inferential_fold_change_barplot_data(
            hit_tbl,
            mark_treatment_significant = spec$mark_treatment_significant,
            significance_label = "*",
            significance_definition = spec$barplot_label
        )

        if (nrow(barplot_data) > 0) {
            dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
            plot_pages <- plot_inferential_fold_change_barplot_pages(
                barplot_data = barplot_data,
                title = sprintf(
                    "%s\n%s\n%s",
                    spec$title_prefix,
                    unique(result_tbl$comparison_label),
                    unique(result_tbl$analysis_method_label)
                ),
                groups_per_page = groups_per_page,
                ncol = ncol,
                nrow = nrow
            )

            plot_paths <- map_chr(seq_along(plot_pages), function(page_i) {
                plot_path <- file.path(
                    plot_output_dir,
                    sprintf("%s_barplot_%s_page_%s.png", filename_prefix, spec$file_stub, page_i)
                )
                temp_plot_path <- tempfile(fileext = ".png")
                ggsave(
                    filename = temp_plot_path,
                    plot = plot_pages[[page_i]],
                    height = 3 + (nrow * 3.1),
                    width = ncol * 3.4,
                    dpi = 150
                )
                dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
                copied <- file.copy(temp_plot_path, plot_path, overwrite = TRUE)
                unlink(temp_plot_path, force = TRUE)
                if (!isTRUE(copied)) {
                    stop(sprintf("Could not copy rendered barplot to '%s'.", plot_path))
                }
                plot_path
            })
            output_paths[[spec$file_stub]] <<- paste(plot_paths, collapse = ";")
        } else {
            output_paths[[spec$file_stub]] <<- NA_character_
        }
    }))

    output_paths
}

#' Plot inferential effect estimates as a waterfall chart
#'
#' @param results_tbl One per-comparison inferential results tibble.
#' @param title Character plot title.
#'
#' @return ggplot object.
plot_inferential_waterfall <- function(results_tbl, title) {
    if (!"effect_se_log2" %in% names(results_tbl)) {
        results_tbl <- results_tbl %>%
            mutate(effect_se_log2 = NA_real_)
    }

    plot_data <- results_tbl %>%
        filter(test_status == "tested", is.finite(effect_estimate_log2)) %>%
        arrange(desc(effect_estimate_log2)) %>%
        mutate(
            Name = factor(Name, levels = Name),
            effect_se_log2 = if_else(is.finite(effect_se_log2) & effect_se_log2 >= 0, effect_se_log2, NA_real_)
        )

    if (nrow(plot_data) == 0) {
        stop("Cannot build inferential waterfall plot because there are no tested analytes with finite effect estimates.")
    }

    interval_bounds <- c(
        plot_data$effect_estimate_log2,
        plot_data$effect_estimate_log2 - plot_data$effect_se_log2,
        plot_data$effect_estimate_log2 + plot_data$effect_se_log2
    )
    fill_limit <- max(abs(interval_bounds[is.finite(interval_bounds)]), na.rm = TRUE)
    has_effect_se <- any(is.finite(plot_data$effect_se_log2))
    subtitle_text <- unique(results_tbl$waterfall_subtitle)
    if (length(subtitle_text) == 0 || is.na(subtitle_text[[1]]) || identical(subtitle_text[[1]], "")) {
        subtitle_text <- "Effect estimate is treatment minus control on the log2 signal scale"
    } else {
        subtitle_text <- subtitle_text[[1]]
    }

    ggplot(plot_data, aes(x = Name, y = effect_estimate_log2, fill = effect_estimate_log2)) +
        geom_col(show.legend = FALSE) +
        geom_errorbar(
            data = plot_data %>% filter(is.finite(effect_se_log2)),
            aes(
                ymin = effect_estimate_log2 - effect_se_log2,
                ymax = effect_estimate_log2 + effect_se_log2
            ),
            inherit.aes = TRUE,
            width = 0.25,
            linewidth = 0.35,
            color = "black"
        ) +
        # Negative effects trend blue, positive effects trend red, and
        # near-zero effects fade toward gray.
        scale_fill_gradient2(
            low = "blue",
            mid = "gray",
            high = "red",
            midpoint = 0,
            limits = c(-fill_limit, fill_limit)
        ) +
        labs(
            title = str_wrap(title, width = 55),
            subtitle = str_wrap(subtitle_text, width = 70),
            caption = sprintf(
                "N analytes plotted = %s%s",
                nrow(plot_data),
                ifelse(has_effect_se, "; whiskers = +/- 1 SE of the effect estimate", "")
            )
        ) +
        ylab("Effect estimate (log2 signal)") +
        xlab("") +
        theme_prism(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = element_text(size = rel(1.1)),
            plot.subtitle = element_text(size = rel(0.9)),
            plot.margin = margin(12, 20, 12, 12),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5))
        )
}

#' Write inferential result tables and plots to disk
#'
#' @param inferential_results Output of `run_within_stratum_differential_analysis()`.
#' @param output_dir Directory where inferential outputs should be written.
#'
#' @return Invisibly returns the output paths written.
build_inferential_workbook_filename <- function(run_index) {
    workbook_method <- run_index %>%
        distinct(analysis_method) %>%
        pull(analysis_method) %>%
        as.character()

    if (length(workbook_method) == 1 && !is.na(workbook_method[[1]]) && workbook_method[[1]] != "") {
        sprintf("%s_results.xlsx", workbook_method[[1]])
    } else {
        "results_workbook.xlsx"
    }
}

inferential_comparison_tables_dir <- function(output_dir, comparison_slug) {
    file.path(output_dir, "comparisons", comparison_slug, "tables")
}

inferential_comparison_waterfall_dir <- function(output_dir, comparison_slug, method_name) {
    file.path(output_dir, "comparisons", comparison_slug, "waterfall_plots", method_name)
}

inferential_comparison_barplot_dir <- function(output_dir, comparison_slug, method_name) {
    file.path(output_dir, "comparisons", comparison_slug, "barplots", method_name)
}

clean_inferential_comparison_dir <- function(comparison_dir, remove_organized_dirs = FALSE) {
    stale_flat_files <- list.files(
        comparison_dir,
        pattern = "^(results|raw_log2_lm_results|normalized_t_test_results|waterfall|raw_log2_lm_waterfall|normalized_t_test_waterfall|raw_log2_lm_barplot|normalized_t_test_barplot).*\\.(tsv|png)$",
        full.names = TRUE
    )
    stale_legacy_dirs <- list.files(
        comparison_dir,
        pattern = "^(raw_log2_lm|normalized_t_test)_bargraphs_",
        full.names = TRUE
    )
    stale_paths <- c(stale_flat_files, stale_legacy_dirs)
    if (remove_organized_dirs) {
        stale_paths <- c(stale_paths, file.path(comparison_dir, c("tables", "waterfall_plots", "barplots")))
    }
    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))
}

write_inferential_workbook <- function(inferential_results, output_dir, run_index = inferential_results$run_index) {
    workbook_sheets <- inferential_results$results
    workbook_sheets$run_index <- run_index

    workbook_path <- file.path(output_dir, build_inferential_workbook_filename(run_index))
    writexl::write_xlsx(workbook_sheets, workbook_path)
    workbook_path
}

write_inferential_outputs <- function(inferential_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stale_workbook_paths <- file.path(output_dir, c(
        "combined_results.xlsx",
        "results_workbook.xlsx",
        "raw_log2_lm_results.xlsx",
        "normalized_t_test_results.xlsx"
    ))
    invisible(walk(stale_workbook_paths, function(path) {
        if (file.exists(path)) {
            unlink(path, force = TRUE)
            if (file.exists(path)) {
                stop(sprintf(
                    paste(
                        "Could not remove stale workbook '%s'.",
                        "Close the file in Excel or another viewer and rerun the analysis to avoid mixing old and new outputs."
                    ),
                    path
                ))
            }
        }
    }))

    output_rows <- map_dfr(names(inferential_results$results), function(comparison_slug) {
        comparison_dir <- file.path(output_dir, "comparisons", comparison_slug)
        dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
        clean_inferential_comparison_dir(comparison_dir, remove_organized_dirs = TRUE)

        result_tbl <- inferential_results$results[[comparison_slug]]
        method_name <- unique(result_tbl$analysis_method)
        if (length(method_name) != 1 || is.na(method_name) || identical(method_name, "")) {
            stop(sprintf(
                "Expected exactly one analysis_method in comparison '%s' when writing inferential outputs.",
                comparison_slug
            ))
        }

        tables_dir <- inferential_comparison_tables_dir(output_dir, comparison_slug)
        waterfall_dir <- inferential_comparison_waterfall_dir(output_dir, comparison_slug, method_name)
        barplot_dir <- inferential_comparison_barplot_dir(output_dir, comparison_slug, method_name)
        dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(waterfall_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)

        result_path <- file.path(tables_dir, sprintf("%s_results.tsv", method_name))
        write_tsv(result_tbl, result_path)

        full_waterfall_path <- file.path(waterfall_dir, sprintf("%s_waterfall.png", method_name))
        waterfall_plot <- plot_inferential_waterfall(
            result_tbl,
            title = sprintf("%s\n%s", unique(result_tbl$comparison_label), unique(result_tbl$analysis_method_label))
        )
        ggsave(
            filename = full_waterfall_path,
            plot = waterfall_plot,
            width = max(12, 0.4 * nrow(result_tbl)),
            height = 10
        )

        threshold_paths <- write_threshold_waterfall_set(
            result_tbl = result_tbl,
            output_dir = waterfall_dir,
            filename_prefix = method_name,
            title_suffix = unique(result_tbl$analysis_method_label)
        )
        bargraph_paths <- write_threshold_bargraph_set(
            result_tbl = result_tbl,
            output_dir = barplot_dir,
            filename_prefix = method_name
        )

        tibble(
            comparison_slug = comparison_slug,
            result_path = result_path,
            full_waterfall_path = full_waterfall_path,
            raw_p_alpha_waterfall_path = threshold_paths$raw_p_lt_alpha,
            fdr_0_20_waterfall_path = threshold_paths$fdr_lt_0_20,
            fdr_0_25_waterfall_path = threshold_paths$fdr_lt_0_25,
            all_tested_barplot_path = bargraph_paths$all_tested,
            raw_p_alpha_barplot_path = bargraph_paths$raw_p_lt_alpha,
            fdr_0_20_barplot_path = bargraph_paths$fdr_lt_0_20,
            fdr_0_25_barplot_path = bargraph_paths$fdr_lt_0_25
        )
    })

    run_index <- inferential_results$run_index %>%
        left_join(output_rows, by = "comparison_slug")
    run_index_path <- file.path(output_dir, "run_index.tsv")
    write_tsv(run_index, run_index_path)

    workbook_path <- write_inferential_workbook(
        inferential_results = inferential_results,
        output_dir = output_dir,
        run_index = run_index
    )

    invisible(list(
        run_index = run_index_path,
        workbook = workbook_path,
        results = output_rows$result_path
    ))
}

write_method_specific_waterfalls <- function(method_results, output_dir, primary_method) {
    comparison_slugs <- method_results[[primary_method]]$run_index %>%
        pull(comparison_slug) %>%
        as.character()

    invisible(walk(comparison_slugs, function(comparison_slug) {
        comparison_dir <- file.path(output_dir, "comparisons", comparison_slug)
        dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
        clean_inferential_comparison_dir(comparison_dir)

        invisible(walk(names(method_results), function(method_name) {
            if (identical(method_name, primary_method)) {
                return(invisible(NULL))
            }

            result_tbl <- method_results[[method_name]]$results[[comparison_slug]]
            tables_dir <- inferential_comparison_tables_dir(output_dir, comparison_slug)
            waterfall_dir <- inferential_comparison_waterfall_dir(output_dir, comparison_slug, method_name)
            barplot_dir <- inferential_comparison_barplot_dir(output_dir, comparison_slug, method_name)
            dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(waterfall_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)

            result_path <- file.path(tables_dir, sprintf("%s_results.tsv", method_name))
            write_tsv(result_tbl, result_path)
            method_plot_path <- file.path(waterfall_dir, sprintf("%s_waterfall.png", method_name))
            method_plot <- plot_inferential_waterfall(
                result_tbl,
                title = sprintf("%s\n%s", unique(result_tbl$comparison_label), unique(result_tbl$analysis_method_label))
            )
            ggsave(
                filename = method_plot_path,
                plot = method_plot,
                width = max(12, 0.4 * nrow(result_tbl)),
                height = 10
            )

            write_threshold_waterfall_set(
                result_tbl = result_tbl,
                output_dir = waterfall_dir,
                filename_prefix = method_name,
                title_suffix = unique(result_tbl$analysis_method_label)
            )
            write_threshold_bargraph_set(
                result_tbl = result_tbl,
                output_dir = barplot_dir,
                filename_prefix = method_name
            )
        }))
    }))
}

#' Write a cross-method comparison table for each inferential comparison
#'
#' @param method_results Named list returned under `$methods` from
#'   `run_replicate_analysis_methods()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written file paths.
write_method_comparison_outputs <- function(method_results, output_dir) {
    build_workbook_sheet_names <- function(labels) {
        used_names <- character()

        vapply(seq_along(labels), function(i) {
            base_name <- labels[[i]] %>%
                str_replace_all("[\\\\/?*\\[\\]:]", "_") %>%
                str_replace_all("\\s+", "_") %>%
                str_sub(1, 31)

            if (base_name == "") {
                base_name <- sprintf("sheet_%02d", i)
            }

            candidate <- base_name
            suffix <- 1
            while (candidate %in% used_names) {
                suffix_text <- sprintf("_%02d", suffix)
                candidate <- str_sub(base_name, 1, 31 - nchar(suffix_text))
                candidate <- paste0(candidate, suffix_text)
                suffix <- suffix + 1
            }

            used_names <<- c(used_names, candidate)
            candidate
        }, character(1))
    }

    combined_long <- map_dfr(names(method_results), function(method) {
        imap_dfr(method_results[[method]]$results, function(result_tbl, comparison_slug) {
            result_tbl %>%
                mutate(comparison_slug = comparison_slug)
        })
    })

    comparison_slugs <- combined_long %>%
        distinct(comparison_slug) %>%
        pull(comparison_slug) %>%
        as.character()

    comparison_tables <- map(comparison_slugs, function(comparison_slug) {
        comparison_long <- combined_long %>%
            filter(comparison_slug == .env$comparison_slug)

        comparison_long %>%
            select(
                Name,
                Coordinate,
                analysis_method,
                analysis_method_label,
                raw_p_value,
                raw_p_lt_alpha,
                adjusted_p_value,
                n_p_adjust_hypotheses,
                fdr_lt_0_20,
                fdr_lt_0_25,
                effect_estimate_log2,
                effect_se_log2,
                fold_change_ratio,
                fold_change_status,
                fold_change_note,
                test_status,
                control_n,
                treatment_n,
                low_signal_flag
            ) %>%
            pivot_wider(
                names_from = analysis_method,
                values_from = c(
                    analysis_method_label,
                    raw_p_value,
                    raw_p_lt_alpha,
                    adjusted_p_value,
                    n_p_adjust_hypotheses,
                    fdr_lt_0_20,
                    fdr_lt_0_25,
                    effect_estimate_log2,
                    effect_se_log2,
                    fold_change_ratio,
                    fold_change_status,
                    fold_change_note,
                    test_status,
                    control_n,
                    treatment_n,
                    low_signal_flag
                ),
                names_glue = "{analysis_method}_{.value}"
            ) %>%
            arrange(Name, Coordinate)
    })
    names(comparison_tables) <- comparison_slugs

    comparison_index <- map_dfr(comparison_slugs, function(comparison_slug) {
        comparison_long <- combined_long %>%
            filter(comparison_slug == .env$comparison_slug)

        comparison_long %>%
            group_by(comparison_slug, subgroup, control, treatment, analysis_method, analysis_method_label) %>%
            summarize(
                n_tested = sum(test_status == "tested", na.rm = TRUE),
                n_raw_p_below_0_05 = sum(raw_p_value < 0.05, na.rm = TRUE),
                n_raw_p_lt_alpha = sum(raw_p_lt_alpha, na.rm = TRUE),
                n_fdr_lt_0_20 = sum(fdr_lt_0_20, na.rm = TRUE),
                n_fdr_lt_0_25 = sum(fdr_lt_0_25, na.rm = TRUE),
                best_raw_p_value = suppressWarnings(min(raw_p_value, na.rm = TRUE)),
                .groups = "drop"
            ) %>%
            mutate(best_raw_p_value = if_else(is.infinite(best_raw_p_value), NA_real_, best_raw_p_value))
    })

    comparison_workbook_path <- file.path(output_dir, "comparison_workbook.xlsx")
    workbook_sheets <- c(
        list(summary = comparison_index),
        stats::setNames(comparison_tables, build_workbook_sheet_names(comparison_slugs))
    )
    writexl::write_xlsx(workbook_sheets, comparison_workbook_path)

    invisible(list(
        comparison_workbook = comparison_workbook_path
    ))
}

#' Write a collaborator-facing method overview for dual-method runs
#'
#' @param analysis_methods Character vector of method slugs.
#' @param output_dir Root inferential output directory.
#'
#' @return Absolute path to the written markdown file.
write_inferential_method_overview <- function(analysis_methods, output_dir) {
    overview_path <- file.path(output_dir, "methods_overview.md")

    lines <- c(
        "# Inferential Methods Overview",
        "",
        paste(
            "Each configured replicate-aware method gets its own root workbook",
            "named `<method>_results.xlsx`.",
            "`comparison_workbook.xlsx` places all configured methods side by side."
        ),
        "Open files in this order:",
        "- `comparison_workbook.xlsx` for the side-by-side analyte table across methods.",
        "- `comparisons/<comparison_slug>/tables/<method>_results.tsv` for one method's full comparison table.",
        "- `comparisons/<comparison_slug>/waterfall_plots/<method>/<method>_waterfall*.png` for one method's visual ranking.",
        "- `comparisons/<comparison_slug>/barplots/<method>/all_tested/<method>_barplot_all_tested_page_<n>.png` for faceted fold-change bars over all tested analytes.",
        "- `comparisons/<comparison_slug>/barplots/<method>/significant_hits/<threshold>/<method>_barplot_<threshold>_page_<n>.png` for faceted fold-change bars over one significant hit set.",
        "",
        "Each comparison directory separates tables from plots, and plot folders are separated by method:",
        "- `tables/<method>_results.tsv` is the full results table for that method and comparison.",
        "- `waterfall_plots/<method>/<method>_waterfall.png` includes all tested analytes for that method.",
        "- `waterfall_plots/<method>/<method>_waterfall_raw_p_lt_alpha.png` includes analytes with unadjusted `p < alpha`.",
        "- `waterfall_plots/<method>/<method>_waterfall_fdr_lt_0_20.png` includes analytes with BH-adjusted `p < 0.20`.",
        "- `waterfall_plots/<method>/<method>_waterfall_fdr_lt_0_25.png` includes analytes with BH-adjusted `p < 0.25`.",
        "- `barplots/<method>/all_tested/<method>_barplot_all_tested_page_<n>.png` includes every tested analyte with a finite, ratio-scale-interpretable fold change and does not add threshold annotations.",
        "- `barplots/<method>/significant_hits/raw_p_lt_alpha/<method>_barplot_raw_p_lt_alpha_page_<n>.png`, `barplots/<method>/significant_hits/fdr_lt_0_20/<method>_barplot_fdr_lt_0_20_page_<n>.png`, and `barplots/<method>/significant_hits/fdr_lt_0_25/<method>_barplot_fdr_lt_0_25_page_<n>.png` use the same threshold sets as the matching waterfall files and add bracketed `*` annotations to treatment bars.",
        "- Barplot page captions define the star threshold, for example `* = FDR < 0.20`, next to the page `N`.",
        "- Barplot analyte titles use fixed-size title strips so short labels and wrapped labels render at the same font size.",
        "- Barplot y-axis limits are fixed across every page within a barplot set and use the smallest rounded range that includes all bars and bracket annotations in that set.",
        "",
        "Important:",
        "- Fold-change values are method-specific and should not be compared across methods as if they were the same quantity.",
        "- `raw_log2_lm` reports effects on the raw `log2(signal)` scale: for one analyte at a time it uses one averaged raw-signal value per biological sample, keeps only finite positive signals, fits `lm(log2(signal) ~ treatment)`, and interprets the treatment coefficient as treated minus control on the log2 raw-signal scale.",
        "- Waterfall whiskers, when shown, are +/- 1 standard error for the plotted effect estimate.",
        "- `normalized_t_test` uses raw-data normalized replicate values and reports ratio-scale fold change only when both group means are positive.",
        "- If ratio-scale fold change is not interpretable, the results table marks that row with `fold_change_status` and `fold_change_note`.",
        "- BH correction uses only the analytes that produced a valid p-value within that comparison family.",
        ""
    )

    for (method in analysis_methods) {
        method_spec <- get_inferential_method_spec(method)
        lines <- c(
            lines,
            sprintf("## %s (`%s`)", method_spec$label, method),
            "",
            "Statistical methods:",
            sprintf("- %s", method_spec$statistical_methods),
            "",
            "So what this means:",
            sprintf("- %s", method_spec$so_what),
            "",
            "Statistical strengths:",
            sprintf("- %s", method_spec$pros),
            "",
            "Statistical limitations:",
            sprintf("- %s", method_spec$cons),
            ""
        )
    }

    writeLines(lines, overview_path)
    overview_path
}

clean_multi_method_output_dir <- function(output_dir) {
    stale_paths <- c(
        file.path(output_dir, "methods"),
        file.path(output_dir, "method_comparison"),
        file.path(output_dir, "method_index.tsv"),
        file.path(output_dir, "comparison_index.tsv"),
        file.path(output_dir, "all_methods_long.tsv")
    )

    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))
}

#' Write one or more inferential methods to disk
#'
#' @param inferential_method_results Output of `run_replicate_analysis_methods()`.
#' @param output_dir Root inferential output directory.
#'
#' @return Invisibly returns written output paths.
write_multi_method_inferential_outputs <- function(inferential_method_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    clean_multi_method_output_dir(output_dir)

    primary_method <- inferential_method_results$primary_method
    method_results <- inferential_method_results$methods
    analysis_methods <- names(method_results)

    primary_paths <- write_inferential_outputs(
        inferential_results = method_results[[primary_method]],
        output_dir = output_dir
    )
    method_workbooks <- setNames(
        map_chr(analysis_methods, function(method_name) {
            if (identical(method_name, primary_method)) {
                return(primary_paths$workbook)
            }

            workbook_run_index <- method_results[[method_name]]$run_index %>%
                mutate(
                    result_path = NA_character_,
                    full_waterfall_path = NA_character_,
                    raw_p_alpha_waterfall_path = NA_character_,
                    fdr_0_20_waterfall_path = NA_character_,
                    fdr_0_25_waterfall_path = NA_character_
                )

            write_inferential_workbook(
                inferential_results = method_results[[method_name]],
                output_dir = output_dir,
                run_index = workbook_run_index
            )
        }),
        analysis_methods
    )
    write_method_specific_waterfalls(
        method_results = method_results,
        output_dir = output_dir,
        primary_method = primary_method
    )

    comparison_paths <- write_method_comparison_outputs(
        method_results = method_results,
        output_dir = output_dir
    )
    overview_path <- write_inferential_method_overview(
        analysis_methods = analysis_methods,
        output_dir = output_dir
    )

    invisible(list(
        primary = primary_paths,
        workbooks = method_workbooks,
        comparison = comparison_paths,
        overview = overview_path
    ))
}

clean_replicate_shortlist_comparison_dir <- function(output_dir) {
    stale_paths <- file.path(output_dir, c(
        "shortlist.tsv",
        "shortlist_summary.tsv",
        "shortlist_waterfall.png",
        "full_results.tsv",
        "comparison_results.tsv",
        "bargraph_index.tsv",
        "bargraphs"
    ))

    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))

    invisible(output_dir)
}

clean_replicate_shortlist_comparison_root <- function(output_dir) {
    stale_paths <- file.path(output_dir, c(
        "comparison_results.tsv",
        "shortlist_index.tsv",
        "method_index.tsv",
        "shortlist.tsv",
        "shortlist_summary.tsv",
        "shortlist_waterfall.png",
        "full_results.tsv",
        "bargraph_index.tsv",
        "bargraphs",
        "raw_p_lt_alpha",
        "fdr_lt_0_20",
        "fdr_lt_0_25",
        validate_analysis_methods(c("raw_log2_lm", "normalized_t_test"))
    ))

    invisible(walk(stale_paths, function(path) {
        if (file.exists(path) || dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        }
    }))

    invisible(output_dir)
}
