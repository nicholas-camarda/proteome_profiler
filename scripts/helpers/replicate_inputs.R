# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

#' Collapse duplicate cached normalized values from one technical-replicate pair
#'
#' @param x Numeric-like vector from collaborator workbook normalized columns.
#'
#' @return First finite nonzero normalized value when present, otherwise the
#'   finite mean or `NA`.
collapse_workbook_normalized_signal <- function(x) {
    numeric_x <- suppressWarnings(as.numeric(x))
    finite_nonzero <- numeric_x[is.finite(numeric_x) & numeric_x != 0]
    if (length(finite_nonzero) > 0) {
        return(finite_nonzero[[1]])
    }

    safe_mean_or_na(numeric_x)
}

#' Recompute workbook-style normalized signal from raw paired spots
#'
#' @param sample_pair_tbl Sample-level paired-spot tibble with raw `signal`,
#'   `Name`, `Coordinate`, and `Sname` columns.
#'
#' @return The input tibble with a `normalized_signal` column.
compute_raw_sheet_normalized_signal <- function(sample_pair_tbl) {
    pair_tbl <- sample_pair_tbl %>%
        mutate(
            compact_coordinate = compact_coordinate_text(Coordinate),
            compact_sname = compact_coordinate_text(str_replace_all(Sname, "([A-Z]{1,2})0([0-9])", "\\1\\2")),
            is_reference_spot = str_to_lower(str_trim(Name)) == "reference spots"
        )

    preferred_reference_rows <- pair_tbl %>%
        filter(is_reference_spot, compact_sname %in% c("A1,2", "J1,2"))

    reference_source <- "preferred_raw_spot_A1_2_J1_2"
    reference_rows <- if (nrow(preferred_reference_rows) > 0) {
        preferred_reference_rows
    } else {
        fallback_reference_rows <- pair_tbl %>%
            filter(is_reference_spot) %>%
            filter(str_detect(compact_coordinate, "^[A-Z]{1,2}1,2$"))

        if (nrow(fallback_reference_rows) > 0) {
            reference_source <- "protocol_reference_spots_coordinate_1_2"
            fallback_reference_rows
        } else {
            reference_source <- "all_protocol_reference_spots"
            pair_tbl %>%
                filter(is_reference_spot)
        }
    }

    if (nrow(reference_rows) == 0) {
        stop("Workbook-style normalized t-test requires at least one 'Reference Spots' pair in the protocol table.")
    }

    normalization_denominator <- safe_mean_or_na(reference_rows$signal)
    if (!is.finite(normalization_denominator) || normalization_denominator <= 0) {
        stop("Workbook-style normalized t-test could not compute a finite positive reference normalization denominator from raw signals.")
    }

    reference_snames <- sort(unique(reference_rows$compact_sname))
    reference_status <- case_when(
        identical(reference_source, "preferred_raw_spot_A1_2_J1_2") && all(c("A1,2", "J1,2") %in% reference_snames) ~ "preferred_complete",
        identical(reference_source, "preferred_raw_spot_A1_2_J1_2") ~ "preferred_partial",
        TRUE ~ "protocol_fallback"
    )
    reference_signal_trace <- reference_rows %>%
        transmute(
            reference_trace = str_c(
                compact_sname,
                " [",
                compact_coordinate,
                "]=",
                signif(signal, 8)
            )
        ) %>%
        pull(reference_trace) %>%
        str_c(collapse = "; ")

    pair_tbl %>%
        mutate(
            normalized_signal = signal / normalization_denominator,
            normalization_denominator = normalization_denominator,
            normalization_reference_source = reference_source,
            normalization_reference_status = reference_status,
            normalization_reference_n = nrow(reference_rows),
            normalization_reference_snames = str_c(reference_snames, collapse = "|"),
            normalization_reference_names = str_c(sort(unique(reference_rows$Name)), collapse = "|"),
            normalization_reference_coordinates = str_c(sort(unique(reference_rows$compact_coordinate)), collapse = "|"),
            normalization_reference_signals = str_c(signif(reference_rows$signal, 8), collapse = "|"),
            normalization_reference_trace = reference_signal_trace
        ) %>%
        select(-compact_coordinate, -compact_sname, -is_reference_spot)
}

#' Normalize paired spot coordinates for protocol joins
#'
#' @param x Character vector of coordinate text.
#'
#' @return Compact coordinate key compatible with protocol coordinates.
normalize_pair_coordinate_key <- function(x) {
    compact_coordinate_text(str_replace_all(x, "([A-Z]{1,2})0([0-9])", "\\1\\2"))
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

    resolved_manifest <- manifest_with_paths %>%
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

    duplicate_workbook_sheets <- resolved_manifest %>%
        count(resolved_workbook_path, resolved_sheet_name, name = "n_rows") %>%
        filter(n_rows > 1)
    if (nrow(duplicate_workbook_sheets) > 0) {
        duplicate_descriptions <- duplicate_workbook_sheets %>%
            mutate(description = sprintf(
                "%s [%s]",
                basename(resolved_workbook_path),
                resolved_sheet_name
            )) %>%
            pull(description)
        stop(sprintf(
            paste(
                "Sample manifest maps multiple sample rows to the same workbook sheet.",
                "Each sample must have a unique workbook_path/sheet_name pair.",
                "Duplicate workbook sheet(s): %s"
            ),
            paste(duplicate_descriptions, collapse = ", ")
        ), call. = FALSE)
    }

    resolved_manifest
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
#' @param text_vector Character vector like `"A1, A2"`, `"A1,2"`, or
#'   `"A1,A2"`.
#'
#' @return Character vector like `"A1,2"`.
compact_coordinate_text <- function(text_vector) {
    output <- sapply(as.character(text_vector), function(text_value) {
        if (is.na(text_value)) {
            return(NA_character_)
        }
        tokens <- str_split(str_trim(text_value), "\\s*,\\s*", simplify = FALSE)[[1]]
        tokens <- tokens[tokens != ""]
        if (length(tokens) == 0) {
            return("")
        }

        current_letters <- NA_character_
        letter_order <- character()
        grouped_numbers <- list()
        for (token in tokens) {
            token_letters <- str_extract(token, "^[A-Za-z]+")
            token_numbers <- str_extract_all(token, "\\d+")[[1]]
            if (!is.na(token_letters) && nzchar(token_letters)) {
                current_letters <- toupper(token_letters)
            }
            if (is.na(current_letters) || length(token_numbers) == 0) {
                next
            }
            if (!current_letters %in% letter_order) {
                letter_order <- c(letter_order, current_letters)
                grouped_numbers[[current_letters]] <- character()
            }
            grouped_numbers[[current_letters]] <- c(grouped_numbers[[current_letters]], token_numbers)
        }

        paste(
            vapply(letter_order, function(letters) {
                paste0(letters, paste(unique(grouped_numbers[[letters]]), collapse = ","))
            }, character(1)),
            collapse = ", "
        )
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

#' Read optional technical-artifact exclusions for exploratory or sample-level inputs
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

#' Summarize reference-spot normalization provenance by sample
#'
#' @param sample_data Tibble returned by `build_sample_level_dataset()`.
#'
#' @return Tibble with one row per sample and the denominator/source used to
#'   compute `normalized_signal`.
summarize_reference_spot_qc <- function(sample_data) {
    collapse_unique_text <- function(x) {
        values <- x[!is.na(x) & as.character(x) != ""]
        values <- sort(unique(as.character(values)))
        if (length(values) == 0) {
            return(NA_character_)
        }
        str_c(values, collapse = " | ")
    }

    first_finite_or_na <- function(x) {
        finite_x <- x[is.finite(x)]
        if (length(finite_x) == 0) {
            return(NA_real_)
        }
        finite_x[[1]]
    }

    reference_columns <- c(
        "normalization_denominator",
        "normalization_reference_source",
        "normalization_reference_status",
        "normalization_reference_n",
        "normalization_reference_snames",
        "normalization_reference_names",
        "normalization_reference_coordinates",
        "normalization_reference_signals",
        "normalization_reference_trace"
    )

    missing_reference_columns <- setdiff(reference_columns, names(sample_data))
    if (length(missing_reference_columns) > 0) {
        stop(sprintf(
            "Sample-level data is missing reference normalization QC column(s): %s",
            paste(missing_reference_columns, collapse = ", ")
        ), call. = FALSE)
    }

    sample_qc_data <- sample_data
    if (!"workbook_path" %in% names(sample_qc_data)) {
        sample_qc_data$workbook_path <- NA_character_
    }
    if (!"sheet_name" %in% names(sample_qc_data)) {
        sample_qc_data$sheet_name <- NA_character_
    }
    if (!"subgroup" %in% names(sample_qc_data)) {
        sample_qc_data$subgroup <- NA_character_
    }
    if (!"sex" %in% names(sample_qc_data)) {
        sample_qc_data$sex <- NA_character_
    }

    group_columns <- c("sample_id", "subgroup", "sex", "treatment", "workbook_path", "sheet_name")

    sample_qc_data %>%
        group_by(across(all_of(group_columns))) %>%
        summarize(
            normalization_denominator = first_finite_or_na(normalization_denominator),
            n_distinct_denominators = n_distinct(normalization_denominator[is.finite(normalization_denominator)]),
            normalization_reference_source = collapse_unique_text(normalization_reference_source),
            normalization_reference_status = collapse_unique_text(normalization_reference_status),
            normalization_reference_n = first_finite_or_na(as.numeric(normalization_reference_n)),
            normalization_reference_snames = collapse_unique_text(normalization_reference_snames),
            normalization_reference_names = collapse_unique_text(normalization_reference_names),
            normalization_reference_coordinates = collapse_unique_text(normalization_reference_coordinates),
            normalization_reference_signals = collapse_unique_text(normalization_reference_signals),
            normalization_reference_trace = collapse_unique_text(normalization_reference_trace),
            n_analyte_rows = n(),
            .groups = "drop"
        ) %>%
        mutate(
            reference_qc_issue = case_when(
                n_distinct_denominators > 1 ~ "inconsistent_denominator_within_sample",
                !is.finite(normalization_denominator) ~ "missing_or_nonfinite_denominator",
                normalization_denominator <= 0 ~ "nonpositive_denominator",
                TRUE ~ "ok"
            ),
            reference_qc_note = case_when(
                str_detect(normalization_reference_status, "protocol_fallback") ~ "used_protocol_reference_spot_rows",
                str_detect(normalization_reference_status, "preferred_partial") ~ "used_partial_preferred_reference_spot_rows",
                TRUE ~ "used_complete_preferred_reference_spot_rows"
            )
        ) %>%
        arrange(reference_qc_issue != "ok", reference_qc_issue, sample_id)
}

#' Summarize reference-spot QC counts for a run
#'
#' @param reference_qc Tibble returned by `summarize_reference_spot_qc()`.
#'
#' @return One-row tibble with sample counts by reference-spot status.
summarize_reference_spot_qc_counts <- function(reference_qc) {
    finite_denominators <- reference_qc$normalization_denominator[
        is.finite(reference_qc$normalization_denominator)
    ]

    tibble(
        n_samples = nrow(reference_qc),
        n_reference_qc_ok = sum(reference_qc$reference_qc_issue == "ok", na.rm = TRUE),
        n_preferred_complete = sum(reference_qc$normalization_reference_status == "preferred_complete", na.rm = TRUE),
        n_preferred_partial = sum(reference_qc$normalization_reference_status == "preferred_partial", na.rm = TRUE),
        n_protocol_fallback = sum(reference_qc$normalization_reference_status == "protocol_fallback", na.rm = TRUE),
        n_reference_qc_issues = sum(reference_qc$reference_qc_issue != "ok", na.rm = TRUE),
        min_reference_denominator = if (length(finite_denominators) > 0) min(finite_denominators) else NA_real_,
        median_reference_denominator = if (length(finite_denominators) > 0) stats::median(finite_denominators) else NA_real_,
        max_reference_denominator = if (length(finite_denominators) > 0) max(finite_denominators) else NA_real_
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

    reference_spot_qc_tbl <- summarize_reference_spot_qc(sample_data)
    reference_spot_summary_tbl <- summarize_reference_spot_qc_counts(reference_spot_qc_tbl)

    list(
        summary = summary_tbl,
        sample_summary = sample_summary_tbl,
        analyte_summary = analyte_summary_tbl,
        issue_rows = issue_rows_tbl,
        reference_spot_summary = reference_spot_summary_tbl,
        reference_spot_qc = reference_spot_qc_tbl
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
    reference_spot_summary_path <- file.path(output_dir, "reference_spot_summary.tsv")
    reference_spot_qc_path <- file.path(output_dir, "reference_spot_qc.tsv")

    write_tsv(cleanliness$summary, summary_path)
    write_tsv(cleanliness$sample_summary, sample_summary_path)
    write_tsv(cleanliness$analyte_summary, analyte_summary_path)
    write_tsv(cleanliness$issue_rows, issue_rows_path)
    write_tsv(cleanliness$reference_spot_summary, reference_spot_summary_path)
    write_tsv(cleanliness$reference_spot_qc, reference_spot_qc_path)

    invisible(list(
        summary = summary_path,
        sample_summary = sample_summary_path,
        analyte_summary = analyte_summary_path,
        issue_rows = issue_rows_path,
        reference_spot_summary = reference_spot_summary_path,
        reference_spot_qc = reference_spot_qc_path
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
