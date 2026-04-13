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

    manifest %>%
        mutate(resolved_workbook_path = resolved_paths)
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
read_licor_sample_workbook <- function(workbook_path, analyte_info, sample_metadata, remove_dat = NULL, treatment_var = "treatment") {
    sample_id_value <- sample_metadata$sample_id[[1]]
    treatment_value <- sample_metadata[[treatment_var]][[1]]

    ds_tmp <- read_excel(workbook_path, skip = 3)

    res_temp_pre <- ds_tmp %>%
        rename(Sname = Name) %>%
        dplyr::select(Sname, Signal) %>%
        filter(!startsWith(Sname, "B")) %>%
        mutate(sname_grouping = rep(seq_len(nrow(.) / 2), each = 2))

    expected_snames <- sprintf("S%03d", seq_len(nrow(res_temp_pre)))
    if (nrow(res_temp_pre) %% 2 != 0 || !identical(res_temp_pre$Sname, expected_snames)) {
        stop(sprintf(
            "Unexpected LI-COR spot ordering in %s. Expected sequential S001... labels after removing background.",
            workbook_path
        ))
    }

    if (is.null(remove_dat)) {
        remove_dat <- tibble(Sname = character(), sample_id = character(), group = character())
    }

    sample_remove_dat <- remove_dat %>%
        filter((is.na(sample_id) | .data$sample_id == .env$sample_id_value)) %>%
        filter((is.na(group) | .data$group == .env$treatment_value)) %>%
        distinct(Sname)

    sname_map <- res_temp_pre %>%
        group_by(sname_grouping) %>%
        reframe(init_sname = Sname, Sname = str_c(Sname, collapse = ", "))

    res_temp <- res_temp_pre %>%
        anti_join(sample_remove_dat, by = "Sname")

    sample_metadata_clean <- sample_metadata %>%
        dplyr::select(-any_of("resolved_workbook_path"))

    result_tbl <- res_temp %>%
        group_by(sname_grouping) %>%
        reframe(
            signal = mean(as.numeric(Signal), na.rm = TRUE),
            Sname = str_c(Sname, collapse = ", ")
        ) %>%
        left_join(sname_map, by = c("Sname" = "init_sname", "sname_grouping")) %>%
        rename(full_sname = Sname.y) %>%
        dplyr::select(-Sname) %>%
        rename(Sname = full_sname) %>%
        left_join(analyte_info, by = "sname_grouping")

    bind_cols(
        result_tbl,
        sample_metadata_clean[rep(1, nrow(result_tbl)), , drop = FALSE]
    ) %>%
        mutate(
            sample_id = .env$sample_id_value,
            treatment = .env$treatment_value,
            workbook_path = workbook_path
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
            list(resolved_workbook_path, sample_id, workbook_path),
            function(resolved_workbook_path, sample_id_value, workbook_path) {
                sample_metadata <- manifest %>%
                    filter(.data$sample_id == .env$sample_id_value) %>%
                    slice_head(n = 1)
                read_licor_sample_workbook(
                    workbook_path = resolved_workbook_path,
                    analyte_info = analyte_info,
                    sample_metadata = sample_metadata,
                    remove_dat = remove_dat,
                    treatment_var = treatment_var
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
fit_one_comparison <- function(comparison_df, control_label, treatment_label, subgroup_label, subgroup_var, alpha, min_reps = 2, low_signal_threshold = NULL) {
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
                    log2_signal = if_else(signal > 0 & is.finite(signal), log2(signal), NA_real_)
                )

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
                    effect_estimate_log2 = NA_real_,
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
                effect_estimate_log2 = effect_estimate_log2,
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
run_within_stratum_differential_analysis <- function(sample_data, comparisons, subgroup_var, p_adjust_method = "BH", alpha = 0.05, min_reps = 2, low_signal_threshold = NULL) {
    validate_p_adjust_method(p_adjust_method)
    validate_requested_comparisons(sample_data, comparisons, subgroup_var)

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

                result_tbl <- fit_one_comparison(
                    comparison_df = comparison_df,
                    control_label = control_label,
                    treatment_label = treatment_label,
                    subgroup_label = subgroup_level,
                    subgroup_var = subgroup_var,
                    alpha = alpha,
                    min_reps = min_reps,
                    low_signal_threshold = low_signal_threshold
                ) %>%
                    mutate(
                        adjusted_p_value = if_else(
                            test_status == "tested",
                            p.adjust(raw_p_value, method = p_adjust_method),
                            NA_real_
                        ),
                        p_adjust_method = p_adjust_method,
                        significant = if_else(
                            !is.na(adjusted_p_value),
                            adjusted_p_value <= alpha,
                            FALSE
                        ),
                        comparison_label = sprintf("%s: %s vs %s", subgroup_level, control_label, treatment_label)
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
                    n_analytes = nrow(result_tbl),
                    n_significant = sum(result_tbl$significant, na.rm = TRUE),
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

#' Build a shortlist table from one inferential comparison
#'
#' @param results_tbl Inferential results table for one comparison.
#' @param selection_names Optional character vector of analyte names to force-include.
#' @param top_n Optional integer number of top tested analytes to retain when no
#'   explicit name list is supplied.
#'
#' @return Tibble of shortlisted analytes with a `shortlist_basis` column.
build_inferential_shortlist <- function(results_tbl, selection_names = NULL, top_n = NULL) {
    tested_tbl <- results_tbl %>%
        filter(test_status == "tested")

    if (!is.null(selection_names) && length(selection_names) > 0) {
        return(
            tested_tbl %>%
                filter(Name %in% selection_names) %>%
                arrange(match(Name, selection_names)) %>%
                mutate(shortlist_basis = "explicit_name_list")
        )
    }

    significant_tbl <- tested_tbl %>%
        filter(significant) %>%
        arrange(adjusted_p_value, desc(abs(effect_estimate_log2)), Name) %>%
        mutate(shortlist_basis = "significant")

    if (nrow(significant_tbl) > 0) {
        return(significant_tbl)
    }

    if (!is.null(top_n) && top_n > 0) {
        return(
            tested_tbl %>%
                arrange(adjusted_p_value, desc(abs(effect_estimate_log2)), Name) %>%
                slice_head(n = top_n) %>%
                mutate(shortlist_basis = "top_n_fallback")
        )
    }

    significant_tbl %>%
        mutate(shortlist_basis = character())
}

#' Plot inferential effect estimates as a waterfall chart
#'
#' @param results_tbl One per-comparison inferential results tibble.
#' @param title Character plot title.
#'
#' @return ggplot object.
plot_inferential_waterfall <- function(results_tbl, title) {
    plot_data <- results_tbl %>%
        filter(test_status == "tested") %>%
        arrange(desc(effect_estimate_log2)) %>%
        mutate(
            Name = factor(Name, levels = Name),
            fill_flag = if_else(significant, "significant", "not_significant")
        )

    ggplot(plot_data, aes(x = Name, y = effect_estimate_log2, fill = fill_flag)) +
        geom_col(show.legend = FALSE) +
        scale_fill_manual(values = c(significant = "#B53530", not_significant = "#9BB3D3")) +
        labs(
            title = str_wrap(title, width = 55),
            subtitle = str_wrap("Effect estimate is treatment minus control on the log2 signal scale", width = 70),
            caption = sprintf("N analytes tested = %s", nrow(plot_data))
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
write_inferential_outputs <- function(inferential_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    result_paths <- map_chr(names(inferential_results$results), function(comparison_slug) {
        comparison_dir <- file.path(output_dir, "comparisons", comparison_slug)
        dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

        result_tbl <- inferential_results$results[[comparison_slug]]
        result_path <- file.path(comparison_dir, "results.tsv")
        write_tsv(result_tbl, result_path)

        waterfall_path <- file.path(comparison_dir, "waterfall.png")
        waterfall_plot <- plot_inferential_waterfall(
            result_tbl,
            title = unique(result_tbl$comparison_label)
        )
        ggsave(
            filename = waterfall_path,
            plot = waterfall_plot,
            width = max(12, 0.4 * nrow(result_tbl)),
            height = 10
        )

        result_path
    })

    result_path_map <- tibble(
        comparison_slug = names(inferential_results$results),
        result_path = unname(result_paths)
    )

    run_index <- inferential_results$run_index %>%
        left_join(result_path_map, by = "comparison_slug")
    run_index_path <- file.path(output_dir, "run_index.tsv")
    write_tsv(run_index, run_index_path)

    workbook_sheets <- inferential_results$results
    workbook_sheets$run_index <- run_index
    workbook_path <- file.path(output_dir, "combined_results.xlsx")
    writexl::write_xlsx(workbook_sheets, workbook_path)

    invisible(list(
        run_index = run_index_path,
        workbook = workbook_path,
        results = result_paths
    ))
}

#' Write shortlist outputs derived from one inferential comparison
#'
#' @param shortlist_tbl Tibble returned by `build_inferential_shortlist()`.
#' @param comparison_tbl Full inferential results table for the selected comparison.
#' @param output_dir Directory where shortlist outputs should be written.
#'
#' @return Invisibly returns output paths written for the shortlist workflow.
write_inferential_shortlist_outputs <- function(shortlist_tbl, comparison_tbl, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    shortlist_path <- file.path(output_dir, "shortlist.tsv")
    write_tsv(shortlist_tbl, shortlist_path)

    full_results_path <- file.path(output_dir, "full_results.tsv")
    write_tsv(comparison_tbl, full_results_path)

    summary_tbl <- tibble(
        comparison_label = unique(comparison_tbl$comparison_label),
        n_tested = sum(comparison_tbl$test_status == "tested", na.rm = TRUE),
        n_significant = sum(comparison_tbl$significant, na.rm = TRUE),
        n_shortlisted = nrow(shortlist_tbl),
        shortlist_basis = if (nrow(shortlist_tbl) > 0) shortlist_tbl$shortlist_basis[[1]] else "empty"
    )
    summary_path <- file.path(output_dir, "shortlist_summary.tsv")
    write_tsv(summary_tbl, summary_path)

    waterfall_path <- NA_character_
    if (nrow(shortlist_tbl) > 0) {
        waterfall_path <- file.path(output_dir, "shortlist_waterfall.png")
        shortlist_plot <- plot_inferential_waterfall(
            shortlist_tbl,
            title = sprintf("Shortlist\n%s", unique(comparison_tbl$comparison_label))
        )
        ggsave(
            filename = waterfall_path,
            plot = shortlist_plot,
            width = max(12, 0.6 * nrow(shortlist_tbl)),
            height = 10
        )
    }

    invisible(list(
        shortlist = shortlist_path,
        full_results = full_results_path,
        summary = summary_path,
        waterfall = waterfall_path
    ))
}
