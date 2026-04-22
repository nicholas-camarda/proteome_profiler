# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

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
#' Exploratory analyses do not have a subgroup label, but their
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
            control_log2_se <- safe_se_or_na(dat$log2_signal[dat$treatment == control_label])
            treatment_log2_se <- safe_se_or_na(dat$log2_signal[dat$treatment == treatment_label])

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
                    control_fold_change_ymin = NA_real_,
                    control_fold_change_ymax = NA_real_,
                    treatment_fold_change_ymin = NA_real_,
                    treatment_fold_change_ymax = NA_real_,
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
            fold_change_ratio <- if_else(is.finite(effect_estimate_log2), 2^effect_estimate_log2, NA_real_)

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
                fold_change_ratio = fold_change_ratio,
                control_fold_change_ymin = if_else(is.finite(control_log2_se), 2^(-control_log2_se), NA_real_),
                control_fold_change_ymax = if_else(is.finite(control_log2_se), 2^control_log2_se, NA_real_),
                treatment_fold_change_ymin = if_else(
                    is.finite(effect_estimate_log2) & is.finite(treatment_log2_se),
                    2^(effect_estimate_log2 - treatment_log2_se),
                    NA_real_
                ),
                treatment_fold_change_ymax = if_else(
                    is.finite(effect_estimate_log2) & is.finite(treatment_log2_se),
                    2^(effect_estimate_log2 + treatment_log2_se),
                    NA_real_
                ),
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
            control_normalized_se <- safe_se_or_na(dat$normalized_signal[dat$treatment == control_label & valid_rows])
            treatment_normalized_se <- safe_se_or_na(dat$normalized_signal[dat$treatment == treatment_label & valid_rows])

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
                    control_fold_change_ymin = NA_real_,
                    control_fold_change_ymax = NA_real_,
                    treatment_fold_change_ymin = NA_real_,
                    treatment_fold_change_ymax = NA_real_,
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
                    control_fold_change_ymin = NA_real_,
                    control_fold_change_ymax = NA_real_,
                    treatment_fold_change_ymin = NA_real_,
                    treatment_fold_change_ymax = NA_real_,
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
            control_fold_change_se <- if (
                is.finite(control_normalized_se) &&
                is.finite(control_mean_normalized) &&
                control_mean_normalized > 0
            ) {
                control_normalized_se / control_mean_normalized
            } else {
                NA_real_
            }
            treatment_fold_change_se <- if (
                is.finite(treatment_normalized_se) &&
                is.finite(control_mean_normalized) &&
                control_mean_normalized > 0
            ) {
                treatment_normalized_se / control_mean_normalized
            } else {
                NA_real_
            }

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
                control_fold_change_ymin = if_else(
                    is.finite(control_fold_change_se),
                    pmax(1 - control_fold_change_se, 0),
                    NA_real_
                ),
                control_fold_change_ymax = if_else(
                    is.finite(control_fold_change_se),
                    1 + control_fold_change_se,
                    NA_real_
                ),
                treatment_fold_change_ymin = if_else(
                    is.finite(fold_change_ratio) & is.finite(treatment_fold_change_se),
                    pmax(fold_change_ratio - treatment_fold_change_se, 0),
                    NA_real_
                ),
                treatment_fold_change_ymax = if_else(
                    is.finite(fold_change_ratio) & is.finite(treatment_fold_change_se),
                    fold_change_ratio + treatment_fold_change_se,
                    NA_real_
                ),
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
