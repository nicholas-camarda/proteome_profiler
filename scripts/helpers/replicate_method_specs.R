# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

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
                "Bars show log2(mean treatment-arm normalized_signal / mean control-arm normalized_signal);",
                "normalized_signal is each sample's averaged raw analyte signal divided by that sample's reference-spot denominator."
            ),
            statistical_methods = c(
                "Recomputes one `normalized_signal` value for each sample and analyte from the raw sample sheets.",
                "`normalized_signal` equals averaged duplicate raw analyte signal divided by that sample's reference-spot denominator; the denominator uses preferred Reference Spots pairs `A1,2` and `J1,2` when present, otherwise available Reference Spots rows from that sample.",
                "The denominator, reference rows, and raw reference signals used for each sample are reported in `input_qc/reference_spot_qc.tsv`.",
                "`low_signal_flag` is separate from reference-spot normalization; it is computed from the configured raw-signal low-signal threshold and is reported for every method.",
                "Runs a two-sided equal-variance two-sample t-test on those per-sample `normalized_signal` values within each subgroup.",
                "Tests equality of mean `normalized_signal` between treatment and control; when both group means are positive, reports fold change as `mean(treatment-arm normalized_signal values) / mean(control-arm normalized_signal values)`."
            ),
            so_what = c(
                "Use this when you want each analyte signal to be divided by the sample's reference-spot intensity first, and then you want to compare treatment and control on that normalized scale.",
                "In practice, this is the better choice if you want the question to be 'do the groups differ after reference-spot normalization?' rather than 'do the raw signals differ on a log2 scale?'"
            ),
            pros = c(
                "Targets differences in mean `normalized_signal` on the same scale used by the normalized replicate workflow.",
                "The test is simple and transparent: per analyte, it compares the two group means of the per-sample `normalized_signal` values directly.",
                "Can still produce a p-value when `normalized_signal` values are finite even if a ratio-scale fold change is not interpretable because one mean is nonpositive."
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

#' Variance helper that returns `NA` when fewer than two finite values exist
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar.
safe_var_or_na <- function(x) {
    finite_x <- as.numeric(x)[is.finite(as.numeric(x))]
    if (length(finite_x) < 2) {
        return(NA_real_)
    }

    stats::var(finite_x)
}

#' Standard-error helper that returns `NA` when fewer than two values exist
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar.
safe_se_or_na <- function(x) {
    finite_x <- as.numeric(x)[is.finite(as.numeric(x))]
    if (length(finite_x) < 2) {
        return(NA_real_)
    }

    stats::sd(finite_x) / sqrt(length(finite_x))
}

#' Approximate SE for a log2 ratio from two arithmetic means
#'
#' This is used for normalized t-test fold-change whiskers. It applies a delta
#' method approximation on the log2 ratio scale using the pooled within-arm
#' variance from the normalized values.
#'
#' @param control_values Numeric control-arm values.
#' @param treatment_values Numeric treatment-arm values.
#' @param control_mean Numeric control mean.
#' @param treatment_mean Numeric treatment mean.
#'
#' @return Numeric standard error on the log2 ratio scale, or `NA`.
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
