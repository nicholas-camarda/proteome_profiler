#' Return the package set required by the active analysis scripts
#'
#' Keeps the install script and runtime package loading in sync so the package
#' list only needs to be updated in one place.
#'
#' @param include_parallel Logical. When `TRUE`, include the packages used by
#'   the parallelized `scripts/main.R` workflow.
#'
#' @return Character vector of package names.
required_analysis_packages <- function(include_parallel = TRUE) {
    packages <- c(
        "tidyverse",
        "readxl",
        "GetoptLong",
        "RColorBrewer",
        "latex2exp",
        "patchwork",
        "ggprism",
        "ggh4x",
        "ggforce",
        "readr"
    )

    if (include_parallel) {
        packages <- c(packages, "progressr", "future", "furrr")
    }

    unique(packages)
}

#' Load the packages required by the active analysis scripts
#'
#' This is a thin wrapper around `require()` that reuses the canonical package
#' list from `required_analysis_packages()`.
#'
#' @param include_parallel Logical. When `TRUE`, also load the packages needed
#'   by the parallel main workflow.
#'
#' @return Invisibly returns the logical results from `require()` calls.
load_analysis_packages <- function(include_parallel = TRUE) {
    invisible(lapply(
        required_analysis_packages(include_parallel = include_parallel),
        require,
        character.only = TRUE
    ))
}

#' Configure progressr handlers for long-running analysis scripts
#'
#' The main workflow can take a while when writing many thresholded plots. This
#' standardizes the progress display so runs are easier to monitor.
#'
#' @return `NULL`, called for side effects.
configure_progress_handlers <- function() {
    handlers(handler_progress(
        format = ":spin [:bar] :current/:total :percent in :elapsed ETA: :eta (:message) ",
        width = 120,
        complete = "+"
    ))
}

#' Build the canonical plotting style for treatment groups
#'
#' Centralizes group ordering and palette assignment so the threshold diagnostic
#' plots and the main analysis plots use consistent labels and colors.
#'
#' @param group_levels Character vector giving the desired treatment-group order.
#' @param scheme One of `"main"` or `"threshold"`, selecting the palette tuned
#'   for the main analysis or the threshold-diagnostic plot.
#'
#' @return Named list with `group_levels`, `fill`, and `outline` entries.
build_group_style <- function(group_levels, scheme = c("main", "threshold")) {
    scheme <- match.arg(scheme)
    group_levels <- factor(group_levels, levels = group_levels)

    if (scheme == "main") {
        fill_values <- c("#ffffff", "#737171", "#CCD8E8", "#E5938A")
        outline_values <- c("black", "black", "#22456F", "#B53530")
    } else {
        fill_values <- c("#ffffff", "#737171", "#9BB3D3", "#B53530")
        outline_values <- c("black", "black", "black", "black")
    }

    # Truncate the palette to the active groups so shorter analyses can reuse
    # the same helper without carrying unused colors around.
    list(
        group_levels = group_levels,
        fill = stats::setNames(fill_values[seq_len(nlevels(group_levels))], group_levels),
        outline = stats::setNames(outline_values[seq_len(nlevels(group_levels))], group_levels)
    )
}
