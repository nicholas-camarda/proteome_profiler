#' Find the repository root from the current working directory
#'
#' Walks upward until it finds the helper path that anchors this repository.
#' This keeps the active scripts runnable from the repo root or from a nested
#' subdirectory during development.
#'
#' @param start Character path to start searching from.
#'
#' @return Absolute path to the repository root.
find_repo_root <- function(start = getwd()) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)

    repeat {
        if (file.exists(file.path(current, "scripts", "helpers", "project_paths.R"))) {
            return(current)
        }

        parent <- dirname(current)
        if (identical(parent, current)) {
            break
        }
        current <- parent
    }

    stop("Could not determine the proteome_profiler repository root.")
}

repo_root <- find_repo_root()

if (!exists("proteome_profiler_config")) {
    stop(paste(
        "Expected `proteome_profiler_config` to already be defined.",
        "Source scripts/config/analysis_config.R before sourcing scripts/helpers/project_paths.R."
    ))
}

#' Build candidate filesystem locations for a project-relative path
#'
#' The repository uses a split layout where code, runtime outputs, and cloud
#' copies may live in different roots. This helper returns the search order used
#' throughout the project: repo checkout, runtime workspace, then cloud backup.
#'
#' @param path Project-relative path to resolve.
#'
#' @return Character vector of candidate absolute paths.
path_candidates <- function(path) {
    runtime_root <- path.expand(proteome_profiler_config$runtime_root)
    cloud_parent <- path.expand(proteome_profiler_config$cloud_parent)
    cloud_project_root <- file.path(cloud_parent, "proteome_profiler")

    unique(c(
        file.path(repo_root, path),
        file.path(runtime_root, path),
        file.path(cloud_project_root, path)
    ))
}

#' Resolve a project-relative path against the known storage roots
#'
#' Looks for an existing match across the repo checkout, runtime workspace, and
#' cloud copy. When `must_exist = FALSE`, it returns the canonical runtime path
#' where new outputs should be written.
#'
#' @param path Project-relative path to resolve.
#' @param must_exist Logical. When `TRUE`, error if the path cannot be found.
#'
#' @return Absolute normalized path.
resolve_project_path <- function(path, must_exist = TRUE) {
    candidates <- path_candidates(path)
    existing <- candidates[file.exists(candidates) | dir.exists(candidates)]

    if (length(existing) > 0) {
        return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
    }

    if (!must_exist) {
        # New outputs always materialize under the runtime tree, even if older
        # copies also exist in the repo or cloud backup locations.
        return(normalizePath(file.path(path.expand(proteome_profiler_config$runtime_root), path), winslash = "/", mustWork = FALSE))
    }

    stop(sprintf(
        "Could not resolve '%s'. Checked:\n%s",
        path,
        paste(sprintf("- %s", candidates), collapse = "\n")
    ))
}

#' Return one named analysis configuration
#'
#' The returned object is metadata only: where the raw data lives, who owns the
#' output subtree, which comparisons to run, and which thresholds to use.
#'
#' @param example_name Name of the entry in
#'   `proteome_profiler_config$examples`.
#'
#' @return Named list describing one analysis run.
get_analysis_config <- function(example_name = "vegfri_dox_cytokine_xl") {
    config <- proteome_profiler_config$examples[[example_name]]
    if (is.null(config)) {
        stop(sprintf("Unknown example config: %s", example_name))
    }
    config
}

#' Validate that an analysis configuration contains the required fields
#'
#' This catches incomplete config entries before the scripts start writing
#' output or looking for data in ambiguous locations.
#'
#' @param example_config Named list returned by `get_analysis_config()`.
#'
#' @return `NULL`, called for side effects.
validate_analysis_config <- function(example_config) {
    required_fields <- c("user", "analysis_slug", "data_dir", "info_fn", "group_levels")
    missing_fields <- required_fields[vapply(required_fields, function(field) {
        is.null(example_config[[field]]) || identical(example_config[[field]], "")
    }, logical(1))]

    if (length(missing_fields) > 0) {
        stop(sprintf(
            "Analysis config is missing required fields: %s",
            paste(missing_fields, collapse = ", ")
        ))
    }
}

#' Return the root output directory for one analysis configuration
#'
#' Every analysis writes into `output/plots/<user>/<analysis_slug>/` so outputs
#' remain grouped by analyst and by analysis.
#'
#' @param example_config Named list returned by `get_analysis_config()`.
#'
#' @return Absolute path to the analysis output root.
get_analysis_output_root <- function(example_config) {
    validate_analysis_config(example_config)

    # This keeps collaborator outputs separated even when they share the same
    # runtime workspace or cloud-backed project tree.
    resolve_project_path(
        file.path("output", "plots", example_config$user, example_config$analysis_slug),
        must_exist = FALSE
    )
}

#' Resolve the generated protocol workbook required by the R workflows
#'
#' The workbook is produced by `scripts/setup/extract_analyte_table.py`. If it
#' is missing, this helper raises an explicit error that points the user at the
#' required Python preprocessing step.
#'
#' @param example_config Named list returned by `get_analysis_config()`.
#'
#' @return Absolute path to the protocol workbook.
get_protocol_workbook_path <- function(example_config) {
    candidates <- path_candidates(example_config$info_fn)
    existing <- candidates[file.exists(candidates)]

    if (length(existing) > 0) {
        return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
    }

    protocol_pdf_hint <- file.path("protocols", basename(example_config$protocol_pdf))
    stop(sprintf(
        paste(
            "Required protocol workbook is missing: %s",
            "Run this first:",
            "python3 scripts/setup/extract_analyte_table.py --preset %s",
            "Expected PDF example: %s"
        ),
        example_config$info_fn,
        example_config$protocol_preset,
        protocol_pdf_hint
    ))
}
