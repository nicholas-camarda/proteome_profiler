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

config_uses_sample_manifest <- function(example_config) {
    !is.null(example_config$sample_manifest) &&
        !identical(example_config$sample_manifest, "") &&
        !is.na(example_config$sample_manifest)
}

#' Return the analysis mode requested by one config entry
#'
#' The repository now separates how inputs are declared from how they are
#' analyzed. Both exploratory and inferential analyses may be manifest-driven;
#' `mode` tells the entry scripts whether to run exploratory plots or
#' replicate-aware inference.
#'
#' @param example_config Named list describing one analysis run.
#'
#' @return Character scalar, either `"legacy"` or `"replicate"`.
get_analysis_mode <- function(example_config) {
    explicit_mode <- example_config[["mode"]]
    if (!is.null(explicit_mode) && !identical(explicit_mode, "") && !is.na(explicit_mode)) {
        if (!explicit_mode %in% c("legacy", "replicate")) {
            stop("Analysis config `mode` must be either 'legacy' or 'replicate'.")
        }
        return(explicit_mode)
    }

    if (config_uses_sample_manifest(example_config)) {
        return("replicate")
    }

    "legacy"
}

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
#' @param analysis_name Name of the entry in
#'   `proteome_profiler_config$analyses`.
#'
#' @return Named list describing one analysis run.
get_analysis_entries <- function() {
    analyses <- proteome_profiler_config$analyses
    if (!is.null(analyses)) {
        return(analyses)
    }

    legacy_examples <- proteome_profiler_config$examples
    if (!is.null(legacy_examples)) {
        return(legacy_examples)
    }

    stop("Expected `proteome_profiler_config$analyses` to be defined.")
}

#' Normalize one user-facing analysis config into the internal flat shape
#'
#' Users edit `scripts/config/analysis_config.R`, which now supports a more
#' compact grouped layout (`protocol`, `input`, `thresholds`, `shortlist`,
#' `stats`). The rest of the pipeline still expects flat field names, so this
#' helper expands the grouped form into the canonical internal structure.
#'
#' @param analysis_config Named list describing one analysis.
#'
#' @return Named list in the flat internal format expected by the scripts.
normalize_analysis_config <- function(analysis_config) {
    config <- analysis_config

    slug_value <- config[["slug"]]
    if (!is.null(slug_value) && is.null(config[["analysis_slug"]])) {
        config$analysis_slug <- slug_value
    }

    protocol <- config[["protocol"]]
    if (!is.null(protocol)) {
        if (!is.null(protocol[["preset"]]) && is.null(config[["protocol_preset"]])) {
            config$protocol_preset <- protocol[["preset"]]
        }
        if (!is.null(protocol[["workbook"]]) && is.null(config[["info_fn"]])) {
            config$info_fn <- protocol[["workbook"]]
        }
        if (!is.null(protocol[["pdf"]]) && is.null(config[["protocol_pdf"]])) {
            config$protocol_pdf <- protocol[["pdf"]]
        }
        if (!is.null(protocol[["pages"]]) && is.null(config[["protocol_pages"]])) {
            config$protocol_pages <- protocol[["pages"]]
        }
    }

    input <- config[["input"]]
    if (!is.null(input)) {
        if (!is.null(input[["data_dir"]]) && is.null(config[["data_dir"]])) {
            config$data_dir <- input[["data_dir"]]
        }
        if (!is.null(input[["groups"]]) && is.null(config[["group_levels"]])) {
            config$group_levels <- input[["groups"]]
        }
        if (!is.null(input[["manifest"]]) && is.null(config[["sample_manifest"]])) {
            config$sample_manifest <- input[["manifest"]]
        }
        if (!is.null(input[["subgroup"]]) && is.null(config[["subgroup_var"]])) {
            config$subgroup_var <- input[["subgroup"]]
        }
        if (!is.null(input[["treatment"]]) && is.null(config[["treatment_var"]])) {
            config$treatment_var <- input[["treatment"]]
        }
    }

    thresholds <- config[["thresholds"]]
    if (!is.null(thresholds)) {
        if (!is.null(thresholds[["ref_coords"]]) && is.null(config[["ref_coords_to_make_filter"]])) {
            config$ref_coords_to_make_filter <- thresholds[["ref_coords"]]
        }
        if (!is.null(thresholds[["ref_signal"]]) && is.null(config[["ref_thresh_to_filter"]])) {
            config$ref_thresh_to_filter <- thresholds[["ref_signal"]]
        }
        if (!is.null(thresholds[["fold_change"]]) && is.null(config[["main_threshold"]])) {
            config$main_threshold <- thresholds[["fold_change"]]
        }
        if (!is.null(thresholds[["groups_per_page"]]) && is.null(config[["groups_per_page"]])) {
            config$groups_per_page <- thresholds[["groups_per_page"]]
        }
    }

    shortlist <- config[["shortlist"]]
    if (!is.null(shortlist)) {
        if (!is.null(shortlist[["control"]]) && is.null(config[["selection_control"]])) {
            config$selection_control <- shortlist[["control"]]
        }
        if (!is.null(shortlist[["treatment"]]) && is.null(config[["selection_group"]])) {
            config$selection_group <- shortlist[["treatment"]]
        }
        if (!is.null(shortlist[["fold_change"]]) && is.null(config[["selection_threshold"]])) {
            config$selection_threshold <- shortlist[["fold_change"]]
        }
        shortlist_comparisons <- shortlist[["comparisons"]]
        if (is.null(shortlist_comparisons)) {
            shortlist_comparisons <- shortlist[["comparison"]]
        }
        if (!is.null(shortlist_comparisons) && is.null(config[["selection_comparison_slugs"]])) {
            config$selection_comparison_slugs <- unname(as.character(shortlist_comparisons))
        }
        if (!is.null(config[["selection_comparison_slugs"]]) &&
            length(config[["selection_comparison_slugs"]]) == 1 &&
            is.null(config[["selection_comparison_slug"]])) {
            config$selection_comparison_slug <- config[["selection_comparison_slugs"]][[1]]
        }
        if (!is.null(shortlist[["top_n"]]) && is.null(config[["selection_top_n"]])) {
            config$selection_top_n <- shortlist[["top_n"]]
        }
    }

    stats <- config[["stats"]]
    if (!is.null(stats)) {
        if (!is.null(stats[["min_reps_per_arm"]]) && is.null(config[["min_reps_per_arm"]])) {
            config$min_reps_per_arm <- stats[["min_reps_per_arm"]]
        }
        if (!is.null(stats[["p_adjust_method"]]) && is.null(config[["p_adjust_method"]])) {
            config$p_adjust_method <- stats[["p_adjust_method"]]
        }
        if (!is.null(stats[["alpha"]]) && is.null(config[["alpha"]])) {
            config$alpha <- stats[["alpha"]]
        }
    }

    config
}

#' Return one named analysis configuration
#'
#' The returned object is metadata only: where the raw data lives, who owns the
#' output subtree, which comparisons to run, and which thresholds to use.
#'
#' @param analysis_name Name of the entry in
#'   `proteome_profiler_config$analyses`.
#'
#' @return Named list describing one analysis run.
get_analysis_config <- function(analysis_name) {
    config <- get_analysis_entries()[[analysis_name]]
    if (is.null(config)) {
        stop(sprintf("Unknown analysis config: %s", analysis_name))
    }
    config <- normalize_analysis_config(config)
    if (config_uses_sample_manifest(config) && is.null(config$p_adjust_method)) {
        config$p_adjust_method <- "BH"
    }
    config
}

#' Resolve which analysis entry should be used for the current run
#'
#' The active analysis should come from the explicit
#' `PROTEOME_PROFILER_ANALYSIS` environment variable when set. Otherwise, the
#' user-editable config may declare `default_analysis`. This keeps example names
#' out of the script bodies while still allowing one local default per machine.
#'
#' @param requested_name Optional explicit analysis name, typically from
#'   `Sys.getenv("PROTEOME_PROFILER_ANALYSIS")`.
#'
#' @return Character scalar naming the selected config entry.
get_selected_analysis_name <- function(requested_name = NULL) {
    if (!is.null(requested_name) && !identical(requested_name, "") && !is.na(requested_name)) {
        return(requested_name)
    }

    default_name <- proteome_profiler_config$default_analysis
    if (!is.null(default_name) && !identical(default_name, "") && !is.na(default_name)) {
        return(default_name)
    }

    analysis_names <- names(get_analysis_entries())
    if (length(analysis_names) == 1) {
        return(analysis_names[[1]])
    }

    stop(paste(
        "No analysis selected.",
        "Set PROTEOME_PROFILER_ANALYSIS or define proteome_profiler_config$default_analysis",
        "in scripts/config/analysis_config.R."
    ))
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
    analysis_mode <- get_analysis_mode(example_config)

    required_fields <- c("user", "analysis_slug", "info_fn", "protocol_preset", "comparisons")
    if (config_uses_sample_manifest(example_config)) {
        required_fields <- c(required_fields, "sample_manifest", "treatment_var")
    } else {
        required_fields <- c(required_fields, "data_dir", "group_levels")
    }

    if (analysis_mode == "replicate") {
        required_fields <- c(
            required_fields,
            "subgroup_var",
            "alpha"
        )
    }

    missing_fields <- required_fields[vapply(required_fields, function(field) {
        is.null(example_config[[field]]) || identical(example_config[[field]], "")
    }, logical(1))]

    if (length(missing_fields) > 0) {
        stop(sprintf(
            "Analysis config is missing required fields: %s",
            paste(missing_fields, collapse = ", ")
        ))
    }

    if (analysis_mode == "replicate" &&
        !example_config$p_adjust_method %in% p.adjust.methods) {
        stop(sprintf(
            "Unsupported p-value adjustment method '%s'. Supported methods: %s",
            example_config$p_adjust_method,
            paste(p.adjust.methods, collapse = ", ")
        ))
    }
}

#' Return the root output directory for one analysis configuration
#'
#' Every analysis writes into `output/<user>/<analysis_slug>/` so outputs
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
        file.path("output", example_config$user, example_config$analysis_slug),
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
