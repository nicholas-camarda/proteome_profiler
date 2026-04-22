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

is_blank_value <- function(value) {
    is.null(value) || length(value) == 0 || all(is.na(value)) || all(str_trim(as.character(value)) == "")
}

env_value <- function(name, default = NULL) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || identical(value, "")) {
        return(default)
    }
    value
}

parse_env_vector <- function(value, required = FALSE, field = "value", delimiter = "\\|") {
    if (is_blank_value(value)) {
        if (required) {
            stop(sprintf("Missing required .env field: %s", field), call. = FALSE)
        }
        return(NULL)
    }

    pieces <- str_split(as.character(value), delimiter, simplify = FALSE)[[1]]
    pieces <- str_trim(pieces)
    pieces <- pieces[pieces != ""]
    if (length(pieces) == 0 && required) {
        stop(sprintf("Missing required .env field: %s", field), call. = FALSE)
    }
    pieces
}

parse_env_numeric_vector <- function(value, required = FALSE, field = "value") {
    pieces <- parse_env_vector(value, required = required, field = field)
    if (is.null(pieces)) {
        return(NULL)
    }

    parsed <- suppressWarnings(as.numeric(pieces))
    if (any(is.na(parsed))) {
        stop(sprintf(
            ".env field %s must contain numeric value(s); got: %s",
            field,
            paste(pieces, collapse = "|")
        ), call. = FALSE)
    }
    parsed
}

parse_env_integer_vector <- function(value, required = FALSE, field = "value") {
    parsed <- parse_env_numeric_vector(value, required = required, field = field)
    if (is.null(parsed)) {
        return(NULL)
    }
    as.integer(parsed)
}

parse_env_comparisons <- function(value, required = TRUE, field = "PROTEOME_PROFILER_COMPARISONS") {
    entries <- parse_env_vector(value, required = required, field = field, delimiter = ";")
    if (is.null(entries)) {
        return(NULL)
    }

    comparisons <- list()
    for (entry in entries) {
        parts <- str_split(entry, "=", n = 2, simplify = TRUE)
        if (ncol(parts) != 2 || is_blank_value(parts[[1]]) || is_blank_value(parts[[2]])) {
            stop(sprintf(
                ".env field %s must use control=treatment1|treatment2 entries; got: %s",
                field,
                entry
            ), call. = FALSE)
        }
        control_label <- str_trim(parts[[1]])
        treatment_labels <- parse_env_vector(parts[[2]], required = TRUE, field = field)
        comparisons[[control_label]] <- treatment_labels
    }

    comparisons
}

env_file_variable_names <- function(env_path) {
    if (!file.exists(env_path)) {
        return(character())
    }

    lines <- readLines(env_path, warn = FALSE)
    lines <- str_trim(lines)
    lines <- lines[lines != "" & !startsWith(lines, "#")]
    matches <- str_match(lines, "^([A-Za-z_][A-Za-z0-9_]*)\\s*=")
    unique(stats::na.omit(matches[, 2]))
}

proteome_profiler_env_names <- function() {
    c(
        "PROTEOME_PROFILER_ANALYSIS",
        "PROTEOME_PROFILER_MODE",
        "PROTEOME_PROFILER_USER",
        "PROTEOME_PROFILER_SLUG",
        "PROTEOME_PROFILER_RUNTIME_ROOT",
        "PROTEOME_PROFILER_CLOUD_PARENT",
        "PROTEOME_PROFILER_PROTOCOL_PRESET",
        "PROTEOME_PROFILER_PROTOCOL_WORKBOOK",
        "PROTEOME_PROFILER_PROTOCOL_PDF",
        "PROTEOME_PROFILER_PROTOCOL_PAGES",
        "PROTEOME_PROFILER_INPUT_MANIFEST",
        "PROTEOME_PROFILER_INPUT_DATA_DIR",
        "PROTEOME_PROFILER_TREATMENT_COLUMN",
        "PROTEOME_PROFILER_SUBGROUP_COLUMN",
        "PROTEOME_PROFILER_GROUP_LEVELS",
        "PROTEOME_PROFILER_COMPARISONS",
        "PROTEOME_PROFILER_REF_COORDS",
        "PROTEOME_PROFILER_REF_SIGNAL",
        "PROTEOME_PROFILER_FOLD_CHANGE",
        "PROTEOME_PROFILER_GROUPS_PER_PAGE",
        "PROTEOME_PROFILER_MIN_REPS_PER_ARM",
        "PROTEOME_PROFILER_P_ADJUST_METHOD",
        "PROTEOME_PROFILER_ALPHA",
        "PROTEOME_PROFILER_ANALYSIS_METHODS",
        "PROTEOME_PROFILER_SHORTLIST_CONTROL",
        "PROTEOME_PROFILER_SHORTLIST_TREATMENT",
        "PROTEOME_PROFILER_SHORTLIST_FOLD_CHANGE",
        "PROTEOME_PROFILER_SHORTLIST_COMPARISONS",
        "PROTEOME_PROFILER_SHORTLIST_METHOD",
        "PROTEOME_PROFILER_SHORTLIST_METHODS",
        "PROTEOME_PROFILER_SHORTLIST_ANALYTES"
    )
}

load_proteome_profiler_env <- function(required = FALSE) {
    env_path <- Sys.getenv("PROTEOME_PROFILER_ENV_FILE", unset = file.path(repo_root, ".env"))
    env_path <- path.expand(env_path)
    if (!file.exists(env_path)) {
        if (required) {
            stop(sprintf(
                paste(
                    "Missing .env file: %s",
                    "Create it with: cp .env.example .env"
                ),
                env_path
            ), call. = FALSE)
        }
        return(invisible(FALSE))
    }

    variable_names <- unique(c(env_file_variable_names(env_path), proteome_profiler_env_names()))
    existing_values <- Sys.getenv(variable_names, unset = NA_character_)

    if (!requireNamespace("dotenv", quietly = TRUE)) {
        stop(paste(
            "Missing required R package: dotenv",
            "Install required packages with: Rscript scripts/install_packages.R"
        ), call. = FALSE)
    }

    Sys.unsetenv(variable_names)
    dotenv::load_dot_env(env_path)

    # `dotenv` overwrites existing environment variables. Restore explicit
    # process-level values so CLI/test overrides remain authoritative.
    for (idx in seq_along(variable_names)) {
        if (!is.na(existing_values[[idx]])) {
            do.call(Sys.setenv, as.list(stats::setNames(existing_values[[idx]], variable_names[[idx]])))
        }
    }

    invisible(TRUE)
}

build_proteome_profiler_config_from_env <- function() {
    analysis_name <- env_value("PROTEOME_PROFILER_ANALYSIS")
    if (is_blank_value(analysis_name)) {
        stop("Missing required .env field: PROTEOME_PROFILER_ANALYSIS", call. = FALSE)
    }

    analysis_mode <- env_value("PROTEOME_PROFILER_MODE")
    if (is_blank_value(analysis_mode)) {
        stop("Missing required .env field: PROTEOME_PROFILER_MODE", call. = FALSE)
    }

    input_manifest <- env_value("PROTEOME_PROFILER_INPUT_MANIFEST")
    input_data_dir <- env_value("PROTEOME_PROFILER_INPUT_DATA_DIR")
    treatment_column <- env_value("PROTEOME_PROFILER_TREATMENT_COLUMN", "treatment")
    subgroup_column <- env_value("PROTEOME_PROFILER_SUBGROUP_COLUMN")
    group_levels <- parse_env_vector(env_value("PROTEOME_PROFILER_GROUP_LEVELS"), field = "PROTEOME_PROFILER_GROUP_LEVELS")

    input <- list(
        manifest = input_manifest,
        data_dir = input_data_dir,
        treatment = treatment_column,
        subgroup = subgroup_column,
        groups = group_levels
    )
    input <- input[!vapply(input, is.null, logical(1))]

    shortlist <- list(
        control = env_value("PROTEOME_PROFILER_SHORTLIST_CONTROL"),
        treatment = env_value("PROTEOME_PROFILER_SHORTLIST_TREATMENT"),
        fold_change = parse_env_numeric_vector(env_value("PROTEOME_PROFILER_SHORTLIST_FOLD_CHANGE"), field = "PROTEOME_PROFILER_SHORTLIST_FOLD_CHANGE"),
        comparisons = parse_env_vector(env_value("PROTEOME_PROFILER_SHORTLIST_COMPARISONS"), field = "PROTEOME_PROFILER_SHORTLIST_COMPARISONS"),
        method = parse_env_vector(
            env_value("PROTEOME_PROFILER_SHORTLIST_METHODS", env_value("PROTEOME_PROFILER_SHORTLIST_METHOD")),
            field = "PROTEOME_PROFILER_SHORTLIST_METHODS"
        ),
        analytes = parse_env_vector(env_value("PROTEOME_PROFILER_SHORTLIST_ANALYTES"), field = "PROTEOME_PROFILER_SHORTLIST_ANALYTES")
    )
    shortlist <- shortlist[!vapply(shortlist, is.null, logical(1))]

    analysis_config <- list(
        mode = analysis_mode,
        user = env_value("PROTEOME_PROFILER_USER"),
        slug = env_value("PROTEOME_PROFILER_SLUG", analysis_name),
        protocol = list(
            preset = env_value("PROTEOME_PROFILER_PROTOCOL_PRESET", "cytokine_xl"),
            workbook = env_value("PROTEOME_PROFILER_PROTOCOL_WORKBOOK"),
            pdf = env_value("PROTEOME_PROFILER_PROTOCOL_PDF"),
            pages = parse_env_integer_vector(env_value("PROTEOME_PROFILER_PROTOCOL_PAGES"), field = "PROTEOME_PROFILER_PROTOCOL_PAGES")
        ),
        input = input,
        comparisons = parse_env_comparisons(env_value("PROTEOME_PROFILER_COMPARISONS"), required = TRUE),
        thresholds = list(
            ref_coords = parse_env_vector(env_value("PROTEOME_PROFILER_REF_COORDS"), field = "PROTEOME_PROFILER_REF_COORDS"),
            ref_signal = parse_env_numeric_vector(env_value("PROTEOME_PROFILER_REF_SIGNAL"), field = "PROTEOME_PROFILER_REF_SIGNAL"),
            fold_change = parse_env_numeric_vector(env_value("PROTEOME_PROFILER_FOLD_CHANGE"), field = "PROTEOME_PROFILER_FOLD_CHANGE"),
            groups_per_page = parse_env_integer_vector(env_value("PROTEOME_PROFILER_GROUPS_PER_PAGE"), field = "PROTEOME_PROFILER_GROUPS_PER_PAGE")
        ),
        stats = list(
            min_reps_per_arm = parse_env_integer_vector(env_value("PROTEOME_PROFILER_MIN_REPS_PER_ARM"), field = "PROTEOME_PROFILER_MIN_REPS_PER_ARM"),
            p_adjust_method = env_value("PROTEOME_PROFILER_P_ADJUST_METHOD", "BH"),
            alpha = parse_env_numeric_vector(env_value("PROTEOME_PROFILER_ALPHA", "0.05"), field = "PROTEOME_PROFILER_ALPHA"),
            methods = parse_env_vector(env_value("PROTEOME_PROFILER_ANALYSIS_METHODS"), field = "PROTEOME_PROFILER_ANALYSIS_METHODS")
        ),
        shortlist = shortlist
    )

    analysis_config$protocol <- analysis_config$protocol[!vapply(analysis_config$protocol, is.null, logical(1))]
    analysis_config$thresholds <- analysis_config$thresholds[!vapply(analysis_config$thresholds, is.null, logical(1))]
    analysis_config$stats <- analysis_config$stats[!vapply(analysis_config$stats, is.null, logical(1))]

    list(
        default_analysis = analysis_name,
        runtime_root = env_value("PROTEOME_PROFILER_RUNTIME_ROOT", "~/ProjectsRuntime/proteome_profiler"),
        cloud_parent = env_value("PROTEOME_PROFILER_CLOUD_PARENT", ""),
        analyses = stats::setNames(list(analysis_config), analysis_name)
    )
}

initialize_runtime_config_from_env <- function(required_env_file = FALSE) {
    load_proteome_profiler_env(required = required_env_file)
    proteome_profiler_config <<- build_proteome_profiler_config_from_env()
    invisible(proteome_profiler_config)
}

ensure_runtime_config <- function() {
    if (!exists("proteome_profiler_config", inherits = TRUE)) {
        initialize_runtime_config_from_env(required_env_file = TRUE)
    }
    invisible(proteome_profiler_config)
}

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
#' @return Character scalar, either `"exploratory"` or `"replicate"`.
get_analysis_mode <- function(example_config) {
    explicit_mode <- example_config[["mode"]]
    if (!is.null(explicit_mode) && !identical(explicit_mode, "") && !is.na(explicit_mode)) {
        if (!explicit_mode %in% c("exploratory", "replicate")) {
            stop("Analysis config `mode` must be either 'exploratory' or 'replicate'.")
        }
        return(explicit_mode)
    }

    if (config_uses_sample_manifest(example_config)) {
        return("replicate")
    }

    "exploratory"
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
    ensure_runtime_config()
    runtime_root <- path.expand(proteome_profiler_config$runtime_root)
    cloud_parent <- proteome_profiler_config$cloud_parent
    cloud_project_root <- if (!is_blank_value(cloud_parent)) {
        file.path(path.expand(cloud_parent), "proteome_profiler")
    } else {
        NULL
    }

    unique(c(
        file.path(repo_root, path),
        file.path(runtime_root, path),
        if (!is.null(cloud_project_root)) file.path(cloud_project_root, path)
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
    if (grepl("^~|^/", path)) {
        expanded <- path.expand(path)
        if (file.exists(expanded) || dir.exists(expanded)) {
            return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
        }
        if (!must_exist) {
            return(normalizePath(expanded, winslash = "/", mustWork = FALSE))
        }
        stop(sprintf("Could not resolve absolute path '%s'.", path), call. = FALSE)
    }

    candidates <- path_candidates(path)
    existing <- candidates[file.exists(candidates) | dir.exists(candidates)]

    if (length(existing) > 0) {
        return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
    }

    if (!must_exist) {
        # Writable outputs materialize under the runtime tree; repo and cloud
        # roots are treated as read/search locations.
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
    ensure_runtime_config()
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

#' Normalize one `.env`-derived analysis config into the internal flat shape
#'
#' The `.env` parser creates a compact grouped layout (`protocol`, `input`,
#' `thresholds`, `shortlist`, `stats`). The rest of the pipeline still expects
#' flat field names, so this helper expands the grouped form into the canonical
#' internal structure.
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
        if (!is.null(shortlist[["basis"]]) && is.null(config[["selection_basis"]])) {
            config$selection_basis <- shortlist[["basis"]]
        }
        if (!is.null(shortlist[["method"]]) && is.null(config[["selection_method"]])) {
            config$selection_method <- shortlist[["method"]]
        }
        if (!is.null(shortlist[["analytes"]]) && is.null(config[["selection_analytes"]])) {
            config$selection_analytes <- unname(as.character(shortlist[["analytes"]]))
        }
        if (!is.null(shortlist[["write_bargraphs"]]) && is.null(config[["selection_write_bargraphs"]])) {
            config$selection_write_bargraphs <- isTRUE(shortlist[["write_bargraphs"]])
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
        if (!is.null(stats[["methods"]]) && is.null(config[["analysis_methods"]])) {
            config$analysis_methods <- unname(as.character(stats[["methods"]]))
        }
    }

    if (is_replicate_aware_config(config)) {
        unsupported_shortlist_fields <- character()
        if (!is.null(config[["selection_basis"]])) {
            unsupported_shortlist_fields <- c(unsupported_shortlist_fields, "`shortlist$basis` / `selection_basis`")
        }
        if (!is.null(config[["selection_top_n"]])) {
            unsupported_shortlist_fields <- c(unsupported_shortlist_fields, "`shortlist$top_n` / `selection_top_n`")
        }
        if (!is.null(config[["selection_write_bargraphs"]])) {
            unsupported_shortlist_fields <- c(unsupported_shortlist_fields, "`shortlist$write_bargraphs` / `selection_write_bargraphs`")
        }
        if (length(unsupported_shortlist_fields) > 0) {
            stop(sprintf(
                paste(
                    "Replicate-aware `select-analytes-analysis.R` uses explicit `shortlist$analytes`",
                    "and writes selected outputs under `select_analytes/<comparison_slug>/<method>/`.",
                    "Remove these unsupported shortlist fields from the analysis config: %s"
                ),
                paste(unsupported_shortlist_fields, collapse = ", ")
            ))
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
    if (get_analysis_mode(config) == "replicate" && is.null(config$analysis_methods)) {
        config$analysis_methods <- "raw_log2_lm"
    }
    config
}

#' Resolve which analysis entry should be used for the current run
#'
#' The active analysis comes from the `.env` run sheet or an explicit process
#' environment override.
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
        "Set PROTEOME_PROFILER_ANALYSIS in .env."
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
