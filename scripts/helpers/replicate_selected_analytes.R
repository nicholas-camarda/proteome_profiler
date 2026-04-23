# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

#' Return configured explicit analyte coordinates for the select-analytes workflow
#'
#' The select-analytes workflow is intentionally explicit: users provide the
#' array coordinates they want to follow up, and the script writes focused plots
#' for those coordinates rather than deriving a shortlist from significance or
#' fold-change thresholds. Coordinates avoid long vendor analyte labels with
#' slashes, aliases, and non-ASCII characters.
#'
#' @param example_config Named analysis config returned by `get_analysis_config()`.
#'
#' @return Character vector of normalized configured coordinates.
get_selected_analyte_coordinates <- function(example_config) {
    selected_coords <- example_config$selection_coords
    if (is.null(selected_coords) || length(selected_coords) == 0 || all(is.na(selected_coords)) || all(selected_coords == "")) {
        stop(paste(
            "Selected analyte coordinates are optional for setup validation and main analysis,",
            "but scripts/select-analytes-analysis.R requires PROTEOME_PROFILER_SHORTLIST_COORDS",
            "in .env to write selected-analyte outputs."
        ), call. = FALSE)
    }

    unique(compact_coordinate_text(as.character(selected_coords)))
}

#' Suggest close protocol/result coordinates for missing selected analytes
#'
#' @param missing_coord Character scalar requested by the user.
#' @param available_coords Character vector of available coordinates.
#' @param n Integer number of suggestions to return.
#'
#' @return Character vector of nearby available coordinates.
suggest_selected_analyte_coordinates <- function(missing_coord, available_coords, n = 5) {
    available_coords <- unique(compact_coordinate_text(as.character(available_coords)))
    if (length(available_coords) == 0) {
        return(character())
    }

    missing_coord <- compact_coordinate_text(missing_coord)
    distances <- utils::adist(tolower(missing_coord), tolower(available_coords))
    available_coords[order(as.numeric(distances), available_coords)][seq_len(min(n, length(available_coords)))]
}

#' Build a lookup table for configured selected-analyte coordinates
#'
#' @param selected_coords Character vector requested in the config.
#' @param available_tbl Tibble containing at least `Name` and `Coordinate`.
#'
#' @return Tibble with selected coordinates, resolved names, and order.
build_selected_coordinate_lookup <- function(selected_coords, available_tbl) {
    required_columns <- c("Name", "Coordinate")
    missing_columns <- setdiff(required_columns, names(available_tbl))
    if (length(missing_columns) > 0) {
        stop(sprintf(
            "Cannot resolve selected-analyte coordinates because required column(s) are missing: %s",
            paste(missing_columns, collapse = ", ")
        ), call. = FALSE)
    }

    selected_lookup <- tibble(
        selection_order = seq_along(selected_coords),
        requested_coordinate = as.character(selected_coords),
        selected_coordinate_key = compact_coordinate_text(selected_coords)
    ) %>%
        distinct(selected_coordinate_key, .keep_all = TRUE)

    available_lookup <- available_tbl %>%
        distinct(Name, Coordinate) %>%
        mutate(selected_coordinate_key = compact_coordinate_text(Coordinate))

    missing_coords <- setdiff(selected_lookup$selected_coordinate_key, available_lookup$selected_coordinate_key)
    if (length(missing_coords) > 0) {
        suggestion_text <- map_chr(missing_coords, function(missing_coord) {
            suggestions <- suggest_selected_analyte_coordinates(missing_coord, available_lookup$Coordinate)
            sprintf(
                "%s (closest available: %s)",
                missing_coord,
                paste(suggestions, collapse = ", ")
            )
        })
        stop(sprintf(
            "Configured selected analyte coordinate(s) were not found: %s",
            paste(suggestion_text, collapse = "; ")
        ), call. = FALSE)
    }

    ambiguous_coords <- available_lookup %>%
        filter(selected_coordinate_key %in% selected_lookup$selected_coordinate_key) %>%
        count(selected_coordinate_key, name = "n_matches") %>%
        filter(n_matches > 1)
    if (nrow(ambiguous_coords) > 0) {
        ambiguity_text <- map_chr(ambiguous_coords$selected_coordinate_key, function(coord_key) {
            matches <- available_lookup %>%
                filter(selected_coordinate_key == coord_key) %>%
                transmute(match_text = sprintf("%s (%s)", Name, Coordinate)) %>%
                pull(match_text)
            sprintf("%s matched: %s", coord_key, paste(matches, collapse = "; "))
        })
        stop(sprintf(
            "Configured selected analyte coordinate(s) are ambiguous: %s",
            paste(ambiguity_text, collapse = "; ")
        ), call. = FALSE)
    }

    selected_lookup %>%
        left_join(available_lookup, by = "selected_coordinate_key") %>%
        arrange(selection_order)
}

#' Validate selected analyte coordinates against an available result table
#'
#' @param selected_coords Character vector requested in the config.
#' @param available_tbl Tibble containing at least `Name` and `Coordinate`.
#'
#' @return Invisibly returns the selected-coordinate lookup.
validate_selected_analyte_coordinates <- function(selected_coords, available_tbl) {
    invisible(build_selected_coordinate_lookup(selected_coords, available_tbl))
}

#' Filter a result table to configured selected-analyte coordinates
#'
#' @param data_tbl Tibble containing at least `Name` and `Coordinate`.
#' @param selected_coords Character vector requested in the config.
#'
#' @return `data_tbl` filtered and ordered by configured coordinate order.
filter_selected_analyte_coordinates <- function(data_tbl, selected_coords) {
    selected_lookup <- build_selected_coordinate_lookup(selected_coords, data_tbl)

    data_tbl %>%
        mutate(selected_coordinate_key = compact_coordinate_text(Coordinate)) %>%
        inner_join(
            selected_lookup %>%
                select(selection_order, selected_coordinate_key),
            by = "selected_coordinate_key"
        ) %>%
        arrange(selection_order) %>%
        select(-selection_order, -selected_coordinate_key)
}

#' Resolve exploratory select-analytes comparisons into control/treatment pairs
#'
#' Exploratory configs can define many comparisons, but select-analytes
#' writes one focused output folder per selected pair. If explicit comparison
#' slugs are configured, they are matched against all configured pairwise
#' comparison slugs; otherwise the exploratory `selection_control` and
#' `selection_group` fields identify the selected pair(s).
#'
#' @param example_config Named analysis config returned by `get_analysis_config()`.
#'
#' @return Tibble with `comparison_slug`, `control`, and `treatment` columns.
resolve_exploratory_select_comparisons <- function(example_config) {
    configured_comparisons <- example_config$comparisons
    if (is.null(configured_comparisons) || length(configured_comparisons) == 0) {
        stop("Exploratory select-analytes requires at least one configured comparison.")
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
                "Configured exploratory select comparison slug(s) were not found: %s. Available slugs: %s",
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
                "Exploratory select comparison '%s' vs '%s' was not found. Available slugs: %s",
                example_config$selection_control,
                paste(example_config$selection_group, collapse = ", "),
                paste(all_pairs$comparison_slug, collapse = ", ")
            ))
        }
        return(selected_pairs)
    }

    all_pairs
}

#' Resolve which inferential comparison(s) should feed the shortlist workflow
#'
#' @param run_index Tibble describing available inferential comparisons.
#' @param example_config Analysis config list.
#'
#' @return Tibble with one row per selected comparison.
resolve_shortlist_comparisons <- function(run_index, example_config) {
    available_comparisons <- run_index %>%
        distinct(comparison_slug, subgroup, control, treatment)
    available_slugs <- available_comparisons$comparison_slug

    configured_slugs <- example_config$selection_comparison_slugs
    if (is.null(configured_slugs) && !is.null(example_config$selection_comparison_slug)) {
        configured_slugs <- example_config$selection_comparison_slug
    }
    if (!is.null(configured_slugs)) {
        configured_slugs <- unique(unname(as.character(configured_slugs)))
        selected <- available_comparisons %>%
            filter(comparison_slug %in% configured_slugs)

        missing_slugs <- setdiff(configured_slugs, selected$comparison_slug)
        if (length(missing_slugs) > 0) {
            stop(sprintf(
                paste(
                    "Configured shortlist comparison slug(s) were not found: %s",
                    "Available slugs: %s",
                    "For replicate-aware .env files, write PROTEOME_PROFILER_SHORTLIST_COMPARISONS as <subgroup value>_<control label>_vs_<treatment label> and separate multiple slugs with |."
                ),
                paste(missing_slugs, collapse = ", "),
                paste(available_slugs, collapse = ", ")
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
        selected <- available_comparisons %>%
            filter(comparison_slug == derived_slug)
        if (nrow(selected) == 0) {
            stop(sprintf(
                "Derived shortlist comparison '%s' was not found. Available slugs: %s",
                derived_slug,
                paste(available_slugs, collapse = ", ")
            ))
        }
        return(slice_head(selected, n = 1))
    }

    if (nrow(available_comparisons) == 1) {
        return(slice_head(available_comparisons, n = 1))
    }

    stop(sprintf(
        paste(
            "Shortlist selection is ambiguous because multiple inferential comparisons are available.",
            "Set PROTEOME_PROFILER_SHORTLIST_COMPARISONS in .env.",
            "Use <subgroup value>_<control label>_vs_<treatment label> and separate multiple slugs with |.",
            "Available slugs: %s"
        ),
        paste(available_slugs, collapse = ", ")
    ))
}

#' Resolve selected-analyte methods for replicate-aware follow-up
#'
#' @param example_config Named analysis config.
#'
#' @return Validated method slugs.
resolve_inferential_shortlist_methods <- function(example_config) {
    configured_method <- example_config$selection_method
    if (is.null(configured_method) || length(configured_method) == 0 || all(is.na(configured_method)) || all(configured_method == "")) {
        configured_method <- example_config$analysis_methods
    }
    validate_analysis_methods(configured_method)
}

#' Read one comparison sheet from a method-specific inferential workbook
#'
#' @param inferential_dir Root `inferential_results/` directory.
#' @param comparison_slug Comparison slug to read.
#' @param analysis_method Method slug whose workbook should be read.
#'
#' @return Tibble of inferential results for the selected comparison and method.
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

#' Remove stale selected-analyte artifacts from one method folder
#'
#' @param output_dir Method-specific selected-analyte output directory.
#'
#' @return Invisibly returns `output_dir`.
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

#' Remove stale selected-analyte artifacts from one comparison folder
#'
#' @param output_dir Comparison-scoped selected-analyte output directory.
#'
#' @return Invisibly returns `output_dir`.
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
