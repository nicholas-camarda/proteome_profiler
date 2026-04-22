# This module is loaded by scripts/helpers/replicate_analysis.R.
# Source replicate_analysis.R rather than this file directly.

#' Define significance tiers used for inferential plots
#'
#' @param alpha Numeric raw p-value threshold.
#'
#' @return List of plot-tier specifications.
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

#' Write raw-p and FDR-specific inferential waterfall plots
#'
#' @param result_tbl One comparison/method inferential results tibble.
#' @param output_dir Directory where waterfall PNGs should be written.
#' @param filename_prefix Optional method prefix for filenames.
#' @param title_suffix Optional title suffix, usually the method label.
#'
#' @return Named list of written plot paths or `NA` for empty tiers.
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

#' Build barplot-ready fold-change rows from inferential results
#'
#' @param hit_tbl Inferential results filtered to the analyte set to plot.
#' @param mark_treatment_significant Logical; mark treatment bars as significant.
#' @param significance_label Label rendered above significant treatment bars.
#' @param significance_definition Caption definition for the significance label.
#'
#' @return Tibble with control-first fold-change plotting rows.
build_inferential_fold_change_barplot_data <- function(hit_tbl, mark_treatment_significant = TRUE, significance_label = "*", significance_definition = "threshold hit") {
    hit_tbl <- hit_tbl %>%
        filter(is.finite(fold_change_ratio))

    if (nrow(hit_tbl) == 0) {
        return(tibble())
    }
    required_columns <- c(
        "control_fold_change_ymin",
        "control_fold_change_ymax",
        "treatment_fold_change_ymin",
        "treatment_fold_change_ymax"
    )
    missing_columns <- setdiff(required_columns, names(hit_tbl))
    if (length(missing_columns) > 0) {
        stop(sprintf(
            "Cannot build inferential fold-change barplots because required column(s) are missing: %s",
            paste(missing_columns, collapse = ", ")
        ), call. = FALSE)
    }

    control_label <- unique(as.character(hit_tbl$control))
    treatment_label <- unique(as.character(hit_tbl$treatment))
    if (length(control_label) != 1 || length(treatment_label) != 1) {
        stop("Expected one control label and one treatment label when building inferential fold-change barplots.")
    }
    group_levels <- c(control_label, treatment_label)
    analysis_method <- if ("analysis_method" %in% names(hit_tbl)) {
        unique(as.character(hit_tbl$analysis_method))
    } else {
        NA_character_
    }
    if (length(analysis_method) != 1) {
        analysis_method <- NA_character_
    }
    value_column_used <- if ("value_column_used" %in% names(hit_tbl)) {
        unique(as.character(hit_tbl$value_column_used))
    } else {
        NA_character_
    }
    if (length(value_column_used) != 1) {
        value_column_used <- NA_character_
    }

    bind_rows(
        hit_tbl %>%
            transmute(
                Name,
                Coordinate,
                group = as.character(control),
                relative_signal = 1,
                relative_signal_ymin = control_fold_change_ymin,
                relative_signal_ymax = control_fold_change_ymax,
                significant = FALSE,
                significance_label = NA_character_,
                significance_definition = NA_character_,
                analysis_method = analysis_method,
                value_column_used = value_column_used
            ),
        hit_tbl %>%
            transmute(
                Name,
                Coordinate,
                group = as.character(treatment),
                relative_signal = fold_change_ratio,
                relative_signal_ymin = treatment_fold_change_ymin,
                relative_signal_ymax = treatment_fold_change_ymax,
                significant = mark_treatment_significant,
                significance_label = if_else(mark_treatment_significant, significance_label, NA_character_),
                significance_definition = if_else(mark_treatment_significant, significance_definition, NA_character_),
                analysis_method = analysis_method,
                value_column_used = value_column_used
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

#' Format the fold-change y-axis description used on barplot pages
#'
#' @param control_label Character scalar naming the control group.
#' @param analysis_method Character scalar method slug, when available.
#' @param value_column_used Character scalar source value column, when available.
#'
#' @return Character scalar suitable for a page-level subtitle.
format_fold_change_y_axis_title <- function(control_label, analysis_method = NA_character_, value_column_used = NA_character_) {
    method_detail <- case_when(
        identical(analysis_method, "raw_log2_lm") || identical(value_column_used, "signal") ~
            "2^log2 raw-signal model effect",
        identical(analysis_method, "normalized_t_test") || identical(value_column_used, "normalized_signal") ~
            "mean per-sample reference-normalized signal ratio",
        TRUE ~ "group fold-change ratio"
    )

    sprintf("Y-axis: linear fold-change ratio relative to %s (%s)", control_label, method_detail)
}

#' Render paginated 25-per-page fold-change barplots
#'
#' @param barplot_data Tibble from a fold-change barplot-data builder.
#' @param title Character plot title.
#' @param groups_per_page Number of analytes per page.
#' @param ncol Optional facet-column count.
#' @param nrow Optional facet-row count.
#' @param common_y_max Optional fixed y-axis maximum shared across pages.
#' @param y_axis_title Optional page-level subtitle describing the y-axis.
#'
#' @return List of ggplot/patchwork page objects.
plot_inferential_fold_change_barplot_pages <- function(barplot_data, title, groups_per_page = 25, ncol = NULL, nrow = NULL, common_y_max = NULL, y_axis_title = NULL) {
    if (nrow(barplot_data) == 0) {
        return(list())
    }
    if (!"relative_signal_ymin" %in% names(barplot_data)) {
        barplot_data$relative_signal_ymin <- NA_real_
    }
    if (!"relative_signal_ymax" %in% names(barplot_data)) {
        barplot_data$relative_signal_ymax <- NA_real_
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

    if (is.null(y_axis_title)) {
        analysis_method <- if ("analysis_method" %in% names(barplot_data)) {
            first(na.omit(as.character(barplot_data$analysis_method)), default = NA_character_)
        } else {
            NA_character_
        }
        value_column_used <- if ("value_column_used" %in% names(barplot_data)) {
            first(na.omit(as.character(barplot_data$value_column_used)), default = NA_character_)
        } else {
            NA_character_
        }
        y_axis_title <- format_fold_change_y_axis_title(
            control_label = control_label,
            analysis_method = analysis_method,
            value_column_used = value_column_used
        )
    }
    significance_definitions <- barplot_data %>%
        filter(significant, !is.na(significance_label), !is.na(significance_definition)) %>%
        distinct(significance_label, significance_definition) %>%
        arrange(significance_label, significance_definition) %>%
        mutate(caption_text = sprintf("%s = %s", significance_label, significance_definition)) %>%
        pull(caption_text)
    finite_set_signals <- c(
        barplot_data$relative_signal,
        barplot_data$relative_signal_ymax
    )
    finite_set_signals <- finite_set_signals[is.finite(finite_set_signals)]
    set_signal_max <- max(finite_set_signals, 1, na.rm = TRUE)
    if (!is.finite(set_signal_max) || set_signal_max <= 0) {
        set_signal_max <- 1
    }
    required_common_y_max <- max(pretty(c(0, set_signal_max * 1.35), n = 5), na.rm = TRUE)
    if (is.null(common_y_max)) {
        common_y_max <- required_common_y_max
    } else {
        common_y_max <- max(common_y_max, required_common_y_max, na.rm = TRUE)
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
            finite_relative_signals <- c(panel_data$relative_signal, panel_data$relative_signal_ymax)
            finite_relative_signals <- finite_relative_signals[is.finite(finite_relative_signals)]
            panel_max <- max(finite_relative_signals, 1, na.rm = TRUE)
            significance_df <- panel_data %>%
                filter(significant, !is.na(significance_label)) %>%
                summarize(
                    group1 = as.character(levels(short_group)[[1]]),
                    group2 = as.character(levels(short_group)[[2]]),
                    y.position = panel_max * 1.18,
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

            errorbar_df <- panel_data %>%
                filter(
                    is.finite(relative_signal),
                    is.finite(relative_signal_ymin),
                    is.finite(relative_signal_ymax)
                )
            panel_plot <- ggplot(panel_data, aes(x = short_group, y = relative_signal, fill = group, color = group)) +
                geom_col(linewidth = rel(1.1)) +
                geom_errorbar(
                    data = errorbar_df,
                    aes(
                        ymin = pmax(relative_signal_ymin, 0),
                        ymax = relative_signal_ymax
                    ),
                    width = 0.25,
                    linewidth = 0.8,
                    color = "black",
                    inherit.aes = TRUE
                ) +
                scale_y_continuous(
                    limits = c(0, common_y_max),
                    expand = expansion(mult = c(0, 0.04))
                ) +
                coord_cartesian(clip = "off") +
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
                caption = paste(c(significance_definitions, sprintf("Analytes shown = %s / %s", end_idx, length(ordered_names))), collapse = "; ")
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

#' Build barplot-ready fold-change rows from exploratory selected-analyte data
#'
#' @param selected_dat Tibble with exploratory pairwise rows from `make_wf_data()`.
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
        stop("Expected one control label when building exploratory selected-analyte barplots.")
    }

    selected_dat %>%
        mutate(
            group = factor(as.character(group), levels = unique(as.character(group))),
            short_group = factor(as.character(group), levels = levels(group)),
            relative_signal_ymin = NA_real_,
            relative_signal_ymax = NA_real_,
            display_name = str_split(as.character(Name), pattern = "/", simplify = TRUE)[, 1],
            display_name = ifelse(display_name == "Complement Component C5", "C5a", display_name),
            Name_Coordinate = str_c(display_name, "\n(", Coordinate, ")"),
            analysis_method = NA_character_,
            value_column_used = NA_character_,
            significant = if (isTRUE(mark_threshold_hits) && !is.null(threshold)) {
                as.character(group) != control_label & is.finite(relative_signal) &
                    (relative_signal > threshold | relative_signal < 1 / threshold)
            } else {
                FALSE
            },
            significance_label = NA_character_,
            significance_definition = NA_character_
        ) %>%
        select(Name, Coordinate, group, short_group, relative_signal, relative_signal_ymin, relative_signal_ymax, significant, significance_label, significance_definition, display_name, Name_Coordinate)
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

    finite_set_signals <- c(plot_data$relative_signal, plot_data$relative_signal_ymax)
    finite_set_signals <- finite_set_signals[is.finite(finite_set_signals)]
    set_signal_max <- max(finite_set_signals, 1, na.rm = TRUE)
    if (!is.finite(set_signal_max) || set_signal_max <= 0) {
        set_signal_max <- 1
    }
    selected_common_y_max <- max(pretty(c(0, set_signal_max * 1.35), n = 5), na.rm = TRUE)

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

#' Write all-tested and significant-hit inferential barplot sets
#'
#' @param result_tbl One comparison/method inferential results tibble.
#' @param output_dir Method-specific barplot directory.
#' @param filename_prefix Method slug used in filenames.
#' @param groups_per_page Number of analytes per page.
#'
#' @return Named list of semicolon-delimited written plot paths or `NA`.
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
        stop("Cannot build inferential waterfall plot because results are missing effect_se_log2.", call. = FALSE)
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
    if (any(!is.finite(plot_data$effect_se_log2))) {
        missing_se_names <- plot_data %>%
            filter(!is.finite(effect_se_log2)) %>%
            pull(Name) %>%
            as.character()
        stop(sprintf(
            paste(
                "Cannot build inferential waterfall plot because plotted effect estimates are missing finite standard errors.",
                "Missing effect_se_log2 for: %s"
            ),
            paste(missing_se_names, collapse = ", ")
        ), call. = FALSE)
    }

    interval_bounds <- c(
        plot_data$effect_estimate_log2,
        plot_data$effect_estimate_log2 - plot_data$effect_se_log2,
        plot_data$effect_estimate_log2 + plot_data$effect_se_log2
    )
    fill_limit <- max(abs(interval_bounds[is.finite(interval_bounds)]), na.rm = TRUE)
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
                "N analytes plotted = %s; whiskers = +/- 1 SE of the effect estimate",
                nrow(plot_data)
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
