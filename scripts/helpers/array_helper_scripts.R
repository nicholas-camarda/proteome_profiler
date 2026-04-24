#' Build the analyte-level dataset from LI-COR workbook exports
#'
#' Each input workbook is assumed to represent one biological sample or one
#' already-collapsed analysis group. Within a workbook, the duplicate membrane
#' spots for each analyte are averaged into one signal value. This function does
#' not model biological replicates across multiple workbooks that share the same
#' group label.
#'
#' @param data_dir Directory containing the LI-COR Excel exports.
#' @param analyte_info Data frame extracted from the protocol appendix, mapping
#'   array positions to analyte names.
#' @param preview Logical. When `TRUE`, print the first part of the assembled
#'   long-format dataset before widening.
#' @param my_group_lvls Character vector giving the treatment groups, in order.
#'
#' @return Wide data frame with one row per analyte and one signal column per
#'   group.
make_plot_ready_dataset <- function(data_dir, analyte_info, preview = FALSE, my_group_lvls = NA) {
    group_levels <- as.character(my_group_lvls)

    #' @note this function is overkill but basically changes "A1, A2" into "A1,2"
    #' just more concise
    convert_text <- function(text_vector) {
        text_list <- str_split(text_vector, ", ")

        output <- sapply(text_list, function(words) {
            letters <- gsub("\\d+", "", words)
            numbers <- gsub("\\D+", "", words)

            unique_letters <- unique(letters)
            output <- character(length(unique_letters))
            for (i in seq_along(unique_letters)) {
                matching_numbers <- numbers[letters == unique_letters[i]]
                output[i] <- paste(unique_letters[i], paste(matching_numbers, collapse = ","), sep = "")
            }

            output <- paste(output, collapse = ", ")

            return(output)
        })

        return(output)
    }

    # read in LICOR data
    # start the data collection with the inital df
    df_temp <- tibble(data_fns = dir(data_dir, pattern = "\\.xlsx$", full.names = TRUE)) %>%
        mutate(
            group_temp = str_split(basename(data_fns), pattern = "\\.", simplify = TRUE)[, 1],
            group_fn = str_split(group_temp, " - ", simplify = TRUE)[, 1], .before = 1
        )

    if (nrow(df_temp) == 0) {
        stop(qq("No .xlsx files were found in @{data_dir}."))
    }

    duplicate_groups <- df_temp %>%
        count(group_fn, name = "n_files") %>%
        filter(n_files > 1)

    if (nrow(duplicate_groups) > 0) {
        stop(qq(
            paste(
                "Found multiple LI-COR workbooks for the same analysis group.",
                "This pipeline currently supports one workbook per group and only averages technical duplicate spots within a membrane.",
                "Biological replicates must be modeled explicitly before this script can be used.",
                "Duplicate groups: @{str_c(duplicate_groups$group_fn, collapse = ', ')}"
            )
        ))
    }

    missing_groups <- df_temp %>%
        filter(!(group_fn %in% group_levels)) %>%
        distinct(group_fn) %>%
        pull(group_fn)

    if (length(missing_groups) > 0) {
        stop(qq(
            "Found workbook prefixes that do not match GROUP_LVLS: @{str_c(missing_groups, collapse = ', ')}"
        ))
    }

    remove_fn <- file.path(data_dir, "..", "bad_analytes.xlsx")
    if (file.exists(remove_fn)) {
        remove_dat <- read_excel(remove_fn) %>%
            mutate(Sname = as.character(.[[1]]), group = as.character(.[[2]]))
    } else {
        remove_dat <- tibble(Sname = "", group = "")
    }

    # grab each individual df and modify to our needs
    init_df <- df_temp %>%
        dplyr::select(-group_temp) %>%
        mutate(group_fn = factor(group_fn, levels = group_levels)) %>%
        arrange(group_fn) %>%
        mutate(lst_df = map2(data_fns, group_fn, .f = function(x, y) {
            # y = init_df$group_fn[1]; x = init_df$data_fns[1]
            ds_tmp <- read_excel(x, skip = 3)
            n_measurements_original <- nrow(ds_tmp) - 1 # remove background

            res_temp_pre <- ds_tmp %>%
                mutate(group = y) %>%
                rename(Sname = Name) %>%
                dplyr::select(Sname, Signal, group) %>% # could use SNR
                filter(!startsWith(Sname, "B")) %>% # remove the background spot because we've already taken it into account
                mutate(sname_grouping = rep(1:(nrow(.) / 2), each = 2))

            expected_snames <- sprintf("S%03d", seq_len(nrow(res_temp_pre)))
            if (nrow(res_temp_pre) %% 2 != 0 || !identical(res_temp_pre$Sname, expected_snames)) {
                stop(qq(
                    "Unexpected LI-COR spot ordering in @{x}. Expected sequential S001... labels after removing background."
                ))
            }

            sname_map <- res_temp_pre %>%
                group_by(sname_grouping) %>%
                reframe(init_sname = Sname, Sname = str_c(Sname, collapse = ", "))

            # remove the bad ones
            res_temp <- res_temp_pre %>%
                anti_join(remove_dat, by = c("Sname", "group"))

            # modify the remove_dat a bit to help with identifying pairs that have had a technical replicate removed
            # doing it this way guarantees no duplicates, since Sname will always be a pair, even if one of the pairs was removed in the data--> hence the logical column
            helper_remove_dat <- remove_dat %>%
                mutate(removed_tech_rep_artifact = TRUE) %>%
                left_join(sname_map, by = c("Sname" = "init_sname")) %>%
                rename(full_sname = Sname.y) %>%
                distinct(group, full_sname, removed_tech_rep_artifact)

            n_measurements_new <- nrow(res_temp)
            message(qq("Removed @{n_measurements_original - n_measurements_new} measurements from @{y} dataset."))

            res <- res_temp %>%
                group_by(sname_grouping, group) %>%
                reframe(
                    Sname = Sname,
                    Sname_signal_pair = str_c(round(as.numeric(Signal), 2), collapse = ", "),
                    signal = mean(as.numeric(Signal), na.rm = TRUE)
                ) %>%
                left_join(sname_map, by = c("Sname" = "init_sname", "sname_grouping")) %>%
                rename(full_sname = Sname.y) %>%
                left_join(helper_remove_dat, by = c("group", "full_sname")) %>%
                dplyr::select(-Sname) %>%
                rename(Sname = full_sname) %>%
                mutate(removed_tech_rep_artifact = ifelse(is.na(removed_tech_rep_artifact), FALSE, removed_tech_rep_artifact)) %>%
                distinct()
            res %>% filter(removed_tech_rep_artifact)
            return(res)
        })) %>%
        dplyr::select(-data_fns, -group_fn) %>%
        unnest(lst_df) %>%
        group_by(group) %>%
        left_join(analyte_info, by = join_by(sname_grouping))

    if (preview) init_df %>% print(n = 50)

    # get into wide format for compatibility with the rest of the script
    my_df <- init_df %>%
        dplyr::select(Name, Sname, Coordinate, group, signal) %>%
        pivot_wider(id_cols = c(Name, Sname, Coordinate), names_from = group, values_from = signal) %>%
        # mutate(Name = make.unique(Name)) %>%
        mutate(Coordinate = convert_text(Coordinate))

    # my_df %>% print(n = Inf)

    return(my_df)
}

# manually identified threshold
#' Derive a heuristic raw-signal cutoff from selected analyte coordinates
#'
#' This helper is intentionally heuristic. When `ref_coords` are supplied, they
#' are treated as a manually chosen low-signal analyte panel whose average
#' signal marks the rough boundary between usable and unusable analytes for this
#' dataset. They are not assumed to be the array's true control spots unless
#' the caller explicitly chooses those coordinates.
#' @param wide_df df produced by make_plot_ready_dataset() function
#' @param ref_coords coordinates to inspect, in the form of a vector of
#'   duplicates e.g. c("A5,6", "A7,8"). In the current VEGFRi/Dox example these
#'   are manually selected low-signal analytes, not true protocol controls.
#' @param filter_threshold signal threshold to filter on,  usually picked by looking at histograms of each of the datasets.
#' this will *keep* analytes with signal strictly greater than this number
#' @param my_colors named vector with groups as names and values as colors for those groups to be displayed
find_filter_thresh <- function(wide_df, filter_threshold = NULL, ref_coords = NULL, my_colors = c()) {
    # Check if at least one of ref_coords or filter_threshold is provided
    if (is.null(ref_coords) && is.null(filter_threshold)) {
        stop("Please provide either 'ref_coords' or 'filter_threshold'")
    }

    cols_to_pivot <- which(!(colnames(wide_df) %in% c("Name", "Sname", "Coordinate")))
    long_df <- wide_df %>%
        pivot_longer(all_of(cols_to_pivot), names_to = "group", values_to = "signal")

    # `Signal Level` <- long_df$signal
    # hist_info <- hist(`Signal Level`, breaks = 100, plot = TRUE)

    # Now, hist_info$breaks contains the break points (edges of the bins)
    # and hist_info$counts contains the counts within each bin.
    # print(hist_info$breaks)
    # print(hist_info$counts)

    if (!is.null(ref_coords) && is.null(filter_threshold)) {
        message("Filtering using ref_coords...")
        # These coordinates define the low-signal reference panel used to set a
        # candidate raw-signal floor for the rest of the analysis.
        ref_coords_long_df <- long_df %>%
            filter(Coordinate %in% ref_coords)
        subtitle <- qq("Region coords: @{str_c(ref_coords, collapse = '; ')}")
    } else if (is.null(ref_coords) && !is.null(filter_threshold)) {
        message("Filtering using threshold...")
        ref_coords_long_df <- long_df %>%
            filter(signal <= filter_threshold) %>%
            mutate(Name = str_c(Name, " (", Coordinate, ")"))
        subtitle <- qq("Signal <= @{filter_threshold}")
    } else {
        message("Filtering using both ref_coords and threshold...")
        # do both if both are defined
        ref_coords_long_df <- long_df %>%
            filter(Coordinate %in% ref_coords, signal <= filter_threshold)
        subtitle <- qq("Region coords: @{str_c(ref_coords, collapse = '; ')}\nSignal <= @{filter_threshold}")
    }

    ref_thresh_to_filter_df <- ref_coords_long_df %>%
        ungroup() %>%
        summarize(
            mean = mean(signal, na.rm = TRUE),
            sd = sd(signal, na.rm = TRUE)
        )

    # print(filtered_long_df, n = Inf)
    boxp <- ggplot(ref_coords_long_df, mapping = aes(x = Name, y = signal)) +
        geom_boxplot(aes(group = Name), outlier.shape = NA) +
        geom_jitter(aes(fill = group), position = position_jitter(0.1), size = rel(2.5), shape = 21) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        scale_fill_manual(values = my_colors) +
        labs(
            title = "Boxplot of analyte signals in select 'low signal' region across groups",
            subtitle = qq("@{subtitle}"),
            caption = qq("mean = @{round(ref_thresh_to_filter_df$mean, 3)}\nsd = @{round(ref_thresh_to_filter_df$sd, 3)}")
        )

    filtered_long_df <- long_df %>% filter(signal > ref_thresh_to_filter_df$mean)
    n_removed <- length(setdiff(long_df$Name, filtered_long_df$Name))
    n_remaining <- length(filtered_long_df$Name %>% unique())

    histop <- ggplot(long_df, mapping = aes(x = signal)) +
        geom_histogram(fill = "black", position = "identity", binwidth = 250, show.legend = FALSE) +
        scale_fill_manual(values = my_colors) +
        geom_histogram(
            data = long_df %>% filter(signal < ref_thresh_to_filter_df$mean),
            mapping = aes(fill = group),
            color = "black", binwidth = 250,
            show.legend = FALSE
        ) +
        geom_vline(xintercept = round(ref_thresh_to_filter_df$mean), color = "#dd0ea9", linewidth = rel(1.25), linetype = 4) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        labs(
            title = "Proportion of analytes removed",
            subtitle = qq("Mean signal intensity threshold = @{round(ref_thresh_to_filter_df$mean)}"),
            caption = qq("Removed = @{n_removed}/@{length(unique(long_df$Name))}")
        )
    histop

    histop2 <- ggplot(filtered_long_df,
        mapping = aes(x = signal)
    ) +
        scale_fill_manual(values = my_colors) +
        geom_histogram(
            mapping = aes(fill = group),
            color = "black", binwidth = 250,
        ) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45), panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        labs(
            title = "Proportion of analytes remaining after gating",
            # subtitle = qq("Mean signal intensity threshold = @{round(ref_thresh_to_filter_df$mean)}"),
            caption = qq("Remaining = @{n_remaining}/@{length(unique(long_df$Name))}")
        )
    histop2

    plt_lst <- boxp / (histop + histop2)

    return(list(ref_thresh_to_filter_df, plt_lst))
}

#' Suggest candidate low-signal analytes for threshold selection
#'
#' This helper provides a starting point for `ref_coords_to_make_filter`. It
#' ranks analytes that are consistently low across the currently configured
#' groups so the user can review plausible candidates instead of choosing every
#' coordinate by hand from scratch. The result is still advisory: users should
#' inspect their own data and confirm that the suggested analytes really sit near
#' the boundary between usable and unusable signal for their assay.
#'
#' @param wide_df Data frame produced by `make_plot_ready_dataset()`.
#' @param candidate_n Integer number of candidates to return.
#' @param exclude_name_patterns Character vector of regex patterns used to drop
#'   rows such as true reference spots or negative controls from the suggestion
#'   list.
#'
#' @return Tibble of suggested low-signal analytes with ranking metadata.
suggest_low_signal_panel <- function(wide_df, candidate_n = 12, exclude_name_patterns = c("^Reference Spots$", "^Negative Control$")) {
    group_cols <- setdiff(colnames(wide_df), c("Name", "Sname", "Coordinate"))

    score_df <- wide_df %>%
        rowwise() %>%
        mutate(
            mean_signal = mean(c_across(all_of(group_cols)), na.rm = TRUE),
            max_signal = max(c_across(all_of(group_cols)), na.rm = TRUE),
            min_signal = min(c_across(all_of(group_cols)), na.rm = TRUE),
            sd_signal = sd(c_across(all_of(group_cols)), na.rm = TRUE),
            cv_signal = if_else(mean_signal > 0, sd_signal / mean_signal, Inf)
        ) %>%
        ungroup()

    if (length(exclude_name_patterns) > 0) {
        exclude_pattern <- str_c(exclude_name_patterns, collapse = "|")
        score_df <- score_df %>%
            filter(!str_detect(Name, exclude_pattern))
    }

    score_df %>%
        mutate(
            mean_rank = min_rank(mean_signal),
            max_rank = min_rank(max_signal),
            cv_rank = min_rank(cv_signal),
            suggestion_score = mean_rank + max_rank + cv_rank
        ) %>%
        arrange(suggestion_score, mean_signal, max_signal, cv_signal, Name, Coordinate) %>%
        transmute(
            suggested_rank = row_number(),
            Name,
            Coordinate,
            mean_signal = round(mean_signal, 3),
            max_signal = round(max_signal, 3),
            min_signal = round(min_signal, 3),
            sd_signal = round(sd_signal, 3),
            cv_signal = round(cv_signal, 3),
            suggestion_score
        ) %>%
        slice_head(n = candidate_n)
}

#' @note make the bar plots
#' @param data looks like dat_filtered
#' @param title title of the barplot
#' @param add_fc whether to add the raw fold changes
#' @param strip_text_rel_size relative sizing of the facet titles for each analyte
#' @param make_facet_plot TRUE = make a large faceted plot (multipage). FALSE = separate each plot into a separate PNG
#' @param groups_per_page if make_facet_plot is true, this is used to determine how many plots get put on each page (make it square)
#' @param caption boolean, add the number of analytes and the filtering procedure or not
#' @return ggplot2 bar chart
make_bar_charts <- function(data, title = "", add_fc = TRUE, strip_text_rel_size = 8, make_facet_plot = TRUE, caption = TRUE, groups_per_page = 25) {
    short_group_data <- data %>%
        # is there a better way to do abbreivations automatically?
        mutate(short_group = group) %>%
        mutate(color = "black") %>%
        arrange(Name) %>%
        group_by(group) %>%
        mutate(
            Name = str_split(as.character(Name), pattern = "/", simplify = TRUE)[, 1],
            Name = ifelse(Name == "Complement Component C5", "C5a", as.character(Name))
        ) %>%
        ungroup() %>%
        mutate(Name_Coordinate = str_c(Name, " (", Coordinate, ")")) %>%
        mutate(
            order = seq_len(n()),
            Name = factor(Name, levels = unique(Name[unique(order)])),
            Name_Coordinate = factor(Name_Coordinate, levels = unique(Name_Coordinate[unique(order)]))
        ) %>%
        ungroup() %>%
        arrange(order, group)

    short_group_data_names <- short_group_data %>%
        distinct(Name_Coordinate, order) %>%
        arrange(order) %>%
        pull("Name_Coordinate")

    n_analytes <- length(unique(short_group_data_names))
    n_groups <- length(unique(short_group_data$group))

    stopifnot(nrow(short_group_data) > 0)

    my_control <- unique(short_group_data$control)
    if (make_facet_plot) {
        # Number of rows to show per page
        # Get unique Name_Coordinate values
        unique_names <- unique(short_group_data$Name_Coordinate)

        # Calculate total pages needed
        facets_per_page <- groups_per_page
        total_pages <- ceiling(length(unique_names) / facets_per_page)

        # lapply
        bar_plot_lst <- lapply(seq_len(total_pages), FUN = function(i) {
            # Calculate start and end indices for unique Name_Coordinate
            start_idx <- (i - 1) * facets_per_page + 1
            end_idx <- min(i * facets_per_page, length(unique_names))

            # Get the subset of unique Name_Coordinate for this page
            subset_names <- unique_names[start_idx:end_idx]

                ordering_group <- setdiff(levels(short_group_data$group), my_control)[1]
                ordered_subset_names <- short_group_data %>%
                    filter(Name_Coordinate %in% subset_names, group == ordering_group) %>%
                    arrange(desc(relative_signal)) %>%
                    mutate(sorafenib_order = seq_len(n())) %>%
                    arrange(sorafenib_order) %>%
                    pull(Name_Coordinate)

                if (length(ordered_subset_names) == 0) {
                    ordered_subset_names <- subset_names
                }

            # Slice the data based on these Name_Coordinate values
            sliced_data <- filter(short_group_data, Name_Coordinate %in% subset_names) %>%
                mutate(Name_Coordinate = factor(Name_Coordinate, levels = ordered_subset_names))

            # Do this to help with the aspect ratio!!
            if (i == max(total_pages)) {
                n_blank <- (groups_per_page - length(subset_names))

                # Create a dummy data frame
                unique_na_names <- paste("IGNORE_SPACER", seq_len(n_blank), sep = "_")
                dummy_data <- expand.grid(
                    Name_Coordinate = factor(unique_na_names),
                    short_group = unique(sliced_data$short_group), # Assuming sliced_data$short_group contains all the groups
                    relative_signal = NA_real_,
                    significant = NA,
                    order = Inf,
                    stringsAsFactors = FALSE
                ) %>%
                    as_tibble() %>%
                    mutate(group = short_group, dummy = TRUE) %>%
                    filter(!is.na(short_group)) %>%
                    arrange(Name_Coordinate)

                # Add the dummy data to the original data
                sliced_data <- bind_rows(sliced_data, dummy_data)
            }

            # real_data <- filter(sliced_data, !dummy)

            bar_plot <- ggplot(sliced_data, aes(x = short_group, y = relative_signal, fill = group, color = group)) +
                geom_col(linewidth = rel(1.1)) +
                # facet_wrap(~Name,
                #     scales = "free_y",
                # ) +
                # ggh4x::facet_wrap2(~Name_Coordinate,
                #     scales = "free",
                #     strip = my_strip_colors,
                #     labeller = labeller(Name_Coordinate = label_wrap_gen(12)),
                #     ncol = 5
                # ) +
                ggforce::facet_wrap_paginate(~Name_Coordinate,
                    scales = "free",
                    # strip = my_strip_colors,
                    labeller = labeller(Name_Coordinate = label_wrap_gen(10)),
                    ncol = ceiling(sqrt(facets_per_page)),
                    page = i
                ) +
                theme_prism(base_size = 16) +
                scale_fill_manual(values = my_colors) +
                scale_color_manual(values = my_outline_colors) +
                labs(x = "Group", y = qq("Fold Change to @{my_control}")) +
                ggtitle(title) +
                theme(
                    strip.text.x = element_text(size = rel(strip_text_rel_size)),
                    axis.text = element_text(size = rel(1.75)),
                    panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
                    panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
                    panel.spacing = unit(1.25, "lines")
                    # aspect.ratio = 1
                ) +
                labs(caption = qq("N = @{end_idx} / @{n_analytes}")) +
                scale_x_discrete(labels = NULL)
            # bar_plot

            if (add_fc) {
                pos_ <- sliced_data %>%
                    group_by(Name_Coordinate, group) %>%
                    summarize(pos = max(!is.infinite(relative_signal), na.rm = TRUE) / 6, .groups = "drop") %>%
                    pull("pos")

                labs_ <- sliced_data %>%
                    mutate(labs = ifelse(significant, str_c(round(relative_signal, 1), " **"), round(relative_signal, 1))) %>%
                    .$labs

                bar_plot <- bar_plot +
                    geom_text(data = sliced_data, aes(
                        label = labs_,
                        y = pos_
                    ), size = 6, hjust = -0.5, angle = 90, color = "black")
            }
            bar_plot
            return(bar_plot)
        })
        return(bar_plot_lst)
    } else {
        all_analytes <- as.character(unique(short_group_data$Name_Coordinate))
        all_single_bar_plots_lst <- lapply(all_analytes, FUN = function(analyte) {
            # message(analyte)
            single_analyte_bar_plot <- ggplot(short_group_data %>% filter(Name_Coordinate == analyte),
                mapping = aes(x = short_group, y = relative_signal, fill = group, color = group)
            ) +
                geom_col(linewidth = rel(1.1), show.legend = TRUE) +
                theme_prism(base_size = 20) +
                scale_fill_manual(values = my_colors) +
                scale_color_manual(values = my_outline_colors) +
                labs(x = "Group", y = qq("Fold Change to @{my_control}")) +
                # ggtitle(title) +
                theme(
                    # strip.text = element_text(size = rel(strip_text_rel_size)),
                    panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
                    panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
                    panel.spacing = unit(0.5, "lines"),
                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                    strip.text = element_text(size = rel(strip_text_rel_size)),
                    axis.title = element_text(size = rel(2)),
                    plot.title = element_text(size = rel(2)),
                    axis.text = element_text(size = rel(1.75)),
                    axis.ticks = element_line(linewidth = rel(1.75)),
                    axis.line = element_line(linewidth = rel(2)),
                    aspect.ratio = 4 / 3
                ) +
                labs(title = analyte, caption = title) +
                scale_x_discrete(labels = NULL)

            if (add_fc) {
                pos_ <- short_group_data %>%
                    filter(Name_Coordinate == analyte) %>%
                    group_by(Name, group) %>%
                    summarize(pos = max(relative_signal) / 6, .groups = "drop") %>%
                    pull("pos")
                labs_ <- round(short_group_data %>% filter(Name_Coordinate == analyte) %>% .$relative_signal, 3)

                single_analyte_bar_plot <- single_analyte_bar_plot +
                    geom_text(aes(
                        label = labs_,
                        y = pos_
                    ), size = rel(15), hjust = -0.5, angle = 90, color = "black")
            }
            single_analyte_bar_plot
            return(single_analyte_bar_plot)
        })
        return(all_single_bar_plots_lst %>% set_names(all_analytes))
    }
}

#' @note convert log2 to raw fc fresh.
#' This is if you want to draw upper and lower bounds, removes anything in between +/- the threshold value in the log2 scale
#' @param log2_fc_thresh log2 change (which will be converted to Raw FC) on which to threshold
#' @return raw fc thresh, num vec
log2_fold_change_to_raw_fc <- function(log2_fc_thresh) {
    # message(qq("Calculating 2^(+/-{log2_fc}))..."))
    fc <- 2^log2_fc_thresh
    fc_recip <- 2^-log2_fc_thresh
    # message("Done!")
    return(list(fc, fc_recip) %>% set_names(c("fc", "fc_recip")))
}

#' @note plot the sorafenib fold change waterfall plot
#' @param data looks like wf_dat
#' @param add_fc boolean for whether to add the raw fold changes
#' @param TITLE NO FILTERING DONE - this is to NAME THE TITLE
#' @param main_cutoff if not NA, include main cutoff line -> this should be in raw FC format (will get converted to log2 automatically)
#' @param line_color the color of the filter line added to waterfall plot, in log2 space
plot_wf_graph <- function(data, add_fc = FALSE, TITLE = "Potential markers of sorafenib toxicity", main_cutoff = NULL, line_color = "purple") {
    final_data_wf_plot <- data %>%
        mutate(color = "black") %>%
        arrange(desc(relative_signal)) %>%
        # need to do this order again
        mutate(new_order = seq_len(n())) %>%
        mutate(Name = factor(Name, levels = Name[new_order])) %>%
        group_by(Name)

    num_analytes <- final_data_wf_plot %>%
        distinct(Name) %>%
        nrow()

    main_comp <- unique(final_data_wf_plot$main_group)
    control <- unique(final_data_wf_plot$control)

    new_title <- str_split(TITLE, "\n", simplify = TRUE)
    main_title <- stringr::str_wrap(new_title[1], width = 55)
    subtitle <- if (ncol(new_title) >= 2) stringr::str_wrap(new_title[2], width = 90) else ""
    helper_caption <- if (ncol(new_title) >= 3) new_title[3] else ""
    caption_text <- if (nzchar(helper_caption)) {
        qq("N = @{num_analytes}\n@{helper_caption}")
    } else {
        qq("N = @{num_analytes}")
    }

    if (nrow(final_data_wf_plot) > 1) {
        my_y_lim <- c(min(c(0, final_data_wf_plot$log_relative_signal, na.rm = TRUE)), max(final_data_wf_plot$log_relative_signal, na.rm = TRUE))
    } else {
        if (final_data_wf_plot$log_relative_signal > 0) {
            my_y_lim <- c(
                0,
                final_data_wf_plot$log_relative_signal + 1
            )
        } else {
            my_y_lim <- c(
                final_data_wf_plot$log_relative_signal - 1,
                0
            )
        }
    }

    main_plot <- ggplot(final_data_wf_plot, aes(
        x = Name, y = log_relative_signal,
        fill = log_relative_signal
    )) +
        geom_bar(
            stat = "identity",
            width = 0.825,
            show.legend = FALSE
        ) +
        scale_fill_gradient2(
            low = "blue", mid = "gray", high = "red",
            midpoint = mean(final_data_wf_plot$log_relative_signal, na.rm = TRUE),
            limits = c(
                min(final_data_wf_plot$log_relative_signal, na.rm = TRUE),
                max(final_data_wf_plot$log_relative_signal, na.rm = TRUE)
            )
        ) +
        ylim(my_y_lim) +
        xlab("") +
        ylab(TeX(str_c("$log_2(\\frac{", main_comp, "}{", control, "})$"))) +
        labs(title = main_title, subtitle = subtitle, caption = caption_text) +
        theme_prism(base_size = 16) +
        theme(
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            # panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
            axis.text.x = element_text(
                angle = 90, vjust = 0.2, hjust = 0.95,
                color = "black", size = rel(0.95)
            ),
            axis.title = element_text(size = rel(2)),
            plot.title = element_text(size = rel(2.1), lineheight = 1),
            plot.subtitle = element_text(size = rel(1.25), lineheight = 1),
            axis.text = element_text(size = rel(1.5)),
            axis.ticks = element_line(linewidth = rel(1.5)),
            axis.line = element_line(linewidth = rel(1.5)),
            plot.margin = margin(15, 15, 15, 15)
            # aspect.ratio = 1
        )

    if (add_fc) {
        sor_bar_text_direction <- sapply(final_data_wf_plot$log_relative_signal, FUN = function(x) ifelse(x < 0, -1, 1))
        main_plot <- main_plot +
            geom_text(aes(
                label = paste0(round(2^(final_data_wf_plot$log_relative_signal), 2)),
                y = 0.05 * sor_bar_text_direction
            ), size = rel(5), hjust = 0, angle = 90 * sor_bar_text_direction, color = "black")
    }

    if (!is.null(main_cutoff)) {
        main_plot <- main_plot +
            geom_hline(aes(yintercept = -log2(main_cutoff)),
                linetype = "dashed", lwd = rel(1.1), color = line_color
            ) +
            geom_hline(aes(yintercept = log2(main_cutoff)),
                linetype = "dashed", lwd = rel(1.1), color = line_color # "orange"
            )
    }
    main_plot
    return(main_plot)
}

#' @note make the wf ready dataset from df
#' @param data df from make_plot_ready_dataset()
#' @param ref_thresh_to_filter_ the signal threshold above which to keep analytes (to avoid analyze shades of nothing)
#' @param my_main_threshold the raw threshold on which to filter for significant fold change
#' @param comparisons A named list where names are control groups and values are vectors of groups to compare against each control.
#' @return lst object containing wf_dat and dat_filtered, for use in plot_wf_graph() and make_bar_charts()
make_wf_data <- function(data, ref_thresh_to_filter_ = 100, my_main_threshold = 2, comparisons = list("Group 1" = c("Group 2", "Group 3"))) {
    #' @note my_group_lvls is globally defined
    group_levels <- as.character(my_group_lvls)

    # Initialize an empty list to hold the results
    result_list <- list()
    # Loop through each control group and its corresponding groups to compare
    for (idx in seq_along(names(comparisons))) {
        # control <- names(comparisons)[4]
        control <- names(comparisons)[idx]
        compare_groups <- comparisons[[idx]]

        # Do all the computations like before, but filter based on the current control group
        dat_init <- data %>%
            mutate(Name = make.unique(Name)) %>%
            pivot_longer(all_of(group_levels), names_to = "group", values_to = "signal") %>%
            mutate(group = factor(group, levels = group_levels)) %>%
            group_by(Coordinate, Sname, Name) %>%
            mutate(
                relative_signal = signal / signal[group == control],
                log_relative_signal = log2(relative_signal)
            )

        # Create a list to hold the wf_dat for this control
        wf_dat_list <- list()
        dat_filtered_lst <- list()
        main_dat_filtered_lst <- list()

        for (compare_group in compare_groups) {
            # compare_group = compare_groups[1]
            comparison_levels <- unique(as.character(c(control, compare_group)))

            pair_dat_filtered <- dat_init %>%
                filter(group %in% c(control, compare_group)) %>%
                mutate(group = factor(group, levels = comparison_levels)) %>%
                group_by(Coordinate, Sname, Name) %>%
                filter(all(signal > ref_thresh_to_filter_)) %>%
                ungroup() %>%
                filter(!(startsWith(Name, "Reference"))) %>%
                filter(!(startsWith(Name, "Negative"))) %>%
                arrange(Name, group)

            pair_sig_names <- pair_dat_filtered %>%
                mutate(significant = abs(log_relative_signal) > log2(my_main_threshold)) %>%
                distinct(Name, group, significant)

            pair_dat_filtered <- pair_dat_filtered %>%
                left_join(pair_sig_names, by = c("Name", "group")) %>%
                arrange(Name, group) %>%
                mutate(control = control)

            wf_dat <- pair_dat_filtered %>%
                filter(group == compare_group) %>%
                arrange(desc(log_relative_signal)) %>%
                ungroup() %>%
                mutate(order = seq_len(n())) %>%
                mutate(Name = factor(Name, levels = Name[order])) %>%
                arrange(Name) %>%
                mutate(main_group = compare_group)

            dat_filtered <- pair_dat_filtered %>%
                filter(group == compare_group | group == control) %>%
                mutate(main_group = compare_group)

            wf_dat_list[[compare_group]] <- wf_dat
            dat_filtered_lst[[compare_group]] <- dat_filtered
            main_dat_filtered_lst[[compare_group]] <- pair_dat_filtered
        }

        main_dat_filtered <- bind_rows(main_dat_filtered_lst) %>%
            distinct() %>%
            mutate(group = factor(group, levels = unique(as.character(c(control, compare_groups))))) %>%
            arrange(Name, group)

        result_name <- str_c(control, " vs ", str_c(compare_groups, collapse = ", "))

        # Add the wf_dat and dat_filtered to the result list, keyed by the control name
        result_list[[result_name]] <- list("wf_dat" = wf_dat_list, "dat_filtered" = dat_filtered_lst, "main_dat_filtered" = main_dat_filtered)
    }
    return(result_list)
}

#' Save Bar Plots to Disk with Optional Faceting
#'
#' This function saves either a list of bar plots or a single faceted bar plot to the specified directory.
#' If a list of bar plots is provided, the function will save each individual plot in parallel.
#' If the `faceted` option is set to TRUE, it saves a faceted bar plot.
#'
#' @param barplots_obj A list of ggplot objects or a single ggplot object. Each ggplot object is a bar plot that needs to be saved.
#' @param barplot_dat A data frame that contains the data used for generating the bar plot(s). It's used for calculating dimensions.
#' @param output_path The directory where the saved plots will be stored.
#' @param toggle Optional; A string that can be set to filter or categorize plots. Default is "all".
#' @param faceted Logical; whether the bar plot is faceted. Default is FALSE.
#' @param groups_per_page Integer; the number of groups to display per page in the faceted plot. Default is 25.
#' @param top_name String; a title or identification string for the major comparison being performed.
#'
#' @return Invisible NULL. The function is called for its side effect of saving plots to disk.
#'
#' @examples
#' # Example usage with a list of ggplot objects
#' save_list_bar_plots(list_of_ggplots, barplot_data, "/path/to/save/plots")
#'
#' # Example usage with a single faceted ggplot object
#' save_list_bar_plots(single_ggplot, barplot_data, "/path/to/save/plots", faceted = TRUE)
#'
#' @note The function uses parallel processing for saving a list of bar plots.
#' The number of cores used is the minimum of one less than the available cores and the length of the list of bar plots.
#'
save_list_bar_plots <- function(barplots_obj, barplot_dat, output_path, toggle = "all", faceted = FALSE, groups_per_page = 25, top_name = "Major comparison string") {
    if (!faceted) {
        # message("Found bar list, saving individual analyte bar plots...")
        idx_bar_lst <- seq_along(barplots_obj)

        available_workers <- suppressWarnings(as.integer(future::availableCores()[[1]]))
        if (!is.finite(available_workers)) {
            available_workers <- parallel::detectCores()
        }
        if (!is.finite(available_workers)) {
            available_workers <- 1
        }
        num_workers <- max(1, min(4, available_workers, length(idx_bar_lst)))
        use_parallel <- num_workers > 1
        if (use_parallel) {
            plan(multisession, workers = num_workers)
            on.exit(plan(sequential), add = TRUE)
        }

        #' @note convenience function for the parallelization
        future_walk_pbar <- function(sub_bar_x, toggle) { #  p_bar,
            # libraries
            options(tidyverse.quiet = TRUE)
            library(tidyverse)
            library(GetoptLong)
            library(latex2exp)
            library(ggprism)
            library(ggh4x)

            # save plot with tryCatch to find errors
            tryCatch(
                {
                    analyte_name <- names(barplots_obj)[sub_bar_x]
                    # message(analyte_name)
                    ggsave(
                        filename = file.path(output_path, GetoptLong::qq("@{analyte_name}-barplot_@{toggle}.png")),
                        plot = barplots_obj[[sub_bar_x]],
                        height = 10,
                        width = 10
                    )
                },
                error = function(e) {
                    message("Error in index: ", sub_bar_x, ". Error message: ", e$message)
                }
            )
            # Update progress
            # p_bar(sprintf("Processing index %s", analyte_name))
        }

        # # hprogress handler info
        # handlers(handler_progress(
        #     format   = ":spin [:bar] :current/:total :percent in :elapsed ETA: :eta (:message) ",
        #     width    = 100,
        #     complete = "+"
        # ))

        # parallelization
        with_progress({
            # Initialize progress
            # p_bar <- progressor(steps = length(idx_bar_lst))
            if (use_parallel) {
                furrr::future_walk(
                    .x = as.list(idx_bar_lst), .f = future_walk_pbar,
                    # p_bar = p_bar,
                    toggle = toggle,
                    .options = furrr::furrr_options(seed = TRUE)
                )
            } else {
                purrr::walk(as.list(idx_bar_lst), future_walk_pbar, toggle = toggle)
            }
        })

        if (use_parallel) {
            plan(sequential)
        }
    } else {
        total_pages <- length(barplots_obj)
        for (i in seq_len(total_pages)) {
            ggsave(
                filename = file.path(output_path, qq("@{top_name}-barplot_page_@{i}__@{toggle}.png")),
                plot = barplots_obj[[i]],
                height = groups_per_page,
                width = groups_per_page
            )
        }
    }
}

#' @note Test whether the n is a perfect square
#' @param n any numeric
#' @return boolean
is_perfect_square <- function(n) {
    if (n < 0) {
        return(FALSE)
    }
    sqrt_n <- as.integer(sqrt(n))
    return(sqrt_n * sqrt_n == n)
}

#' Execute Analysis Pipeline and Generate Plots
#'
#' This function serves as the main engine for executing the entire analysis pipeline.
#' It performs data transformations, generates waterfall plots and bar graphs,
#' and saves these visualizations to disk.
#'
#' @param df Data frame containing the raw data for analysis.
#' @param ref_thresh_to_filter Numeric, reference threshold for filtering the data.
#' @param main_threshold Numeric, main fold-change threshold for significance.
#' @param output_dir_full_path String, full path to the main-analysis output
#'   directory for one user-owned analysis tree.
#' @param groups_per_page Integer, number of groups per page for faceted plots. Default is 25.
#' @param comparisons List, a named list specifying the groups to compare. Default compares "Ctrl Preg" with "sFlt1 Preg" and "Ctrl PP".
#'
#' @return A list containing two named lists: 'wf_plots' for waterfall plots and 'bar_plots' for bar graphs.
#'
#' @examples
#' # Example usage
#' result <- make_graphs(df, ref_thresh_to_filter = 100, main_threshold = 1.5, output_dir_full_path = "path/to/output", comparisons = list("GroupA" = c("GroupB", "GroupC")))
#'
#' @note This function calls other helper functions like 'make_wf_data', 'plot_wf_graph', and 'make_bar_charts'.
#' It also utilizes parallel processing for efficiency.
#'
#' @seealso \code{\link{make_wf_data}}, \code{\link{plot_wf_graph}}, \code{\link{make_bar_charts}}, \code{\link{save_list_bar_plots}}
#'
make_graphs <- function(df, ref_thresh_to_filter = NA, main_threshold = NA, output_dir_full_path = NA, groups_per_page = 25, comparisons = list("Ctrl Preg" = c("sFlt1 Preg", "Ctrl PP"))) {
    # `output_dir_full_path` should be the script-specific output directory for
    # the main analysis. This function then creates one subfolder per
    # reference-threshold / fold-change-threshold combination.
    ref_threshold_dir <- file.path(output_dir_full_path, qq("ref_threshold_@{ref_thresh_to_filter}"))
    threshold_hits_dir <- file.path(ref_threshold_dir, "fold_change_hits", qq("threshold_@{main_threshold}"))

    # make wf data
    lst_obj <- make_wf_data(
        data = df,
        ref_thresh_to_filter_ = ref_thresh_to_filter,
        my_main_threshold = main_threshold,
        comparisons = comparisons
    )
    # str(lst_obj, max.level = 2)

    transposed_lst_obj <- transpose(lst_obj)$wf_dat
    num_comparisons <- seq_along(transposed_lst_obj)
    # message("\nMaking waterfall plots...")

    lst_wf_plots <- map(num_comparisons, .f = function(i) {
        lst <- transposed_lst_obj[[i]]
        top_name <- names(transposed_lst_obj)[i]
        # message(top_name)

        # each comparison might have sub comparisons
        num_sub_comps <- seq_along(lst)
        mini_lst <- map(.x = num_sub_comps, .f = function(xz) {
            bot_name <- names(lst)[xz]
            sub_wf_df <- lst[[xz]]
            # message(str_c("     ", bot_name))

            if (nrow(sub_wf_df) == 0) {
                message("Sub df was empty, skipping...")
                return(list(tibble(), tibble()) %>% set_names("all", "significant"))
            }
            sub_path <- file.path(ref_threshold_dir, "all_comparisons", "waterfalls", top_name)
            dir.create(sub_path, showWarnings = FALSE, recursive = TRUE)

            # all data
            plot_res <- plot_wf_graph(
                data = sub_wf_df,
                add_fc = FALSE,
                TITLE = qq("@{top_name}\n@{bot_name} comparison on All Data\nRaw Signal > @{ref_thresh_to_filter}"),
                main_cutoff = main_threshold,
                line_color = "#d900ff"
            )
            # Calculate the dimensions
            all_height <- with(sub_wf_df, {
                18 + as.integer(abs((max(log_relative_signal, na.rm = TRUE) %/% 1.5))) + as.integer(abs((min(log_relative_signal, na.rm = TRUE) %/% 1.5)))
            })
            all_width <- 20 + nrow(sub_wf_df) %/% 10
            # Save it
            ggsave(
                filename = file.path(sub_path, qq("@{top_name}--@{bot_name}-main_waterfall-all.png")),
                plot = plot_res,
                height = all_height, width = all_width
            )
            write_tsv(sub_wf_df, file.path(sub_path, qq("@{top_name}--@{bot_name}-main_waterfall-all.tsv")))

            # filtered wf
            sig_sub_path <- file.path(threshold_hits_dir, "waterfalls", top_name)
            dir.create(sig_sub_path, showWarnings = FALSE, recursive = TRUE)

            filtered_sub_wf_df <- sub_wf_df %>% filter(significant)
            if (nrow(filtered_sub_wf_df) == 0) {
                message("Sig sub df was empty, skipping...")
                return(list(tibble(), tibble()) %>% set_names("all", "significant"))
            }
            plot_res_significant <- plot_wf_graph(
                data = filtered_sub_wf_df,
                add_fc = FALSE,
                TITLE = qq("@{top_name}\n@{bot_name} comparison on Significant Data\nRaw Signal > @{ref_thresh_to_filter} || (FC > @{main_threshold} and FC < @{round(1 / main_threshold, 3)})"),
                main_cutoff = NULL,
                line_color = NA
            )

            # Calculate the dimensions
            sig_height <- with(filtered_sub_wf_df, {
                15 + as.integer(abs((max(log_relative_signal, na.rm = TRUE) %/% 1.5))) + as.integer(abs((min(log_relative_signal, na.rm = TRUE) %/% 1.5)))
            })
            sig_width <- 25 + nrow(filtered_sub_wf_df) %/% 10
            # Save it
            ggsave(
                filename = file.path(sig_sub_path, qq("@{top_name}--@{bot_name}-main_waterfall-significant.png")),
                plot = plot_res_significant,
                height = sig_height, width = sig_width
            )
            write_tsv(filtered_sub_wf_df, file.path(sig_sub_path, qq("@{top_name}--@{bot_name}-main_waterfall-significant.tsv")))

            return(list(plot_res, plot_res_significant) %>% set_names("all", "significant"))
        }) %>%
            suppressWarnings() %>%
            set_names(names(lst))

        return(mini_lst)
    }) %>% set_names(names(transposed_lst_obj))

    # now filter down the dataframes and make bar graph

    # message("\nMaking bar graphs plots...")
    transposed_lst_dat_obj <- transpose(lst_obj)$dat_filtered
    transposed_lst_dat_obj_MAIN <- transpose(lst_obj)$main_dat_filtered
    # loop through
    lst_of_dat_filtered <- map(seq_along(lst_obj), .f = function(ix) {
        top_name <- names(transposed_lst_dat_obj)[ix]
        # message(top_name)

        # make one plot of the entire thing, all comparisons, only the significant ones
        main_sig_data_plot <- transposed_lst_dat_obj_MAIN[[ix]] %>%
            group_by(Name) %>%
            filter(any(significant))

        if (nrow(main_sig_data_plot) == 0) {
            message("Sub df was empty, skipping...")
            return(list(tibble(), tibble()) %>% set_names("all", "significant"))
        }

        main_path <- file.path(threshold_hits_dir, "barplots", top_name, "combined")
        dir.create(main_path, showWarnings = FALSE, recursive = TRUE)
        # message("Plotting significant analytes in multi-page faceted plot...")
        barplot_main_plot <- make_bar_charts(
            data = main_sig_data_plot,
            add_fc = TRUE,
            title = qq("All passing base threshold raw signal > @{ref_thresh_to_filter}\nMain FC threshold >@{main_threshold} and <@{round(1 / main_threshold, 3)}"),
            strip_text_rel_size = 0.99,
            make_facet_plot = TRUE,
            caption = TRUE,
            groups_per_page = groups_per_page
        )

        save_list_bar_plots(
            barplots_obj = barplot_main_plot, barplot_dat = main_sig_data_plot, output_path = main_path,
            toggle = "significant", faceted = TRUE,
            groups_per_page = groups_per_page,
            top_name = top_name
        )


        # back to the individual comparisons
        lst <- transposed_lst_dat_obj[[ix]]
        # each comparison might have sub comparisons
        # loop through the sub conditions
        num_sub_comps <- seq_along(lst)
        mini_lst <- map(num_sub_comps, .f = function(xz) {
            bot_name <- names(lst)[xz]
            sub_bar_plot_df <- lst[[xz]]
            # message(str_c("      ", bot_name))
            if (nrow(sub_bar_plot_df) == 0) {
                message("Sub df was empty, skipping...")
                return(list(tibble(), tibble()) %>% set_names("all", "significant"))
            }

            # # all --> this is probably not what you want
            # sub_path <- file.path(ref_threshold_dir, "all_comparisons", "barplots", top_name)

            # dir.create(sub_path, showWarnings = FALSE, recursive = TRUE)

            # # message("Plotting all bar plots in individual pngs...")
            # barplot_all_data <- make_bar_charts(
            #     data = sub_bar_plot_df,
            #     add_fc = TRUE,
            #     title = qq("All passing base threshold signal > @{ref_thresh_to_filter}"),
            #     strip_text_rel_size = 0.99,
            #     make_facet_plot = FALSE,
            #     caption = FALSE,
            #     groups_per_page = groups_per_page
            # )

            # save_list_bar_plots(
            #     barplots_obj = barplot_all_data, barplot_dat = sub_bar_plot_df,
            #     output_path = sub_path, toggle = "all",
            #     faceted = FALSE, top_name = top_name
            # )

            # significant
            sub_filtered_bar_plot_df <- sub_bar_plot_df %>%
                group_by(Name) %>%
                filter(any(significant))

            if (nrow(sub_filtered_bar_plot_df) == 0) {
                message("Sub df was empty, skipping...")
                return(list(tibble(), tibble()) %>% set_names("all", "significant"))
            }


            main_path <- file.path(threshold_hits_dir, "barplots", top_name, bot_name, "combined")
            dir.create(main_path, showWarnings = FALSE, recursive = TRUE)
            # message("Plotting significant analytes in multi-page faceted plot...")
            barplot_main_plot <- make_bar_charts(
                data = main_sig_data_plot %>% filter(Name %in% unique(sub_filtered_bar_plot_df$Name)),
                add_fc = TRUE,
                title = qq("All passing base threshold raw signal > @{ref_thresh_to_filter}\nMain FC threshold >@{main_threshold} and <@{round(1 / main_threshold, 3)}"),
                strip_text_rel_size = 3,
                make_facet_plot = TRUE,
                caption = TRUE,
                groups_per_page = groups_per_page
            )

            save_list_bar_plots(
                barplots_obj = barplot_main_plot, barplot_dat = main_sig_data_plot, output_path = main_path,
                toggle = "significant", faceted = TRUE,
                groups_per_page = groups_per_page,
                top_name = str_c(top_name, bot_name, sep = "--")
            )



            # get path
            sub_sig_path <- file.path(threshold_hits_dir, "barplots", top_name, bot_name, "individual")
            dir.create(sub_sig_path, showWarnings = FALSE, recursive = TRUE)

            # message("Plotting all significant plots in individual pngs...")
            barplot_sig_data <- make_bar_charts(
                data = sub_filtered_bar_plot_df,
                add_fc = TRUE,
                title = qq("Passing base threshold raw signal > @{ref_thresh_to_filter}\nMain FC threshold >@{main_threshold} and <@{round(1 / main_threshold, 3)}"),
                strip_text_rel_size = 0.98,
                make_facet_plot = FALSE,
                caption = FALSE,
                groups_per_page = groups_per_page
            )
            save_list_bar_plots(
                barplots_obj = barplot_sig_data,
                barplot_dat = sub_filtered_bar_plot_df,
                output_path = sub_sig_path, toggle = "significant",
                faceted = FALSE, top_name = top_name
            )

            return(list(barplot_main_plot, barplot_sig_data) %>% set_names("all", "significant"))
        })
        return(mini_lst)
    })

    return(list(lst_wf_plots, lst_of_dat_filtered) %>% set_names(c("wf_plots", "bar_plots")))
}
