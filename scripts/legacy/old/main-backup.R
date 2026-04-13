rm(list = ls())
library(tidyverse)
# library(ComplexHeatmap)
library(readxl)
library(GetoptLong)
library(RColorBrewer)
library(latex2exp)
library(patchwork)
library(ggprism)
library(ggh4x)

# protocol data
# info_fn <- "output/Mouse Angiogenesis Array Kit - Protocol.xlsx"
info_fn <- "output/cytoXL array kit - protocol.xlsx"

# Signal column is the background corrected signal
data_dir <- "/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/projects/VEGFRi and Dox/in-vivo projects/vehicle vs sor dox lis/protein arrays/Veh vs Sor Dox Lis - Cytokine XL/data"

# read in analyte info
analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

# output_dir <- "output/plots/nick/angiogenesis_array"
output_dir <- "output/plots/nick/cytokine_xl_array"

vehicle_point_color <- "#ffffff"
sorafenib_point_color <- "#737171"
sor_dox_point_color <- "#9BB3D3"
sor_lis_point_color <- "#B53530"
my_colors <- c(vehicle_point_color, sorafenib_point_color, sor_dox_point_color, sor_lis_point_color) %>%
    set_names(c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis"))

#' @note takes as import the data directory of the LICOR data and magically makes the dataset
make_plot_ready_dataset <- function(data_dir, preview = FALSE) {
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
    df_temp <- tibble(data_fns = dir(data_dir, full.names = TRUE)) %>%
        mutate(
            group_temp = str_split(basename(data_fns), pattern = "\\.", simplify = TRUE)[, 1],
            group_fn = str_to_title(str_split(group_temp, " - ", simplify = TRUE)[, 1]), .before = 1
        )

    # grab each individual df and modify to our needs
    init_df <- df_temp %>%
        dplyr::select(-group_temp) %>%
        mutate(group_fn = factor(group_fn, levels = c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis"))) %>%
        arrange(group_fn) %>%
        mutate(lst_df = map2(data_fns, group_fn, .f = function(x, y) {
            read_excel(x, skip = 3) %>%
                mutate(group = y) %>%
                rename(Sname = Name) %>%
                dplyr::select(Sname, Signal, group) %>% # could use SNR
                filter(!startsWith(Sname, "B")) %>% # remove the background spot because we've already taken it into account
                mutate(sname_grouping = rep(1:(nrow(.) / 2), each = 2)) %>%
                group_by(sname_grouping, group) %>%
                summarize(
                    Sname = str_c(Sname, collapse = ", "),
                    Sname_signal_pair = str_c(round(as.numeric(Signal), 2), collapse = ", "),
                    signal = mean(as.numeric(Signal), na.rm = TRUE), .groups = "drop"
                )
        })) %>%
        dplyr::select(-data_fns, -group_fn) %>%
        unnest(lst_df) %>%
        group_by(group) %>%
        left_join(analyte_info_df, by = join_by(sname_grouping))

    if (preview) init_df %>% print(n = 50)

    # get into wide format for compatibility with the rest of the script
    my_df <- init_df %>%
        dplyr::select(Name, Sname, Coordinate, group, signal) %>%
        pivot_wider(id_cols = c(Name, Sname, Coordinate), names_from = group, values_from = signal) %>%
        mutate(Name = make.unique(Name)) %>%
        mutate(Coordinate = convert_text(Coordinate))
    return(my_df)
}

# manually identified threshold
#' @note takes wide df and calculates the mean intensity of a set of coordinates
find_filter_thresh <- function(wide_df, ref_coords = c("A5,6")) {
    cols_to_pivot <- which(!(colnames(wide_df) %in% c("Name", "Sname", "Coordinate")))
    long_df <- wide_df %>%
        pivot_longer(all_of(cols_to_pivot), names_to = "group", values_to = "signal")

    # `Signal Level` <- long_df$signal
    # hist_info <- hist(`Signal Level`, breaks = 100, plot = TRUE)

    # Now, hist_info$breaks contains the break points (edges of the bins)
    # and hist_info$counts contains the counts within each bin.
    # print(hist_info$breaks)
    # print(hist_info$counts)
    filtered_long_df <- long_df %>%
        filter(Coordinate %in% ref_coords)

    ref_thresh_to_filter_df <- filtered_long_df %>%
        ungroup() %>%
        summarize(
            mean = mean(signal, na.rm = TRUE),
            sd = sd(signal, na.rm = TRUE)
        )

    # print(filtered_long_df, n = Inf)
    boxp <- ggplot(filtered_long_df, mapping = aes(x = Name, y = signal, fill = group, group = Name)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(position = position_jitter(0.1), size = rel(2.5), shape = 21) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45), panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        scale_fill_manual(values = my_colors) +
        labs(
            title = "Boxplot of analyte signals in region across groups",
            subtitle = qq("region coords: @{str_c(ref_coords, collapse = '; ')}"),
            caption = qq("mean = @{round(ref_thresh_to_filter_df$mean, 3)}\nsd = @{round(ref_thresh_to_filter_df$sd, 3)}")
        )


    histop <- ggplot(long_df, mapping = aes(x = signal)) +
        geom_histogram(fill = "black", position = "identity", binwidth = 500, show.legend = FALSE) +
        scale_fill_manual(values = my_colors) +
        geom_histogram(
            data = long_df %>% filter(signal < ref_thresh_to_filter_df$mean),
            mapping = aes(fill = group),
            color = "black", binwidth = 500,
            show.legend = FALSE
        ) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45), panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        labs(
            title = "Proportion of analytes removed",
            subtitle = qq("Mean signal intensity threshold = @{round(ref_thresh_to_filter_df$mean)}")
        )
    # histop

    histop2 <- ggplot(long_df %>% filter(signal > ref_thresh_to_filter_df$mean),
        mapping = aes(x = signal)
    ) +
        scale_fill_manual(values = my_colors) +
        geom_histogram(
            mapping = aes(fill = group),
            color = "black", binwidth = 500,
        ) +
        theme_prism(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45), panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
        ) +
        labs(
            title = "Proportion of analytes remaining after gating",
            subtitle = qq("Mean signal intensity threshold = @{round(ref_thresh_to_filter_df$mean)}")
        )

    # histop2

    plt_lst <- boxp / (histop + histop2)

    return(list(ref_thresh_to_filter_df, plt_lst))
}


## Run it
df <- make_plot_ready_dataset(data_dir, preview = FALSE)
df

## Identify the threshold to filter out analytes
# ref_coords_to_make_filter <- c("A5,6", "B3,4", "B5,6", "C3,4", "C5,6", "D3,4", "D5,6")
ref_coords_to_make_filter <- c("H1,2", "H3,4", "G1,2", "G3,4", "F9,10", "F11,12", "E9,10", "E11,12")
filt_lst <- find_filter_thresh(wide_df = df, ref_coords = ref_coords_to_make_filter)

ggsave(
    filename = file.path(output_dir, "region_stats.png"), plot = filt_lst[[2]],
    width = 10, height = 10
)

# this can be played with
my_sor_threshold <- 0.4
my_sor_threshold2 <- 0.6

# this should be decided above
# # angiogenesis
# my_ref_thresh_to_filter <- 60
# my_ref_thresh_to_filter2 <- 100

# cyto XL
my_ref_thresh_to_filter <- 100
my_ref_thresh_to_filter2 <- 150


# sizes of the graphs...
my_size_tibble <- tribble(
    ~Name, ~width, ~height,
    "sor_plot", 20, 15,
    "bar_plot_big", 20, 20,
    "bar_plot_small", 10, 10
)


###########################################

# will go into the 'done' argument
already_performed_ELISAs <- tibble(Name = c(
    "MIP-1a", "MCP-1", "PD-ECGF",
    "HGF", "Serpin E1", "Prolactin",
    "PDGF-AA"
), color = "#1fc61f")

# ref thresh to filter calculated below
make_graphs <- function(
    df, ref_thresh_to_filter = NA, sor_threshold = NA, done = NA,
    output_dir_full_path = NA, size_tibble = NA) {
    # artifact replicates removed manually before this
    dat_init <- df %>%
        mutate(Name = make.unique(Name)) %>% # shouldn't hurt to do this again...
        pivot_longer(Vehicle:`Sor + Lis`, names_to = "group", values_to = "signal") %>%
        mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis"))) %>%
        group_by(Coordinate, Sname, Name) %>%
        mutate(
            relative_signal = signal / signal[1L],
            log_relative_signal = log2(relative_signal),
            sor_relative_signal = signal / signal[2L],
            sor_log_relative_signal = log2(sor_relative_signal)
        )

    dat_filtered <- dat_init %>%
        filter(all(signal > ref_thresh_to_filter)) %>%
        filter(!(Name %in% c("Reference spots", "Reference spots.1", "Reference spots.2", "Negative control")))

    # create the data to plot, order it
    filtered_sorafenib_wf <- dat_filtered %>%
        filter(group == "Sorafenib") %>%
        arrange(desc(log_relative_signal)) %>%
        ungroup() %>%
        mutate(order = seq_len(n())) %>%
        mutate(Name = factor(Name, levels = Name[order])) %>%
        arrange(Name)
    # filtered_sorafenib_wf %>% print(n = Inf)

    soraf_disp_names_filtered <- subset(filtered_sorafenib_wf, abs(log_relative_signal) > sor_threshold)$Name

    # plot it
    make_sor_plot <- function(data, threshold, add_fc = FALSE, done_df = done) {
        final_data_wf_plot <- data %>%
            left_join(done_df, by = "Name") %>%
            mutate(color = ifelse(is.na(color), "black", color)) %>%
            arrange(order) %>%
            mutate(Name = factor(Name, levels = Name[order]))

        sor_plot <- ggplot(final_data_wf_plot, aes(
            x = Name, y = log_relative_signal,
            fill = log_relative_signal
        )) +
            geom_bar(
                stat = "identity",
                show.legend = FALSE
            ) +
            theme_prism(base_size = 16) +
            scale_fill_gradient2(
                low = "blue", mid = "gray", high = "red",
                midpoint = mean(final_data_wf_plot$log_relative_signal, na.rm = TRUE),
                limits = c(
                    min(final_data_wf_plot$log_relative_signal, na.rm = TRUE),
                    max(final_data_wf_plot$log_relative_signal, na.rm = TRUE)
                )
            ) +
            ylim(c(min(final_data_wf_plot$log_relative_signal, na.rm = TRUE) - 0.5, max(final_data_wf_plot$log_relative_signal, na.rm = TRUE) + 0.5)) +
            xlab("") +
            ylab(TeX("$log_2(\\frac{sor}{vehicle})$")) +
            # geom_hline(aes(yintercept = sor_threshold), linewidth = rel(0.5), linetype = 3) +
            # geom_hline(aes(yintercept = -sor_threshold), linewidth = rel(0.5), linetype = 3) +
            ggtitle(qq("Potential Biomarkers of Sorafenib Treatment Response\n(Filtered with signal > @{ref_thresh_to_filter})")) +
            # labs(caption = qq("Threshold fold change >= +/-log2(@{sor_threshold})\nUpper fold change = @{round(2^sor_threshold, 3)}\nLower fold change = @{round(abs(2^(-sor_threshold)), 3)}")) +
            theme(
                panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
                # panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
                axis.text.x = element_text(
                    angle = 90, vjust = 0.2, hjust = 0.95,
                    color = final_data_wf_plot$color
                ),
                # plot.margin = margin(0.5, 0.75, 0.5, 0.5, "cm"),
                axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
                plot.title = element_text(size = rel(1.25))
            )

        if (add_fc) {
            sor_bar_text_direction <- sapply(final_data_wf_plot$log_relative_signal, FUN = function(x) ifelse(x < 0, -1, 1))
            sor_plot <- sor_plot +
                geom_text(aes(
                    label = paste0(round(2^(final_data_wf_plot$log_relative_signal), 2)),
                    y = 0.05 * sor_bar_text_direction
                ), size = 4, hjust = 0, angle = 90 * sor_bar_text_direction, color = "black")
        }
        sor_plot
        return(sor_plot)
    }

    sor_plot_filtered <- make_sor_plot(data = filtered_sorafenib_wf, threshold = sor_threshold, add_fc = TRUE)
    # sor_plot_filtered
    ggsave(
        filename = file.path(output_dir_full_path, qq("filtered_sor_waterfall-thresh@{ref_thresh_to_filter}.png")),
        plot = sor_plot_filtered,
        height = size_tibble %>%
            filter(Name == "sor_plot") %>%
            pluck("height"),
        width = size_tibble %>%
            filter(Name == "sor_plot") %>%
            pluck("width")
    )

    dat_filtered_by_sor <- dat_filtered %>%
        filter(Name %in% soraf_disp_names_filtered)

    #' @note make the bar plots
    make_bar_charts <- function(data, title = "", add_fc = TRUE, done_df = done, strip_text_rel_size = 1.1) {
        helper_dat <- tibble(
            group = c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis"),
            short_group = factor(c("V", "S", "D", "L"), levels = c("V", "S", "D", "L"))
        )

        short_group_data <- left_join(data, helper_dat, by = "group") %>%
            left_join(done_df, by = "Name") %>%
            mutate(color = ifelse(is.na(color), "black", color)) %>%
            arrange(Name) %>%
            group_by(group) %>%
            mutate(order = seq_len(n()), Name = factor(Name, levels = unique(Name[unique(order)]))) %>%
            mutate(Name_Coordinate = str_c(Name, " (", Coordinate, ")"))

        my_strip_colors_dat <- short_group_data %>%
            ungroup() %>%
            distinct(Name, color) %>%
            pull("color")

        my_strip_colors <- strip_themed(text_x = elem_list_text(color = my_strip_colors_dat))

        stopifnot(nrow(short_group_data) > 0)

        bar_plot <- ggplot(short_group_data,
            mapping = aes(x = short_group, y = relative_signal, fill = group)
        ) +
            geom_col(color = "black") +
            # facet_wrap(~Name,
            #     scales = "free_y",
            # ) +
            ggh4x::facet_wrap2(~Name_Coordinate,
                scales = "free_y",
                strip = my_strip_colors,
                labeller = labeller(Name_Coordinate = label_wrap_gen(12))
            ) +
            theme_prism(base_size = 16) +
            scale_fill_manual(values = my_colors) +
            labs(x = "Group", y = "Fold Change to Vehicle") +
            ggtitle(title) +
            theme(
                strip.text = element_text(size = rel(strip_text_rel_size)),
                panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
                panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25))
            )

        if (add_fc) {
            pos_ <- short_group_data %>%
                group_by(Name, group) %>%
                summarize(pos = max(relative_signal) / 6, .groups = "drop") %>%
                pull("pos")
            labs_ <- round(short_group_data$relative_signal, 3)

            bar_plot <- bar_plot +
                geom_text(aes(
                    label = labs_,
                    y = pos_
                ), size = 4, hjust = -0.5, angle = 90, color = "black")
        }
        bar_plot
    }

    barplot_all_data <- make_bar_charts(
        data = dat_filtered,
        add_fc = TRUE,
        title = qq("All passing base threshold signal > @{ref_thresh_to_filter}"),
        strip_text_rel_size = 1.01
    )
    ggsave(
        filename = file.path(output_dir_full_path, qq("barplot_all_passing-base_thresh@{ref_thresh_to_filter}-sor_thresh_log2(@{sor_threshold}).png")),
        plot = barplot_all_data,
        height = size_tibble %>%
            filter(Name == "bar_plot_big") %>%
            pluck("height"),
        width = size_tibble %>%
            filter(Name == "bar_plot_big") %>%
            pluck("width")
    )

    barplot_filtered_data <- make_bar_charts(
        data = dat_filtered_by_sor,
        add_fc = TRUE,
        title = qq("Passing base threshold signal > @{ref_thresh_to_filter}\nand SOR FC threshold >@{round(2^sor_threshold, 3)} and <@{round(2^-sor_threshold, 3)}"),
        strip_text_rel_size = 0.99
    )
    ggsave(
        filename = file.path(output_dir_full_path, qq("barplot_significant_passing-thresh@{ref_thresh_to_filter}-sor_thresh_log2(@{sor_threshold}).png")),
        plot = barplot_filtered_data,
        height = size_tibble %>%
            filter(Name == "bar_plot_small") %>%
            pluck("height"),
        width = size_tibble %>%
            filter(Name == "bar_plot_small") %>%
            pluck("width")
    )
}

make_graphs(df,
    ref_thresh_to_filter = my_ref_thresh_to_filter,
    sor_threshold = my_sor_threshold,
    done = already_performed_ELISAs,
    output_dir_full_path = output_dir,
    size_tibble = my_size_tibble
)

make_graphs(df,
    ref_thresh_to_filter = my_ref_thresh_to_filter,
    sor_threshold = my_sor_threshold2,
    done = already_performed_ELISAs,
    output_dir_full_path = output_dir,
    size_tibble = my_size_tibble
)

make_graphs(df,
    ref_thresh_to_filter = my_ref_thresh_to_filter2,
    sor_threshold = my_sor_threshold,
    done = already_performed_ELISAs,
    output_dir_full_path = output_dir,
    size_tibble = my_size_tibble
)

make_graphs(df,
    ref_thresh_to_filter = my_ref_thresh_to_filter2,
    sor_threshold = my_sor_threshold2,
    done = already_performed_ELISAs,
    output_dir_full_path = output_dir,
    size_tibble = my_size_tibble
)






# OLD
# angiogenesis project
# home_dir <- "/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/projects/VEGFRi and Dox/in-vivo projects/vehicle vs sor dox lis/protein arrays"
# data_dir <- "Veh vs Sor Dox Lis - Angiogenesis Protein Array/array results/array pptx and excel"
# fn <- file.path(home_dir, data_dir, "angiogenesis_analysis-v2.xlsx")
# df <- read_excel(fn,
#     sheet = 3
# ) %>%
#     # remove blank rows
#     na.omit()

# ## SOR DOX
# interesting_soraf_targets <- soraf_disp_names %>% droplevels()

# dat_sordoxlis_init <- dat %>%
#     filter(group == "Sor + Dox") %>%
#     filter(Name %in% interesting_soraf_targets)

# sor_wf_for_doxlis <- dat_sordoxlis_init %>%
#     arrange(desc(log_relative_signal)) %>%
#     ungroup() %>%
#     mutate(order = seq_len(n())) %>%
#     mutate(Name = factor(Name, levels = Name[order])) %>%
#     arrange(Name)

# sordox_disp_names <- subset(sor_wf_for_doxlis, abs(log_relative_signal) <= thresh_dox)$Name %>% droplevels()
# sordox_disp_names

# ## SOR LIS

# dat_sorlis_init <- dat %>%
#     filter(group == "Sor + Lis") %>%
#     filter(Name %in% interesting_soraf_targets)

# sorlis_disp_names <- subset(sorlis_wf, abs(log_relative_signal) <= thresh_lis)$Name

# sorlis_disp_names


# # plot them
# sordox_plot <- ggplot(sor_wf_for_doxlis, aes(x = Name, y = log_relative_signal, fill = log_relative_signal)) +
#     geom_bar(
#         stat = "identity",
#         show.legend = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         plot.title = element_text(size = rel(1.25))
#     ) +
#     scale_fill_gradient2(
#         low = "blue", mid = "gray", high = "red",
#         midpoint = mean(sordox_wf$log_relative_signal, na.rm = TRUE),
#         limits = c(
#             min(sordox_wf$log_relative_signal, na.rm = TRUE),
#             max(sordox_wf$log_relative_signal, na.rm = TRUE)
#         )
#     ) +
#     ylim(c(
#         min(sordox_wf$log_relative_signal, na.rm = TRUE) - 0.5,
#         max(sordox_wf$log_relative_signal, na.rm = TRUE) + 0.5
#     )) +
#     ylab(TeX("$log_2(\\frac{sor+dox}{vehicle})$")) +
#     xlab("") +
#     geom_hline(aes(yintercept = thresh_dox), linewidth = rel(0.5), linetype = 3) +
#     geom_hline(aes(yintercept = -thresh_dox), linewidth = rel(0.5), linetype = 3) +
#     scale_x_discrete(breaks = sordox_disp_names) +
#     geom_text(aes(label = paste0(round(2^(sordox_wf$log_relative_signal), 2)), y = log_relative_signal + 0.05 * log_relative_signal), size = rel(1.5), color = "black") +
#     ggtitle(str_wrap("Potential Biomarkers of Doxazosin co-Treatment Response", 40)) +
#     labs(caption = qq("Threshold fold change <= +/-log2(@{thresh_dox})\nUpper fold change = @{round(2^thresh_dox, 3)}\nLower fold change = @{round(abs(2^(-thresh_dox)), 3)}"))
# sordox_plot

# sorlis_plot <- ggplot(sor_wf_for_doxlis, aes(x = Name, y = log_relative_signal, fill = log_relative_signal)) +
#     geom_bar(
#         stat = "identity",
#         show.legend = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         plot.title = element_text(size = rel(1.25))
#     ) +
#     scale_fill_gradient2(
#         low = "blue", mid = "gray", high = "red",
#         midpoint = mean(sorlis_wf$log_relative_signal, na.rm = TRUE),
#         limits = c(
#             min(sorlis_wf$log_relative_signal, na.rm = TRUE),
#             max(sorlis_wf$log_relative_signal, na.rm = TRUE)
#         )
#     ) +
#     ylim(c(
#         min(sorlis_wf$log_relative_signal, na.rm = TRUE) - 0.5,
#         max(sorlis_wf$log_relative_signal, na.rm = TRUE) + 0.5
#     )) +
#     ylab(TeX("$log_2(\\frac{sor+lis}{vehicle})$")) +
#     xlab("") +
#     geom_hline(aes(yintercept = thresh_lis), linewidth = rel(0.5), linetype = 3) +
#     geom_hline(aes(yintercept = -thresh_lis), linewidth = rel(0.5), linetype = 3) +
#     scale_x_discrete(breaks = sorlis_disp_names) +
#     ggtitle(str_wrap("Potential Biomarkers of Lisinopril co-Treatment Response", 40)) +
#     geom_text(aes(label = paste0(round(2^(sorlis_wf$log_relative_signal), 2)), y = log_relative_signal + 0.05 * log_relative_signal), size = rel(1.5), color = "black") +
#     labs(caption = qq("Threshold fold change <= +/-log2(@{thresh_lis})\nUpper fold change = @{round(2^thresh_lis, 3)}\nLower fold change = @{round(abs(2^(-thresh_lis)), 3)}"))
# sorlis_plot




# sor_dox_lis_interesting_names <- tibble(Name = c(sordox_disp_names, sorlis_disp_names)) %>%
#     mutate(color = c(
#         rep("dodgerblue", length(sordox_disp_names)),
#         rep("#ef5e5e", length(sorlis_disp_names))
#     )) %>%
#     group_by(Name) %>%
#     mutate(
#         count = n(),
#         color = ifelse(count > 1, "purple", color)
#     ) %>%
#     ungroup() %>%
#     select(-count)

# all_soraf_color_names <- sorafenib_wf %>%
#     dplyr::select(Name) %>%
#     left_join(sor_dox_lis_interesting_names, ,
#         by = join_by(Name)
#     ) %>%
#     replace_na(list(color = "black")) %>%
#     distinct() %>%
#     mutate(new_name = ifelse(Name %in% interesting_soraf_targets, str_c("** ", Name), as.character(Name)))


# sor_plot <- sor_plot +
#     theme(
#         axis.text.x = element_text(
#             angle = 90, vjust = 0.2, hjust = 0.95,
#             size = rel(1.05), color = all_soraf_color_names$color
#         )
#     ) +
#     scale_x_discrete(labels = all_soraf_color_names$new_name)
# sor_plot

# dox_colors <- all_soraf_color_names %>%
#     inner_join(sordox_wf) %>%
#     filter(Name %in% sordox_disp_names) %>%
#     arrange(desc(log_relative_signal))
# dox_colors

# sordox_plot <- sordox_plot +
#     theme(
#         axis.text.x = element_text(
#             angle = 90, vjust = 0.2, hjust = 0.95,
#             size = rel(1.1), color = dox_colors$color
#         ),
#         plot.title = element_text(size = rel(1.25), color = "dodgerblue")
#     )
# sordox_plot

# lis_colors <- all_soraf_color_names %>%
#     inner_join(sorlis_wf) %>%
#     filter(Name %in% sorlis_disp_names) %>%
#     arrange(desc(log_relative_signal))
# lis_colors

# sorlis_plot <- sorlis_plot +
#     theme(
#         axis.text.x = element_text(
#             angle = 90, vjust = 0.2, hjust = 0.95,
#             size = rel(1.1), color = lis_colors$color
#         ),
#         plot.title = element_text(size = rel(1.25), color = "#ef5e5e")
#     )

# # NEITHER INTERESTING
# tibble(Name = as.character(interesting_soraf_targets)) %>% arrange(Name)
# sor_dox_lis_interesting_names %>%
#     mutate(Name = as.character(Name)) %>%
#     arrange(Name)

# neither_dox_interesting <- anti_join(tibble(Name = interesting_soraf_targets), tibble(Name = sordox_disp_names))
# neither_dox_colors <- all_soraf_color_names %>% filter(Name %in% neither_dox_interesting$Name)

# dat_dox_neither_init <- dat %>%
#     filter(group == "Sor + Dox") %>%
#     filter(Name %in% neither_dox_colors$Name)

# dox_neither_wf <- dat_dox_neither_init %>%
#     arrange(desc(log_relative_signal)) %>%
#     ungroup() %>%
#     mutate(order = seq_len(n())) %>%
#     mutate(Name = factor(Name, levels = Name[order])) %>%
#     arrange(Name)

# # plot it
# sordox_neither_plot <- ggplot(dox_neither_wf, aes(x = Name, y = log_relative_signal, fill = log_relative_signal)) +
#     geom_bar(
#         stat = "identity",
#         show.legend = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(
#             angle = 90, vjust = 0.2, hjust = 0.95,
#             size = rel(1.1), color = "black"
#         ),
#         axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         plot.title = element_text(size = rel(1.25))
#     ) +
#     scale_fill_gradient2(
#         low = "blue", mid = "gray", high = "red",
#         midpoint = mean(dox_neither_wf$log_relative_signal, na.rm = TRUE),
#         limits = c(
#             min(dox_neither_wf$log_relative_signal, na.rm = TRUE),
#             max(dox_neither_wf$log_relative_signal, na.rm = TRUE)
#         )
#     ) +
#     ylab(TeX("$log_2(\\frac{sor+dox}{vehicle})$")) +
#     xlab("") +
#     geom_hline(aes(yintercept = thresh_dox), linewidth = rel(0.5), linetype = 3) +
#     geom_hline(aes(yintercept = -thresh_dox), linewidth = rel(0.5), linetype = 3) +
#     ggtitle(str_wrap("Potential Non-responder biomarkers of Doxazosin co-Treatment Response", 40)) +
#     geom_text(aes(label = paste0(round(2^(dox_neither_wf$log_relative_signal), 2)), y = log_relative_signal + 0.05 * log_relative_signal), size = rel(1.5), color = "black") +
#     labs(caption = qq("Threshold fold change <= +/-log2(@{thresh_dox})\nUpper fold change = @{round(2^thresh_dox, 3)}\nLower fold change = @{round(abs(2^(-thresh_dox)), 3)}"))
# sordox_neither_plot

# neither_lis_interesting <- anti_join(tibble(Name = interesting_soraf_targets), tibble(Name = sorlis_disp_names))
# neither_lis_colors <- all_soraf_color_names %>% filter(Name %in% neither_lis_interesting$Name)

# dat_lis_neither_init <- dat %>%
#     filter(group == "Sor + Lis") %>%
#     filter(Name %in% neither_lis_colors$Name)

# lis_neither_wf <- dat_lis_neither_init %>%
#     arrange(desc(log_relative_signal)) %>%
#     ungroup() %>%
#     mutate(order = seq_len(n())) %>%
#     mutate(Name = factor(Name, levels = Name[order])) %>%
#     arrange(Name)

# # plot it
# sorlis_neither_plot <- ggplot(
#     lis_neither_wf,
#     aes(x = Name, y = log_relative_signal, fill = log_relative_signal)
# ) +
#     geom_bar(
#         stat = "identity",
#         show.legend = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(
#             angle = 90, vjust = 0.2, hjust = 0.95,
#             size = rel(1.1), color = "black"
#         ),
#         axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.1)),
#         plot.title = element_text(size = rel(1.25))
#     ) +
#     scale_fill_gradient2(
#         low = "blue", mid = "gray", high = "red",
#         midpoint = mean(lis_neither_wf$log_relative_signal, na.rm = TRUE),
#         limits = c(
#             min(lis_neither_wf$log_relative_signal, na.rm = TRUE),
#             max(lis_neither_wf$log_relative_signal, na.rm = TRUE)
#         )
#     ) +
#     ylab(TeX("$log_2(\\frac{sor+lis}{vehicle})$")) +
#     xlab("") +
#     geom_hline(aes(yintercept = thresh_lis), linewidth = rel(0.5), linetype = 3) +
#     geom_hline(aes(yintercept = -thresh_lis), linewidth = rel(0.5), linetype = 3) +
#     ggtitle(str_wrap("Potential Non-responder biomarkers of Lisinopril co-Treatment Response", 40)) +
#     geom_text(aes(label = paste0(round(2^(lis_neither_wf$log_relative_signal), 2)), y = log_relative_signal + 0.05 * log_relative_signal), size = rel(1.5), color = "black") +
#     labs(caption = qq("Threshold fold change <= +/-log2(@{thresh_lis})\nUpper fold change = @{round(2^thresh_lis, 3)}\nLower fold change = @{round(abs(2^(-thresh_lis)), 3)}"))
# sorlis_neither_plot


# sor_dox_line <- dat_init %>%
#     filter(
#         Name %in% sordox_disp_names,
#         group %in% c("Vehicle", "Sorafenib", "Sor + Dox")
#     ) %>%
#     mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Dox")))
# sor_dox_line
# sordox_bar_plot <- ggplot(sor_dox_line) +
#     # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
#     geom_bar(
#         stat = "identity", mapping = aes(x = group, y = relative_signal, fill = group, group = Name),
#         color = "black",
#         width = 0.7, show.legend = FALSE
#     ) +
#     scale_fill_manual(values = c("#838080", "#a8a8a7", "#ded9d9")) +
#     geom_text(
#         aes(
#             label = paste0(round(sor_dox_line$relative_signal, 2)),
#             x = group,
#             y = mean(sor_dox_line$relative_signal) / 2
#         ),
#         size = rel(2),
#         color = "black"
#     ) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_bw() +
#     facet_wrap(~Name, ncol = 4, scales = "free_y") +
#     ylab(TeX("$Fold Change(\\frac{group}{vehicle})$")) +
#     scale_x_discrete(labels = c("V", "S", "D")) +
#     xlab("") +
#     ggtitle(str_wrap("Biomarker Fold Changes of Doxazosin Co-Treatment Response", 40))
# sordox_bar_plot

# sor_lis_line <- dat_init %>%
#     filter(
#         Name %in% sorlis_disp_names,
#         group %in% c("Vehicle", "Sorafenib", "Sor + Lis")
#     ) %>%
#     mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Lis")))
# sor_lis_line
# sorlis_bar_plot <- ggplot(sor_lis_line) +
#     # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
#     geom_bar(
#         stat = "identity", mapping = aes(x = group, y = relative_signal, group = Name, fill = group),
#         color = "black",
#         width = 0.7, show.legend = FALSE
#     ) +
#     scale_fill_manual(values = c("#838080", "#a8a8a7", "#ded9d9")) +
#     geom_text(
#         aes(
#             label = paste0(round(sor_lis_line$relative_signal, 2)),
#             x = group,
#             y = mean(sor_lis_line$relative_signal) / 2
#         ),
#         size = rel(2),
#         color = "black"
#     ) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_bw() +
#     facet_wrap(~Name, ncol = 4, scales = "free_y") +
#     ylab(TeX("$Fold Change(\\frac{group}{vehicle})$")) +
#     scale_x_discrete(labels = c("V", "S", "L")) +
#     xlab("") +
#     ggtitle(str_wrap("Biomarker Fold Changes of Lisinopril Co-Treatment Response", 40))
# sorlis_bar_plot


# sor_dox_neither <- dat_init %>%
#     filter(
#         Name %in% setdiff(interesting_soraf_targets, sordox_disp_names),
#         group %in% c("Vehicle", "Sorafenib", "Sor + Dox")
#     ) %>%
#     mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Dox")))
# sor_dox_neither
# sordox_neither_bar_plot <- ggplot(sor_dox_neither) +
#     # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
#     geom_bar(
#         stat = "identity", mapping = aes(x = group, y = relative_signal, fill = group, group = Name),
#         color = "black",
#         width = 0.7, show.legend = FALSE
#     ) +
#     scale_fill_manual(values = c("#838080", "#a8a8a7", "#ded9d9")) +
#     geom_text(
#         aes(
#             label = paste0(round(sor_dox_neither$relative_signal, 2)),
#             x = group,
#             y = mean(sor_dox_neither$relative_signal) / 2
#         ),
#         size = rel(2),
#         color = "black"
#     ) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_bw() +
#     facet_wrap(~Name, ncol = 4, scales = "free_y") +
#     ylab(TeX("$Fold Change(\\frac{group}{vehicle})$")) +
#     scale_x_discrete(labels = c("V", "S", "D")) +
#     xlab("") +
#     ggtitle(str_wrap("Non-responder Biomarker Fold Changes of Doxazosin Co-Treatment Response", 40))
# sordox_neither_bar_plot


# sor_lis_neither <- dat_init %>%
#     filter(
#         Name %in% setdiff(interesting_soraf_targets, sorlis_disp_names),
#         group %in% c("Vehicle", "Sorafenib", "Sor + Lis")
#     ) %>%
#     mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Lis")))
# sor_lis_neither
# sorlis_neither_bar_plot <- ggplot(sor_lis_neither) +
#     # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
#     geom_bar(
#         stat = "identity", mapping = aes(x = group, y = relative_signal, fill = group, group = Name),
#         color = "black",
#         width = 0.7, show.legend = FALSE
#     ) +
#     scale_fill_manual(values = c("#838080", "#a8a8a7", "#ded9d9")) +
#     geom_text(
#         aes(
#             label = paste0(round(sor_lis_neither$relative_signal, 2)),
#             x = group,
#             y = mean(sor_lis_neither$relative_signal) / 2
#         ),
#         size = rel(2),
#         color = "black"
#     ) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_bw() +
#     facet_wrap(~Name, ncol = 4, scales = "free_y") +
#     scale_x_discrete(labels = c("V", "S", "L")) +
#     xlab("") +
#     ylab(TeX("$Fold Change(\\frac{group}{vehicle})$")) +
#     ggtitle(str_wrap("Non-responder Biomarker Fold Changes of Doxazosin Co-Treatment Response", 40))
# sorlis_neither_bar_plot

# final_plot <- sor_plot / (sordox_plot | sordox_bar_plot) / (sordox_neither_bar_plot) / (sorlis_plot | sorlis_bar_plot) / sorlis_neither_bar_plot
# final_plot

# # library(autoggsaveR)
# # list()
# # autoggsaveR::auto_save_plot(plot_lst = list(), file_name = "waterfalls.pdf", relative_output_dir = "output/plots/Nick_final")

# ggsave(plot = final_plot, filename = "output/plots/Nick_final/waterfalls.pdf", width = 10, height = 16)
# ggsave(plot = sor_plot, filename = "output/plots/Nick_final/sor_waterfall.pdf", width = 8, height = 6)

# ### BARS
# interest_dat <- dat_init %>%
#     filter(
#         Name %in% interesting_soraf_targets,
#         group %in% c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis")
#     ) %>%
#     mutate(group = factor(group, levels = c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis")))

# interesting_plot <- ggplot(interest_dat) +
#     # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
#     geom_bar(
#         stat = "identity", mapping = aes(x = group, y = relative_signal, fill = group, group = Name),
#         color = "black",
#         width = 0.7, show.legend = FALSE
#     ) +
#     scale_fill_grey() +
#     # scale_fill_manual(values = c("#838080", "#a8a8a7", "#ded9d9")) +
#     geom_text(
#         aes(
#             label = paste0(round(interest_dat$relative_signal, 2)),
#             x = group,
#             y = mean(interest_dat$relative_signal) / 2
#         ),
#         size = rel(2),
#         color = "black"
#     ) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_bw() +
#     facet_wrap(~Name, ncol = 2, scales = "free_y") +
#     scale_x_discrete(labels = c("V", "S", "D", "L")) +
#     xlab("") +
#     ylab(TeX("$Fold Change(\\frac{group}{vehicle})$")) +
#     ggtitle(str_wrap("Potentially Interesting Biomarkers", 40))
# interesting_plot
# ggsave(plot = interesting_plot, filename = "output/plots/Nick_final/interesting_analytes.pdf", width = 4, height = 7)

# # Waterfall plus interesting plot
# wf_plus_bar_interesting_plot <- sor_plot / interesting_plot
# ggsave(plot = wf_plus_bar_interesting_plot, filename = "output/plots/Nick_final/wf_and_interesting_analytes.pdf", width = 7, height = 10)
