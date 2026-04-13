rm(list = ls())
library(tidyverse)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
library(GetoptLong)
library(RColorBrewer)
library(latex2exp)
library(patchwork)
library(ggprism)

library(rstatix)
library(ggpubr)
# for hacking the facet text color title
library(ggh4x)
library(autoggsaveR)
library(openxlsx)


df <- read_excel("arrays/Lauren Inflammation Protein Array/2023-04-27 XL mouse cytokine array analysis.xlsx",
    sheet = 1
) %>%
    # remove blank rows
    na.omit() %>%
    rename(
        Sname = `Licor analysis label`,
        Name = Analyte
    )

dat_init1 <- df %>%
    # group_by(Coordinate, Sname) %>%
    mutate(Name = make.unique(Name)) %>%
    pivot_longer(`Control-5359`:`sFlt-A43`,
        names_to = "group-ID", values_to = "signal"
    ) %>%
    filter(signal >= 0) %>%
    separate(col = `group-ID`, into = c("group", "ID")) %>%
    group_by(Coordinate, Sname, Name) %>%
    mutate(
        mean_ref_signal = mean(signal[group == "Control"], na.rm = TRUE),
        relative_signal = signal / mean_ref_signal, # this is the fold change
        log_relative_signal = log2(relative_signal), # this is the log fold change
        # fc = signal / mean_ref_signal, # this is the fold change
        # logfc = log2(fc) # this is the log fold change
    ) %>%
    group_by(Name, group) %>%
    mutate(
        num_samples = max(seq_len(n()))
    ) %>%
    ungroup()

samples_to_remove <- dat_init1 %>%
    filter(
        num_samples < 4,
        Name != "NEGATIVE CONTROL"
    ) %>%
    distinct(Name, num_samples)

message("Removing these samples because of low data:")
print(samples_to_remove)

my_dat_all_temp <- dat_init1 %>%
    mutate(num_patients = max(length(unique(ID[group == "sFlt"])))) %>%
    # remove any analyte that doesn't have the max number of samples
    filter(num_samples == num_patients) %>%
    dplyr::select(-num_patients)

dat <- my_dat_all_temp %>%
    filter(group != "Control")

wf_temp <- dat %>%
    filter(!(Name %in% c("REF SPOT", "REF SPOT.1", "REF SPOT.2", "NEGATIVE CONTROL"))) %>%
    ungroup() %>%
    filter(!is.na(log_relative_signal))
wf_temp


my_wf_mean <- wf_temp %>%
    group_by(Name, Coordinate, Sname) %>%
    # calculate the sem
    summarize(
        mean_relative_signal = mean(relative_signal, na.rm = TRUE),
        sem_relative_signal = sd(relative_signal, na.rm = TRUE) / sqrt(length(Name)),
        upper_relative_signal = mean_relative_signal + sem_relative_signal,
        lower_relative_signal = mean_relative_signal - sem_relative_signal,
        mean_relative_log_signal = log2(mean_relative_signal),

        # mean_relative_log_signal = mean(log_relative_signal, na.rm = TRUE),
        sem_relative_log_signal = sd(log_relative_signal, na.rm = TRUE) / sqrt(length(Name)),
        upper_relative_log_signal = mean_relative_log_signal + sem_relative_log_signal,
        lower_relative_log_signal = mean_relative_log_signal - sem_relative_log_signal,
        .groups = "drop"
    ) %>%
    arrange(desc(mean_relative_log_signal)) %>%
    mutate(
        order = seq_len(n()),
        Name = factor(Name, levels = Name[order])
    )
my_wf_mean

my_wf <- left_join(wf_temp, my_wf_mean) %>%
    ungroup() %>%
    arrange(order)

# dat_init1 %>%
#     pivot_wider(Name, names_from = c(group, ID), values_from = relative_signal) %>%
#     write_tsv("~/Downloads/test.tsv")

shapiro_test_res <- dat_init1 %>%
    group_by(Name, group) %>%
    filter(!(num_samples < 4)) %>%
    shapiro_test(signal) %>%
    # shapiro_test(relative_signal) %>%
    mutate(normal_temp = ifelse(p < 0.05, FALSE, TRUE)) %>%
    group_by(Name) %>%
    summarize(normal_sum = sum(normal_temp)) %>%
    mutate(normal = normal_sum == 2)

shapiro_test_res %>% filter(Name == "IGFBP-6")
# # A tibble: 2 × 6
#   Name    group   variable statistic      p normal_temp
#   <chr>   <chr>   <chr>        <dbl>  <dbl> <lgl>
# 1 IGFBP-6 Control signal       0.758 0.0453 FALSE
# 2 IGFBP-6 sFlt    signal       0.846 0.214  TRUE

shapiro_test_res

t_test_normals <- dat_init1 %>%
    left_join(shapiro_test_res, by = c("Name")) %>%
    filter(normal) %>%
    group_by(Name) %>%
    # https://github.com/kassambara/rstatix/issues/28#issuecomment-587196783
    # code performs a separate t-test for each value of 'Name'
    # the p-values should be adjusted within each test, not across all tests
    # so this is correct but won't do anything since we're only comparing two groups (1 comparison)
    # this would scale with the number of levels in the comparison
    # so a simple g1 vs g2 would have no correction while a
    # g1 vs g2 | g2 vs g3 | g1 vs g3 would have 3 comparisons and require a correction
    #
    #
    # In the case of the pairwise_t_test() function in the rstatix package, it is performing
    # multiple testing correction separately for each level of the grouping variable. This is
    # appropriate only
    # **if you consider each level of the grouping variable to be a separate family of hypotheses**,
    # and you want to control the false discovery rate within each family separately.
    # However, in your case (and in many gene expression analyses), you might consider all the
    # tests across all genes to be one family of hypotheses. In this case, you would want to apply
    # multiple testing correction across all the tests, not separately for each gene.
    pairwise_t_test(relative_signal ~ group,
        p.adjust.method = "fdr",
        alternative = "two.sided",
    ) %>%
    arrange(p.adj) %>%
    mutate(test_type = "t_test") %>%
    ungroup() %>%
    add_xy_position()
t_test_normals

# Q <- 0.25
# xxx <- t_test_normals$p * (length(t_test_normals$p) / seq_len(length(t_test_normals$p)))*Q

t_test_non_normals <- dat_init1 %>%
    left_join(shapiro_test_res, by = c("Name")) %>%
    filter(!normal) %>%
    group_by(Name) %>%
    # code performs a separate t-test for each value of 'Name'
    # the p-values should be adjusted within each test, not across all tests
    # so this is correct but won't do anything since we're only comparing two groups
    # this would scale with the number of groups in the comparison
    # so a simple g1 vs g2 would have no correction while a
    # g1 vs g2 | g2 vs g3 | g1 vs g3 would have 3 comparison correction
    wilcox_test(relative_signal ~ group,
        p.adjust.method = "fdr",
        alternative = "two.sided",
    ) %>%
    add_significance() %>%
    dplyr::select(-statistic) %>%
    mutate(p.adj = p, p.adj.signif = p.signif) %>%
    arrange(p.adj) %>%
    mutate(test_type = "wilcox") %>%
    ungroup() %>%
    add_xy_position()
t_test_non_normals
t_test <- bind_rows(t_test_normals, t_test_non_normals) %>%
    arrange(p.adj)

my_dat_all <- left_join(my_dat_all_temp, t_test, by = "Name")

#' @note do a quick manual t-test on gene name
#' @param my_gene_name any gene name that exists in the dataset
#' @param var_equal whether or not to assume equal variance
#' @return t.test result
man_calc_t_test <- function(my_gene_name = "VCAM", var_equal = FALSE) {
    vcam1_sig_ctrl <- dat_init1 %>%
        filter(Name == my_gene_name, group == "Control") %>%
        pluck("signal")
    vcam1_sig_sflt <- dat_init1 %>%
        filter(Name == my_gene_name, group == "sFlt") %>%
        pluck("signal")
    t.test(vcam1_sig_ctrl, vcam1_sig_sflt, var.equal = var_equal)
}

make_bar_charts <- function(data, title = "", add_fc = TRUE, strip_text_rel_size = 1.1) {
    helper_dat <- tibble(
        group = c("Control", "sFlt"),
        short_group = factor(c("Ctrl", "sFlt"), levels = c("Ctrl", "sFlt"))
    )

    pos_labs_df <- data %>%
        group_by(Name, group) %>%
        summarize(
            pos = max(relative_signal) / 6,
            labs = round(mean(relative_signal), 3), .groups = "drop"
        )

    short_group_data <- left_join(data, helper_dat, by = "group") %>%
        left_join(pos_labs_df, by = c("group", "Name")) %>%
        arrange(Name) %>%
        group_by(group) %>%
        mutate(order = seq_len(n()), Name = factor(Name, levels = unique(Name[unique(order)]))) %>%
        mutate(Name_Coordinate = str_c(Name, " (", Coordinate, ")"))

    test_data <- short_group_data %>%
        ungroup() %>%
        distinct(
            Name, relative_signal, `.y.`, group, short_group, group1, group2, n1, n2, p, p.signif, p.adj, p.adj.signif,
            test_type, y.position, xmin, xmax, Name_Coordinate
        ) %>%
        # group_by(Name) %>%
        #     rename(old_y_position = y.position) %>%
        #     group_by(Name, group) %>%
        #     mutate(y.position_temp = mean(relative_signal, na.rm = TRUE)) %>%
        #     group_by(Name) %>%
        #     mutate(y.position = max(relative_signal, na.rm = TRUE) - 0.) %>%
        # ungroup() %>%
        dplyr::select(-relative_signal) %>%
        distinct()
    test_data
    # short_group_data <- short_group_data %>% filter(Sname %in% c("S004", "S006", "S010"))

    stopifnot(nrow(short_group_data) > 0)

    ctrl_color <- "#ffffff" # "#737171"
    sFlt1_color <- "#B53530"
    my_colors <- c(ctrl_color, sFlt1_color) %>%
        set_names(c("Control", "sFlt"))

    bar_plot <- ggplot(short_group_data,
        mapping = aes(x = short_group, y = relative_signal, fill = group, group = Name_Coordinate)
    ) +
        stat_summary(geom = "col", fun = "mean", color = "black") +
        geom_jitter(position = position_jitter(width = 0.2), color = "black") +
        facet_wrap(~Name_Coordinate,
            scales = "free_y",
        ) +
        theme_prism(base_size = 16) +
        scale_fill_manual(values = my_colors) +
        labs(x = "Group", y = "Fold Change to Ctrl") +
        ggtitle(title) +
        theme(
            strip.text = element_text(size = rel(strip_text_rel_size)),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = rel(0.5)),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = rel(0.25)),
            aspect.ratio = 1
        ) +
        stat_pvalue_manual(test_data)
    bar_plot

    if (add_fc) {
        bar_plot <- bar_plot +
            geom_text(
                mapping = aes(
                    x = short_group,
                    label = labs,
                    y = pos, group = Name
                ), size = rel(2), color = "black"
            )
    }
    bar_plot
}

# plot it
make_plot <- function(wf_mean1, wf1, dat_all, t_test, threshold, include_stars_name = FALSE, color_x_axis_names = FALSE) {
    disp_names <- subset(wf_mean1, abs(mean_relative_log_signal) > threshold)$Name
    bar_text_direction <- sapply(wf_mean1$mean_relative_log_signal, FUN = function(x) ifelse(x < 0, -1, 1))

    # my_wf %>%
    #     arrange(mean_relative_log_signal) %>%
    #     print(n = 50)
    # my_wf_mean %>%
    #     mutate(NAME = as.character(Name), .before = 1) %>%
    #     arrange(NAME) %>%
    #     dplyr::select(Name, NAME, mean_relative_log_signal) %>%
    #     filter(Name == "Fetuin A")

    # wf_mean %>% filter(Name %in% disp_names)
    upper_fold_change_thresh <- round(2^threshold, 3)
    upper_fold_change_thresh
    lower_fold_change_thresh <- round(2^(-threshold), 3)
    lower_fold_change_thresh

    wf_mean <- wf_mean1 %>%
        filter(Name %in% disp_names)
    wf <- wf1 %>%
        filter(Name %in% disp_names)

    og_order <- wf_mean %>%
        distinct(Name) %>%
        mutate(order = row_number())

    waterfall_plot <- ggplot(wf_mean, aes(
        x = Name, y = mean_relative_log_signal,
        fill = mean_relative_log_signal, group = Name
    )) +
        geom_bar(
            stat = "identity",
            show.legend = FALSE
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95, face = "bold", size = rel(1.2)),
            axis.text.y = element_text(angle = 90, vjust = 0.2, hjust = 0.95, size = rel(1.3)),
            plot.title = element_text(size = rel(1))
        ) +
        # color the bars by taking into account all analytes, not just those that pass filter
        scale_fill_gradient2(
            low = "blue", mid = "gray", high = "red",
            midpoint = mean(wf_mean1$mean_relative_log_signal, na.rm = TRUE),
            limits = c(
                min(wf_mean1$mean_relative_log_signal, na.rm = TRUE),
                max(wf_mean1$mean_relative_log_signal, na.rm = TRUE)
            )
        ) +
        ylim(c(min(wf_mean$mean_relative_log_signal, na.rm = TRUE) - 0.01, max(wf_mean$mean_relative_log_signal, na.rm = TRUE) + 0.01)) +
        xlab("") +
        ylab(TeX("$log_2(\\frac{sFlt1}{\\it{mean}(Control)})$")) +
        geom_hline(aes(yintercept = threshold), linewidth = rel(0.5), linetype = 4) +
        geom_hline(aes(yintercept = -threshold), linewidth = rel(0.5), linetype = 4) +
        geom_hline(aes(yintercept = 0), linewidth = rel(0.25), linetype = 1, color = "black") +
        ggtitle("Potential Biomarkers of sFlt1 Treatment Response") +
        labs(caption = qq("Threshold fold change >= +/-log2(@{threshold})\nUpper fold change = @{upper_fold_change_thresh}\nLower fold change = @{lower_fold_change_thresh}\nRemoved @{str_c(samples_to_remove$Name, collapse = ', ')}"))
    waterfall_plot


    interesting_targets <- disp_names %>% droplevels()

    if (include_stars_name) {
        interesting_names <- wf_mean %>%
            left_join(t_test, by = join_by(Name)) %>%
            mutate(
                new_name = ifelse(Name %in% interesting_targets, str_c("** ", Name), as.character(Name)),
                color = ifelse(Name %in% interesting_targets, "purple", "black"),
                pass_manual_thresh = ifelse(Name %in% interesting_targets, TRUE, FALSE),
                .before = 1
            ) %>%
            filter(pass_manual_thresh)
        interesting_qvalue_targets <- wf_mean %>%
            left_join(t_test, by = join_by(Name)) %>%
            mutate(
                new_name = ifelse(p.adj.signif != "ns", str_c("** ", Name), Name),
                color = ifelse(p.adj.signif != "ns", "#7272e9", "black"),
                pass_pval_thresh = ifelse(p.adj.signif != "ns", TRUE, FALSE),
                .before = 1
            ) %>%
            filter(pass_pval_thresh)
    } else {
        interesting_names <- wf_mean %>%
            left_join(t_test, by = join_by(Name)) %>%
            mutate(
                new_name = as.character(Name),
                color = ifelse(Name %in% interesting_targets, "purple", "black"),
                pass_manual_thresh = ifelse(Name %in% interesting_targets, TRUE, FALSE),
                .before = 1
            ) %>%
            filter(pass_manual_thresh)

        interesting_qvalue_targets <- wf_mean %>%
            left_join(t_test, by = join_by(Name)) %>%
            mutate(
                new_name = Name,
                color = ifelse(p.adj.signif != "ns", "#7272e9", "black"),
                pass_pval_thresh = ifelse(p.adj.signif != "ns", TRUE, FALSE),
                .before = 1
            ) %>%
            filter(pass_pval_thresh)
    }


    interesting_in_both <- tibble(new_name = intersect(interesting_qvalue_targets$new_name, interesting_names$new_name)) %>%
        mutate(
            color = "#00a200",
            pass_both = TRUE,
            pass_pval_thresh = TRUE,
            pass_manual_thresh = TRUE
        ) %>%
        mutate(Name = str_trim(gsub(x = new_name, pattern = "\\** ", replacement = ""))) %>%
        left_join(wf_mean, by = join_by(Name)) %>%
        left_join(t_test, by = join_by(Name))

    not_interesting <- anti_join(wf_mean, bind_rows(interesting_names, interesting_qvalue_targets), by = join_by(
        Name, Coordinate, Sname, mean_relative_signal, sem_relative_signal, upper_relative_signal, lower_relative_signal, mean_relative_log_signal,
        sem_relative_log_signal, upper_relative_log_signal, lower_relative_log_signal, order
    )) %>%
        mutate(new_name = Name) %>%
        left_join(t_test, by = join_by(Name))

    stopifnot(nrow(not_interesting) + nrow(interesting_qvalue_targets) + nrow(interesting_names) - nrow(interesting_in_both) == nrow(wf_mean))

    interesting_qvalue_targets <- interesting_qvalue_targets %>%
        filter(!(new_name %in% interesting_in_both$new_name))
    interesting_names <- interesting_names %>%
        filter(!(new_name %in% interesting_in_both$new_name))

    final_data_with_names <- bind_rows(interesting_qvalue_targets, interesting_names, not_interesting, interesting_in_both) %>%
        dplyr::select(new_name, color, pass_pval_thresh, pass_manual_thresh, pass_both, everything()) %>%
        replace_na(list(
            pass_pval_thresh = FALSE,
            pass_manual_thresh = FALSE,
            pass_both = FALSE,
            color = "black"
        )) %>%
        dplyr::select(-order) %>%
        left_join(og_order, by = c("new_name" = "Name")) %>%
        arrange(order)
    nrow(final_data_with_names)

    if (color_x_axis_names) {
        waterfall_plot_final <- waterfall_plot +
            theme(
                axis.text.x = element_text(
                    angle = 90, vjust = 0.2, hjust = 0.95,
                    color = final_data_with_names$color
                ), plot.title = element_text(size = rel(1.5)),
                plot.margin = margin(0.5, 0.75, 0.5, 0.5, "cm")
            ) +
            labs(subtitle = "Purple = Analytes passing manual threshold\nOrange = Analytes passing t-test significance threshold\nGreen = Analytes passing both t-test significance and manual threshold")
    } else {
        waterfall_plot_final <- waterfall_plot +
            theme(
                axis.text.x = element_text(
                    angle = 90, vjust = 0.2, hjust = 0.95
                ), plot.title = element_text(size = rel(1.5)),
                plot.margin = margin(0.5, 0.75, 0.5, 0.5, "cm")
            ) +
            labs(subtitle = "Purple = Analytes passing manual threshold\nOrange = Analytes passing t-test significance threshold\nGreen = Analytes passing both t-test significance and manual threshold")
    }
    waterfall_plot_final

    # one hypothesis per gene
    t_test_normals_new_one_hypoth <- dat_all %>%
        filter(Name %in% disp_names) %>%
        left_join(shapiro_test_res, by = c("Name")) %>%
        filter(normal) %>%
        group_by(Name) %>%
        pairwise_t_test(relative_signal ~ group,
            p.adjust.method = "fdr",
            alternative = "two.sided",
        ) %>%
        arrange(p.adj) %>%
        mutate(test_type = "t_test") %>%
        ungroup() %>%
        add_xy_position()
    t_test_normals_new_one_hypoth

    # Q <- 0.25
    # xxx <- t_test_normals$p * (length(t_test_normals$p) / seq_len(length(t_test_normals$p)))*Q

    t_test_non_normals_new_one_hypoth <- dat_all %>%
        filter(Name %in% disp_names) %>%
        left_join(shapiro_test_res, by = c("Name")) %>%
        filter(!normal) %>%
        group_by(Name) %>%
        wilcox_test(relative_signal ~ group,
            p.adjust.method = "fdr",
            alternative = "two.sided",
        ) %>%
        add_significance() %>%
        dplyr::select(-statistic) %>%
        mutate(p.adj = p, p.adj.signif = p.signif) %>%
        arrange(p.adj) %>%
        mutate(test_type = "wilcox") %>%
        ungroup() %>%
        add_xy_position()
    t_test_non_normals_new_one_hypoth
    t_test_new_one_hypothesis <- bind_rows(t_test_normals_new_one_hypoth, t_test_non_normals_new_one_hypoth) %>%
        arrange(p.adj) %>%
        dplyr::select(-c(y.position, groups, xmin, xmax))
    t_test_new_one_hypothesis


    # each gene is a part of one large hypothesis
    t_test_normals_new_all <- dat_all %>%
        filter(Name %in% disp_names) %>%
        left_join(shapiro_test_res, by = c("Name")) %>%
        filter(normal) %>%
        group_by(Name) %>%
        pairwise_t_test(relative_signal ~ group,
            p.adjust.method = "none",
            alternative = "two.sided",
        ) %>%
        ungroup() %>%
        mutate(p.adj = p.adjust(p, method = "fdr")) %>%
        add_significance(p.col = "p.adj", output.col = "p.adj.signif", cutpoints = c(0.1, 0.25, 1), symbols = c("*", "ns")) %>%
        arrange(p.adj) %>%
        mutate(test_type = "t_test") %>%
        ungroup() %>%
        add_xy_position()
    t_test_normals_new_all

    # Q <- 0.25
    # xxx <- t_test_normals$p * (length(t_test_normals$p) / seq_len(length(t_test_normals$p)))*Q

    t_test_non_normals_new_all <- dat_all %>%
        filter(Name %in% disp_names) %>%
        left_join(shapiro_test_res, by = c("Name")) %>%
        filter(!normal) %>%
        group_by(Name) %>%
        wilcox_test(relative_signal ~ group,
            p.adjust.method = "none",
            alternative = "two.sided",
        ) %>%
        dplyr::select(-statistic) %>%
        ungroup() %>%
        add_significance() %>%
        mutate(p.adj = p.adjust(p, method = "fdr")) %>%
        add_significance(p.col = "p.adj", output.col = "p.adj.signif", cutpoints = c(0.1, 0.25, 1), symbols = c("*", "ns")) %>%
        arrange(p.adj) %>%
        mutate(test_type = "wilcox") %>%
        ungroup() %>%
        add_xy_position()
    t_test_non_normals_new_all
    t_test_new_all <- bind_rows(t_test_normals_new_all, t_test_non_normals_new_all) %>%
        arrange(p.adj) %>%
        dplyr::select(-c(y.position, groups, xmin, xmax))
    t_test_new_all
    # bar_charts <- make_bar_charts(data = dat_all %>% filter(Name %in% disp_names), title = "", add_fc = TRUE, strip_text_rel_size = 1.1)

    return(list(waterfall_plot_final, final_data_with_names, t_test_new_one_hypothesis, t_test_new_all) %>%
        set_names(c("waterfall_plot", "final_data", "stats_one_hypothesis_per_gene", "stats_one_hypothesis_for_all_genes")))
}

# manual
threshold <- 0.25
# plots
my_include_stars <- FALSE
my_plot_and_info_lst <- make_plot(
    wf_mean1 = my_wf_mean, wf1 = my_wf, dat_all = my_dat_all, t_test, threshold,
    include_stars_name = my_include_stars, color_x_axis_names = FALSE
)
waterfall_plot1 <- my_plot_and_info_lst[[1]]
waterfall_plot1
final_data_with_names <- my_plot_and_info_lst[[2]]
stats_one_hypothesis_per_gene <- my_plot_and_info_lst[[3]]
stats_one_hypothesis_per_gene
stats_one_hypothesis_for_all_genes <- my_plot_and_info_lst[[4]]
stats_one_hypothesis_for_all_genes

stars_dir <- ifelse(my_include_stars, "with_stars", "no_stars")
dir.create(qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}"), showWarnings = FALSE, recursive = TRUE)


{
    write.xlsx(x = stats_one_hypothesis_per_gene, file = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/stats_one_hypothesis_per_gene.xlsx"))
    write.xlsx(x = stats_one_hypothesis_for_all_genes, file = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/stats_one_hypothesis_for_all_genes.xlsx"))

    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_8x8.png"),
        width = 8, height = 8
    )
    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_9x9.png"),
        width = 9, height = 9
    )
    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_11x8.png"),
        width = 11, height = 8
    )

    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_15x7.png"),
        width = 15, height = 7
    )
    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_14x7.png"),
        width = 14, height = 7
    )
    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_13x7.png"),
        width = 13, height = 7
    )
    ggsave(
        plot = waterfall_plot1,
        filename = qq("output/plots/Lauren_final/@{stars_dir}/thresh_@{threshold}/Lauren-waterfall-threshold_@{threshold}_13x8.png"),
        width = 13, height = 8
    )
}

# auto_save_plot(
#     plot_lst = list(waterfall_plot1),
#     relative_output_dir = "output/plots/Lauren_final",
#     file_name = qq("Lauren-waterfall-threshold_@{threshold}.pdf"),
#     ncol = 1
# )
# ggsave(filename = qq("output/plots/Lauren-waterfalls-unfiltered-no_ebars-no_fc_1-1_threshold_@{threshold}.pdf"), unfiltered_plot, width = 13, height = 15)
# ggsave(filename = qq("output/plots/Lauren-waterfall-pval_thresh-unfiltered-no_ebars-no_fc_1-1_threshold_@{threshold}.pdf"), waterfall_plot2, width = 15, height = 8)


# interesting_names_chr <- str_trim(gsub(
#     pattern = "\\**", replacement = "",
#     x = final_data_with_names$new_name[startsWith(as.character(final_data_with_names$new_name), "**")]
# ))

# bar_wf <- dat_all %>%
#     ungroup() %>%
#     filter(Name %in% interesting_names_chr)

# mean_bar_wf <- bar_wf %>%
#     group_by(Name, group) %>%
#     summarize(
#         mean_relative_signal = mean(relative_signal, na.rm = TRUE),
#         my_sd = sd(relative_signal, na.rm = TRUE),
#         len_name = length(Name),
#         sem_relative_signal = my_sd / sqrt(len_name),
#         upper_relative_signal = mean_relative_signal + sem_relative_signal,
#         lower_relative_signal = mean_relative_signal - sem_relative_signal,
#         max_val = max(relative_signal, na.rm = TRUE), .groups = "keep"
#     ) %>%
#     left_join(final_data_with_names %>% distinct(Name, color), by = join_by(Name))
# group_names <- c("Control", "sFlt")
# bar_colors <- c("black", "#FB070A")
# names(bar_colors) <- group_names
# shapes <- c(16, 15)
# names(shapes) <- group_names

# ttest_dat <- t_test %>%
#     filter(Name %in% (bar_wf %>% pluck("Name") %>% unique())) %>%
#     arrange(p.adj) %>%
#     group_by(Name) %>%
#     add_xy_position() %>%
#     ungroup()

# facet_title_colors <- mean_bar_wf %>%
#     ungroup() %>%
#     distinct(Name, color) %>%
#     .$color

# my_strip_colors <- strip_themed(text_x = elem_list_text(color = facet_title_colors))

# interesting_marks_plot <- ggplot(mean_bar_wf) +
#     geom_bar(
#         stat = "identity",
#         mapping = aes(x = group, y = mean_relative_signal, color = group, group = Name),
#         fill = "white",
#         linewidth = rel(1.25),
#         width = 0.7, show.legend = FALSE
#     ) +
#     geom_errorbar(
#         mapping = aes(
#             x = group,
#             ymin = lower_relative_signal,
#             ymax = upper_relative_signal,
#             color = group
#         ),
#         width = 0.25
#     ) +
#     geom_jitter(
#         data = bar_wf,
#         mapping = aes(
#             x = group, y = relative_signal, shape = group,
#             color = group, fill = group, group = Name
#         ), width = 0.1,
#         size = rel(2.5)
#     ) +
#     ggh4x::facet_wrap2(~Name, scales = "free_y", strip = my_strip_colors, ncol = 8) +
#     scale_color_manual(values = bar_colors) +
#     scale_fill_manual(values = bar_colors) +
#     scale_shape_manual(values = shapes) +
#     geom_text(
#         data = mean_bar_wf,
#         mapping = aes(
#             label = paste0(round(mean_relative_signal, 2)),
#             x = group,
#             y = mean(mean_relative_signal) / 2
#         ),
#         size = rel(4),
#         color = "black"
#     ) +
#     # geom_hline(aes(yintercept = 0), color = "black", linetype = 1, linewidth = rel(1)) +
#     theme_prism(base_size = 16) +
#     scale_x_discrete(labels = c("Control", "sFlt")) +
#     xlab("") +
#     ylab(TeX("$FoldChange(\\frac{sFlt1}{\\it{mean}(Control)})$")) +
#     ggtitle(str_wrap("Potentially Interesting Biomarkers of sFlt Treatment", 450)) +
#     stat_pvalue_manual(
#         inherit.aes = FALSE,
#         ttest_dat,
#         # color = "anti_htn_group",
#         show.legend = FALSE,
#         label.size = rel(5), # rel(5)
#         bracket.size = 1,
#         # hide.ns = TRUE,
#         label = "p = {p.adj}" # "p = {round(p.adj, 6)}"
#     ) +
#     labs(subtitle = "Purple = Analytes passing manual threshold\nOrange = Analytes passing t-test significance threshold\nGreen = Analytes passing both t-test significance and manual threshold") +
#     theme(
#         plot.title = element_text(size = rel(2), vjust = 1, hjust = 0),
#         plot.subtitle = element_text(size = rel(1), vjust = 5, hjust = 0),
#         panel.spacing = unit(1.5, "lines"),
#         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
#     )
# interesting_marks_plot


# ggsave(filename = qq("output/plots/Lauren_final/Lauren-interesting_analytes-threshold_@{threshold}.pdf"), plot = interesting_marks_plot, dpi = 300, width = 15, height = 20)
# ggsave(
#     filename = qq("output/plots/Lauren_final/Lauren-interesting_analytes.pdf"),
#     plot = interesting_marks_plot,
#     dpi = 300, width = 20, height = 15
# )

# # auto_save_plot(
# #     plot_lst = list(interesting_marks_plot),
# #     relative_output_dir = "output/plots/Lauren_final",
# #     file_name = qq("Lauren-interesting_analytes.pdf"),
# #     ncol = 1
# # )


# wf_plus_marks_plot <- waterfall_plot1 / interesting_marks_plot
# ggsave(filename = qq("output/plots/Lauren_final/Lauren-analytes_and_wf-threshold_@{threshold}.pdf"), plot = wf_plus_marks_plot, dpi = 300, width = 20, height = 40)
