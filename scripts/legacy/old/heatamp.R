library(tidyverse)
library(ComplexHeatmap)
library(readxl)
library(RColorBrewer)

df <- read_excel("/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/projects/VEGFRi and Dox/in-vivo projects/vehicle vs sor dox lis/protein arrays/Veh vs Sor Dox Lis - Angiogenesis Protein Array/angiogenesis_analysis-v2.xlsx",
    sheet = 3
) %>%
    # mutate(across(everything(), .fns = function(x) ifelse(x < 0, 0, x))) %>%
    na.omit()
# mutate(across(Vehicle:`Sor + Lis`, .fns = function(x) log(x + 1))) # unnormalized

lis_path <- "output/lis"
dox_path <- "output/dox"
both_path <- "output/both"
dir.create(lis_path, showWarnings = FALSE, recursive = TRUE)
dir.create(dox_path, showWarnings = FALSE, recursive = TRUE)
dir.create(both_path, showWarnings = FALSE, recursive = TRUE)

df_mat_temp <- df %>%
    dplyr::select(Vehicle:`Sor + Lis`)


df_mat <- as.matrix(df_mat_temp)
rownames(df_mat) <- df$Name

# dev.off()
# create the heatmap
ht_map <- Heatmap(df_mat,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = "Treatment", # add a column title
    row_title = "Analyte",
    name = "Raw Signal", # add a title
    width = 6, # set the plot width
    height = 6
)
ht_map

# create a new PDF device
pdf("heatmap.pdf", width = 8, height = 10)

# draw the heatmap to the PDF device
draw(ht_map)

# close the PDF device
dev.off()

## RELATIVE HEATMAP TO VEHICLE
rel_df_long <- df %>%
    pivot_longer(
        cols = Vehicle:`Sor + Lis`,
        names_to = "group",
        values_to = "value"
    ) %>%
    group_by(group) %>%
    mutate(Name = make.unique(Name)) %>%
    group_by(Name) %>%
    mutate(rel_val = value / value[1L]) %>%
    ungroup() %>%
    mutate(
        rel_val = map_dbl(rel_val, function(x) log(x + 1))
    )

rel_df <- rel_df_long %>%
    dplyr::select(-value) %>%
    pivot_wider(names_from = group, values_from = rel_val) %>%
    ungroup()
rel_df

rel_df_mat <- as.matrix(rel_df %>% dplyr::select(Vehicle:`Sor + Lis`))
rel_df_mat_scaled <- t(scale(t(rel_df_mat), center = FALSE, scale = apply(rel_df_mat, 1, median)))
rownames(rel_df_mat_scaled) <- rel_df$Name


# dev.off()
# create the heatmap
rel_ht_map <- Heatmap(rel_df_mat_scaled,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "Treatment", # add a column title
    row_title = "Analyte",
    name = "Relative Signal to Vehicle", # add a title
    width = 6, # set the plot width
    height = 6
)
rel_ht_map

# create a new PDF device
pdf("relative-heatmap.pdf", width = 8, height = 10)

# draw the heatmap to the PDF device
draw(rel_ht_map)

# close the PDF device
dev.off()

raw_df <- rel_df_long %>%
    dplyr::select(-rel_val) %>%
    pivot_wider(names_from = group, values_from = value) %>%
    ungroup()

raw_df %>%
    mutate(check = case_when(
        Vehicle <= Sorafenib & Vehicle >= `Sor + Dox` & Vehicle <= `Sor + Lis` ~ TRUE,
        Vehicle > Sorafenib & Vehicle < `Sor + Dox` & Vehicle > `Sor + Lis` ~ FALSE
    ))

bar_plot_rel_val_temp <- rel_df_long %>%
    mutate(group = factor(group,
        levels = c("Vehicle", "Sorafenib", "Sor + Dox", "Sor + Lis")
    ))

sum_bar_plot_rel_val <- bar_plot_rel_val_temp %>%
    group_by(Name) %>%
    nest(data = c(group, rel_val, value)) %>%
    reframe(lims = map(data, function(x) c(min(x$value, na.rm = TRUE) - 0.5, max(x$value, na.rm = TRUE) + 0.5)))

bar_plot_rel_val <- bar_plot_rel_val_temp %>%
    left_join(sum_bar_plot_rel_val, by = join_by(Name))


# Set the y-axis limits for each facet
y_limits <- function(values) {
    c(min(values[values != 0]), max(values))
}

barplots <- ggplot(bar_plot_rel_val) +
    # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
    geom_line(mapping = aes(x = group, y = value, group = Name)) +
    geom_point(mapping = aes(x = group, y = value, color = group)) +
    theme_bw() +
    facet_wrap(~Name, ncol = 4, scales = "free_y") +
    # scale_y_continuous(limits = y_limits) +
    # coord_cartesian(ylim = c(min(bar_plot_rel_val$value), NA)) +
    # scale_y_continuous(limits = sum_bar_plot_rel_val$lims) +
    theme(axis.text.x = element_blank())
ggsave(filename = "lineplots.pdf", plot = barplots, width = 10, height = 22)



# my chosen ones
interesting_analytes_dox <- c(
    "Angiopoietin-3", "Aniogenin", "DLL4", "Endostatin",
    "FGF basic", "Fractalkine", "HGF", "Prolactin",
    "Thrombospondin-2", "Timp-1"
)
chosen_dox <- bar_plot_rel_val %>%
    filter(Name %in% interesting_analytes_dox)
wide_chosen_dox <- chosen_dox %>%
    pivot_wider(id_cols = Name, names_from = group, values_from = value)
write_tsv(wide_chosen_dox, file.path(dox_path, "doxazosin_analytes.tsv"))

lineplot_dox <- ggplot(chosen_dox) +
    # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
    geom_line(mapping = aes(x = group, y = value, group = Name)) +
    geom_point(mapping = aes(x = group, y = value, color = group)) +
    theme_bw() +
    facet_wrap(~Name, ncol = 4, scales = "free_y") +
    # scale_y_continuous(limits = y_limits) +
    # coord_cartesian(ylim = c(min(bar_plot_rel_val$value), NA)) +
    # scale_y_continuous(limits = sum_bar_plot_rel_val$lims) +
    theme(axis.text.x = element_blank()) +
    ggtitle("Interesting analytes with respect to doxazosin tx")
ggsave(filename = file.path(dox_path, "lineplots-dox.pdf"), plot = lineplot_dox, width = 8, height = 5.5)

interesting_analytes_lis <- c(
    "IL-10", "IP-10", "KC", "DPPIV", "IGFBP-3",
    "MIP-1a", "Timp-4", "Cyr61", "PlGF-2", "PDGF-AB/PDGF-BB", "Pentraxin-3"
)
chosen_lis <- bar_plot_rel_val %>%
    filter(Name %in% interesting_analytes_lis)
wide_chosen_lis <- chosen_lis %>%
    pivot_wider(id_cols = Name, names_from = group, values_from = value)
write_tsv(wide_chosen_lis, file.path(lis_path, "lisinopril_analytes.tsv"))

lineplot_lis <- ggplot(chosen_lis) +
    # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
    geom_line(mapping = aes(x = group, y = value, group = Name)) +
    geom_point(mapping = aes(x = group, y = value, color = group)) +
    theme_bw() +
    facet_wrap(~Name, ncol = 4, scales = "free_y") +
    # scale_y_continuous(limits = y_limits) +
    # coord_cartesian(ylim = c(min(bar_plot_rel_val$value), NA)) +
    # scale_y_continuous(limits = sum_bar_plot_rel_val$lims) +
    theme(axis.text.x = element_blank()) +
    ggtitle("Interesting analytes with respect to lisinopril tx")
ggsave(filename = file.path(lis_path, "lineplots-lis.pdf"), plot = lineplot_lis, width = 8, height = 4)


## FOR NOW DON"T WORRY
interesting_analytes_both <- c(
    "Cyr61", "MCP-1", "PD-ECGF", "PDGF-AA"
)
chosen_both <- bar_plot_rel_val %>%
    filter(Name %in% interesting_analytes_both)
wide_chosen_both <- chosen_both %>%
    pivot_wider(id_cols = Name, names_from = group, values_from = value)
write_tsv(wide_chosen_both, file.path(both_path, "both_doxlis_analytes.tsv"))

lineplot_both <- ggplot(chosen_both) +
    # geom_col(mapping = aes(x = group, y = value, fill = group), width = 0.5) +
    geom_line(mapping = aes(x = group, y = value, group = Name)) +
    geom_point(mapping = aes(x = group, y = value, color = group)) +
    theme_bw() +
    facet_wrap(~Name, ncol = 4, scales = "free_y") +
    # scale_y_continuous(limits = y_limits) +
    # coord_cartesian(ylim = c(min(bar_plot_rel_val$value), NA)) +
    # scale_y_continuous(limits = sum_bar_plot_rel_val$lims) +
    theme(axis.text.x = element_blank()) +
    ggtitle("Interesting analytes with respect to both lis and dox tx")
ggsave(filename = file.path(both_path, "lineplots-both.pdf"), plot = lineplot_both, width = 10, height = 4)
