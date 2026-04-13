rm(list = ls())
library(tidyverse)
library(readxl)
library(GetoptLong)
library(RColorBrewer)
library(latex2exp)
library(patchwork)
library(ggprism)
library(ggh4x)

source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number1 protocol data
#' @note run scripts/setup/extract_analyte_table.py first if the protocol workbook does not exist yet.

example_config <- get_analysis_config("vegfri_dox_cytokine_xl")
info_fn <- example_config$info_fn
data_dir <- example_config$data_dir
output_dir <- example_config$output_dir
ref_coords_to_make_filter <- example_config$ref_coords_to_make_filter
my_group_lvls <- example_config$group_levels
info_fn <- resolve_project_path(info_fn, must_exist = TRUE)
data_dir <- resolve_project_path(data_dir, must_exist = TRUE)
output_dir <- resolve_project_path(output_dir, must_exist = FALSE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)



# read in analyte info
analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

## Run it
df <- make_plot_ready_dataset(data_dir, analyte_info = analyte_info_df, preview = TRUE, my_group_lvls = my_group_lvls)

long_df <- df %>%
    pivot_longer(all_of(my_group_lvls), names_to = "group", values_to = "signal")

# CHeck this first if you don't want to look through the blots themselves and select regions
# ggplot(long_df, mapping = aes(x = signal)) +
#     geom_histogram(binwidth = 1000) +
#     facet_wrap(~group)


control_point_color <- "#ffffff"
tx1_point_color <- "#737171"
tx2_point_color <- "#9BB3D3"
tx3_point_color <- "#B53530"
my_group_lvls <- factor(my_group_lvls)
my_colors <- c(control_point_color, tx1_point_color, tx2_point_color, tx3_point_color)[seq_len(nlevels(my_group_lvls))] %>%
    set_names(my_group_lvls)


## Identify the threshold to filter out analytes
filt_lst <- find_filter_thresh(wide_df = df, ref_coords = ref_coords_to_make_filter, my_colors = my_colors)
# filt_lst <- find_filter_thresh(wide_df = df, filter_threshold = 100, my_colors = my_colors)

ggsave(
    filename = file.path(output_dir, "region_stats.png"), plot = filt_lst[[2]],
    width = 10, height = 10
)
