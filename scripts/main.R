rm(list = ls())
library(tidyverse)
library(readxl)
library(GetoptLong)
library(RColorBrewer)
library(latex2exp)
library(patchwork)
library(ggprism)
library(ggh4x)
library(progressr)
library(future)
library(furrr)

# progress
handlers(handler_progress(
    format   = ":spin [:bar] :current/:total :percent in :elapsed ETA: :eta (:message) ",
    width    = 120,
    complete = "+"
))

# helper scripts
source(file.path("scripts", "project_paths.R"))
source(file.path("scripts", "array_helper_scripts.R"))

#' @Number0 This code requires that the analytes are measured in LICOR software starting in order of the protocol datasheet.
#' So the S001 corresponds to the very first item on the protocol coordinate system. Background should be measured last, after negative controls

#' @Number1 Parse the protocol data
#' Note: You must run scripts/extract_analyte_table.py first!!

#' @Number2 Make an .env file with the appropriate arguments.
#' Note: Regarding the .env file:
#' commas to separate the elements of the vectors
#' pipe (|) and colon (:) to represent the list of comparisons,
#' e.g "Group1:Group2" group1 is the control and group2 is the comparatore.g.

#' @Number3 Load the appropriate .env file
# dotenv::load("arrays/Lauren Pregnancy Cytokine XL Array/pregnancy_cytokine_xl_array.env")
# dotenv::load("arrays/Veh vs Sor Dox Lis - Cytokine XL/dox_cytokine_xl_array.env")
# dotenv::load_dot_env("arrays/Veh vs Sor Dox Lis - Angiogenesis Protein Array/angiogenesis_array.env")

#

example_config <- get_analysis_config("vegfri_dox_cytokine_xl")

info_fn <- resolve_project_path(example_config$info_fn, must_exist = TRUE)
data_dir <- resolve_project_path(example_config$data_dir, must_exist = TRUE)
output_dir <- resolve_project_path(example_config$output_dir, must_exist = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

my_main_threshold <- example_config$main_threshold
my_ref_thresh_to_filter <- example_config$ref_thresh_to_filter
my_comparisons <- example_config$comparisons
my_groups_per_page <- example_config$groups_per_page
my_group_lvls <- factor(example_config$group_levels, levels = example_config$group_levels)

if (!is_perfect_square(n = my_groups_per_page)) {
    message(qq("Warning: @{my_groups_per_page} is not a perfect square. Facet plot results may be ugly.\n"))
}

control_point_color <- "#ffffff"
tx1_point_color <- "#737171"
tx2_point_color <- "#CCD8E8"
tx2_outline_color <- "#22456F"
tx3_point_color <- "#E5938A"
tx3_outline_color <- "#B53530"
my_colors <- c(control_point_color, tx1_point_color, tx2_point_color, tx3_point_color)[seq_len(nlevels(my_group_lvls))] %>%
    set_names(my_group_lvls)
my_outline_colors <- c("black", "black", tx2_outline_color, tx3_outline_color)[seq_len(nlevels(my_group_lvls))] %>%
    set_names(my_group_lvls)

message("Using analysis config: vegfri_dox_cytokine_xl")
message(qq("Resolved data dir: @{data_dir}"))
message(qq("Resolved protocol table: @{info_fn}"))
message(qq("Resolved output dir: @{output_dir}"))

analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

my_initial_ready_df <- make_plot_ready_dataset(
    data_dir = data_dir,
    analyte_info = analyte_info_df,
    preview = FALSE,
    my_group_lvls = my_group_lvls
)

num_cores <- max(1, min(parallel::detectCores() - 1, length(my_main_threshold)))
message(qq("Switching to multisession mode. Using @{num_cores} cores..."))
plan(multisession, workers = num_cores)

total_iterations <- length(my_ref_thresh_to_filter) * length(my_main_threshold)

with_progress({
    message(qq("\nStarting processing of @{data_dir}"))
    message(qq("Output directory: @{output_dir}"))

    p_inner <- progressor(steps = total_iterations)
    for (i in seq_along(my_ref_thresh_to_filter)) {
        furrr::future_walk(.x = seq_along(my_main_threshold), .f = function(j, p_inner) {
            p_inner(qq("Processing ref thresh @{my_ref_thresh_to_filter[i]} and main thresh @{my_main_threshold[j]}"))

            make_graphs(
                df = my_initial_ready_df,
                ref_thresh_to_filter = my_ref_thresh_to_filter[i],
                main_threshold = my_main_threshold[j],
                output_dir_full_path = output_dir,
                comparisons = my_comparisons,
                groups_per_page = my_groups_per_page
            )
        }, p_inner = p_inner, .options = furrr::furrr_options(seed = TRUE))
    }
})
plan(sequential)
