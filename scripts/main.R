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
library(dotenv)

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

default_env_file <- resolve_example_env_file(
    project_name = "Veh vs Sor Dox Lis - Cytokine XL",
    env_filename = "dox_cytokine_xl_array.env"
)

dot_string_vec <- c(Sys.getenv("PROTEOME_PROFILER_ENV_FILE", unset = default_env_file))

for (dot_string in dot_string_vec) {
    dot_string <- resolve_project_path(dot_string, must_exist = TRUE)
    dotenv::load_dot_env(dot_string)

    # Load variables from .env
    info_fn <- resolve_project_path(Sys.getenv("INFO_FN"), must_exist = TRUE)
    data_dir <- resolve_project_path(Sys.getenv("DATA_DIR"), must_exist = TRUE)
    output_dir <- resolve_project_path(Sys.getenv("OUTPUT_DIR"), must_exist = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Convert comma-separated strings to vectors
    my_main_threshold <- as.numeric(unlist(strsplit(Sys.getenv("MAIN_THRESHOLD"), ",")))
    my_ref_thresh_to_filter <- as.numeric(unlist(strsplit(Sys.getenv("REF_THRESH_TO_FILTER"), ",")))

    # Convert comma-separated strings to vectors for group levels
    my_group_lvls <- unlist(strsplit(Sys.getenv("GROUP_LVLS"), ","))

    # Parse comparisons as a list
    comparisons_str <- str_split(Sys.getenv("COMPARISONS"), "\\|", simplify = TRUE)
    my_comparisons <- lapply(comparisons_str, function(x) {
        key_val <- str_split(x, "\\:", simplify = TRUE)
        key <- key_val[1]
        val <- str_split(key_val[2], ",", simplify = TRUE)[1, ]
        return(list(key = key, val = val))
    })
    names(my_comparisons) <- sapply(my_comparisons, function(x) x$key)
    my_comparisons <- lapply(my_comparisons, function(x) x$val)

    # Groups per page
    my_groups_per_page <- as.integer(Sys.getenv("GROUPS_PER_PAGE"))


    ################################################################################################################################# s
    #################################################################################################################################
    #################################################################################################################################
    if (!is_perfect_square(n = my_groups_per_page)) {
        message(qq("Warning: @{my_groups_per_page} is not a perfect square. Facet plot results may be ugly.\n"))
    }

    control_point_color <- "#ffffff"
    tx1_point_color <- "#737171"
    tx2_point_color <- "#CCD8E8"
    tx2_outline_color <- "#22456F"
    tx3_point_color <- "#E5938A"
    tx3_outline_color <- "#B53530"
    my_group_lvls <- factor(my_group_lvls, levels = my_group_lvls)
    my_colors <- c(control_point_color, tx1_point_color, tx2_point_color, tx3_point_color)[seq_len(nlevels(my_group_lvls))] %>%
        set_names(my_group_lvls)
    my_outline_colors <- c("black", "black", tx2_outline_color, tx3_outline_color)[seq_len(nlevels(my_group_lvls))] %>%
        set_names(my_group_lvls)

    message(qq("Using env file: @{dot_string}"))
    message(qq("Resolved data dir: @{data_dir}"))
    message(qq("Resolved protocol table: @{info_fn}"))
    message(qq("Resolved output dir: @{output_dir}"))


    # read in analyte info
    analyte_info_df <- read_excel(info_fn) %>%
        mutate(sname_grouping = row_number()) %>%
        rename(Name = `Analyte/Control`)

    ## Run it
    my_initial_ready_df <- make_plot_ready_dataset(
        data_dir = data_dir,
        analyte_info = analyte_info_df,
        preview = FALSE,
        my_group_lvls = my_group_lvls
    )

    # RUN
    # DEBUG:
    # comparisons <- my_comparisons
    # i <- 4
    # j <- 4
    # ref_thresh_to_filter <- my_ref_thresh_to_filter[i]
    # main_threshold <- my_main_threshold[j]
    # output_dir_full_path <- output_dir
    # df <- my_initial_ready_df
    # data <- df
    # groups_per_page = my_groups_per_page
    # comparisons = my_comparisons

    num_cores <- max(1, min(parallel::detectCores() - 1, length(my_main_threshold)))
    message(qq("Switching to multisession mode. Using @{num_cores} cores..."))
    plan(multisession, workers = num_cores)

    total_iterations <- length(my_ref_thresh_to_filter) * length(my_main_threshold)

    with_progress({
        # Create a progressor function
        message(qq("\nStarting processing of @{data_dir}"))
        message(qq("Output directory: @{output_dir}"))

        p_inner <- progressor(steps = total_iterations)
        for (i in seq_along(my_ref_thresh_to_filter)) {
            # Update outer loop progress
            furrr::future_walk(.x = seq_along(my_main_threshold), .f = function(j, p_inner) {
                # Update inner loop progress
                p_inner(qq("Processing ref thresh @{my_ref_thresh_to_filter[i]} and main thresh @{my_main_threshold[j]}"))

                sink_res <- make_graphs(
                    df = my_initial_ready_df,
                    ref_thresh_to_filter = my_ref_thresh_to_filter[i],
                    main_threshold = my_main_threshold[j],
                    output_dir_full_path = output_dir,
                    comparisons = my_comparisons,
                    groups_per_page = my_groups_per_page
                )
                # tick
                # p(sprintf("Processing index %s", analyte_name))
            }, p_inner = p_inner, .options = furrr::furrr_options(seed = TRUE))

            # }
        }
    })
    plan(sequential)
}
