rm(list = ls())
source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = FALSE)

source(file.path("scripts", "config", "analysis_config.R"))
source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number1 protocol data
#' @note run scripts/setup/extract_analyte_table.py before this script if the protocol workbook has not been created yet.

# Load the shared analysis metadata, then write this script's outputs into the
# per-user analysis tree under `threshold_diagnostics/`.
# The configured coordinates are a manually chosen low-signal analyte panel
# used to estimate a practical raw-signal floor for this dataset.
example_config <- get_analysis_config("vegfri_dox_cytokine_xl")
info_fn <- get_protocol_workbook_path(example_config)
data_dir <- example_config$data_dir
analysis_output_root <- get_analysis_output_root(example_config)
output_dir <- file.path(analysis_output_root, "threshold_diagnostics")
ref_coords_to_make_filter <- example_config$ref_coords_to_make_filter
my_group_lvls <- example_config$group_levels
data_dir <- resolve_project_path(data_dir, must_exist = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)



# read in analyte info
analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

## Build the analyte-level dataset once, then inspect the user-selected
## background-control coordinates to choose a raw-signal floor.
df <- make_plot_ready_dataset(data_dir, analyte_info = analyte_info_df, preview = TRUE, my_group_lvls = my_group_lvls)

# Suggest a reviewable starting panel of low-signal analytes so the user can
# choose `ref_coords_to_make_filter` with actual data in hand rather than by
# trial and error alone.
candidate_low_signal_df <- suggest_low_signal_panel(df)
write_tsv(candidate_low_signal_df, file.path(output_dir, "candidate_low_signal_analytes.tsv"))
message("Suggested low-signal analytes written to:")
message(file.path(output_dir, "candidate_low_signal_analytes.tsv"))

long_df <- df %>%
    pivot_longer(all_of(my_group_lvls), names_to = "group", values_to = "signal")

# CHeck this first if you don't want to look through the blots themselves and select regions
# ggplot(long_df, mapping = aes(x = signal)) +
#     geom_histogram(binwidth = 1000) +
#     facet_wrap(~group)


style_config <- build_group_style(my_group_lvls, scheme = "threshold")
my_group_lvls <- style_config$group_levels
my_colors <- style_config$fill


## This writes the diagnostic figure used to pick a reference threshold before
## running the full `main.R` workflow.
filt_lst <- find_filter_thresh(wide_df = df, ref_coords = ref_coords_to_make_filter, my_colors = my_colors)
# filt_lst <- find_filter_thresh(wide_df = df, filter_threshold = 100, my_colors = my_colors)

ggsave(
    filename = file.path(output_dir, "region_stats.png"), plot = filt_lst[[2]],
    width = 10, height = 10
)
