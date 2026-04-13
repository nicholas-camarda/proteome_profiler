rm(list = ls())
source(file.path("scripts", "helpers", "runtime_setup.R"))
load_analysis_packages(include_parallel = FALSE)

source(file.path("scripts", "config", "analysis_config.R"))
source(file.path("scripts", "helpers", "project_paths.R"))
source(file.path("scripts", "helpers", "array_helper_scripts.R"))

#' @Number1
#' This script writes the sorafenib shortlist views under the per-user analysis
#' tree in `select_analytes/`.
example_config <- get_analysis_config("vegfri_dox_cytokine_xl")
info_fn <- get_protocol_workbook_path(example_config)
data_dir <- example_config$data_dir
analysis_output_root <- get_analysis_output_root(example_config)
output_dir <- file.path(analysis_output_root, "select_analytes")
data_dir <- resolve_project_path(data_dir, must_exist = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

my_group_lvls <- example_config$group_levels

# # CytoXL labels
# already_performed_ELISAs <- tibble(Name = c(
#     "MIP-1a", "MCP-1", "PD-ECGF",
#     "HGF", "Serpin E1/PAI-1", "Prolactin",
#     "PDGF-AA", "Angiopoietin-2", "CCL21/6Ckine", "Complement Component C5/C5a", "CXCL13/BLC/BCA-1",
#     "Pentraxin 2/SAP", "E-Selectin/CD62E", "P-Selectin/CD62P", "MMP-9", "IGFBP-1"
# ), color = "#1fc61f") # "Pentraxin 3/TSG-14",
already_performed_ELISAs <- tibble(Name = c(""), color = "")


# these had dox/lis ~> 0.9
names_to_check <- c(
    "CCL21/6Ckine", "Angiopoietin-2", "Angiopoietin-1", "IGFBP-2", "Angiopoietin-like 3",
    "Complement Component C5/C5a", "CXCL13/BLC/BCA-1",
    "E-Selectin/CD62E", "P-Selectin/CD62P",
    "Serpin E1/PAI-1", # already done
    "Pentraxin 2/SAP", "MMP-9",
    "IGFBP-1" # honorable mentions
)
# names_to_check <- c("Angiopoietin-2", "Angiopoietin-like 3", "VCAM-1/CD106", "P-Selectin/CD62P", "Pentraxin 2/SAP", "BAFF/BLyS/TNFSF13B", "C-Reactive Protein/CRP", "CCL11/Eotaxin")

style_config <- build_group_style(my_group_lvls, scheme = "main")
my_group_lvls <- style_config$group_levels
my_colors <- style_config$fill
my_outline_colors <- style_config$outline



# read in analyte info
analyte_info_df <- read_excel(info_fn) %>%
    mutate(sname_grouping = row_number()) %>%
    rename(Name = `Analyte/Control`)

df <- make_plot_ready_dataset(
    data_dir = data_dir,
    analyte_info = analyte_info_df,
    preview = TRUE,
    my_group_lvls = my_group_lvls
)

# CHANGE THIS FOR REF THRESH FILTERING
my_ref_thresh <- example_config$ref_thresh_to_filter[[1]]
sor_thresh <- example_config$selection_threshold

# Reuse the pairwise comparison helper, but narrow it to the specific
# vehicle-vs-sorafenib contrast used for the shortlist view.
special_lst_obj <- make_wf_data(df,
    my_main_threshold = sor_thresh, ref_thresh_to_filter_ = my_ref_thresh,
    comparisons = setNames(list(c(example_config$selection_group)), example_config$selection_control)
)
result_key <- qq("@{example_config$selection_control} vs @{example_config$selection_group}")
special_wf <- special_lst_obj[[result_key]]$wf_dat[[example_config$selection_group]]
special_dat <- special_lst_obj[[result_key]]$dat_filtered[[example_config$selection_group]]
# fc_lst <- log2_fold_change_to_raw_fc(sor_log2_thresh)
# upper_thresh <- fc_lst[[1]]
# lower_thresh <- fc_lst[[2]]
special_wf_final <- special_wf %>%
    filter(relative_signal > sor_thresh | relative_signal < 1 / sor_thresh) %>%
    mutate(order = factor(seq_along(order))) %>%
    mutate(Name = if_else(as.character(Name) == "Complement Component C5/C5a", "C5a", Name))


gg_wf_plot <- plot_wf_graph(
    data = special_wf_final,
    add_fc = FALSE,
    # done_df = already_performed_ELISAs,
    TITLE = qq("Potential biomarkers of sorafenib treatment\nRef thresh > @{my_ref_thresh} || Sor thresh < @{round(1 / sor_thresh, 3)} or Sor thresh > @{round(sor_thresh, 3)}"),
    main_cutoff = NULL,
    line_color = "#f20f75"
)
gg_wf_plot


output_subdir <- output_dir
dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(
    output_subdir,
    qq("barplot_select_REF-@{gsub(x = my_ref_thresh, pattern = '\\\\.', replacement = '_')}_SOR-@{gsub(x = round(sor_thresh,3),  pattern = '\\\\.', replacement = '_')}.png")
)

ggsave(
    plot = gg_wf_plot,
    filename = file.path(output_subdir, qq("waterfall_plot_REF-@{gsub(x = my_ref_thresh, pattern = '\\\\.', replacement = '_')}_SOR-@{gsub(x = round(sor_thresh, 3),  pattern = '\\\\.', replacement = '_')}.png")),
    width = 14, height = 11, dpi = 300
)

select_faceted_bargraphs_lst <- make_bar_charts(
    data = special_dat %>%
        filter(Name %in% unique(special_wf_final$Name)),
    # Name %in% names_to_check |
    title = qq("Potential biomarkers of sorafenib treatment\nRef thresh > @{my_ref_thresh} || Sor thresh < @{round(1 / sor_thresh, 3)} or Sor thresh > @{round(sor_thresh, 3)}"),
    add_fc = TRUE,
    caption = FALSE,
    make_facet_plot = FALSE # returns a list of single bar plots if FALSE
)
filtered_names <- str_remove_all(names(select_faceted_bargraphs_lst), pattern = " \\([A-Z]*[0-9]+,[0-9]+\\)$")
# The plot titles still include coordinate suffixes from the protocol table.
# Strip those off so exported filenames are readable for collaborators.
dir.create(file.path(output_subdir, "select_bargraphs"), showWarnings = FALSE, recursive = TRUE)
message("\n")
for (i in seq_along(filtered_names)) {
    message(filtered_names[i])
    ggsave(
        plot = select_faceted_bargraphs_lst[[i]],
        filename = file.path(output_subdir, "select_bargraphs", str_c(filtered_names[i], ".png")),
        width = 10, height = 10
    )
}



message(qq("Check @{output_path}"))
