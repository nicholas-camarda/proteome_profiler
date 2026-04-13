# This is the one file users are expected to edit.
#
# What to change for a new project:
# 1. Set `runtime_root`
# 2. Set `cloud_parent`
# 3. Add one entry under `analyses`
# 4. Run with: export PROTEOME_PROFILER_ANALYSIS=your_analysis_name
#
# The VEGFRi/Dox block below is a worked reference only.

proteome_profiler_config <- list(
    default_analysis = NULL,
    runtime_root = "~/ProjectsRuntime/proteome_profiler",
    cloud_parent = "~/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects",
    analyses = list(
        vegfri_dox_cytokine_xl = list(
            mode = "legacy",
            user = "nick",
            slug = "cytokine_xl_array",
            protocol = list(
                preset = "cytokine_xl",
                workbook = "output/cytoXL array kit - protocol.xlsx",
                pdf = "protocols/cytoXL array kit - protocol.pdf",
                pages = c(17, 18, 19)
            ),
            input = list(
                manifest = "manifests/vegfri_dox_cytokine_xl_samples.csv",
                treatment = "treatment"
            ),
            comparisons = list(
                vehicle = c("sorafenib", "sor + dox", "sor + lis")
            ),
            thresholds = list(
                ref_coords = c("H1,2", "H3,4", "G1,2", "G3,4", "F9,10", "F11,12", "E9,10", "E11,12"),
                ref_signal = c(150, 200),
                fold_change = c(1.46),
                groups_per_page = 25
            ),
            shortlist = list(
                control = "vehicle",
                treatment = "sorafenib",
                fold_change = 1.49
            )
        ),
        replicate_aware_template = list(
            mode = "replicate",
            user = "nicole",
            slug = "mf_veh_vs_aldo",
            protocol = list(
                preset = "cytokine_xl",
                workbook = "output/cytoXL array kit - protocol.xlsx",
                pdf = "protocols/cytoXL array kit - protocol.pdf",
                pages = c(17, 18, 19)
            ),
            input = list(
                manifest = "manifests/example_samples.csv", ## change
                subgroup = "sex",
                treatment = "treatment"
            ),
            comparisons = list(
                control = c("treated")
            ),
            thresholds = list(
                ref_coords = c("A3,4"), ## change
                ref_signal = c(150) ## change
            ),
            stats = list(
                min_reps_per_arm = 2,
                p_adjust_method = "BH",
                alpha = 0.05
            ),
            shortlist = list(
                comparisons = c("male_control_vs_treated", "female_control_vs_treated"),
                top_n = 10
            )
        )
    )
)
