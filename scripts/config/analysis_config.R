# This is the user-editable project configuration file.
# Change the paths and analysis entries here when running the project on a new
# machine or for a new dataset.
#
# The VEGFRi/Dox entry below is an example configuration, not a required
# directory structure.

proteome_profiler_config <- list(
    runtime_root = "~/ProjectsRuntime/proteome_profiler",
    cloud_parent = "~/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects",
    examples = list(
        vegfri_dox_cytokine_xl = list(
            # Analyst name used to namespace outputs under output/plots/<user>/.
            user = "nick",
            # Stable folder name for this analysis under the user subtree.
            analysis_slug = "cytokine_xl_array",
            # Preset that tells the setup step which array layout/protocol schema
            # to expect when extracting the analyte workbook.
            protocol_preset = "cytokine_xl",
            # Generated workbook of analyte names / coordinates extracted from
            # the vendor protocol. Created by scripts/setup/extract_analyte_table.py.
            info_fn = "output/cytoXL array kit - protocol.xlsx",
            # Source protocol PDF used to generate / verify the workbook above.
            protocol_pdf = "protocols/cytoXL array kit - protocol.pdf",
            # Pages in the protocol PDF that contain the analyte spot map.
            protocol_pages = c(17, 18, 19),
            # Directory containing the raw LI-COR spot-export workbooks for the
            # experimental groups in this analysis.
            data_dir = "projects/Veh vs Sor Dox Lis - Cytokine XL/data",
            # Expected group names and plotting order used throughout the run.
            group_levels = c("vehicle", "sorafenib", "sor + dox", "sor + lis"),
            # Pairwise comparisons to generate. Names are control/reference
            # groups; values are treatment groups to compare against them.
            comparisons = list("vehicle" = c("sorafenib", "sor + dox", "sor + lis")),
            # Minimum reference-signal cutoffs used to drop low-intensity
            # analytes before fold-change plots are generated.
            ref_thresh_to_filter = c(150, 200),
            # Fold-change cutoff(s) used to define "hits" in the main analysis.
            main_threshold = c(1.46),
            # Number of analytes per faceted barplot page.
            groups_per_page = 25,
            # Manually selected low-signal analyte coordinates used as a
            # heuristic reference panel when choosing a raw-signal cutoff.
            # These are not the protocol's true reference or negative-control
            # spots; they are dataset-specific analytes chosen to sit near the
            # boundary between usable and unusable signal.
            ref_coords_to_make_filter = c("H1,2", "H3,4", "G1,2", "G3,4", "F9,10", "F11,12", "E9,10", "E11,12"),
            # Default control arm used by the analyte-selection helper script.
            selection_control = "vehicle",
            # Default treatment arm used by the analyte-selection helper script.
            selection_group = "sorafenib",
            # Fold-change cutoff used by the analyte-selection helper script.
            selection_threshold = 1.49
        )
    )
)
