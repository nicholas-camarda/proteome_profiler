find_repo_root <- function(start = getwd()) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)

    repeat {
        if (file.exists(file.path(current, "scripts", "project_paths.R"))) {
            return(current)
        }

        parent <- dirname(current)
        if (identical(parent, current)) {
            break
        }
        current <- parent
    }

    stop("Could not determine the proteome_profiler repository root.")
}

repo_root <- find_repo_root()

# Edit this block when running the project on a different machine or with a
# different storage layout. The VEGFRi/Dox entry is an example configuration.
proteome_profiler_config <- list(
    runtime_root = "~/ProjectsRuntime/proteome_profiler",
    cloud_parent = "~/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects",
    examples = list(
        vegfri_dox_cytokine_xl = list(
            info_fn = "output/cytoXL array kit - protocol.xlsx",
            data_dir = "projects/Veh vs Sor Dox Lis - Cytokine XL/data",
            output_dir = "output/plots/nick/cytokine_xl_array",
            group_levels = c("vehicle", "sorafenib", "sor + dox", "sor + lis"),
            comparisons = list("vehicle" = c("sorafenib", "sor + dox", "sor + lis")),
            ref_thresh_to_filter = c(150, 200),
            main_threshold = c(1.46),
            groups_per_page = 25,
            ref_coords_to_make_filter = c("H1,2", "H3,4", "G1,2", "G3,4", "F9,10", "F11,12", "E9,10", "E11,12"),
            selection_control = "vehicle",
            selection_group = "sorafenib",
            selection_threshold = 1.49
        )
    )
)

path_candidates <- function(path) {
    runtime_root <- path.expand(proteome_profiler_config$runtime_root)
    cloud_parent <- path.expand(proteome_profiler_config$cloud_parent)
    cloud_project_root <- file.path(cloud_parent, "proteome_profiler")

    unique(c(
        file.path(repo_root, path),
        file.path(runtime_root, path),
        file.path(cloud_project_root, path)
    ))
}

resolve_project_path <- function(path, must_exist = TRUE) {
    candidates <- path_candidates(path)
    existing <- candidates[file.exists(candidates) | dir.exists(candidates)]

    if (length(existing) > 0) {
        return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
    }

    if (!must_exist) {
        return(normalizePath(file.path(path.expand(proteome_profiler_config$runtime_root), path), winslash = "/", mustWork = FALSE))
    }

    stop(sprintf(
        "Could not resolve '%s'. Checked:\n%s",
        path,
        paste(sprintf("- %s", candidates), collapse = "\n")
    ))
}

get_analysis_config <- function(example_name = "vegfri_dox_cytokine_xl") {
    config <- proteome_profiler_config$examples[[example_name]]
    if (is.null(config)) {
        stop(sprintf("Unknown example config: %s", example_name))
    }
    config
}
