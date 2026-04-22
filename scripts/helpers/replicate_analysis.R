# Public loader for replicate-aware helper modules.
# Entry scripts and tests should source this file rather than individual modules.

replicate_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
replicate_helper_candidates <- unique(c(
    if (!is.null(replicate_source_file) && !is.na(replicate_source_file) && nzchar(replicate_source_file)) {
        dirname(normalizePath(replicate_source_file, winslash = "/", mustWork = FALSE))
    },
    file.path("scripts", "helpers"),
    file.path("..", "scripts", "helpers"),
    file.path("..", "..", "scripts", "helpers")
))
replicate_helper_dir <- replicate_helper_candidates[
    file.exists(file.path(replicate_helper_candidates, "replicate_method_specs.R"))
][[1]]

replicate_target_env <- environment()

invisible(lapply(
    file.path(
        replicate_helper_dir,
        c(
            "replicate_method_specs.R",
            "replicate_inputs.R",
            "replicate_models.R",
            "replicate_plotting.R",
            "replicate_outputs.R",
            "replicate_selected_analytes.R"
        )
    ),
    source,
    local = replicate_target_env
))
