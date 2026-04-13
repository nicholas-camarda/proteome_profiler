is_nonempty_string <- function(x) {
    is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
}

is_absolute_path <- function(path) {
    grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

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

    stop("Could not determine the proteome_profiler repository root from the current working directory.")
}

get_repo_root <- local({
    cached <- NULL

    function() {
        if (is.null(cached)) {
            cached <<- find_repo_root()
        }
        cached
    }
})

first_existing_dir <- function(paths) {
    for (path in paths) {
        if (!is_nonempty_string(path)) {
            next
        }

        expanded <- path.expand(path)
        if (dir.exists(expanded)) {
            return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
        }
    }

    NA_character_
}

get_runtime_root <- function() {
    runtime_root <- first_existing_dir(c(
        Sys.getenv("PROTEOME_PROFILER_RUNTIME_ROOT", unset = ""),
        "~/ProjectsRuntime/proteome_profiler"
    ))

    if (is.na(runtime_root)) {
        get_repo_root()
    } else {
        runtime_root
    }
}

get_cloud_root <- function() {
    first_existing_dir(c(
        Sys.getenv("PROTEOME_PROFILER_CLOUD_ROOT", unset = ""),
        "~/Library/CloudStorage/OneDrive-Personal/phd/projects/VEGFRi and Dox/in-vivo mouse projects/proteome_profiler"
    ))
}

legacy_path_variants <- function(path) {
    variants <- c(path)

    if (startsWith(path, "arrays/")) {
        variants <- c(variants, sub("^arrays/", "projects/", path))
    }

    unique(variants)
}

resolve_project_path <- function(path, must_exist = TRUE) {
    if (!is_nonempty_string(path)) {
        stop("Expected a non-empty path string.")
    }

    if (is_absolute_path(path)) {
        expanded <- path.expand(path)
        if (!must_exist || file.exists(expanded) || dir.exists(expanded)) {
            return(normalizePath(expanded, winslash = "/", mustWork = must_exist))
        }

        stop(sprintf("Path does not exist: %s", expanded))
    }

    search_roots <- unique(c(
        get_repo_root(),
        get_runtime_root(),
        get_cloud_root()
    ))
    search_roots <- search_roots[!is.na(search_roots)]

    variants <- legacy_path_variants(path)
    candidates <- unlist(lapply(search_roots, function(root) {
        file.path(root, variants)
    }))

    existing <- candidates[file.exists(candidates) | dir.exists(candidates)]
    if (length(existing) > 0) {
        return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
    }

    if (!must_exist) {
        return(normalizePath(file.path(get_runtime_root(), variants[[1]]), winslash = "/", mustWork = FALSE))
    }

    stop(sprintf(
        "Could not resolve '%s'. Checked:\n%s",
        path,
        paste(sprintf("- %s", candidates), collapse = "\n")
    ))
}

resolve_example_env_file <- function(project_name, env_filename) {
    resolve_project_path(file.path("projects", project_name, env_filename), must_exist = TRUE)
}
