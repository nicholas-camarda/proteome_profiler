source(file.path("scripts", "helpers", "runtime_setup.R"))

required_packages <- required_analysis_packages(include_parallel = TRUE)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) == 0) {
    message("All required packages are already installed.")
} else {
    if (!requireNamespace("pak", quietly = TRUE)) {
        install.packages("pak", repos = "https://cloud.r-project.org")
    }

    pak::pkg_install(missing_packages, ask = FALSE)
}
