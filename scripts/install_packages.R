source(file.path("scripts", "helpers", "runtime_setup.R"))

required_packages <- unique(c(
    required_analysis_packages(include_parallel = TRUE),
    required_test_packages()
))

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) == 0) {
    message("All required packages are already installed.")
} else {
    if (!requireNamespace("pak", quietly = TRUE)) {
        install.packages("pak", repos = "https://cloud.r-project.org")
    }

    pak::pkg_install(missing_packages, ask = FALSE)
}
