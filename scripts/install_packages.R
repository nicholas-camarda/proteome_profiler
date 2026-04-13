required_packages <- c(
    "tidyverse",
    "readxl",
    "GetoptLong",
    "RColorBrewer",
    "latex2exp",
    "patchwork",
    "ggprism",
    "ggh4x",
    "progressr",
    "future",
    "furrr",
    "ggforce",
    "readr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) == 0) {
    message("All required packages are already installed.")
} else {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
