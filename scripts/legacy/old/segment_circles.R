# Load required libraries
# library(imager)
# library(EBImage)

library(ROpenCVLite)
library(Rvision)
library(GetoptLong)
library(tidyverse)
# devtools::install_github("swarm-lab/Rvision")
# if (!require("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install("EBImage")

# Load the image
# image_path <- "images/VEHICLE flipped-001.tif"
# image_path <- "images/SORAFENIB flipped-002.tif"

images_paths <- c("images/VEHICLE flipped-001.tif", "images/SORAFENIB flipped-002.tif", "images/SOR + DOX flipped-004.tif", "images/SOR + LIS flipped-003.tif")
for (image_path in images_paths) {
    # image_path <- images_paths[1]
    print(image_path)
    nuc <- EBImage::readImage(image_path) 
    nmask <- EBImage::thresh(nuc, w = 10, h = 10, offset = 0.05)
    nmask <- EBImage::opening(nmask, EBImage::makeBrush(5, shape = "disc"))
    nmask <- EBImage::fillHull(nmask)
    nmask <- EBImage::bwlabel(nmask)
    ctmask <- EBImage::opening(nmask, EBImage::makeBrush(5, shape = "disc"))
    cmask <- EBImage::propagate(nuc, seeds = nmask, mask = ctmask)
    res <- EBImage::paintObjects(nuc, nuc, col = "#ff00ff")
    plot(res)

    # str(image)
    image <- EBImage::channel(image, "gray")
    # str(image)
    # plot(image)
    # Background subtraction
    print("Subtracting background...")
    # hist(image)
    image_bg_subtracted <- EBImage::thresh(image, w = 30, h = 30, offset = 0.01)
    # plot(image_bg_subtracted)

    # Blob segmentation
    print("Segmenting circles...")
    image_smooth <- EBImage::gblur(image_bg_subtracted, sigma = 1)
    image_threshold <- EBImage::thresh(image_smooth, w = 30, h = 30, offset = 0.01)
    image_threshold <- EBImage::closing(image_threshold, EBImage::makeBrush(19, shape = "Gaussian"))
    image_labeled <- EBImage::bwlabel(image_threshold)
    plot(image_labeled)

    # Quantify and label blobs
    print("Quantifying and labeling...")
    stats <- EBImage::computeFeatures.moment(image_labeled) %>% as_tibble()
    shape <- EBImage::computeFeatures.shape(image_labeled) %>% as_tibble()
    print(nrow(stats))
    centroids <- stats %>% dplyr::select(m.cx, m.cy)
    areas <- shape$`s.area`
    pixel_densities <- areas / sum(areas)

    # Show the image
    image_name <- str_c(str_split(str_split(image_path, pattern = "/", simplify = TRUE)[, 2], "\\.", simplify = TRUE)[, 1], ".pdf")
    # dir.create("output/seg_res", showWarnings = FALSE)
    # pdf(qq("output/seg_res/@{image_name}"))
    # par(mar = c(0, 0, 0, 0))
    plot(EBImage::combine(image_bg_subtracted, image_labeled))
    # Label the blobs
    text(centroids$m.cx, centroids$m.cy, labels = sprintf("%d\n%.2f", seq_along(pixel_densities), pixel_densities * 100), col = "red", cex = 0.8)
    dev.off()
    print("Done")
}
