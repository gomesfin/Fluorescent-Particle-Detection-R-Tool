# workflow_updated.R
# This script performs particle detection on a confocal microscopy image.
# It has been updated to use modern R packages for image analysis,
# primarily 'imager', making the workflow more efficient and robust.

# --- 1. Setup ---
# Load necessary libraries
library(imager)
library(dplyr)
library(ggplot2)

# Source the updated helper functions
source("source_updated.R")

# Define the input image file.
# IMPORTANT: Place your image in the same directory as this script,
# or provide the full path to the image.
image_file <- "c1_green_C0_T4.ome.tif"

# --- 2. Load and Prepare Image ---
# Load the image using the imager package
if (!file.exists(image_file)) {
  stop("Image file not found: ", image_file, ". Please check the path.")
}
image_original <- load.image(image_file)

# The original script crops the image. We will do the same.
# Note: imager coordinates start at 1.
image_cropped <- imsub(image_original, x %in% 600:725, y %in% 450:575)

# Standardize pixel values to be between 0 and 1
image_std <- standardize(image_cropped)

# --- 3. Image Filtering ---
# The goal of filtering is to enhance the features (particles) we want to detect.

# A) Low-pass (box average) filter to smooth the image
# This is equivalent to a mean filter.
image_box_filtered <- boxblur(image_std, boxsize = 3)

# B) High-pass (Gaussian) filter
# The original script used a custom kernel. We will apply it using imager's convolve function.
gaussian_kernel <- matrix(c( 0.0585, 0.0965, 0.0585, 
                             0.0965, 0.1592, 0.0965,
                             0.0585, 0.0965, 0.0585),
                          nrow = 3, byrow = TRUE)
image_gauss_filtered <- convolve(image_box_filtered, as.cimg(gaussian_kernel))

# C) Laplacian of Gaussian (LoG) filter to enhance edges/blobs
log_kernel <- matrix(c( 0, -1,  0, 
                       -1,  4, -1,
                        0, -1,  0),
                     nrow = 3, byrow = TRUE)
# Note: The original LoG kernel was slightly different. A standard LoG is often more effective.
# We will convolve the Gaussian-filtered image with the LoG kernel.
image_log_filtered <- convolve(image_gauss_filtered, as.cimg(log_kernel))

# The original script added the LoG result back to the original. This enhances the particles.
image_enhanced <- image_log_filtered + image_std


# --- 4. Detect Particles (Local Maxima) ---
# Use our updated function to find local maxima on the enhanced image.
# These local maxima correspond to potential particles.
detection_summary <- get_local_max(image_enhanced, kernel_window = 3)

# The key metric for classifying a spot is 'mean_x_DHess'
spot_response <- detection_summary$mean_x_DHess


# --- 5. Analyze and Threshold Results ---
# Plot the distribution of the spot classification response to find a good threshold.
png("Spot_Classification_Response_Updated.png", res = 150, height = 800, width = 1000)
plot(sort(spot_response), 
     1:length(spot_response),
     main = "Spot Classification Response",
     xlab = "Spot Response Value (mean_x_DHess)",
     ylab = "Cumulative Count of Local Maxima",
     col = "#8b0000",
     type = "l",
     lwd = 2,
     bty = "l")
grid()

# Set a threshold to identify true particles.
# A common method is to use a high quantile, like the 99th percentile.
# This means we select the top 1% of spots as the most likely particles.
cutoff_value <- quantile(spot_response, 0.99, na.rm = TRUE)

# Add lines to the plot to show the chosen threshold
abline(v = cutoff_value, lty = 2, col = "blue", lwd = 2)
legend("bottomright", legend = paste("99% Quantile Threshold =", round(cutoff_value, 4)), 
       col = "blue", lty = 2, bty = "n")
dev.off()

cat("Chosen threshold value:", cutoff_value, "\n")

# Filter the detected spots to get the final list of particles
final_particles <- detection_summary %>%
  filter(mean_x_DHess > cutoff_value)

cat("Detected", nrow(final_particles), "particles above the threshold.\n")


# --- 6. Visualize Final Results ---
# Save an image of the detected particles overlaid on the original cropped image.
png("Detected_Particles.png", res = 150, height = 800, width = 800)
plot_image_with_points(image_cropped, 
                       points_df = final_particles, 
                       point_color = "red",
                       main = paste("Detected", nrow(final_particles), "Particles"))
dev.off()

# Display the final image in the RStudio plot pane
plot_image_with_points(image_cropped, 
                       points_df = final_particles, 
                       point_color = "red",
                       main = paste("Detected", nrow(final_particles), "Particles"))

# --- 7. Save Particle Data ---
# Save the coordinates and properties of the final particles to a CSV file.
write.csv(final_particles, "detected_particles_summary.csv", row.names = FALSE)

cat("Workflow complete. Results saved to 'detected_particles_summary.csv' and plots saved to the directory.\n")
