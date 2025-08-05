# source_updated.R
# This source file includes modernized functions for the particle detection workflow.
# Key changes:
# - Replaced manual loops with optimized functions from the 'imager' package.
# - Improved efficiency and readability.
# - Removed dependency on the outdated 'raster' package for image operations.

# Ensure the 'imager' package is installed and loaded.
if (!requireNamespace("imager", quietly = TRUE)) {
  install.packages("imager")
}
library(imager)

############################### 1. Standardize Pixel Intensity ##############################
# This function normalizes image pixel values to a range of [0, 1].
# The original function was efficient and remains largely unchanged.
standardize <- function(image_obj) {
  # imager objects can be standardized directly with imager::renorm()
  # but to keep the logic identical to the original, we'll use the same math.
  i_min <- min(image_obj, na.rm = TRUE)
  i_max <- max(image_obj, na.rm = TRUE)
  i_range <- i_max - i_min

  if (i_range == 0) {
    return(image_obj - i_min) # Avoid division by zero if the image is flat
  }
  
  return((image_obj - i_min) / i_range)
}


########################### 2. Get Image Determinant of Hessian #############################
# This function calculates the determinant of the Hessian matrix for each pixel,
# which is useful for blob detection (identifying round objects).
# It correctly uses the imager package.
Hess_det <- function(image_obj) {
  # Ensure the input is a 'cimg' object for imager functions.
  if (!is.cimg(image_obj)) {
    image_obj <- as.cimg(image_obj)
  }
  
  hess <- imhessian(image_obj)
  # The determinant of the Hessian (DoH) is calculated as Ixx*Iyy - Ixy^2
  det_hess <- with(hess, (xx * yy - xy^2))
  
  # Return as a standard R matrix for compatibility with other functions.
  return(as.matrix(det_hess))
}


############################# 3. Get Image Local Maxima Summary #############################
# This function detects local maxima in an image and computes summary statistics.
# It has been significantly updated to be more efficient by removing slow loops
# and using vectorized operations.
get_local_max <- function(image_obj, kernel_window = 3, thresh_coef = 1.3) {
  
  if (!is.cimg(image_obj)) {
    img <- as.cimg(image_obj)
  } else {
    img <- image_obj
  }

  # Calculate the determinant of the Hessian
  det_hess <- Hess_det(img)
  
  # Find local maxima using imager's built-in function.
  # This is much faster than manual iteration.
  local_maxima_df <- as.data.frame(localMaxima(img, w = kernel_window))
  
  # Get coordinates of all pixels for background calculation
  all_pixels_df <- as.data.frame(img, "wide")
  
  # Identify background pixels (those that are not local maxima)
  background_pixels <- dplyr::anti_join(all_pixels_df, local_maxima_df, by = c("x", "y"))
  
  # Calculate background statistics
  background_mean <- mean(background_pixels$value, na.rm = TRUE)
  background_sd <- sd(background_pixels$value, na.rm = TRUE)
  
  # --- Build the summary table for all pixels ---
  # This part is for full analysis but can be slow for large images.
  # For this workflow, we focus on the local maxima only.
  
  # --- Build the summary table for LOCAL MAXIMA only ---
  lm_summary <- local_maxima_df
  names(lm_summary)[names(lm_summary) == 'value'] <- 'intensity'
  
  # Get the Determinant of Hessian value for each local maximum
  lm_summary$DHess <- mapply(function(x, y) det_hess[x, y], lm_summary$x, lm_summary$y)
  
  # Calculate the mean intensity in the neighborhood of each local maximum
  # Using a boxcar filter (mean filter) is equivalent to this.
  mean_filtered_img <- boxblur(img, boxsize = kernel_window)
  lm_summary$kernel_mean_intensity <- mapply(function(x, y) as.matrix(mean_filtered_img)[x, y], lm_summary$x, lm_summary$y)
  
  # Calculate the final metric from the original script
  lm_summary$mean_x_DHess <- lm_summary$kernel_mean_intensity * lm_summary$DHess
  
  # Calculate significance threshold
  lm_summary$significance <- thresh_coef * (lm_summary$intensity - background_mean) / background_sd
  lm_summary$is_sig <- lm_summary$intensity > lm_summary$significance
  
  # The original script returned a list with two tables.
  # We will return just the summary for local maxima as it's the one used.
  return(lm_summary)
}


################################# 4. Plot Image with Points #################################
# This function plots an image and optionally overlays points on it.
# It has been simplified to use imager's plotting functions.
plot_image_with_points <- function(image_obj, points_df = NULL, point_color = "red", ...) {
  
  if (!is.cimg(image_obj)) {
    image_obj <- as.cimg(image_obj)
  }
  
  plot(image_obj, ...)
  
  if (!is.null(points_df)) {
    # Ensure the dataframe has 'x' and 'y' columns
    if (all(c("x", "y") %in% names(points_df))) {
      points(points_df$x, points_df$y, col = point_color, pch = 20, cex = 0.8)
    } else {
      warning("points_df must contain 'x' and 'y' columns.")
    }
  }
}
