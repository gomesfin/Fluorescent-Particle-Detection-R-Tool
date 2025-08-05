# Project Overview

This project provides an R-based workflow for automatically identifying and analyzing fluorescent probes (particles) in confocal microscopy images. The scripts load an image, apply a series of filters to enhance particle features, detect potential particles by finding local intensity maxima, and then apply a statistical threshold to isolate the most likely candidates.

The code has been updated from its original version to use modern, efficient R packages, primarily `imager`, which significantly speeds up image processing tasks.

# Key Files

* `workflow_updated.R`: The main script that executes the entire analysis pipeline from start to finish.
* `source_updated.R`: A supplementary script containing the helper functions used by the main workflow.
* `c1_green_C0_T4.ome.tif`: An example image file (you should replace this with your own).

# Dependencies

You will need a current version of R and RStudio. The following R packages are required to run the analysis. You can install them by running this code chunk in your R console.

```{r, eval=FALSE, message=FALSE}
# eval=FALSE prevents this chunk from running automatically when the document is "knitted".
# Run this code manually in your console if you don't have these packages.
install.packages(c("imager", "dplyr", "ggplot2"))
```

* **imager**: The core package for image loading, processing, and analysis.
* **dplyr**: Used for efficient data manipulation and filtering.
* **ggplot2**: Used for creating plots (though the current script uses base R plotting, it's good practice to have it).

# How to Run the Analysis

1.  **Set Up Your Project:**
    * Create a new folder for your project.
    * Save `workflow_updated.R` and `source_updated.R` into this folder.
    * Place the confocal image you want to analyze (e.g., in `.tif`, `.png`, `.jpg` format) into the **same folder**.

2.  **Configure the Workflow:**
    * Open `workflow_updated.R` in RStudio.
    * Find the line `image_file <- "c1_green_C0_T4.ome.tif"` and change the filename to match your image file.

3.  **Execute the Script:**
    * With `workflow_updated.R` open, you can run the entire script by clicking the "Source" button in RStudio or by pressing `Ctrl+Shift+Enter`.
    * The script will run automatically, performing all steps and printing progress messages to the console.

# The Analysis Pipeline

The workflow script performs the following steps:

1.  **Load and Prepare Image**: The script loads your specified image file, crops it to a region of interest, and standardizes the pixel intensity to a range of \[0, 1\].

2.  **Image Filtering**: A series of filters are applied to make the particles easier to detect:
    * **Box Blur**: Smooths the image to reduce noise.
    * **Gaussian Filter**: A more advanced smoothing filter.
    * **Laplacian of Gaussian (LoG)**: An edge/blob detection filter that highlights circular features like fluorescent probes.

3.  **Particle Detection**: The script identifies all "local maxima" on the filtered image. These are pixels that are brighter than all their immediate neighbors and are considered potential particle candidates.

4.  **Thresholding**: For each candidate, a "spot response" score is calculated. A threshold is automatically determined by finding the 99th percentile of these scores. Candidates with a score above this threshold are classified as true particles.

5.  **Visualization and Output**:
    * The script generates and saves two plots:
        * `Spot_Classification_Response_Updated.png`: Shows the distribution of spot scores and the chosen threshold.
        * `Detected_Particles.png`: Displays the original cropped image with the final detected particles marked in red.
    * A CSV file named `detected_particles_summary.csv` is created, containing the coordinates and detailed properties of each detected particle.

# License

This project is under the MIT License, as per the original repository.
