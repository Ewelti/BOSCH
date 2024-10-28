# Function to add silhouette to existing plot
add_silhouette <- function (img = NULL, name = NULL, uuid = NULL, upload_img = NULL, filter = NULL, 
                          x = NULL, y = NULL, ysize = NULL, height = NULL, 
                          width = NULL, alpha = 1, color = NA, fill = "black", 
                          horizontal = FALSE, vertical = FALSE, angle = 0, 
                          hjust = 0.5, vjust = 0.5, remove_background = TRUE, 
                          verbose = FALSE) {
  # Install and load required packages
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    install.packages("rsvg")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    install.packages("grid")
  }
  if (!requireNamespace("lifecycle", quietly = TRUE)) {
    install.packages("lifecycle")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    install.packages("png")
  }
  if (!requireNamespace("rphylopic", quietly = TRUE)) {
    install.packages("rphylopic")
  }
  library(rsvg)
  library(grid)
  library(lifecycle)
  library(png)
  library(rphylopic)
  
  # Check for required arguments
  message("Checking required arguments...")
  if (all(sapply(list(img, name, uuid, upload_img), is.null))) {
    stop("One of `img`, `name`, `uuid`, or `upload_img` is required.")
  }
  if (sum(sapply(list(img, name, uuid, upload_img), is.null)) < 3) {
    stop("Only one of `img`, `name`, `uuid`, or `upload_img` may be specified")
  }
  
  # Load image from file if `upload_img` is specified and convert to raster
  if (!is.null(upload_img)) {
    message("Loading image from file...")
    # Convert SVG to PNG and read it as a raster array
    temp_file <- tempfile(fileext = ".png")
    rsvg::rsvg_png(upload_img, file = temp_file)
    img <- png::readPNG(temp_file)
    unlink(temp_file)  # Remove the temporary file after reading
    message("Image successfully loaded and converted.")
  }
  
  # Validation checks
  if (any(alpha > 1 | alpha < 0)) {
    stop("`alpha` must be between 0 and 1.")
  }
  if (any(hjust > 1 | hjust < 0)) {
    stop("`hjust` must be between 0 and 1.")
  }
  if (any(vjust > 1 | vjust < 0)) {
    stop("`vjust` must be between 0 and 1.")
  }
  if (!is.logical(verbose)) {
    stop("`verbose` should be a logical value.")
  }
  
  # Deprecation warning for `ysize`
  if (!is.null(ysize)) {
    lifecycle::deprecate_warn("1.5.0", "add_phylopic(ysize)", 
                              "add_phylopic(height)")
    if (is.null(height)) 
      height <- ysize
  }
  if (!is.null(height) & !is.null(width)) {
    stop("At least one of `height` or `width` must be NULL.")
  }
  
  # Define plot boundaries
  message("Defining plot boundaries...")
  usr <- par()$usr
  usr_x <- if (par()$xlog) 10^usr[1:2] else usr[1:2]
  usr_y <- if (par()$ylog) 10^usr[3:4] else usr[3:4]
  
  # Set default x and y positions if not specified
  if (is.null(x)) {
    mn <- mean(usr[1:2])
    x <- if (par()$xlog) 10^mn else mn
  }
  if (is.null(y)) {
    mn <- mean(usr[3:4])
    y <- if (par()$ylog) 10^mn else mn
  }
  
  # Set default height and width if not specified
  if (is.null(height) && is.null(width)) {
    height <- abs(diff(usr_y)) * 0.1  # Default to 10% of y-axis range if not set
    width <- NA
  } else {
    if (!is.null(height)) {
      base_y <- grconvertY(ifelse(par()$ylog, 1, 0), to = "ndc")
      height <- grconvertY(height, to = "ndc") - base_y
    }
    if (!is.null(width)) {
      base_x <- grconvertX(ifelse(par()$xlog, 1, 0), to = "ndc")
      width <- grconvertX(width, to = "ndc") - base_x
    }
  }
  
  message(paste("Calculated width:", width))
  message(paste("Calculated height:", height))
  
  # Convert x and y coordinates
  message("Converting coordinates...")
  x <- grconvertX(x, to = "ndc")
  y <- grconvertY(y, to = "ndc")
  
  # Render the image
  message("Rendering the image...")
  if (is.null(img)) {
    stop("No image data available for rendering.")
  }
  if (horizontal || vertical) img <- flip_phylopic(img, horizontal, vertical)
  if (angle != 0) img <- rotate_phylopic(img, angle)
  if (is.na(color) || color == "original") color <- NULL
  if (is.na(fill)) {
    fill <- color
    color <- NULL
  }
  if (fill == "original") fill <- NULL
  img <- recolor_phylopic(img, alpha, color, fill, remove_background)
  
  # Display the image as a raster
  grid.raster(img, x = x, y = y, width = width, height = height, 
              hjust = hjust, vjust = vjust)
  message("Image rendering completed.")
}