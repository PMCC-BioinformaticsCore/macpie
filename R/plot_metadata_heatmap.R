utils::globalVariables(c(".data", "Plate_ID"))

#' Generate Heatmaps of Metadata Function
#'
#' This function generates heatmaps from metadata, either from an existing
#' object or from a CSV file.
#' @param metadata Metadata object of the class data frame, matrix or tibble.
#' @param metadata_file Path to the metadata CSV file path.
#' @param legend A character value ("show" or "none") to control whether to display legends.
#' @param output_file A file path and name to save the heatmaps as a graph (png, pdf, or jpg).
#'
#' @return Displays the plot as a ggplot object or saves it as a file (jpg, png, or pdf).
#'
#' @importFrom readr read_csv
#' @importFrom dplyr %>% mutate filter
#' @importFrom stringr str_extract
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom grid unit
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 ggplot aes geom_tile labs coord_fixed theme_minimal theme
#' @importFrom ggplot2 element_text element_rect scale_fill_viridis_c scale_fill_viridis_d

#' @examples
#' #Example
#' metadata_file_path <- system.file("extdata", "PMMSq033/PMMSq033_metadata.csv", package = "macpie")
#' metadata<-read_metadata(metadata_file_path)
#' metadata_heatmap(metadata=metadata)
#' @export


plot_metadata_heatmap <- function(metadata = NULL, metadata_file = NULL, legend = TRUE, output_file = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, metadata_file, legend, output_file) {

    #If neither present
    if (is.null(data) && is.null(metadata_file)) {
      stop("No metadata found and no file provided. Please provide a CSV file.")
    }
    #If metadata object present but wrong format
    if (!inherits(as.data.frame(data), "data.frame") && !is.null(data)) {
      stop("Metdata object cannot be read into a data frame.")
    } else if (nrow(as.data.frame(data)) > 0 && !is.null(data)) {
      #If metadata object present but right format
      metadata <- data
    } else if (!is.null(metadata_file)) {
      #If metadata file present
      if (inherits(metadata_file, "character") && file.exists(metadata_file) && length(metadata_file) == 1) {
        tryCatch({
          metadata <- read_csv(metadata_file, show_col_types = FALSE)
        }, error = function(e) {
          stop("Could not load the metadata file, or multiple metadat files provided: ", e$message)
        })
      } else {
        stop("Could not load the metadata file, or multiple metadata files provided.")
      }
    } else {
      stop("Check the format of inputs.")
    }

    if (!(inherits(legend, "logical") && length(legend) == 1)) {
      stop("The format of `legend` should be a single logical value.")
    }

    if (!(inherits(output_file, "NULL") || inherits(output_file, "character"))) {
      stop("The `output_file` should be either empty or provide the name of the output file.")
    }
    return(list(metadata = metadata, legend = legend, output_file = output_file))
  }

  # Validate inputs
  validated <- validate_inputs(metadata, metadata_file, legend, output_file)
  metadata <- validated$metadata
  legend <- validated$legend
  output_file <- validated$output_file

  metadata <- metadata %>%
    mutate(
      # keep row as A-P instead of numeric 1-16
      row = factor(str_extract(metadata$Well_ID, "^[A-P]")),
      col = as.numeric(str_extract(metadata$Well_ID, "\\d+"))
    )

  if (!"Plate_ID" %in% colnames(metadata)) {
    stop("The metadata is missing a 'Plate_ID' column.")
  }

  annotation_cols <- setdiff(names(metadata), c("Well_ID", "row", "col", "Row", "Column", "Plate_ID"))
  if (length(annotation_cols) == 0) {
    stop("No annotation columns found to generate heatmaps. Please check your input data.")
  }

  hide_legends <- sapply(annotation_cols, function(col) {
    length(unique(metadata[[col]])) > 20
  })

  unique_plate_ids <- unique(metadata$Plate_ID)
  for (plate_id in unique_plate_ids) {
    plate_data <- metadata %>% filter(Plate_ID == plate_id)

    heatmaps <- lapply(seq_along(annotation_cols), function(i) {
      annotation <- annotation_cols[i]
      scale_func <- if (is.numeric(plate_data[[annotation]])) {
        scale_fill_viridis_c(option = "C", na.value = "grey80")
      } else {
        scale_fill_viridis_d(option = "C", na.value = "grey80")
      }

      legend_setting <- if (legend && !hide_legends[i]) "right" else "none"

      ggplot(plate_data, aes(x = as_factor(col), y = fct_rev(as_factor(row)), fill = .data[[annotation]])) +
        geom_tile(color = "white") +
        scale_func +
        labs(title = annotation, x = "Column", y = "Row") +
        coord_fixed() +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 4),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.position = legend_setting,
          panel.border = element_rect(color = "black", fill = NA),
          plot.background = element_rect(color = "black")
        )
    })

    if (length(heatmaps) == 0) {
      stop("No heatmaps were generated. Please check the input metadata and annotation columns.")
    }

    combined_plot <- wrap_plots(heatmaps, ncol = 3) +
      plot_annotation(
        title = paste("Metadata Heatmap Grid - Plate ID:", plate_id),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold")
        )
      )

    if (!is.null(output_file)) {
      dir_name <- dirname(output_file)
      if (!dir.exists(dir_name)) {
        stop("The directory for the output file does not exist: ", dir_name)
      }

      file_name <- normalizePath(paste0(
        sub(".(png|jpg|pdf)$", "", output_file),
        "_Plate_", plate_id, ".", tools::file_ext(output_file)
      ), mustWork = FALSE)

      message("Saving plot to: ", file_name)
      ggsave(file_name, combined_plot, width = 10, height = 6 * ceiling(length(heatmaps) / 3), units = "in")
      message("Saved plot for Plate ID: ", plate_id, " to ", file_name)
    } else {
      return(combined_plot)
    }
  }
}
