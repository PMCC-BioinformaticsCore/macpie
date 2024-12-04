#' Generate Heatmaps of Metadata Function
#'
#' This function generates heatmaps from metadata, either from an existing
#' `cleaned_metadata` object in the global environment or from a CSV file.
#'
#' @param metadata_file Either the existing 'cleaned_metadata' file or a metadata CSV file path.
#' @param legend A character value ("show" or "none") to control whether to display legends.
#' @param output_file A file path and name to save the heatmaps as a graph (png, pdf, or jpg).
#'
#' @return Displays the plot as a ggplot object or saves it as a file (jpg, png, or pdf).
#' @export
#'
#' @examples
#' generate_heatmaps("~/Documents/metadata_test.csv")
#' generate_heatmaps("~/Documents/metadata_test.csv", legend = "none")
#' generate_heatmaps("~/Documents/metadata_test.csv",output_file = "~/Documents/metadata_graph.png")
#'
generate_heatmaps <- function(metadata_file = NULL, legend = "show", output_file = NULL) {
  # Check if 'cleaned_metadata' exists in the global environment
  if (exists("cleaned_metadata", where = .GlobalEnv)) {
    message("Using existing 'cleaned_metadata' from the global environment.")
    metadata <- get("cleaned_metadata", envir = .GlobalEnv)
  } else {
    # If not, check if a file was provided
    if (is.null(metadata_file)) {
      stop("No 'cleaned_metadata' found and no file provided. Please provide a CSV file.")
    }
    # Read the metadata CSV
    metadata <- read_csv(metadata_file, show_col_types = FALSE)
  }

  # Validate the Well_ID column
  if (!"Well_ID" %in% colnames(metadata)) {
    stop("The metadata is missing a 'Well_ID' column.")
  }

  # Check for valid Well_ID format
  valid_well_id <- grepl("^[A-P](0[1-9]|1[0-9]|2[0-4])$", metadata$Well_ID)
  if (!all(valid_well_id)) {
    stop("The 'Well_ID' column contains invalid entries. It must consist of a letter A-P followed by a number between 01 and 24.")
  }

  # Extract row and column numbers from Well_ID
  metadata <- metadata %>%
    mutate(
      row = as.numeric(factor(str_extract(.data$Well_ID, "^[A-P]"))),
      col = as.numeric(str_extract(.data$Well_ID, "\\d+"))
    )

  # Validate the Plate_ID column
  if (!"Plate_ID" %in% colnames(metadata)) {
    stop("The metadata is missing a 'Plate_ID' column.")
  }

  # List of all annotation columns excluding 'Well_ID', 'row', 'Row', 'Column', and 'col'
  annotation_cols <- setdiff(names(metadata), c("Well_ID", "row", "col", "Row", "Column", "Plate_ID"))

  # Check for columns with more than 30 unique values
  hide_legends <- sapply(annotation_cols, function(col) {
    length(unique(metadata[[col]])) > 20
  })

  # Iterate over unique Plate_IDs
  unique_plate_ids <- unique(metadata$Plate_ID)
  for (plate_id in unique_plate_ids) {
    # Filter data for the current Plate_ID
    plate_data <- metadata %>% filter(Plate_ID == plate_id)

    # Generate heatmaps for each annotation column
    heatmaps <- lapply(seq_along(annotation_cols), function(i) {
      annotation <- annotation_cols[i]

      # Determine scale type based on annotation data
      if (is.numeric(plate_data[[annotation]])) {
        scale_func <- scale_fill_viridis_c(option = "C", na.value = "grey80")
      } else {
        scale_func <- scale_fill_viridis_d(option = "C", na.value = "grey80")
      }

      # Determine if legend should be hidden
      legend_setting <- if (legend == "show" && !hide_legends[i]) "right" else "none"

      # Create the heatmap
      ggplot(plate_data, aes(x = col, y = row, fill = .data[[annotation]])) +
        geom_tile(color = "white") +
        scale_func +
        labs(title = annotation, x = "Column", y = "Row") +
        coord_fixed() +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 5),
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

    # Combine heatmaps for this Plate_ID using patchwork
    combined_plot <- wrap_plots(heatmaps, ncol = 3) +
      plot_annotation(
        title = paste("macpie Metadata Heatmap Grid - Plate ID:", plate_id),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold")
        )
      )

    # Save or display the plot
    if (!is.null(output_file)) {
      file_name <- paste0(sub(".(png|jpg|pdf)$", "", output_file), "_Plate_", plate_id, ".", tools::file_ext(output_file))
      ggsave(file_name, combined_plot, width = 10, height = 6 * ceiling(length(heatmaps) / 3), units = "in")
      message("Saved plot for Plate ID: ", plate_id, " to ", file_name)
    } else {
      print(combined_plot)
    }
  }
}
