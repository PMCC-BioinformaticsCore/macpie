utils::globalVariables(c(".data", "Plate_ID"))

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
#' @importFrom readr read_csv
#' @importFrom dplyr %>% mutate filter
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot aes geom_tile labs coord_fixed theme_minimal theme element_text element_rect scale_fill_viridis_c scale_fill_viridis_d
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom grid unit
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' temp_file <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(Well_ID = c("A01", "A02", "B01", "B02"),
#'   Plate_ID = c("P1", "P1","P1", "P1"),
#'   Compound_ID = c("C1", "C2", "C3", "C4")),
#'   temp_file, row.names = FALSE
#' )
#' metadata_heatmap(temp_file)


metadata_heatmap <- function(metadata_file = NULL, legend = "show", output_file = NULL) {
  if (exists("cleaned_metadata", where = .GlobalEnv)) {
    message("Using existing 'cleaned_metadata' from the global environment.")
    metadata <- get("cleaned_metadata", envir = .GlobalEnv)
  } else {
    if (is.null(metadata_file)) {
      stop("No 'cleaned_metadata' found and no file provided. Please provide a CSV file.")
    }
    metadata <- read_csv(metadata_file, show_col_types = FALSE)
  }

  if (!"Well_ID" %in% colnames(metadata)) {
    stop("The metadata is missing a 'Well_ID' column.")
  }

  valid_well_id <- grepl("^[A-P](0[1-9]|1[0-9]|2[0-4])$", metadata$Well_ID)
  if (!all(valid_well_id)) {
    stop("The 'Well_ID' column contains invalid entries.")
  }

  metadata <- metadata %>%
    mutate(
      row = as.numeric(factor(str_extract(metadata$Well_ID, "^[A-P]"))),
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

      legend_setting <- if (legend == "show" && !hide_legends[i]) "right" else "none"

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

    if (length(heatmaps) == 0) {
      stop("No heatmaps were generated. Please check the input metadata and annotation columns.")
    }

    combined_plot <- wrap_plots(heatmaps, ncol = 3) +
      plot_annotation(
        title = paste("macpie Metadata Heatmap Grid - Plate ID:", plate_id),
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
