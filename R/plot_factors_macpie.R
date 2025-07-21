plot_factors_macpie<-function (object, factors = c(1, 2), groups = "all", show_missing = TRUE, 
                               scale = FALSE, color_by = NULL, shape_by = NULL, size_by = NULL,
                               color_name = NULL, shape_name = NULL, dot_size = 2, alpha = 1, 
                               legend = TRUE, stroke = NULL, return_data = FALSE) 
{
  env <- environment(MOFA2::plot_factors)
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    .args <- .args[names(.args) != "factors"]
    return(do.call(plot_factor, c(.args, list(factors = unique(factors)))))
  } else if (length(factors) > 2) {
    .args <- as.list(match.call()[-1])
    p <- do.call(.plot_multiple_factors, .args)
    return(p)
  }
  
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name)) 
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name)) 
    shape_name <- shape_by
  
  factors <- MOFA2::.check_and_get_factors(object, factors)
  Z <- MOFA2::get_factors(object, factors = factors, groups = groups, as.data.frame = TRUE)
  
  #color_by <- .set_colorby(object, color_by)
  #shape_by <- .set_shapeby(object, shape_by)
  #size_by <- .set_colorby(object, size_by)
  
  Z <- Z[complete.cases(Z), ]
  df <- merge(Z, color_by, by = "sample")
  df <- merge(df, shape_by, by = "sample")
  df$shape_by <- as.character(df$shape_by)
  
  # Optional size_by merge
  if (!is.null(size_by)) {
    size_df <- data.frame(sample = size_by$sample, size_val = as.numeric(size_by$color_by))
    df <- merge(df, size_df, by = "sample", all.x = TRUE)
  } else {
    df$size_val <- dot_size  # fallback to default
  }
  
  if (isFALSE(show_missing)) 
    df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  df <- spread(df, key = "factor", value = "value")
  df <- df[, c(colnames(df)[seq_len(4)], "size_val", factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "size", "x", "y"))
  
  if (scale) {
    df$x <- df$x / max(abs(df$x))
    df$y <- df$y / max(abs(df$y))
  }
  
  if (return_data) 
    return(df)
  
  if (is.null(stroke)) {
    stroke <- .select_stroke(N = length(unique(df$sample)))
  }
  
  p <- ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$color_by, 
                      shape = .data$shape_by, size = .data$size)) +
    geom_point(alpha = alpha, stroke = stroke) +
    labs(x = factors[1], y = factors[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.8), color = "black"),
      axis.title = element_text(size = rel(1.1), color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5)
    )
  
  p <- .add_legend(p, df, legend, color_name, shape_name)
  if (!is.null(color_name)) p <- p + labs(fill = color_name)
  if (!is.null(shape_name)) p <- p + labs(shape = shape_name)
  p <- p + guides(size = guide_legend(title = "Dot Size"))
  
  return(p)
}