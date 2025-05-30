#' Color palettes used by macpie
#'
#' A named list of discrete and continuous color palettes for consistent styling in macpie plots.
#'
#' @format A named list with the following elements:
#' * `discrete`: Preferred unikn palette colors.
#' * `discrete_40`: 60-color extension sampled for high cardinality.
#' * `discrete_400`: 400-color diverging Zissou1 palette (reversed).
#' * `high`: Single high-emphasis color.
#' * `low`: Single low-emphasis color.
#' * `divergent`: 100-step divergent palette.
#' * `continuous`: 100-step continuous palette (high to white to low).
#' * `continuous_rev`: Reverse of `continuous`.
#' * `scale_3`: Three-color scale for small categorical data.
#' @import unikn
#' @import ggplot2
#' @import colorspace
#' @examples
#' macpie_colours$discrete
#' barplot(rep(1, 10), col = macpie_colours$discrete_40, border = NA)
#' @export

# n_colors <- 60  # Adjust as needed
# 
# # Get base colors
# unikn_colors <- c(usecol(pal_unikn_pref))
# 
# # Select well-separated colors using farthest-point sampling
# unikn_colors_extended <- colorRampPalette(unikn_colors)(n_colors)
# set.seed(60)
# distinct_colors <- unikn_colors_extended[sample(1:n_colors)]  # Pick spaced-out colors
#barplot(rep(1, 10), col = distinct_colors, border = NA,space = 0, axes = FALSE, main = "Extended unikn palette")

# pal1 <- c(rev(pal_seeblau), "white", pal_bordeaux)

# Replace colors with unikn palettes
macpie_colours <- list(
  'discrete' = usecol(pal_unikn_pref),
  'discrete_40' = {
    n_colors <- 60
    base_cols <- usecol(pal_unikn_pref)
    extended <- colorRampPalette(base_cols)(n_colors)
    set.seed(60)
    extended[sample(seq_len(n_colors))]
  },
  'discrete_400' = hcl.colors(400, palette = "Zissou1", rev = TRUE),
  'high' = pal_signal[[1]],
  'low' = Karpfenblau[[1]],
  'divergent' = colorRampPalette(c(pal_signal[[1]], Karpfenblau[[1]]))(100),
  'continuous' = colorRampPalette(c(pal_signal[[1]], "white", Karpfenblau[[1]]))(100),  # Signal yellow to green
  'continuous_rev' = rev(colorRampPalette(c(pal_signal[[1]], "white", Karpfenblau[[1]]))(100)),  # Signal yellow to green
  'scale_3' = {
    pal1 <- c(rev(pal_seeblau), "white", pal_bordeaux)
    c(pal1[[11]], "white", Karpfenblau[[1]])
    }  # Same 3-color scale
   # Using a diverging palette
)


#' colour theme for macpie plots
#'
#' A custom `ggplot2` theme optimized for macpie figures, with clean backgrounds,
#' consistent text sizing, and configurable axis/legend elements.
#'
#' @param show_x_title Logical; if `TRUE`, displays the x-axis title.
#' @param show_y_title Logical; if `TRUE`, displays the y-axis title.
#' @param legend_position Character; position of the legend (e.g., `'bottom'`, `'none'`).
#' @param x_labels_angle Numeric; rotation angle (in degrees) for x-axis text labels.
#' @return A `ggplot2` theme object that can be added to a plot.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point(color = macpie_colours$high) +
#'   macpie_theme()
#' @export

macpie_theme <- function(show_x_title = TRUE, show_y_title = TRUE, legend_position_ = 'bottom', x_labels_angle = 0) {
  
  # Set axis label alignment based on the x_labels_angle
  if(x_labels_angle != 0){
    hjust_ <- 1
    vjust_ <- 1
  } else {
    hjust_ <- 0.5
    vjust_ <- 0.5
  }
  
  # Apply the theme with the customizations
  theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(face = 'plain', size = rel(1.2), hjust = 0.5),  # Slightly larger title
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank(),
      #plot.background = element_rect(fill = NULL, colour = "white"),
      #panel.background = element_rect(fill = "white"),
      
      # Axis settings
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.text.x = element_text(size = rel(0.8), angle = x_labels_angle, hjust = hjust_, vjust = vjust_),  # Keep x-axis ticks
      axis.text.y = element_text(size = rel(0.8)),  # Keep y-axis ticks
      axis.title.x = element_blank(),  # Remove x-axis label
      axis.title.y = element_blank(),  # Remove y-axis label
      axis.ticks = element_line(colour = "black", linewidth = 0.5),  # Slightly thinner ticks
      
      # Heatmap Text Labels (Make this large)
      plot.subtitle = element_text(size = rel(1.5)),  # Larger text inside heatmap
      
      # Legend settings (Same size as axis ticks)
      legend.title = element_text(face = 'plain', size = rel(0.8)),  # Match tick size
      legend.text = element_text(size = rel(0.8)),  # Match tick size
      legend.background = element_blank(),
      
      # Facet strip background
      strip.background = element_rect(fill = "white")
    )
}


