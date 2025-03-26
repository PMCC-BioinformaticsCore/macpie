library(unikn)
library(ggplot2)
library(colorspace)

n_colors <- 40  # Adjust as needed

# Convert unikn colors to HCL space for better spacing
unikn_colors <- usecol(pal_unikn_pref)  # Get base colors
unikn_hcl <- as(hex2RGB(unikn_colors), "polarLUV")  # Convert to HCL color space

# Select well-separated colors using farthest-point sampling
selected_indices <- seq(1, length(unikn_colors), length.out = n_colors)  # Even spacing
distinct_colors <- unikn_colors[selected_indices]  # Pick spaced-out colors

# Shuffle slightly to break clusters (optional)
set.seed(2)
distinct_colors <- sample(distinct_colors)

pal1 <- c(rev(pal_seeblau), "white", pal_bordeaux)

# Replace colors with unikn palettes
macpie_colours <- list(
  'discrete' = distinct_colors,
  'continuous' = colorRampPalette(c(pal_signal[[1]], "white", Karpfenblau[[1]]))(100),  # Signal yellow to green
  'scale_3' = c(pal1[[11]], "white", Karpfenblau[[1]]),  # Same 3-color scale
  'discrete_400' = hcl.colors(400, palette = "Zissou1", rev = TRUE)  # Using a diverging palette
)


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
      plot.background = element_rect(fill = NULL, colour = "white"),
      panel.background = element_rect(fill = "white"),
      
      # Axis settings
      axis.line = element_line(colour = "black", linewidth = 1),
      axis.text.x = element_text(size = rel(0.8), angle = x_labels_angle, hjust = hjust_, vjust = vjust_),  # Keep x-axis ticks
      axis.text.y = element_text(size = rel(0.8)),  # Keep y-axis ticks
      #axis.title.x = element_blank(),  # Remove x-axis label
      #axis.title.y = element_blank(),  # Remove y-axis label
      axis.ticks = element_line(colour = "black", linewidth = 0.8),  # Slightly thinner ticks
      
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

