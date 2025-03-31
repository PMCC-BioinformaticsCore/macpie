library(unikn)
library(ggplot2)
library(colorspace)

n_colors <- 60  # Adjust as needed

# Get base colors
unikn_colors <- c(usecol(pal_unikn_pref))

# Select well-separated colors using farthest-point sampling
unikn_colors_extended <- colorRampPalette(unikn_colors)(n_colors)
set.seed(60)
distinct_colors <- unikn_colors_extended[sample(1:n_colors)]  # Pick spaced-out colors
#barplot(rep(1, 10), col = distinct_colors, border = NA,space = 0, axes = FALSE, main = "Extended unikn palette")

pal1 <- c(rev(pal_seeblau), "white", pal_bordeaux)

# Replace colors with unikn palettes
macpie_colours <- list(
  'discrete' = usecol(pal_unikn_pref),
  'discrete_40' = distinct_colors,
  'discrete_400' = hcl.colors(400, palette = "Zissou1", rev = TRUE),
  'high' = pal_signal[[1]],
  'low' = Karpfenblau[[1]],
  'continuous' = colorRampPalette(c(pal_signal[[1]], "white", Karpfenblau[[1]]))(100),  # Signal yellow to green
  'continuous_rev' = rev(colorRampPalette(c(pal_signal[[1]], "white", Karpfenblau[[1]]))(100)),  # Signal yellow to green
  'scale_3' = c(pal1[[11]], "white", Karpfenblau[[1]])  # Same 3-color scale
   # Using a diverging palette
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

# For multi_plates, when combining plates, genes are different
# Define a helper to fill missing genes with 0s
pad_sparse_matrix <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    # Create a zero matrix for missing genes
    zero_mat <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(mat)
    
    # Combine existing and missing gene matrices
    mat <- rbind(mat, zero_mat)
  }
  
  # Reorder rows to match all_genes
  mat <- mat[all_genes, , drop = FALSE]
  return(mat)
}
