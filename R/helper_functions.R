macpie_colours <- list(
  'discrete' = c(
    "#3CB371",  # medium sea green
    "#FFD700",  # gold
    "#8A2BE2",  # blue violet
    "#00BFFF",  # deep sky blue
    "#DA70D6",  # orchid
    "#FF8C00",  # dark orange
    "#ADFF2F",  # green yellow
    "#20B2AA",  # light sea green
    "#FF4500",  # orange red
    "#003366",  # Prussian Blue
    "#32CD32",  # lime green
    "#FF1493",  # deep pink
    "#00CED1",  # dark turquoise
    "#FF6347",  # tomato
    "#FF00FF",  # fuchsia
    'black', #black
    "#00FF7F",  # spring green
    "#800080",  # purple
    "#DC143C",  # crimson
    "#8B4513",  # saddle brown
    'white', # white
    "#B22222",   # firebrick
    "#FF6347",  # tomato
    "#B0E0E6",  # powder blue
    "#C71585",  # medium violet red
    "#808000",  # olive
    "#D2691E",  # chocolate
    "#A52A2A",  # brown
    "#FFB6C1",  # light pink
    "#4B0082",  # indigo
    "#F0E68C",  # khaki
    "#E6E6FA",  # lavender
    "#FFFACD",  # lemon chiffon
    "#D3D3D3",  # light gray
    "#FF7F50",  # coral
    "#F5DEB3",  # wheat
    "#98FB98",  # pale green
    "#B0C4DE",  # light steel blue
    "#FFDEAD",  # navajo white
    "#E0FFFF"   # light cyan
  ),
  'continuous' = colorRampPalette(c("darkorange", "white", "darkgreen"))(100),
  'scale_3' = c("darkorange", "white", "darkgreen"),
  'discrete_400' = hcl.colors(400, palette = "Spectral", rev = TRUE)
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
      plot.title = element_text(face = 'bold', size = rel(0.9), hjust = 0.5),
      panel.grid.major.y = element_blank(), #element_line(colour = 'gray'),
      panel.grid.minor.y = element_blank(), #element_line(colour = 'gray'),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.background = element_rect(fill = NULL, colour = 'white'),
      panel.background = element_rect(fill = 'white'),
      # Axis stuff
      axis.line = element_line(colour = 'black', linewidth = 1),
      axis.text = element_text(colour = "black", face = 'bold'),
      axis.text.x = element_text(size = rel(0.75), angle = x_labels_angle, hjust = hjust_, vjust = vjust_),
      axis.text.y = element_text(size = rel(0.75)),
      axis.title.x = element_text(size = rel(ifelse(show_x_title, 1, 0))),
      axis.title.y = element_text(size = rel(ifelse(show_x_title, 1, 0))),
      axis.ticks = element_line(colour = 'black', linewidth = 1.1),
      # Legend stuff
      legend.position = legend_position_,
      legend.margin = margin(6, 6, 6, 6),
      legend.title = element_text(face = 'bold', size = rel(0.6)),
      legend.text = element_text(size = rel(0.5)),
      legend.background = element_blank(),
      # legend.box.background = element_rect(colour = "black")
      strip.background = element_rect(fill = "white")
    )
}
