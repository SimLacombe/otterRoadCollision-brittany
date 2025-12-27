### Plots ----------------------------------------------------------------------

simsim_map <- function(background,
                       backgroundfield = NULL,
                       urban = NULL,
                       cultures = NULL,
                       vegetation = NULL,
                       rivers = NULL,
                       roads = NULL,
                       waterbodies = NULL,
                       map_palette = simsim_map_palette()) {
  mabbox <- sf::st_bbox(background)
  
  if (methods::is(background, "SpatRaster")) {
    if (!any(names(background) == backgroundfield))
      stop("Provide appropriate value for backgroundfield")
    out <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = as.data.frame(background, xy = TRUE),
        aes(x = .data$x, y = .data$y, fill = .data[[backgroundfield]]),
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_gradient(low = map_palette$background_low,
                                   high = map_palette$background_high)
  } else {
    out <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = background, fill = map_palette$background_low)
  }
  
  out +
    ggplot2::geom_sf(
      data = cultures,
      fill = map_palette$cultures,
      color = NA,
      alpha = 0.5
    ) +
    ggplot2::geom_sf(
      data = urban,
      fill = map_palette$urban,
      color = NA,
      alpha = 0.5
    ) +
    ggplot2::geom_sf(
      data = vegetation,
      fill =   map_palette$vegetation,
      color = NA,
      alpha = 0.5
    ) +
    ggplot2::geom_sf(data = rivers,
                     color = map_palette$rivers,
                     linewidth = 1.33) +
    ggplot2::geom_sf(data = waterbodies,
                     fill = map_palette$waterbodies,
                     color = NA) +
    ggplot2::geom_sf(data = roads,
                     color = map_palette$roads,
                     linewidth = 1.33) +
    ggplot2::scale_x_continuous(
      limits = c(mabbox["xmin"], mabbox["xmax"]),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(mabbox["ymin"], mabbox["ymax"]),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed()+
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.border = element_rect(fill = NA),
      panel.background = element_rect(fill = map_palette$waterbodies)
    )
}

simsim_scatter <- function(data, xVar, yVar, groupVar = NULL, method = "gam",
                           palette_points = simsim_palette_d_1(factor = 0.25),
                           palette_smooth = simsim_palette_d_1(factor = -0.25), 
                           mode = "light"){
  
  if(is.null(groupVar)) {
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_point(
        aes(
          x = .data[[xVar]],
          y = .data[[yVar]]
        ),
        pch = 21,
        fill = NA,
        show.legend = FALSE,
        color  = "black"
      ) +
      ggplot2::geom_smooth(
        aes(
          x = .data[[xVar]],
          y = .data[[yVar]]
        ),
        method = method
      )
  } else {
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_point(
        aes(
          x = .data[[xVar]],
          y = .data[[yVar]],
          col = .data[[groupVar]],
          fill = .data[[groupVar]]
        ),
        pch = 21,
        fill = NA,
        show.legend = FALSE
      ) +
      scale_color_manual(values = palette_points) +
      ggnewscale::new_scale_color() +
      ggplot2::scale_fill_manual(values = palette_smooth) +
      ggplot2::scale_color_manual(values = palette_smooth) +
      ggplot2::geom_smooth(
        aes(
          x = .data[[xVar]],
          y = .data[[yVar]],
          col = .data[[groupVar]],
          fill = .data[[groupVar]]
        ),
        method = method
      )
    
  }
  p +
    simsim_theme_scatter(mode = mode) +
    ggplot2::xlab("") + ggplot2::ylab("")
}

### Themes ---------------------------------------------------------------------


simsim_theme_scatter <- function(mode = "light"){ 
  mycol <- case_when(
    mode == "light" ~ "black",
    mode == "dark" ~ "#747474"
  )
  myfill <- case_when(
    mode == "light" ~ NA,
    mode == "dark" ~ "#3A3A3A"
  )

  mytext <- case_when(
    mode == "light" ~ "black",
    mode == "dark" ~ "#D9D9D9"
  )
  
  theme_minimal() %+replace%    
    theme(
      panel.background = element_rect(colour = mycol, fill = myfill),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = mycol, linetype = "dotted"),
      text = element_text(color = mytext, face = "bold"),
      axis.text = element_text(color = mytext, face = "bold"),
      legend.position = "bottom"
    )
}

simsim_theme_map <- function(mode = "light"){ 
  mycol <- case_when(
    mode == "light" ~ "lightgrey",
    mode == "dark" ~ "#747474"
  )
  myfill <- case_when(
    mode == "light" ~ NA,
    mode == "dark" ~ "#3A3A3A"
  )
  mytext <- case_when(
    mode == "light" ~ "black",
    mode == "dark" ~ "#D9D9D9"
  )
  
  theme_minimal() %+replace%  
    theme(panel.background = element_rect(color = mycol, fill = myfill),
          panel.grid.major = element_line(color = mycol),
          text = element_text(color = mytext, face = "bold"),
          axis.text = element_text(color = mytext, face = "bold"),
          axis.ticks = element_blank())
}
### Utilities ------------------------------------------------------------------

adjust_colors <- function(pal, factor = 0) {
  rgb <- grDevices::col2rgb(pal) / 255
  if (factor >= 0) {
    adjusted_rgb <- rgb + (1 - rgb) * factor
  } else {
    adjusted_rgb <- rgb * (1 + factor)
  }
  grDevices::rgb(adjusted_rgb[1,], adjusted_rgb[2,], adjusted_rgb[3,], maxColorValue = 1)
}

strip_y <- theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank())

simsim_title <- function(plot, title, mycol){
  plot +
    ggplot2::facet_wrap(~ "...",
                        labeller = ggplot2::labeller(.default = function(x) title)) +
    ggplot2::theme(strip.background = element_blank(),
                   strip.text = element_text(size = 10, face = "bold", color = mycol))
}

### Palettes -------------------------------------------------------------------

simsim_palette_d_1 <- function(n = 4, factor = 0){
  pal <- paletteer::paletteer_d("ggthemes::Tableau_10")[1:n]
  adjust_colors(pal, factor = factor)
}

simsim_map_palette <- function(background_low = "#f5f5dc",
                               background_high = "#c07b3f",
                               river = "lightskyblue",
                               waterbodies = "lightblue",
                               urban = "grey70",
                               roads = "grey30",
                               vegetation = "darkolivegreen2",
                               cultures = "yellow") {
  list(
    background_low = background_low,
    background_high = background_high,
    rivers = river,
    roads = roads,
    waterbodies = waterbodies,
    urban = urban,
    vegetation = vegetation,
    cultures = cultures
  )
}