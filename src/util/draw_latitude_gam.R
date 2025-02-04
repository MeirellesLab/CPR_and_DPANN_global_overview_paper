#' @title draw_latitude_gam
#' @description Draw a latitude plot with a gam model.
#' @param data the data frame with the data to plot.
#' @param breaks the breaks of the y axis.
#' @param title_y the title of the y axis.
#' @param main_title the main title of the plot.
#' @param title_x the title of the x axis.
#' @return the plot.
#' @author Leonardo Brait
#' @date 2023
library(ggplot2)
library(ggpubr)
draw_latitude_gam <- function(data, breaks, title_y, main_title, title_x) {

  plot <-
    ggplot(data = data) +
    theme_pubr() +
    theme(
      text = element_text(size = unit(8, "points"), family = "Arial"),
      axis.title.x = element_text(
        size = unit(8, "points"), face = "bold", family = "Arial"
      ),
      axis.title.y = element_text(
        size = unit(8, "points"), face = "bold", family = "Arial"
      ),
      strip.text.x = element_blank()
    ) +
    geom_smooth(
      aes(x = latitude, y = richness),
      method = "gam",
      formula = y ~ s(x, bs = "tp", k = 20)
    ) +
    scale_x_continuous(breaks = c(-80, -40, 0, 40, 80)) +
    scale_y_continuous(breaks = breaks) +
    scale_color_manual(values = ecosystem_colors) +
    ggtitle(main_title) +
    theme(plot.title = element_text(
      size = unit(8, "points"), face = "bold", hjust = 0.5, family = "Arial"
    )) +
    labs(
      y = title_y,
      x = title_x,
      color = "ecosystem",
      size = "Relative abundance",
      shape = "Microbial group"
    ) +
    theme(
      strip.background = element_blank(),
      legend.position = ""
    )

  return(plot)
}