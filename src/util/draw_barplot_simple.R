#' @title draw_barplot_simple
#' @description Draw a barplot with error bars.
#' @param data the data frame with the data to plot.
#' @param title the title of the plot.
#' @param x_var the variable to plot in the x axis.
#' @param y_var the variable to plot in the y axis.
#' @param error_var the variable to plot as error bars.
#' @param title_y the title of the y axis.
#' @param title_x the title of the x axis.
#' @param legend_title the title of the legend.
#' @param legend_position the position of the legend.
#' @param breaks the breaks of the x axis.
#' @param break_labels the labels of the breaks in the y axis.
#' @param colors the colors of the bars.
#' @author Leonardo Brait
#' @date 2023
library(ggplot2)
library(ggpubr)

draw_barplot_simple <- function(
  data,
  x_var,
  y_var,
  error_var = NULL,
  title = "",
  title_x = "",
  title_y = "",
  legend_title = "",
  legend_position = "none",
  breaks = NULL,
  break_labels = NULL,
  limits = NULL,
  colors) {
  
  barplot <- ggplot(data) +
    
    ## Bars
    geom_bar(
      aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]]),
      stat = "identity",
      width = .5) +
    
    ## Error bars
    geom_errorbar(
      aes(
        x = .data[[x_var]],
        ymin = .data[[y_var]] - .data[[error_var]],
        ymax = .data[[y_var]] + .data[[error_var]]
      ),
      width = 0.3,
      colour = "black",
      alpha = 0.5,
      linewidth = 0.5) +
    
    ## Titles
    ggtitle(title) +
    labs(x = title_x, y = gsub(" ", "\n", title_y)) +
    
    ## Theme and legend
    theme_pubr() +
    theme(
      text = element_text(size = unit(12, "points"), family = "Arial"),
      plot.title = element_text(
        size = unit(15, "points"),
        face = "bold",
        family = "Arial",
        hjust = 0.5),
      strip.background = element_blank(),
      legend.position = legend_position,
      legend.text = element_text(size = unit(12, "points"), family = "Arial"),
      legend.title = element_text(
        face = "bold", size = unit(15, "points"), family = "Arial"
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(
        face = "bold", size = unit(15, "points"), family = "Arial"
      ),
      axis.title.y = element_text(
        face = "bold", size = unit(15, "points"), family = "Arial"
      ),
      legend.spacing.x = unit(0.1, "points"),
      legend.spacing.y = unit(0.1, "points")) +
    
    ## Legend
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(
      title = legend_title, size = unit(15, "points"), shape = 16, nrow = 2
    )) +
    
    ## Breaks
    scale_y_continuous(breaks = breaks, labels = break_labels) +
    coord_cartesian(ylim = limits)
  
  return(barplot)
}
