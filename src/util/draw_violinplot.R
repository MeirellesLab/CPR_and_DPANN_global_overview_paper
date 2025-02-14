#' @title draw_violinplot
#' @description Draw a violin plot with optional boxplots and jittered points.
#' @param data the data frame with the data to plot.
#' @param title the title of the plot.
#' @param x_var the variable to plot on the x-axis.
#' @param y_var the variable to plot on the y-axis.
#' @param title_y the title of the y-axis.
#' @param title_x the title of the x-axis.
#' @param legend_title the title of the legend.
#' @param legend_position the position of the legend.
#' @param breaks the breaks of the y-axis.
#' @param break_labels the labels of the breaks on the y-axis.
#' @param colors the colors of the violin plots.
#' @param add_boxplot logical, whether to add a boxplot inside the violin plot.
#' @param add_jitter logical, whether to add jittered points.
#' @author Diogo Burgos
#' @date 2025
library(ggplot2)
library(ggpubr)

draw_violinplot <- function(
  data,
  title = "",
  x_var,
  y_var,
  title_y = "",
  title_x = "",
  legend_title = "",
  legend_position = "none",
  breaks = NULL,
  break_labels = NULL,
  colors,
  add_boxplot = TRUE,
  add_jitter = TRUE
) {
  violinplot <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    
    ## Violin plot
    geom_violin(trim = TRUE, alpha = 0.7) +
    
    ## Optional boxplot
    if (add_boxplot) geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) else NULL +
    
    ## Optional jittered points
    if (add_jitter) geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) else NULL +
    
    ## Colors
    scale_fill_manual(values = colors) +
    
    ## Theme and legend
    theme_pubr() +
    theme(
      text = element_text(size = unit(8, "points"), family = "Arial"),
      legend.position = legend_position,
      legend.text = element_text(size = unit(8, "points"), family = "Arial"),
      legend.title = element_text(face = "bold", size = unit(8, "points"), family = "Arial"),
      axis.title.x = element_text(face = "bold", size = unit(8, "points"), family = "Arial"),
      axis.title.y = element_text(face = "bold", size = unit(8, "points"), family = "Arial"),
      plot.title = element_text(size = unit(8, "points"), face = "bold", family = "Arial", hjust = 0.5)
    ) +
    
    ## Titles
    ggtitle(title) +
    labs(x = title_x, y = title_y, fill = legend_title) +
    
    ## Y-axis breaks
    scale_y_continuous(breaks = breaks, labels = break_labels)
  
  return(violinplot)
}
