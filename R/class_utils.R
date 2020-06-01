expand_fpop_intervals <- function(df, PLOT_MIN, PLOT_MAX, ni) {
  n_segments <- dim(df)[[1]]
  df_plot <- NULL
  
  for (seg_i in 1:n_segments) {
    if (df$max_mean[seg_i] >= PLOT_MIN && df$min_mean[seg_i] <= PLOT_MAX) {
      xs <- seq(from = max(PLOT_MIN, df$min_mean[seg_i]), 
                to = min(df$max_mean[seg_i], PLOT_MAX), 
                length.out = ni)
      ys <- df$square[seg_i] * xs ^2 + df$linear[seg_i] * xs + df$constant[seg_i]
      if (!is.null(df$contained)) {
        df_tmp <- data.frame(x = xs, y = ys, contained = df$contained[seg_i])  
      } else { 
        df_tmp <- data.frame(x = xs, y = ys, data_i = df$data_i[seg_i])
      }
      df_plot <- rbind(df_plot, df_tmp)
    }
  }
  return(df_plot)
}

#' Plot optimal cost and/or S as a function of \eqn{\phi}
#' 
#' @param  x ChangepointInference_L0_changepoints_pvals
#' @param thj changepoint 
#' @param plot_cost plot the optimal cost (TRUE), or plot S (FALSE)
#' @param PLOT_MIN minimum phi 
#' @param PLOT_MAX maximum phi 
#' @param ni number of values to calculate the optimal cost at in each segment 
#'
#' @importFrom magrittr %>%
#' @export
plot.SpikeInference <- function(x, thj, plot_cost = TRUE, PLOT_MIN = -10, PLOT_MAX = 10, ni = 1000, ...) {
  if (is.null(x$conditioning_sets)) { 
    stop("Re-run changepoint_inference with parameter 'return_conditioning_sets' set to true")
  }
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  
  col_red <- "#d95f02"
  col_blue <- "#1f78b4"
  
  if (plot_cost) {
    df <- x$conditioning_sets[[ind]]
    df_plot <- expand_fpop_intervals(df, PLOT_MIN, PLOT_MAX, ni)
    par(mfrow = c(1, 1))
    plot(df_plot$x, df_plot$y, pch = 20, cex = 0.1,
         col = ifelse(df_plot$contained == 1, col_blue, col_red),
         xlab = latex2exp::TeX("$\\phi$"), 
         ylab = latex2exp::TeX("Cost$(\\phi)$"))
    legend("bottomright", 
           pch = 20,
           col = c(col_blue, col_red),
           c(latex2exp::TeX("$\\phi$ in  S"), latex2exp::TeX("$\\phi$ not in S")))
    
  } else {
    x$conditioning_sets[[ind]] %>% dplyr::mutate(y = 1, contained = factor(contained, levels = c(0, 1), labels = c("Not in S", "In S"))) %>% 
      ggplot2::ggplot() + 
      ggplot2::geom_rect(ggplot2::aes(xmin = min_mean, xmax = max_mean, ymin = -10, ymax = 10, fill = contained)) +
      ggplot2::coord_cartesian(xlim = c(PLOT_MIN, PLOT_MAX)) + 
      ggplot2::xlab(latex2exp::TeX("$\\phi$")) + 
      ggplot2::ylab('') + 
      ggplot2::scale_fill_manual(values=c(col_red, col_blue)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
                     panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), 
                     legend.position="bottom", legend.title = ggplot2::element_blank())
  }
}



#' Evaluate Cost(phi) for any phi 
#' 
#' @param x ChangepointInference_L0_changepoints_pvals
#' @param thj changepoint 
#' @param phi evaluate cost at pertubation phi 
#' 
#' @export
eval_cost_at.SpikeInference <- function(x, thj, phi) {
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  df <- x$conditioning_sets[[ind]]
  return(eval_cost_at_helper(df, phi))
}