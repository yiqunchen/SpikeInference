#' Plot the solution to an L0 segmentation problem
#' @param x output from running estimate_spikes
#' @param xlims optional parameter to specify the x-axis limits
#' @param  arguments to be passed to methods
#' @param ... arguments to be passed to methods
#' @export
#' @import graphics
#'
plot.spike_estimates <- function(x, xlims = NULL, ...) {
  if (sum(is.na(x$estimated_calcium))) {
    warning("Calcium concentration must be estimated before plotting. Automatically running estimate_calcium(fit), however estimated_calicum is not saved.)")
    x <- estimate_calcium(x)
  }
  ind <- 1:length(x$dat)
  rng <- range(c(x$dat, x$estimated_calcium))
  ylims <- rng 
  if (is.null(xlims)){
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlab = "Index")
  } else {
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims, xlab = "Time")
  }
  lines(ind, x$estimated_calcium, col = "blue", lwd = 2)
  
  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
}

MAD_var_estimator <- function(y, decay_rate){
  lag_1_diff <- (y[2:length(y)]-decay_rate*y[1:(length(y)-1)])/sqrt(2)
  sigma_hat <- stats::mad(lag_1_diff)
  return(sigma_hat)
}

JNFL_var_estimator <- function(y, decay_rate){
  lag_1_diff <- (y[2:length(y)]-decay_rate*y[1:(length(y)-1)])/sqrt(2)
  lag_2_diff <- (y[3:length(y)]-decay_rate*decay_rate*y[1:(length(y)-2)])/sqrt(2)
  sigma_hat <- sqrt(max(0, 2*var(lag_1_diff)-var(lag_2_diff)))
  return(sigma_hat)
}

#' Generate the contrast vector for testing the presence of spikes
#' @param n length of the vector
#' @param thj location of the spike of interest
#' @param window_size parameter h
#' @param gam AR-1 decaying parameter
#' @export
construct_v <- function(n, thj, window_size, gam) {
  tL <- max(1, thj-window_size+1)
  tR <- min(n, thj+window_size)
  v <- 0 * numeric(n)
  thj <- max(thj,1)
  if(length(-(thj - tL):0)!=(thj-tL+1)){
    stop(paste0("wrong size!",n," ", thj," ", window_size," ", gam,"\n"))
  }
  v[tL:thj] <- -gam * (gam ^ 2 - 1) / (gam ^ 2 - gam ^ (2 * (tL - thj) ) ) * gam ^ (-(thj - tL):0)
  v[(thj + 1):tR] <- (gam ^ 2 - 1) / (gam ^ (2 * (tR - thj)) - 1) * gam ^ (0:(tR - (thj + 1)))
  return(v)
}

nearest_changepoint <- function(pt, cps) {
  M <- length(pt)
  out <- numeric(M)
  for (i in 1:M) {
    out[i] <- cps[which.min(abs(pt[i] - cps))]
  }
  return(out)
}

loss_function_n_spikes <- function(lam, gcamp, decay_rate, num_spikes_target){
  fit <- spike_estimates(dat = gcamp, decay_rate = decay_rate, tuning_parameter = lam, 
                         functional_pruning_out = FALSE)
  num_spiked_spike <- length(fit$spikes)
  spike_num_diff <- abs(num_spiked_spike-num_spikes_target)
  return(spike_num_diff)
}

###
### 
### 
one_d_binary_search <- function(gcamp, decay_rate, lam_min, lam_max, 
                                num_spikes_target, max_iters=50, tolerance=5){
  iter_i <- 0
  loss_i <- Inf
  verbose <- 0
  
  while(iter_i <= max_iters & loss_i > tolerance){
    
    lam_1 <- (3 * lam_min + lam_max) / 4
    lam_2 <- (3 * lam_max + lam_min) / 4
    
    loss_lam_1 = loss_function_n_spikes(lam_1, gcamp, decay_rate, num_spikes_target)
    loss_lam_2 = loss_function_n_spikes(lam_2, gcamp, decay_rate, num_spikes_target)
    
    if(verbose){
      cat('at iteration i = ', iter_i, "\n")
      cat("loss at lam 1 (" , lam_1 , ") = " , loss_lam_1, "\n")
      cat("loss at lam 2 (" , lam_2 , ") = " , loss_lam_2, "\n")
    }
    
    if (loss_lam_1 < loss_lam_2){
      if (loss_lam_1 < tolerance){
        return(lam_1)
      }
      lam_max = (lam_min + lam_max) / 2
    }else{
      if (loss_lam_2 < tolerance){
        return(lam_2)
      }
      lam_min = (lam_min + lam_max) / 2
    }
      
    iter_i = iter_i + 1
  }
  
  lam_star = (lam_min + lam_max) / 2
  
  return(lam_star)
}


#' search for the right lambda to meet the targeting firing rate
#' @param dat Numeric vector; observed data.
#' @param decay_rate Numeric; specified AR(1) decay rate \eqn{\gamma}, a 
#' number between 0 and 1 (non-inclusive).
#' @param target_firing_rate Numeric; a number between 0 to 1 
#' indicating the average probability of firing
#' @param lam_min Numeric; minimal lambda to consider
#' @param lam_max Numeric; maximal lambda to consider
#' @param max_iters Numeric; maximal iterations to search
#' @param tolerance Numeric; tolerance level for differences in firing rate
#' @export 


estimate_spike_by_spike_number <- function(dat, decay_rate, target_firing_rate, 
                                           lam_min = 1e-6, lam_max = 1, max_iters=50, tolerance=5){
  
  fps = 1 # sampling ratio - default to one
  gcamp = dat$fl
  # transform target firing rate into a number of spikes
  n = length(gcamp)
  nbin_spike_target = floor(target_firing_rate * n / fps)
  lam_star = one_d_binary_search(gcamp, decay_rate, lam_min, lam_max,
                                  nbin_spike_target, max_iters, tolerance)
  fit = spike_estimates(gcamp, decay_rate, lam_star, 
                        functional_pruning_out = FALSE)
  return(fit)
  
}

### 
### Formatting utilities
###
fpop_inference_intervals_formatter <- function(phi_list) {
  phi_outs <- list()
  for (i in 1:length(phi_list)) {
    phi_i <- phi_list[[i]]
    phi_out <- data.frame(phi_i)
    colnames(phi_out) <- c("square", "linear", "constant", "min_mean", "max_mean", "prev_last_mean", "contained", "s")
    phi_out$contained <- phi_out$contained - 1
    phi_outs[[i]] <- phi_out[,c("square", "linear", "constant", "min_mean", "max_mean", "contained")]
  }
  return(phi_outs) 
}


#' Plot simulated data
#' @param x output data from simulate_ar1
#' @param xlims optional parameter to specify the x-axis limits
#' @param  ... to be passed to methods
#' @return Plot with simulated fluorescence (dark grey circles), calcium concentration (dark green line) and spikes (dark green tick marks on x-axis)
#' @examples
#'
#' sim <- simulate_ar1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' @export
#'
plot.simdata <- function(x, xlims = NULL, ...) {
  rng <- range(x$fl)
  ylims <- c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims)
  } else {
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims)
  }
  lines(x$conc, col = "darkgreen", lwd = 2)
  
  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = spike, x1 = spike, y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "darkgreen", lwd = 1)
  }
}

#' Print simulated data
#' @param x simulated data
#' @param ... arguments to be passed to methods
print.simdata <- function(x, ...){
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")
  
  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$fl), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Pois mean \t\t", x$poisMean, "\n")
  cat("SD \t\t\t", x$sd, "\n")
}

#' Print estimated spikes
#'
#' @param x estimated spikes
#' @param ... arguments to be passed to methods
#' @export
print.estimated_spikes <- function(x, ...)
{
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")
  
  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$dat), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Lambda \t\t\t", x$lambda, "\n")
}





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
#' @param ... arguments to be passed to methods
#' @importFrom magrittr %>%
#' 
#' 
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
#'
eval_cost_at.SpikeInference <- function(x, thj, phi) {
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  df <- x$conditioning_sets[[ind]]
  return(eval_cost_at_helper(df, phi))
}


