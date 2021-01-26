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
#' @param ... arguments to be passed to methods
#' @seealso
#' \code{\link{estimate_spikes}},
#' \code{\link{estimate_calcium}},
#' @return Plot with simulated fluorescence (dark grey circles), calcium concentration (dark green line) and spikes (dark green tick marks on x-axis)
#'
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

#' Plot number of spikes vs. tuning parameter
#' @param x output from running estimate_spike_paths
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
#'
#' @seealso
#' \code{\link{estimate_spike_paths}},
#' \code{\link{estimate_spikes}},
#' \code{\link{estimate_calcium}},
#' @export
#' @import graphics
#'
plot.estimated_spike_paths <- function(x, xlims = NULL, ...) {
  ind <- sort.int(x$path_stats[, 1], index.return = T)$ix
  
  xx <- x$path_stats[ind, 1]
  y <- x$path_stats[ind, 2]
  
  if (x$approximate_path) {
    title_str = "Approximate # spikes vs tuning parameter"
  } else {
    title_str = "# spikes vs tuning parameter"
  }
  
  if (is.null(xlims)) {
    plot(xx,y,type="n", ylab = "Number of spikes",
         xlab = "Tuning parameter (lambda)", main = title_str)
  } else {
    plot(xx,y,type="n", xlim = xlims, ylab = "Number of spikes",
         xlab = "Tuning parameter (lambda)", main = title_str)  
  }
  
  segments(xx[-length(xx)],y[-length(xx)],xx[-1],y[-length(xx)])
  points(xx[-length(xx)],y[-length(xx)],pch=16)
  points(xx[-1],y[-length(xx)],pch=1)
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




