#' @export
construct_v <- function(n, thj, window_size, gam) {
  tL <- max(1, thj-window_size+1)
  tR <- min(n, thj+window_size)
  v <- 0 * numeric(n)
  v[tL:thj] <- -gam * (gam ^ 2 - 1) / (gam ^ 2 - gam ^ (2 * (tL - thj) ) ) * gam ^ (-(thj - tL):0)
  v[(thj + 1):tR] <- (gam ^ 2 - 1) / (gam ^ (2 * (tR - thj)) - 1) * gam ^ (0:(tR - (thj + 1)))
  return(v)
}




#' @export
nearest_changepoint <- function(pt, cps) {
  M <- length(pt)
  out <- numeric(M)
  for (i in 1:M) {
    out[i] <- cps[which.min(abs(pt[i] - cps))]
  }
  return(out)
}

#' @export
loss_function_n_spikes <- function(lam, gcamp, decay_rate, num_spikes_target){
  fit <- spike_estimates(dat = gcamp, decay_rate = decay_rate, tuning_parameter = lam, 
                         functional_pruning_out = FALSE)
  num_spiked_spike <- length(fit$spikes)
  spike_num_diff <- abs(num_spiked_spike-num_spikes_target)
  return(spike_num_diff)
}

#' @export
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


#' @param u - vector of spike times [n]
#' @param v - vector of spike times [m]
#' @param cost - cost parameter that parameterizes distance
#' @return VP distance, with cost parameter cost, between u and v 
VictorPurpuraDist <- function(u, v, cost){
  lenU <- length(u)
  lenV <- length(v)
  if (cost == 0){
    return(abs(lenU-lenV))
  }else if(cost == Inf){
    return(lenU+lenV)
  }
  # non-degenerate case
  cost_mat <- matrix(0,nrow=lenU+1,ncol=lenV+1)
  cost_mat[,1] <- c(0:lenU)
  cost_mat[1,] <- c(0:lenV)
  
  if (lenU & lenV){
    for (i in c(2:(lenU+1))){
      for (j in c(2:(lenV+1))){
        choice_vec = c(cost_mat[i-1,j]+1, cost_mat[i,j-1]+1, 
                       cost_mat[i-1,j-1]+cost*abs(u[i-1]-v[j-1]))
        cost_mat[i,j] = min(choice_vec)
      }
    }
  }
  
 result_d = cost_mat[lenU+1,lenV+1]
 return(result_d)
  
}


#' @param u - vector of spike times [n]
#' @param v - vector of spike times [n]
#' @return correlation distance between u and v (1-cor(u,v))
#' @export
corrDistance <- function(u,v){
  if(length(u)!=length(v)){
    stop("Incompatible dimensions")
  }
  result = 1-cor(u,v)
  return(result)
}  
  
#' @param u - vector of spike times [n]
#' @param v - vector of spike times [m]
#' @param tau - timescale parameter that parameterizes distance
#' @return VR distance, with timescale tau, between u and v 
#' @export
vanRossumDist <- function(u, v, tau) {
  lenU <- length(u)
  lenV <- length(v)
  d2 <- 0
  s1 <- 0
  
  for (i in 1:lenU) {
    for (j in 1:lenU) {
      s1 <- s1 + exp(-abs(u[i] - u[j]) / tau)
    }
  }
  s3 <- 0
  for (i in 1:lenU) {
    for (j in 1:lenV) {
      s3 <- s3 - 2 * exp(-abs(u[i] - v[j]) / tau)
    }
  }
  s2 <- 0
  for (i in 1:lenV) {
    for (j in 1:lenV) {
      s2 <- s2 + exp(-abs(v[i] - v[j]) / tau)
    }
  }
  return(sum(c(s1, s2, s3)))
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

#' Plot the solution to an L0 segmentation problem
#' @param x output from running estimate_spikes
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
#'
#' @seealso
#' \code{\link{estimate_spikes}},
#' \code{\link{estimate_calcium}},
#' @export
#' @import graphics
#'
plot.SpikeInference_estimated_changes <- function(x, xlims = NULL, ...) {
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
  # estimated calcium curve
  lines(ind, x$estimated_calcium, col = "blue", lwd = 2)

  hh <- 0.01 * diff(ylims)
  for (spike in x$change_pts)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
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




