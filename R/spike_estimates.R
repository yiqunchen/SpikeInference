#' Estimate spikes with an L0 penalty
#'
#' @param dat Numeric vector; observed data.
#' @param decay_rate Numeric; specified AR(1) decay rate \eqn{\gamma}, a 
#' number between 0 and 1 (non-inclusive).
#' @param tuning_parameter Numeric; tuning parameter \eqn{\lambda} for 
#' L0 spike estimation, a non-negative number.
#' @param functional_pruning_out Logical; if TRUE, return
#' cost functions for L0 spike estimation. Defaults to FALSE.
#'
#' @return For L0 spike estimation, returns a list with elements:
#' \itemize{
#' \item \code{estimated_calcium} Estimated calcium levels
#' \item \code{spikes} The set of estimated spikes
#' \item \code{cost} The cost at each time point
#' \item \code{n_intervals} The number of piecewise quadratics used at each point
#' \item \code{piecewise_square_losses} A data frame of optimal cost functions
#'  Cost_s*(mu) for s = 1,..., T.
#'}
#' @details
#'
#' Estimation:
#'  This function estimates spikes via an L0 penalty based on the following 
#'  optimization problem:
#' \deqn{minimize_{c_1,...,c_T\geq 0} \frac{1}{2} \sum_{t=1}^T ( y_t - c_t )^2 + \lambda \sum_{t=2}^{T} 
#' 1(c_t \neq \gamma c_t-1) ,}
#' where \eqn{y_t} is the observed fluorescence at the t-th timestep.
#' 
#' @examples
#' 
#' ### Generate sample data
#' sim <- simulate_ar1(n = 500, gam = 0.998, poisMean = 0.01, sd = 0.05, seed = 1)
#' ### Fit the spike
#' fit_spike <- spike_estimates(sim$fl, decay_rate = 0.998, tuning_parameter = 0.01)
#' ### Plot estimated spikes
#' plot(fit_spike)
#' ### summarize estimated spike times
#' summary_fit_spike <- summary(fit_spike)
#' @references
#' Jewell, S. W., Hocking, T. D., Fearnhead, P., & Witten, D. M. (2019). 
#' Fast nonconvex deconvolution of calcium imaging data. Biostatistics. 
#'
#' Maidstone, R., Hocking, T., Rigaill, G., & Fearnhead, P. (2017). On optimal
#' multiple changepoint algorithms for large data. Statistics and Computing,
#' 27(2), 519-533.
#'
#' Rigaill, G. (2015). A pruned dynamic programming algorithm to recover the
#' best segmentations with 1 to K_max change-points. Journal de la Societe
#' Francaise de Statistique, 156(4), 180-205.
#'
#' @export
spike_estimates <- structure(function(dat, decay_rate, tuning_parameter, functional_pruning_out = FALSE){
  n.data <- length(dat)
  stopifnot(3 <= n.data)
  stopifnot(is.numeric(tuning_parameter))
  stopifnot(0 < decay_rate)
  stopifnot(1 > decay_rate)

  stopifnot(length(tuning_parameter) == 1)
  stopifnot(0 <= tuning_parameter)
  cost_mat_r <- double(n.data)
  end_vec_r <- integer(n.data)
  mean_vec_r <- double(n.data)
  intervals_mat_r <- integer(n.data)
  min_mean <- 0
  max_mean <- Inf

  fpop_out <- .fpop(dat, decay_rate, tuning_parameter, min_mean, max_mean, cost_mat_r, end_vec_r, mean_vec_r, intervals_mat_r)
  end_vec_r <- end_vec_r +  1

  piecewise_square_losses <- NULL
  if (functional_pruning_out) {
    for (d in fpop_out) {
      piecewise_square_losses <- rbind(piecewise_square_losses, data.frame(d))
    }
    colnames(piecewise_square_losses) <- c("square", "linear", "constant", "min_mean", "max_mean", 
                                           "prev_last_mean", "data_i", "s")  
  }
  
  estimated_means = rev(mean_vec_r)
  change_pts = rev(unique(end_vec_r[end_vec_r > 0]))
  out <-list(estimated_means = estimated_means,
        estimated_calcium = estimated_means,
        dat = dat,
        decay_rate = decay_rate, 
        change_pts = change_pts,
        spikes = change_pts,
        call = match.call(),
        tuning_parameter = tuning_parameter,
        cost = cost_mat_r,
        n_intervals = intervals_mat_r,
        end_vec = end_vec_r,
        piecewise_square_losses = piecewise_square_losses)
  class(out) <- "spike_estimates"
  return(out)
})




