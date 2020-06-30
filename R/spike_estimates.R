#' Estimate spikes via L0 segmentation
#'
#' @param dat observed data
#' @param decay_rate specifies AR(1) decay rate
#' @param tuning_parameter L0 segmentation: tuning parameter \eqn{\lambda};
#'   Binary segmentation: number of segments K.
#' @param functional_pruning_out boolean specifying whether, in the case of L0
#'   segmentation, cost functions are returned. Defaults to FALSE.
#'
#'
#' @return For L0 segmentation, returns a list with elements:
#' @return \code{estimated_means} estimated means
#' @return \code{change_pts} the set of changepoints
#' @return \code{cost} the cost at each time point
#' @return \code{n_intervals} the number of piecewise quadratics used at each
#'   point
#' @return \code{piecewise_square_losses} dataframe of optimal cost functions
#'   Cost_s*(mu) for s = 1,..., T.
#'
#' @details
#'
#' Estimation:
#'
#' This function estimates changepoints via L0 or binary segmentation based on
#' noisy observations \eqn{y_t, t = \ldots, T}. In the case of L0 segmentation,
#' changepoints are estimated by solving the optimization problem
#'
#' @examples
#'
#' ### Generate sample data
#'
#' @references
#'
#' Jewell, S., Fearnhead, P., and Witten, D. (2019+). Testing for a change in
#' mean after changepoint detection. Technical report.
#'
#' Maidstone, R., Hocking, T., Rigaill, G., & Fearnhead, P. (2017). On optimal
#' multiple changepoint algorithms for large data. Statistics and Computing,
#' 27(2), 519-533.
#'
#' Rigaill, G. (2015). A pruned dynamic programming algorithm to recover the
#' best segmentations with 1 to K_max change-points. Journal de la Societe
#' Francaise de Statistique, 156(4), 180-205.
#'
#'
#' @export
spike_estimates <- structure(function(dat, decay_rate, tuning_parameter, functional_pruning_out = FALSE) {
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
    colnames(piecewise_square_losses) <- c("square", "linear", "constant", "min_mean", "max_mean", "prev_last_mean", "data_i", "s")  
  }
  
  # edit to encode the positive/negative spike information
  estimated_means = rev(mean_vec_r)
  change_pts = rev(unique(end_vec_r[end_vec_r > 0]))
  est_diff_mean <- estimated_means[2:length(estimated_means)]-decay_rate*estimated_means[1:(length(estimated_means)-1)]
  est_spike_size <- est_diff_mean[change_pts]
  
  spike_sign <- rep("Positive", times = length(change_pts))
  spike_sign[est_spike_size<0] <- "Negative"
  
  out <-
      list(
        estimated_means = estimated_means,
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
        piecewise_square_losses = piecewise_square_losses,
        spike_sign = spike_sign
      )
    class(out) <- "SpikeInference_estimated_changes"
    return(out)
})

