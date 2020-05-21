#' Estimation and inference for AR1 spike problem 
#'
#' @param dat observed data
#' @param decay_rate specifies AR(1) decay rate
#' @param tuning_parameter L0 segmentation: tuning parameter \eqn{\lambda};
#'   Binary segmentation: number of segments K.
#' @param window_size window size for fixed window hypothesis testing.
#' @param sig noise variance. If unknown, the sample variance is used for to
#'   compute p-values.
#'
#' @return Returns a list with elements:
#' @return \code{change_pts} the set of changepoints
#' @return \code{pvals}  p-values associated with each changepoint
#'
#' @details
#'
#'
#' @examples
#'
#'
#' @seealso \code{\link{changepoint_estimates}},
#'   \code{\link{changepoint_inference}},
#'
#' @references
#'
#' @export
spike_inference <- structure(function(dat, decay_rate, tuning_parameter, window_size = NULL, sig = NULL, return_conditioning_sets = FALSE) {
  stopifnot(decay_rate > 0)
  stopifnot(decay_rate < 1)
  stopifnot(tuning_parameter > 0)
  ## L0 segmentation
  if (is.null(sig)) {
    estimated = TRUE
  } else {
    estimated = FALSE
  }
  
  if (is.null(sig)) {
      stop("not implemented")
      # fit <- changepoint_estimates(dat, "L0", tuning_parameter)
      # sig <- var(dat - fit$estimated_means)
    }
    
    out_fpop_inference <- .fpop_inference(dat, decay_rate, tuning_parameter, window_size, sig, return_conditioning_sets)
    
    if (return_conditioning_sets) { 
      conditioning_sets = fpop_inference_intervals_formatter(out_fpop_inference[[2]])
    } else { 
      conditioning_sets = NA
    }
    
    out <- list(
      dat = dat,
      decay_rate = decay_rate,
      call = match.call(),
      tuning_parameter = tuning_parameter,
      window_size = window_size,
      sig = sig,
      change_pts = out_fpop_inference[[1]][, 1],
      pvals = out_fpop_inference[[1]][, 2],
      conditioning_sets = conditioning_sets, 
      estimated_variance = estimated
    )
    class(out) <- "SpikeInference"
    return(out)
  })