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
spike_inference <- structure(function(dat, decay_rate, tuning_parameter, window_size = NULL, sig = NULL, 
                                      sig_estimation = NULL,
                                      return_conditioning_sets = FALSE, return_ci = TRUE, alpha = 0.05
                                    ) {
  stopifnot(decay_rate > 0)
  stopifnot(decay_rate < 1)
  stopifnot(tuning_parameter > 0)
  ## L0 segmentation
  if (is.null(sig)) {
    estimated = TRUE
  } else {
    estimated = FALSE
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
  
  
  if (is.null(sig)) {
       if(is.null(sig_estimation)){
      stop("not implemented")
       }
  }else{
      if(sig_estimation=='MAD'){
        sig = MAD_var_estimator(dat, decay_rate)^2 # get varianxe
      }else if(sig_estimation=='JNFL'){
        sig = JNFL_var_estimator(dat, decay_rate)^2
        if(sig==0){
          sig = (MAD_var_estimator(dat, decay_rate)^2)
        }
      }else{
        stop("Only MAD or JNFL is implemented")
      }
    }
    
    out_fpop_inference <- .fpop_inference(dat, decay_rate, tuning_parameter, window_size, sig,
                                          return_conditioning_sets, return_ci, alpha)
    
    if (return_conditioning_sets) { 
      conditioning_sets = fpop_inference_intervals_formatter(out_fpop_inference[[2]])
    } else { 
      conditioning_sets = NA
    }
    
    if (return_ci){
      out <- list(
        dat = dat,
        decay_rate = decay_rate,
        call = match.call(),
        tuning_parameter = tuning_parameter,
        window_size = window_size,
        sig = sig, 
        change_pts = out_fpop_inference[[1]][, 1], # thj+1 - thj \neq 0
        pvals = out_fpop_inference[[1]][, 2],
        LCB = out_fpop_inference[[1]][, 4],
        UCB = out_fpop_inference[[1]][, 5],
        conditioning_sets = conditioning_sets, 
        estimated_variance = estimated,
        alpha = alpha
      )
    }else{
    out <- list(
      dat = dat,
      decay_rate = decay_rate,
      call = match.call(),
      tuning_parameter = tuning_parameter,
      window_size = window_size,
      sig = sig, 
      change_pts = out_fpop_inference[[1]][, 1], # thj+1 - thj \neq 0
      pvals = out_fpop_inference[[1]][, 2],
      conditioning_sets = conditioning_sets, 
      estimated_variance = estimated
    )
    }
    class(out) <- "SpikeInference"
    return(out)
  })
