#' Estimation and inference for AR-1 spike problem using an L0 penalty
#' @param dat Numeric vector; observed data.
#' @param decay_rate Numeric; specified AR-1 decay rate \eqn{\gamma}, a 
#' number between 0 and 1 (non-inclusive).
#' @param tuning_parameter Numeric; tuning parameter \eqn{\lambda} for 
#' L0 spike estimation, a non-negative number.
#' @param window_size Numeric; window size for fixed window hypothesis testing,
#' a non-negative integer.
#' @param sig2 Numeric; noise variance for the observed data, a non-negative number.
#'  If unknown (NULL), sample variance of residuals is used instead.
#' @param return_conditioning_sets Logical; Should the conditioning set S be returned?
#' @param return_ci Logical; if TRUE, the confidence interval for the change in calcium 
#' is computed and returned.
#' @param two_sided Logical; if TRUE, a 2-sided p-value is computed and returned.
#' @param alpha Numeric; significance level for the hypothesis test, 
#' a number between 0 and 1 (non-inclusive).
#' @param mu Numeric; parameter for the test  \eqn{\nu^{T}c\leq mu}.

#' @return Returns a list with elements:
#' \itemize{
#' \item \code{spikes} the set of spikes,
#' \item \code{pvals}  p-values associated with each spike,
#' \item \code{LCB} lower confidence band for each spike,
#' \item \code{UCB} upper confidence band for each spike.
#' }

#' @details
#' Consider the AR-1 generative model \deqn{Y_t = c_t + \epsilon_t, \epsilon_t \sim N(0, \sigma^2),}
#' where \eqn{c_t = \gamma  c_{t-1} + z_t} and \eqn{z_t \sim Poisson(poisMean)}. In words,
#' this says between spikes (when \eqn{z_t=0}), calcium decays exponentially at a known rate \eqn{\gamma\in(0,1)}, 
#' which is taken to be known. Further denote the locations of true spikes, \eqn{\{t:z_t\geq 0\}} as 
#' \eqn{\{0 = \tau_0 < \tau_1 < \ldots < \tau_K < \tau_{K+1} = T\}}.
#' 
#' This function first estimates spikes via L0 penalty based on
#' noisy observations \eqn{y_t, t = 1,  \ldots, T} by solving the following optimization problem
#' \deqn{ minimize_{c_1,...,c_T\geq 0} \frac{1}{2} \sum_{t=1}^{T} ( y_t - c_t )^2 + \lambda \sum_{t=2}^{T} 
#' 1(c_t \neq \gamma c_t-1).}
#' Estimated spikes correspond to the time t such that estimated calcium does not decay exponentially, i.e.,
#' \eqn{\{\cdots,\hat{\tau}_j,\cdots\}  = \{t: \hat{c}_{t+1} -\gamma c_{t} \neq 0\} }.
#' 
#' Now suppose that we want to test whether the calcium is exponentially decaying near an estimated spike
#' \eqn{\hat{\tau}_j}; or equivalently, the null hypothesis of the form \eqn{H_{0}:  \nu^T c = 0} versus 
#' \eqn{H_{1}:  \nu^T c> 0} for suitably chosen \eqn{\nu} (see Section 2 in Chen et al. (2021+) for details). 
#' This function computes the following p-value
#' \deqn{P(\nu^T Y \geq \nu^T y |  \hat{\tau}_j \in M(Y), \nu^T Y> 0, \Pi_\nu^\perp Y  =  \Pi_\nu^\perp y)}
#' where \eqn{M(Y)} is the set of spikes estimated from Y via the L0 method, \eqn{\Pi_\nu^\perp} is 
#' the orthogonal projection to the orthogonal complement of \eqn{\nu}. In particular, this p-value controls the 
#' selective Type I error (see Section 3 in  Chen et al. (2021+) for details). 
#' 
#' In addition, we implement a \eqn{1-\alpha} confidence interval for the parameter \eqn{ \nu^T c}, the increase
#' in calcium associated with an estimated spike \eqn{\hat{\tau}_j} (see Section 4 in Chen et al. (2021+) for details).
#'
#'
#' @examples
#'
#' gam <- 0.98
#' LAMBDA <- 0.7
#' sigma <- 0.3
#' n_length <- 1000
#' curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.01, sd = sigma, seed = 2)
#' curr_fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA)
#' curr_inference_spike <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
#'  tuning_parameter = LAMBDA, window_size = 2, sig2 = sigma*sigma, return_ci = TRUE)

#' @references
#'
#' Chen YT, Jewell SW, Witten DM. (2021+) Quantifying uncertainty in spikes estimated from calcium imaging data.
#'  arXiv:2103.0781 [statME].
#'  
#' Jewell, S. W., Hocking, T. D., Fearnhead, P., & Witten, D. M. (2019). 
#' Fast nonconvex deconvolution of calcium imaging data. Biostatistics. 
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
spike_inference <- structure(function(dat, decay_rate, tuning_parameter, 
                                      window_size, sig2 = NULL, 
                                      return_conditioning_sets = FALSE, return_ci = FALSE, 
                                      two_sided = FALSE, alpha = 0.05, mu = 0){
  stopifnot(decay_rate > 0)
  stopifnot(decay_rate < 1)
  stopifnot(tuning_parameter > 0)
  stopifnot(sum(is.na(dat))==0)
  
  ## L0 
  if (is.null(sig2)){
    estimated = TRUE
  } else {
    estimated = FALSE
    stopifnot(sig2>0)
  }
  
  if (is.null(sig2)){
    
    fit_spike <- SpikeInference::spike_estimates(dat, decay_rate,
                                                 tuning_parameter,
                                                 functional_pruning_out = FALSE)
    sig2 <- var(dat-fit_spike$estimated_calcium)
    }
    
    #sig <- max(sig,1e-5) #truncate sigma2 for numerical stability
    
    out_fpop_inference <- .fpop_inference(dat, decay_rate, tuning_parameter, window_size, sig2,
                                          return_conditioning_sets, return_ci, two_sided, alpha,
                                          mu,0)
    
    if (return_conditioning_sets) { 
      conditioning_sets = fpop_inference_intervals_formatter(out_fpop_inference[[2]])
    } else { 
      conditioning_sets = NA
    }
    
    n_cps <- as.numeric(out_fpop_inference[[3]])
    
    if(n_cps==0){
      out <- list(
        dat = dat,
        decay_rate = decay_rate,
        tuning_parameter = tuning_parameter,
        window_size = window_size,
        sig2 = sig2, 
        change_pts = NA, # thj+1 - thj \neq 0
        spikes = NA,
        pvals = NA,
        LCB = NA,
        UCB = NA,
        conditioning_sets = NA, 
        estimated_variance = estimated,
        two_sided = two_sided,
        alpha = alpha,
        mu = mu
      )
    }else{
      if (return_ci){
        out <- list(
          dat = dat,
          decay_rate = decay_rate,
          tuning_parameter = tuning_parameter,
          window_size = window_size,
          sig2 = sig2, 
          change_pts = pmax(out_fpop_inference[[1]][1:n_cps, 1],1), # thj+1 - thj \neq 0
          spikes = pmax(out_fpop_inference[[1]][1:n_cps, 1],1),
          pvals = out_fpop_inference[[1]][1:n_cps, 2],
          LCB = out_fpop_inference[[1]][1:n_cps, 4],
          UCB = out_fpop_inference[[1]][1:n_cps, 5],
          conditioning_sets = conditioning_sets, 
          estimated_variance = estimated,
          two_sided = two_sided,
          alpha = alpha,
          mu = mu
        )
      }else{
      out <- list(
        dat = dat,
        decay_rate = decay_rate,
        tuning_parameter = tuning_parameter,
        window_size = window_size,
        sig2 = sig2, 
        change_pts = pmax(out_fpop_inference[[1]][1:n_cps, 1],1), # thj+1 - thj \neq 0
        spikes = pmax(out_fpop_inference[[1]][1:n_cps, 1],1),
        pvals = out_fpop_inference[[1]][1:n_cps, 2],
        conditioning_sets = conditioning_sets, 
        estimated_variance = estimated,
        two_sided = two_sided,
        alpha = alpha,
        mu = mu
      )
      }
    }
    class(out) <- "spike_inference"
    return(out)
  })



