#' Estimation and inference for AR1 spike problem using an L0 penalty
#'
#' @param dat Numeric vector; observed data.
#' @param decay_rate Numeric; specified AR(1) decay rate \eqn{\gamma}, a 
#' number between 0 and 1 (non-inclusive).
#' @param tuning_parameter Numeric; tuning parameter \eqn{\lambda} for 
#' L0 spike estimation, a nonnegative number.
#' @param window_size Numeric; window size for fixed window hypothesis testing,
#' a nonnegative integer.
#' @param sig Numeric; noise variance for the observed data, a nonnegative number.
#'  If unknown (NULL), sample variance of residuals is used instead.
#' @param return_conditioning_sets. Logical; Should the conditioning set S be returned?
#' @param return_ci Logical; if TRUE, the confidence interval for the change in calcium 
#' is computed and returned.
#' @param two_sided Logical; if TRUE, a 2-sided p-value is computed and returned.
#' @param alpha Numeric; significance level for the hypothesis test, 
#' a number between 0 and 1 (non-inclusive).
#' @param mu Numeric; parameter for the test  \eqn{\nu^{T}c\leq 0}.

#' @return Returns a list with elements:
#' @return \code{spikes} the set of spikes,
#' @return \code{pvals}  p-values associated with each spike,
#' @return \code{LCB} lower confidence band for each spike,
#' @return \code{UCB} upper confidence band for each spike.

#' @details
#' Consider the AR(1) generative model \deqn{Y_t = c_t + \epsilon_t, \epsilon_t \sim N(0, \sigma^2)}
#' where \eqn{c_t = \gamma * c_{t-1} + z_t} and \eqn{z_t ~ Pois(poisMean)}. In words,
#' this says between spikes (when \eqn{z_t=0}), calcium decays exponentially at a known rate \eqn{\gamma\in(0,1)}, 
#' which is taken to be known. Further denote the locations of true spikes, \eqn{t:z_t\geq 0} as 
#' \eqn{0 = \tau_0 < \tau_1 < \ldots < \tau_K < \tau_{K+1} = T}.
#'
#' Estimation:
#'
#' This function first estimates spikes via L0 penalty based on
#' noisy observations \eqn{y_t, t = 1,  \ldots, T} by solving the optimization problem
#'
#' \deqn{\hat{c}_1, \cdots, \hat{c}_T \in minimize{c_1,...,c_T\geq 0} 0.5 \sum_{t=1}^T ( y_t - c_t )^2 + \lambda \sum_t=2^T 
#' 1(c_t != \gamma c_t-1) }.
#' 
#' Estimated spikes correspond to the time t such that estimated calcium does not decay exponentially, i.e.,
#' \eqn{\{\cdots,\hat{\tau}_j,\cdots\}  = \{t: \hat{c}_{t+1} -\gamma c_{t} \neq 0\} }.
#'
#' See Jewell et al. (2019), Rigaill, G. (2015), 
#' and Maidstone, R. et al. (2017) for more information and discussion of
#' the detailed algorithm.
#'
#'
#' Inference:
#' 
#' For each estimated changepoint \eqn{\hat{\tau}_j}, we test the null
#' hypothesis that the calcium is decaying exponentially *near* \eqn{\hat{\tau}_j}.
#' In particular, near is defined based on a fixed window around
#' \eqn{\hat{\tau}_j} with left endpoint \deqn{tL = max(1, \hat{\tau}_j - window_size + 1),} and
#' right endpoint \deqn{tR = min(T, \hat{\tau}_j + window_size).}
#' 
#' We then test the null hypothesis \eqn{H_0: \sum_t v_t * \mu_t = 0}. 
#' for a linear contrast vector \eqn{v in R^T} defined as \eqn{v_t = 0}, if
#' \eqn{t < tL} or \eqn{t > tR}; 
#' \eqn{v_t = \gamma \cdot (\gamma^2-1)/(\gamma^{2(tL-\hat{\tau}_j)}-\gamma^2) \cdot \gamma^{t-\hat{\tau}_j}}, if
#' \eqn{tL \le t \le \hat{\tau}_j}; 
#' and \eqn{v_t = (\gamma^2-1)/(\gamma^{2(tR-\hat{\tau}_j)}-1) \cdot \gamma^{t-(\hat{\tau}_j}+1)}, if
#' \eqn{\hat{\tau}_j + 1 \le t \le tR}. 
#' 
#' Our framework allows us to efficiently compute p-values based on the 
#' test statistic \eqn{t(v) * y = \sum_t v_t \cdot y_t} for conditioning set \eqn{\{ \hat{\tau}_j in
#' M(y'(\phi)) \} }, where \eqn{\phi = \sum_t v_t * Y_t}. 
#' 
#' In Chen et al. (2020+), we show that the p-value corresponding
#' to the test \eqn{H0: \sum_t v_t * c_t = \mu} can be written as
#' \deqn{Pr(|\phi| \ge |t(v) * y| | \phi in S)} for a conditioning set S.
#'
#' This function computes p-values for the following test statistic and
#' conditioning set S. In what follows, \eqn{y'(\phi)} is a perturbation of the
#' observed data y. See Proposition 1 of Chen et al. (2020+) for additional
#' details.
#'
#' \deqn{Pr(\phi >= t(v) * y |  \hat{\tau}_j in M(y'(\phi))}
#'
#' Since the observation \eqn{y_t} is normally distributed, \eqn{\phi | S} follows a truncated 
#' normal distribution, and we obtain the (one-sided or two-sided) p-value associated with  \eqn{ \hat{\tau}_j}
#' by computing the set S and find the CDF of the resulting truncated normal distribution.
#'
#' By the duality of hypothesis testing and confidence interval, we can construct a \eqn{1-\alpha} confidence interval for
#' the parameter \eqn{t(v) * c = \sum_t v_t \cdot c_t}. This function computes the lower confidence band and upper confidence
#' band using a bisection method. See Proposition 8 of Chen et al. (2020+) for additional details.
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
#'  tuning_parameter = LAMBDA, window_size = 2,
#'   sig = sigma*sigma, return_ci = TRUE)

#' @references
#'
#' Chen, Y., Jewell, S., and Witten, D. (2020+). Uncertainty quantification for
#' spikes estimated from calcium imaging data. Technical report.
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
                                      window_size, sig = NULL, 
                                      return_conditioning_sets = FALSE, return_ci = FALSE, 
                                      two_sided = FALSE, alpha = 0.05, mu = 0, lower_trunc = -Inf
                                    ){
  stopifnot(decay_rate > 0)
  stopifnot(decay_rate < 1)
  stopifnot(tuning_parameter > 0)
  stopifnot(sum(is.na(dat))==0)
  
  ## L0 
  if (is.null(sig)){
    estimated = TRUE
  } else {
    estimated = FALSE
    stopifnot(sig>0)
  }
  
  if (is.null(sig)){
    
    fit_spike <- SpikeInference::spike_estimates(dat, decay_rate,
                                                 tuning_parameter,
                                                 functional_pruning_out = FALSE)
    sig <- var(dat-fit_spike$estimated_calcium)
    }
    
    #sig <- max(sig,1e-5) #truncate sigma2 for numerical stability
    
    out_fpop_inference <- .fpop_inference(dat, decay_rate, tuning_parameter, window_size, sig,
                                          return_conditioning_sets, return_ci, two_sided, alpha,
                                          mu,lower_trunc)
    
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
        call = match.call(),
        tuning_parameter = tuning_parameter,
        window_size = window_size,
        sig = sig, 
        change_pts = NA, # thj+1 - thj \neq 0
        spikes = NA,
        pvals = NA,
        LCB = NA,
        UCB = NA,
        conditioning_sets = NA, 
        estimated_variance = estimated,
        two_sided = two_sided,
        alpha = alpha,
        mu = mu,
        lower_trunc= lower_trunc
      )
    }else{
      if (return_ci){
        out <- list(
          dat = dat,
          decay_rate = decay_rate,
          call = match.call(),
          tuning_parameter = tuning_parameter,
          window_size = window_size,
          sig = sig, 
          change_pts = pmax(out_fpop_inference[[1]][1:n_cps, 1],1), # thj+1 - thj \neq 0
          spikes = pmax(out_fpop_inference[[1]][1:n_cps, 1],1),
          pvals = out_fpop_inference[[1]][1:n_cps, 2],
          LCB = out_fpop_inference[[1]][1:n_cps, 4],
          UCB = out_fpop_inference[[1]][1:n_cps, 5],
          conditioning_sets = conditioning_sets, 
          estimated_variance = estimated,
          two_sided = two_sided,
          alpha = alpha,
          mu = mu,
          lower_trunc= lower_trunc
        )
      }else{
      out <- list(
        dat = dat,
        decay_rate = decay_rate,
        call = match.call(),
        tuning_parameter = tuning_parameter,
        window_size = window_size,
        sig = sig, 
        change_pts = pmax(out_fpop_inference[[1]][1:n_cps, 1],1), # thj+1 - thj \neq 0
        spikes = pmax(out_fpop_inference[[1]][1:n_cps, 1],1),
        pvals = out_fpop_inference[[1]][1:n_cps, 2],
        conditioning_sets = conditioning_sets, 
        estimated_variance = estimated,
        two_sided = two_sided,
        alpha = alpha,
        mu = mu,
        lower_trunc=lower_trunc
      )
      }
    }
    class(out) <- "SpikeInference"
    return(out)
  })



