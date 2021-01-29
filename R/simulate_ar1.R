#' Simulate fluorescence trace based on simple AR(1) generative model
#'
#' @details
#' Simulate fluorescence trace based on simple AR(1) generative model
#'
#' \deqn{y_t = c_t + \epsilon_t, \epsilon_t \sim N(0, \sigma^2)}
#'
#' \deqn{c_t = gam * c_{t-1} + s_t}
#'
#' \deqn{s_t ~ Pois(poisMean)}
#'
#' @param n Numeric; number of timesteps
#' @param seed Numeric; random seed
#' @param gam Numeric; AR(1) decay rate
#' @param poisMean Numeric; mean for Poisson distributed spikes
#' @param sd Numeric; standard deviation
#'
#' @return spikes, fluorescence, and calcium concentration
#'
#' @examples
#' sim <- simulate_ar1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#' @import stats
#' 
#' @export
simulate_ar1 <- function(n, gam, poisMean, sd, seed, c0 = 0)
{
  set.seed(seed)
  stopifnot(poisMean>=0)
  stopifnot(sd>=0)
  eta <- numeric(n)
  c <- numeric(n)
  f <- numeric(n)
  for (i in 1:n)
  {
    eta[i] <- rpois(1, poisMean)
    if (i > 1){
      c[i] <- gam * c[i - 1] + eta[i]
    }
    else{
      c[i] <- c0+eta[i]
    }
    f[i] <- c[i] + rnorm(n = 1, mean = 0, sd = sd)
  }
  
  spikesOut <- unique((eta > 0) * c(1:n))
  out <- list(spikes = spikesOut[-1]-1, fl = f, conc = c, call = match.call(),
              gam = gam, poisMean = poisMean, type = "ar1",
              sd = sd, seed = seed)
  class(out) <- "simdata"
  return(out)
}