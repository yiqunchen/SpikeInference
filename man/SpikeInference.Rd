\name{SpikeInference}
\alias{SpikeInference}
\docType{package}
\title{Quantifying uncertainty for spikes estimated from calcium imaging data}
\description{
SpikeInference is an R package for quantifying uncertainty (i.e., obtaining valid p-values and confidence intervals) for spikes estimated from calcium imaging data via the  penalized algorithm described in Jewell and Witten (2018) and Jewell et al. (2020). Details can be found in our manuscript (Chen et al. 2021). In addition to the manual in R, a detailed tutorial can be found at https://yiqunchen.github.io/SpikeInference/index.html.
}

\details{
Maintainer: Yiqun Chen <yiqunc@uw.edu>

The main functions included in the package are: \code{\link{simulate_ar1}}, \code{\link{spike_estimates}}, and \code{\link{spike_inference}}.
}

\examples{
gam <- 0.98
LAMBDA <- 0.7
sigma <- 0.3
n_length <- 1000
curr_sim <- simulate_ar1(n = n_length, gam = gam, poisMean = 0.01, sd = sigma, seed = 2)
curr_fit_spike <- spike_estimates(dat = curr_sim$fl, decay_rate = gam, tuning_parameter = LAMBDA)
curr_inference_spike <- spike_inference(dat = curr_sim$fl, decay_rate = gam,
 tuning_parameter = LAMBDA, window_size = 2, sig2 = sigma*sigma, return_ci = TRUE)
### a summary of the inferential results
summary(curr_inference_spike)
}

\references{
Chen YT, Jewell SW, Witten DM. (2021) Quantifying uncertainty in spikes estimated from calcium imaging data. arXiv:2103.0781 [statME].

Jewell, S. and Witten, D. (2018). Exact spike train inference via l0 optimization. Ann. Appl. Stat., 12(4):2457–2482.

Jewell SW, Hocking TD, Fearnhead P, Witten DM. Fast nonconvex deconvolution of calcium imaging data. Biostatistics. 2020;21(4):709-726. doi:10.1093/biostatistics/kxy083
}