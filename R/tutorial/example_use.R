

sim <- simulate_ar1(n = 10000, gam = 0.95, poisMean = 0.01, sd = 0.15, seed = 1)
plot(sim)

fit_spike <- spike_estimates(dat = sim$fl, decay_rate = 0.95, 
                             tuning_parameter = 1, functional_pruning_out = FALSE)

inference_spike <- spike_inference(dat = sim$fl, decay_rate = 0.95, tuning_parameter = 1, window_size = 10, 
                sig = 0.15, return_conditioning_sets = TRUE)

plot(inference_spike, thj = inference_spike$change_pts[1])

# 1-sided: 
# 3.920091e-04 3.510919e-07 4.734645e-07 1.297679e-08 5.090982e-07 1.239837e-08 7.686889e-07
# 2-sided:
# 7.840182e-04 7.021838e-07 9.469291e-07 2.595358e-08 1.018196e-06 2.479674e-08 1.537378e-06


sim <- simulate_ar1(n = 100000, gam = 0.95, poisMean = 0.01, sd = 0.15, seed = 1)
system.time(fit_spike <- spike_estimates(dat = sim$fl, decay_rate = 0.95, tuning_parameter = 1, functional_pruning_out = FALSE))

system.time(inference_spike <- spike_inference(dat = sim$fl, decay_rate = 0.95, tuning_parameter = 1, window_size = 10, 
                                               sig = 0.15, return_conditioning_sets = TRUE))

system.time(inference_spike <- spike_inference(dat = sim$fl, decay_rate = 0.95, tuning_parameter = 1, window_size = 10, 
                                               sig = 0.15, return_conditioning_sets = FALSE))

system.time(FastLZeroSpikeInference::estimate_spikes(sim$fl, 0.95, 1, FALSE))




sim_one_changepoint <- function(n, gam, seed = 1, sig = 0.015) {
  stopifnot(gam > 0)
  stopifnot(gam < 1)
  set.seed(seed)
  thj <- floor(n / 2)
  calc <- 3 * gam ^ (0:thj)
  underlying_conc <- c(calc, calc)
  y <- underlying_conc + rnorm(length(underlying_conc), sd = sig)
}

n <- 5
gam <- 0.9
h <- 2
sig <- 0
y <- sim_one_changepoint(n, gam, seed = 1, sig = sig)
spike_inference(dat = y, decay_rate = 0.9, tuning_parameter = 0.1, window_size = 2, 
                sig = 0, return_conditioning_sets = FALSE)
