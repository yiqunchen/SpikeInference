

sim <- simulate_ar1(n = 10000, gam = 0.95, poisMean = 0.01, sd = 0.15, seed = 1)

spike_estimates(dat, decay_rate, tuning_parameter, functional_pruning_out = FALSE)
spike_inference(dat, decay_rate, tuning_parameter, window_size = NULL, sig = NULL, return_conditioning_sets = FALSE)