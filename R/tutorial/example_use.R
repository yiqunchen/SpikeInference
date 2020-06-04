

sim <- simulate_ar1(n = 10000, gam = 0.95, poisMean = 0.01, sd = 0.15, seed = 1)
plot(sim)

fit_spike <- spike_estimates(dat = sim$fl, decay_rate = 0.95, 
                             tuning_parameter = 1, functional_pruning_out = FALSE)

inference_spike <- spike_inference(dat = sim$fl, decay_rate = 0.95, tuning_parameter = 1, window_size = 5, 
                sig = 0.15, return_conditioning_sets = TRUE)


plot(inference_spike, thj = inference_spike$change_pts[1])
