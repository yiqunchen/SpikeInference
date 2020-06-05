

sim <- simulate_ar1(n = 10000, gam = 0.95, poisMean = 0.001, sd = 0.15, seed = 1)
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

