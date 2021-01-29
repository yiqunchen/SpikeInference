context("Spike estimation")
library(SpikeInference)

test_that("Spike estimation on simple examples should agree with the analytical solution", {
  test_dat <- c(8,4,2,5)
  decay_rate <- 0.5
  tuning_parameter <- 1
  expect_equal(3, spike_estimates(test_dat, decay_rate, tuning_parameter)$spikes)
  expect_equal(4, spike_estimates(c(3, 2.7, 2.43, 2.18, 2.7, 2.43), 
                                  0.9, 0.1)$spikes)
  
})

