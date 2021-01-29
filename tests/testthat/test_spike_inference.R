context("Spike inference")
library(SpikeInference)

test_that("Spike estimation on simple examples should agree with the analytical solution", {
  expect_equal(4, spike_inference(c(3, 2.7, 2.43, 2.18, 2.7, 2.43), 
                                0.9, 0.1)$spikes)
  expect_equal(TRUE, spike_inference(c(3, 2.7, 2.43, 2.18, 2.7, 2.43), 
                                  0.9, 0.1)$pvals<=0.05)
})

