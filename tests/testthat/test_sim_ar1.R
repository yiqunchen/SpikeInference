context("AR 1 Simulation")
library(SpikeInference)

test_that("Testing behavior of ar1 model without noise", {
  expect_equal(simulate_ar1(3,0.9,0,0,1,c0=1)$fl, c(1.00, 0.90, 0.81))
})
