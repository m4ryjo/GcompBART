test_that("sim_lfried returns correct structure and values", {
  set.seed(123)
  sim <- sim_lfried(N = 10, n_time = 4, lags = c(0.4, 0.2, 0.1),
                    sigma = 10, lambda = c(0.25, 0.5, 0.75, 1))

  # Check class and dimensions
  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 10)
  expect_equal(ncol(sim), 22) # 20 predictors + Y + mu

  # Check column names
  expect_true(all(c("Y", "mu") %in% names(sim)))

  # Check Y and mu are numeric
  expect_type(sim$Y, "double")
  expect_type(sim$mu, "double")

  # Check reproducibility
  set.seed(123)
  sim2 <- sim_lfried(N = 10, n_time = 4, lags = c(0.4, 0.2, 0.1),
                     sigma = 10, lambda = c(0.25, 0.5, 0.75, 1))
  expect_equal(sim, sim2)
})
