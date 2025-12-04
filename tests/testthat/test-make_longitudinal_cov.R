test_that("make_longitudinal_cov builds correct covariance matrix", {
  Sig <- make_longitudinal_cov(n_time = 4, lags = c(0.4, 0.2, 0.1))

  # Check dimensions
  expect_equal(dim(Sig), c(20, 20))

  # Diagonal should be 1
  expect_equal(diag(Sig), rep(1, 20))

})
