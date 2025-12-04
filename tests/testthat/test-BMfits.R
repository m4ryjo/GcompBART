test_that("BMfits runs and returns expected structure", {
  set.seed(123)
  sim <- sim_lfried(N = 50, n_time = 4,lags = c(0.4, 0.2, 0.1))

  BM <- BMfits(sim[, 1:21],
               var.type = c(rep("X0", 5), rep("X", 15), "Y"),
               tgroup = c(rep(1:4, each = 5), 4),
               opts = Opts(num_burn = 10, num_thin = 1, num_save = 10),
               base_hypers = BaseHypers(num_tree = 5))

  # Check it's a list
  expect_type(BM, "list")

  # Check for expected components
  expect_true(all(c("BModels", "n_Y", "n_S", "n_Reg", "n_save", "iterations") %in% names(BM)))

  # Check that BModels is not NULL
  expect_false(is.null(BM$BModels))
})
