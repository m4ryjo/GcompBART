test_that("gcompbart runs and returns expected output", {
  set.seed(123)
  sim <- sim_lfried(N = 50, n_time = 4, lags = c(0.4, 0.2, 0.1))

  BM <- BMfits(sim[, 1:21],
               var.type = c(rep("X0", 5), rep("X", 15), "Y"),
               tgroup = c(rep(1:4, each = 5), 4),
               opts = Opts(num_burn = 10, num_thin = 1, num_save = 10),
               base_hypers = BaseHypers(num_tree = 5))

  out <- gcompbart(sim[, 1:21],
                   var.type = c(rep("X0", 5), rep("X", 15), "Y"),
                   J = 100,
                   tgroup = c(rep(1:4, each = 5), 4),
                   BModels = BM)

  expect_type(out, "list")
  expect_true("summary_out" %in% names(out))
  expect_true("y_hat" %in% names(out))
})
