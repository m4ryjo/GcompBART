#' @title SoftBART Regression
#' @description Fits a Bayesian Additive Regression Trees (SoftBART) model for continuous outcomes.
#' @param X A numeric matrix of predictors for training.
#' @param Y A numeric vector of outcomes for training.
#' @param X_test Optional numeric matrix of predictors for testing.
#' @param num_tree Number of trees in the ensemble (default = 20).
#' @param k Prior parameter controlling shrinkage (default = 2).
#' @param hypers Optional list of hyperparameters (created by `Hypers()`).
#' @param opts Optional list of options (created by `Opts()`).
#' @param verbose Logical; if TRUE, prints progress information.
#' @return A list containing:
#' \describe{
#'   \item{mu_train_mean}{Posterior mean predictions for training data.}
#'   \item{mu_test_mean}{Posterior mean predictions for test data (if provided).}
#'   \item{mu_Y}{Training outcome mean.}
#'   \item{sd_Y}{Training outcome standard deviation.}
#'   \item{ecdfs}{List of normalization functions for predictors.}
#'   \item{opts}{Options used in the model.}
#'   \item{forest}{The fitted Forest object.}
#' }
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(200), nrow = 20)
#' Y <- rnorm(20)
#' fit <- gc_softbart_regression(X, Y, num_tree = 50)
#' str(fit)
#' }
#' @export
#'
gc_softbart_regression <- function(X, Y, X_test = NULL, num_tree = 20, k = 2,
                                   hypers = NULL, opts = NULL, verbose = TRUE) {
  # Ensure numeric outcome
  stopifnot(is.numeric(Y))

  # Normalize Y
  mu_Y <- mean(Y)
  sd_Y <- sd(Y)
  Y_train <- (Y - mu_Y) / sd_Y
  X_train <- X

  # Set up hypers
  if (is.null(hypers)) {
    hypers <- Hypers(X = X_train, Y = Y_train, normalize_Y = FALSE)
  } else {
    hypers$X <- X_train
    hypers$Y <- Y_train
  }

  hypers$sigma_mu <- 3 / k / sqrt(num_tree)
  hypers$num_tree <- num_tree
  hypers$group <- (1:ncol(X_train) - 1)

  # Set up opts
  if (is.null(opts)) {
    opts <- Opts()
  }
  opts$num_print <- .Machine$integer.max

  # Normalize X
  make_01_norm <- function(x) {
    a <- min(x)
    b <- max(x)
    function(y) (y - a) / (b - a)
  }

  ecdfs <- vector("list", ncol(X_train))
  for (i in seq_along(ecdfs)) {
    xi <- X_train[, i]
    if (length(unique(xi)) == 1) {
      ecdfs[[i]] <- identity
    } else if (length(unique(xi)) == 2) {
      ecdfs[[i]] <- make_01_norm(xi)
    } else {
      ecdfs[[i]] <- ecdf(xi)
    }
  }

  for (i in seq_along(ecdfs)) {
    X_train[, i]
    if (!is.null(X_test)) {
      X_test[, i]
    }
  }

  # Make forest
  reg_forest <- MakeForest(hypers, opts)

  # Warmup
  for (i in 1:opts$num_burn) {
    mu <- reg_forest$do_gibbs(X_train, Y_train, X_train, 1)
  }

  # Accumulate mean predictions
  mu_train_sum <- rep(0, length(Y_train))
  mu_test_sum <- if (!is.null(X_test)) rep(0, nrow(X_test)) else NULL
  sigma <- numeric(opts$num_save)

  for (i in 1:opts$num_save) {
    for (j in 1:opts$num_thin) {
      mu <- reg_forest$do_gibbs(X_train, Y_train, X_train, 1)
    }
    sigma[i] <- reg_forest$get_sigma() * sd_Y

    mu_train_sum <- mu_train_sum + mu
    if (!is.null(X_test)) {
      mu_test_sum <- mu_test_sum + reg_forest$do_predict(X_test)
    }
  }

  mu_train_mean <- mu_train_sum / opts$num_save
  mu_test_mean <- if (!is.null(mu_test_sum)) mu_test_sum / opts$num_save else NULL

  # Return minimal output
  out <- list(
    mu_train_mean = mu_train_mean * sd_Y + mu_Y,
    mu_test_mean = if (!is.null(mu_test_mean)) mu_test_mean * sd_Y + mu_Y else NULL,
    mu_Y = mu_Y,
    sd_Y = sd_Y,
    sigma = mean(sigma),
    ecdfs = ecdfs,
    opts = opts,
    forest = reg_forest
  )

  class(out) <- "softbart_regression"
  return(out)
}
