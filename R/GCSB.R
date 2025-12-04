#' @title Create Base Hyperparameters for Soft BART models.
#'
#' @description
#' Defines a set of default hyperparameters for soft BART models, including options
#' tailored for longitudinal data. These parameters control aspects of the Bayesian
#' tree ensemble such as prior distributions, shrinkage, and the number of trees.
#' Users can override any defaults by specifying arguments explicitly.
#'
#' @details
#' For longitudinal settings, the `alpha_vec` parameter implements the ordered
#' hyperprior structure described in [Paper Citation]. This structure introduces
#' sparsity based on temporal proximity to the response: predictors measured earlier
#' in time receive stronger shrinkage (higher sparsity), while those closer to the
#' response are penalized less, reflecting the assumption that recent measurements
#' are more informative.
#'
#' The ordering is achieved by setting weights:
#' \deqn{c_j = 1 - \frac{t - j}{t - 1} \times 0.5}
#' where \eqn{t} is the total number of time points and \eqn{j} indexes the predictor's
#' time point. This results in \eqn{c_1 = 0.5} for the earliest predictors and
#' \eqn{0.5 < c_{j-1} < c_j < 1} for later predictors.
#'
#' The returned list is intended to be combined with model-specific arguments
#' (e.g., `X`, `Y`, `tgroup`) when calling the internal `Hypers()` constructor.
#'
#' @param alpha Numeric. Shape parameter for the tree prior. Default = 1.
#' @param eta Numeric. Controls variance shrinkage for time-varying effects. Default = 1.
#' @param phi Numeric or vector. Controls time-varying effect scaling. Default = 1.
#' @param alpha_vec Numeric or vector. Specifies ordered hyperpriors for predictor-specific
#' shrinkage parameters in longitudinal models. By default, all predictors share the same
#' shrinkage level (1). When using the longitudinal prior, `alpha_vec` can be set to reflect
#' temporal ordering, where earlier time points receive stronger shrinkage and later time
#' points less shrinkage, as described above.
#' @param alpha_shape_1 Numeric. First shape parameter for alpha prior. Default = 0.5.
#' @param alpha_shape_2 Numeric. Second shape parameter for alpha prior. Default = 1.
#' @param num_tree Integer. Number of trees in the ensemble. Default = 20.
#' @param ... Additional arguments to override defaults or pass to `Hypers()`.
#'
#' @return A list of hyperparameters.
#'
#' @examples
#' # Create a base hyperparameter list with defaults
#' base_hypers <- BaseHypers()
#'
#' # Override some defaults
#' base_hypers <- BaseHypers(num_tree = 50)
#'
#' # Combine with model-specific arguments
#' X_train <- matrix(rnorm(100), ncol = 5)
#' Y_train <- rnorm(20)
#' tmp_tg <- rep(1:2, each = 10)
#' hypers_args <- c(base_hypers, list(X = X_train, Y = Y_train, tgroup = tmp_tg))
#' hypers <- do.call(Hypers, hypers_args)
#'
#' @export
BaseHypers <- function(alpha = 1, eta = 1, phi = 1,
                       alpha_vec = 1, alpha_shape_1 = 0.5,
                       alpha_shape_2 = 1, num_tree = 20, ...)
{  c(list(alpha = alpha, eta = eta, phi = phi,
          alpha_vec = alpha_vec, alpha_shape_1 = alpha_shape_1,
          alpha_shape_2 = alpha_shape_2, num_tree = num_tree), list(...))
}

#' @title Extended Hypers for Soft BART Models
#'
#' @description
#' Constructs a hyperparameter list by calling \code{SoftBart::Hypers()} for base fields
#' and then extends it with GcompBART-specific components: \code{alpha_vec}, \code{eta}, \code{phi},
#' and time-grouping information (\code{tgroup}, \code{tgroup_size}). This approach preserves
#' SoftBart defaults while adding longitudinal extensions.
#'
#' @details
#' For longitudinal settings, the \code{alpha_vec} parameter implements the ordered hyperprior
#' structure described in [Paper Citation]. This structure introduces sparsity based on
#' temporal proximity to the response: predictors measured earlier in time receive stronger
#' shrinkage (higher sparsity), while those closer to the response are penalized less, reflecting
#' the assumption that recent measurements are more informative.
#'
#' The ordering is achieved by setting weights:
#' \deqn{c_j = 1 - \frac{t - j}{t - 1} \times 0.5}
#' where \eqn{t} is the total number of time points and \eqn{j} indexes the predictor's time point.
#' This results in \eqn{c_1 = 0.5} for the earliest predictors and
#' \eqn{0.5 < c_{j-1} < c_j < 1} for later predictors.
#'
#' @param X Matrix of predictors (same semantics as \code{SoftBart::Hypers}).
#' @param Y Response vector (same semantics as \code{SoftBart::Hypers}).
#' @param group Optional grouping vector for SoftBart trees.
#' @param tgroup Integer vector of length \code{ncol(X)} specifying time-group assignments for predictors.
#'   If \code{NULL}, defaults to all predictors in one group.
#' @param alpha_vec Numeric or vector. Specifies ordered hyperpriors for predictor-specific
#' shrinkage parameters in longitudinal models. By default, all predictors share the same shrinkage level (1).
#' When using the longitudinal prior, \code{alpha_vec} can be set to reflect temporal ordering.
#' @param eta Numeric; GcompBART-specific prior parameter for variance shrinkage. Default = 1.
#' @param phi Numeric or vector; GcompBART-specific prior parameter for time-varying effect scaling. Default = 1.
#' @inheritParams SoftBart::Hypers
#'
#' @return A list containing all fields from \code{SoftBart::Hypers()} plus:
#' \itemize{
#'   \item \code{alpha_vec} – longitudinal alpha prior
#'   \item \code{eta}, \code{phi} – additional prior parameters
#'   \item \code{tgroup} – integer vector of predictor group assignments
#'   \item \code{tgroup_size} – sizes of each time group
#' }
#'
#' @examples
#' # Example: Construct hypers for longitudinal data
#' X_train <- matrix(rnorm(100), ncol = 5)
#' Y_train <- rnorm(20)
#' tmp_tg <- rep(1:2, each = 10)
#'
#' hypers <- Hypers(X = X_train, Y = Y_train, tgroup = tmp_tg,
#'                  alpha_vec = c(0.5, 0.75, 0.9, 0.95, 1),
#'                  eta = 1, phi = 1)
#'
#' @export
Hypers <- function(
    X, Y, group = NULL, tgroup = NULL,
    alpha = 1, eta = 1, phi = 1, alpha_vec = 1,
    beta = 2, gamma = 0.95, k = 2,
    sigma_hat = NULL, shape = 1, width = 0.1, num_tree = 20,
    alpha_scale = NULL, alpha_shape_1 = 0.5,
    alpha_shape_2 = 1, tau_rate = 10, num_tree_prob = NULL,
    temperature = 1.0, weights = NULL, normalize_Y = TRUE
) {
  # 1) Build the base structure using SoftBart (all defaults are preserved)
  sb <- SoftBart::Hypers(
    X = X, Y = Y, group = group,
    alpha = alpha, beta = beta, gamma = gamma, k = k,
    sigma_hat = sigma_hat, shape = shape, width = width, num_tree = num_tree,
    alpha_scale = alpha_scale, alpha_shape_1 = alpha_shape_1,
    alpha_shape_2 = alpha_shape_2, tau_rate = tau_rate,
    num_tree_prob = num_tree_prob, temperature = temperature,
    weights = weights, normalize_Y = normalize_Y
  )

  # 2) Inject GcompBART-specific fields
  sb$alpha_vec <- alpha_vec
  sb$eta       <- eta
  sb$phi       <- phi

  # 3) tgroup handling
  if (is.null(tgroup)) {
    sb$tgroup <- rep(1L, ncol(X)) - 1L
  } else {
    tg <- as.integer(tgroup)
    sb$tgroup <- tg + ifelse(tg == 0L, 1L, 0L) * max(tg)
  }
  sb$tgroup_size <- as.vector(table(sb$tgroup))

  return(sb)
}



#' @title Extended Opts for GcompBART
#' @description Builds the base options via SoftBart::Opts() and then extends with GcompBART-specific update flags for time-varying pieces and alpha_vec.
#'
#' @param num_burn Integer; number of burn-in iterations (SoftBart semantics).
#' @param num_thin Integer; thinning interval (SoftBart semantics).
#' @param num_save Integer; number of saved iterations (SoftBart semantics).
#' @param num_print Integer; frequency of printing progress (SoftBart semantics).
#' @param update_sigma_mu Logical; update sigma_mu flag (SoftBart semantics).
#' @param update_s Logical; update s flag (SoftBart semantics).
#' @param update_alpha Logical; update alpha flag (SoftBart semantics).
#' @param update_beta Logical; update beta flag (SoftBart semantics).
#' @param update_gamma Logical; update gamma flag (SoftBart semantics).
#' @param update_tau Logical; update tau flag (SoftBart semantics).
#' @param update_tau_mean Logical; update tau_mean flag (SoftBart semantics).
#' @param update_sigma Logical; update sigma flag (SoftBart semantics).
#' @param cache_trees Logical; cache trees for efficiency (SoftBart semantics).
#' @param update_tvp Logical; GcompBART-specific extension for time-varying pieces.
#' @param update_eta Logical; GcompBART-specific extension for eta updates.
#' @param update_phi Logical; GcompBART-specific extension for phi updates.
#' @param update_alpha_vec Logical; GcompBART-specific extension for alpha_vec updates.
#'
#' @return A named list identical to SoftBart::Opts() plus GcompBART-specific flags.
#' @examples
#' opts <- Opts(num_burn = 1000, update_tvp = TRUE)
#' @export
Opts <- function(
    num_burn = 250, num_thin = 1, num_save = 250, num_print = 100,
    update_sigma_mu = TRUE, update_s = TRUE, update_alpha = TRUE,
    update_beta = FALSE, update_gamma = FALSE, update_tau = TRUE,
    update_tau_mean = FALSE, update_sigma = TRUE,
    cache_trees = TRUE,
    # --- GcompBART extensions (default FALSE to preserve SoftBart behavior) ---
    update_tvp = FALSE, update_eta = FALSE, update_phi = FALSE,
    update_alpha_vec = FALSE
) {
  # 1) Build the base list with SoftBart (covers all shared fields + defaults)
  sb <- SoftBart::Opts(
    num_burn = num_burn, num_thin = num_thin, num_save = num_save,
    num_print = num_print,
    update_sigma_mu = update_sigma_mu, update_s = update_s,
    update_alpha = update_alpha,
    update_beta = update_beta, update_gamma = update_gamma,
    update_tau = update_tau, update_tau_mean = update_tau_mean,
    update_sigma = update_sigma,
    cache_trees = cache_trees
  )

  # 2) Inject GcompBART-only flags
  sb$update_tvp        <- isTRUE(update_tvp)
  sb$update_eta        <- isTRUE(update_eta)
  sb$update_phi        <- isTRUE(update_phi)
  sb$update_alpha_vec  <- isTRUE(update_alpha_vec)

  # 3) Keep explicit false for update_num_tree to mirror your current shape
  if (is.null(sb$update_num_tree)) sb$update_num_tree <- FALSE

  return(sb)
}


#' @title MakeForest
#' @description Creates a Forest object from Rcpp module.
#' @param hypers List of hyperparameters.
#' @param opts List of options.
#' @return An object of class Forest.
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10)
#' Y <- rnorm(10)
#' hypers <- Hypers(X = X, Y = Y)
#' opts <- Opts()
#' forest <- MakeForest(hypers, opts)
#' }
#' @export
MakeForest <- function(hypers, opts) {
  mf <- Rcpp::Module(module = "mod_forest", PACKAGE = "GcompBART")
  return(new(mf$Forest, hypers, opts))
}



#' @title #' Create a Longitudinal Covariance Matrix
#' @description Constructs a block-structured covariance matrix for longitudinal data,
#' where correlations between time points decay according to specified lag values.
#'
#' @param n_time Integer. Number of time points.
#' @param lags Numeric vector. Correlation values for each lag (e.g., \code{c(a1, a2, a3)}).
#'
#' @return A \code{(n_time * 5) x (n_time * 5)} covariance matrix.
#' @examples
#' # Example: 4 time points, 3 lags
#' make_longitudinal_cov(n_time = 4, lags = c(0.4, 0.2, 0.1))
#'
#' @export
# Helper function for covariance
make_longitudinal_cov <- function(n_time, lags) {
  p <- n_time * 5 # 5 variables per time point.
  Sig <- matrix(0, p, p)
  diag(Sig) <- 1

  for (lag in seq_along(lags)) {
    a <- lags[lag]
    for (t in 1:(p - lag * 5)) {
      Sig[t, t + lag * 5] <- a
      Sig[t + lag * 5, t] <- a
    }
  }
  Sig
}


#' @title Simulate Longitudinal Data Using Modified Friedman Function
#'
#' @description
#' This function extends the classic Friedman function to multiple blocks of predictors,
#' each contributing nonlinearly to the outcome. Predictors are correlated according
#' to a user-specified covariance matrix and transformed to uniform via the normal CDF.
#'
#' @details
#' Data are generated according to a modified version of Friedman's five-dimensional test function
#' extended to a longitudinal setting with multiple time points.
#'
#' For each subject \eqn{i}, the true mean function is:
#' \deqn{
#' f(\mathbf{x}_i) = \sum_{t=1}^T c_t \left[ 10 \sin(\pi x_{1,it} x_{2,it}) +
#' 20 (x_{3,it} - 0.5)^2 + 10 x_{4,it} + 5 x_{5,it} \right]
#' }
#'
#' where \eqn{T} is the number of time points, and \eqn{c_t} determines the time-specific
#' association between predictors and the outcome. The outcome at the final time point is:
#' \deqn{
#' Y_{iT} \sim \mathrm{Normal}\left( f(\mathbf{X}_i), \sigma_T^2 \right)
#' }
#'
#' with default values \eqn{T = 4}, \eqn{\sigma_T = 10}, and decreasing weights \eqn{c_t}
#' to reflect stronger associations closer to the outcome time.
#'
#' Correlated predictors are generated using a multivariate normal distribution with covariance
#' matrix \code{Sigma}, then transformed to uniform via \code{pnorm()}.
#'
#' @param N Integer. Number of observations to simulate.
#' @param n_time Integer. Number of time points (blocks).
#' @param lags Numeric vector. Correlation values for each lag between time points.
#' @param sigma Numeric. Standard deviation of the noise term added to the outcome.
#' @param lambda Numeric vector of weights for different predictor blocks
#'   (length should match \code{n_time}).
#'
#' @details
#' The predictors are simulated from a multivariate normal distribution with
#' a block-structured covariance matrix. The outcome \code{Y} is generated as:
#' \deqn{
#' Y = \sum_{b=1}^{n\_time} \lambda_b \left[ 10 \sin(\pi X_{b1} X_{b2}) +
#' 20 (X_{b3} - 0.5)^2 + 10 X_{b4} + 5 X_{b5} \right] + \epsilon
#' }
#' where \eqn{\epsilon \sim N(0, \sigma^2)}.
#'
#' @return A data frame with simulated predictors (\code{X}), outcome (\code{Y}),
#' and true mean (\code{mu}).
#'
#' @examples
#' set.seed(123)
#' sim <- sim_lfried(N = 100, n_time = 4, lags = c(0.4, 0.2, 0.1),
#'                   sigma = 1, lambda = c(0.25, 0.5, 0.75, 1))
#'
#' @export
#'
sim_lfried <- function(N, n_time = 4, lags = c(0.4, 0.2, 0.1),
                       sigma = 1, lambda =  c(0.25,0.5,0.75,1)) {

  P <- n_time * 5
  Sigma <- make_longitudinal_cov(n_time, lags)

  # Simulate predictors
  rawvars <- MASS::mvrnorm(N, rep(0, P), Sigma)
  X <- pnorm(rawvars)

  # Compute mu based on blocks
  mu <- numeric(N)
  for (block in seq_len(n_time)) {
    start <- (block - 1) * 5 + 1
    end <- block * 5
    mu <- mu + lambda[block] * (
      10 * sin(pi * X[, start] * X[, start + 1]) +
        20 * (X[, start + 2] - 0.5)^2 +
        10 * X[, start + 3] +
        5 * X[, start + 4]
    )
  }

  # Outcome
  Y <- mu + sigma * rnorm(N)

  data.frame(X, Y = Y, mu = mu)
}


#' @title rmvnorm
#' @description Internal helper wrapping C++ for generating random multivariate normal samples.
#' @param mean Numeric vector of means.
#' @param Precision Precision matrix (inverse covariance).
#' @return A numeric vector or matrix of samples.
#' @keywords internal
#' @name rmvnorm
NULL

#' @title choll
#' @description Internal helper wrapping C++ for Cholesky decomposition.
#' @param Sigma Covariance matrix.
#' @return Cholesky factor.
#' @keywords internal
#' @name choll
NULL

#' @title update_sigma
#' @description Internal helper wrapping C++ for updating sigma in MCMC.
#' @param r Numeric value.
#' @param sigma_hat Estimated sigma.
#' @param sigma_old Previous sigma value.
#' @param temperature Temperature parameter for MH step.
#' @return Updated sigma value.
#' @keywords internal
#' @name update_sigma
NULL

#' @title rlgam
#' @description Internal helper wrapping C++ for sampling from a gamma distribution.
#' @param shape Shape parameter.
#' @return A random gamma value.
#' @keywords internal
#' @name rlgam
NULL

#' @title do_mh
#' @description Internal helper wrapping C++ for Metropolis-Hastings acceptance step.
#' @param loglik_new New log-likelihood.
#' @param loglik_old Old log-likelihood.
#' @param new_to_old Proposal ratio new→old.
#' @param old_to_new Proposal ratio old→new.
#' @return Logical indicating acceptance.
#' @keywords internal
#' @name do_mh
NULL


utils::globalVariables(c("opts2", "continuous"))


