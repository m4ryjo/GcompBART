#' @title BMfits
#' @description Fits Bayesian models (e.g., SoftBart) for each variable in the data according to its type and temporal structure. The fitted models are used in subsequent g-computation via `gcompbart()`.
#'
#' @param data A matrix or data frame containing variables in temporal order: confounders, exposures, outcomes, survival, and dropout indicators.
#' @param var.type A character vector specifying the type of each variable in `data`. Options: `"X"` (confounder), `"Fi"` (fixed exposure), `"Ra"` (random exposure), `"Y"` (outcome), `"S"` (survival), `"D"` (dropout).
#' @param fixed.regime Optional specification for fixed treatment regime.
#' @param random.regime Optional specification for random treatment regime.
#' @param param Optional parameter vector for regime specification.
#' @param above Logical; indicates if regime applies above a threshold.
#' @param threshold Numeric threshold for regime application.
#' @param nat_value Logical; if TRUE, uses natural value regime.
#' @param incremental Logical; if TRUE, applies incremental regime.
#' @param drop_param Optional parameter for dropout modeling.
#' @param J Integer; number of MCMC iterations (default = 2000).
#' @param opts List of options passed to SoftBart.
#' @param Suppress Logical; if TRUE, suppresses console output (default = TRUE).
#' @param By Integer; if `Suppress = FALSE`, output is printed every `By` iterations.
#' @param tgroup Integer vector indicating time grouping for each variable. Use `0` for static variables and positive integers for time-varying ones.
#' @param Lag1 Logical; if TRUE, applies lag structure to time-varying variables.
#' @param base_hypers A list of base hyperparameters.
#' @param ... Additional arguments passed to internal modeling functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{BModels}{List of fitted Bayesian models for each variable.}
#'   \item{n_Y}{Number of outcome variables.}
#'   \item{n_S}{Number of survival variables.}
#'   \item{n_Reg}{Number of regimes (if applicable).}
#'   \item{num_save}{Number of saved iterations.}
#'   \item{continuous}{Continuous variables.}
#'   \item{iterations}{Number of iterations.}
#'   \item{...}{Other metadata needed for `gcompbart()`.}
#' }
#'
#' @note Future versions will include an optional argument `return_full = TRUE`
#' to return additional objects such as model options (`opts`), hyperparameters
#' (`base_hypers`), and diagnostics for model fit.
#'
#'
#' @examples
#' \dontrun{
#' n_burn <- 10
#' n_thin <- 1
#' n_save <- 10
#' n_J <- 100
#'
#' set.seed(1234)
#'
#' # Simulate dataset
#' n <- 100
#' times_sim <- rep(1:4, each = 5)
#' lambda <- c(0.25, 0.5, 0.75, 1)
#' a1 <- 0.4; a2 <- 0.2; a3 <- 0.1; sig <- 10
#'
#' Sig <- matrix(c(...), 20, 20, byrow = TRUE) # truncated for brevity
#'
#' vartype_bl <- c(rep("X0", 5), rep("X", 15), "Y")
#' tgroup <- c(rep(1:4, each = 5), 4)
#'
#' training_data <- sim_lfried(n, 20, Sig, sig, lambda)
#'
#' # Fit the model
#' BM <- BMfits(training_data[, 1:21],
#'              var.type = vartype_bl,
#'              opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
#'                          update_s = TRUE, update_alpha = TRUE, update_tau = TRUE,
#'                          update_sigma_mu = TRUE),
#'              tgroup = tgroup)
#'
#' # Perform g-computation
#' out_gcomp <- gcompbart(training_data[, 1:21],
#'                        BModels = BM,
#'                        var.type = vartype_bl,
#'                        J = n_J,
#'                        opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
#'                                    update_s = TRUE, update_alpha = TRUE, update_tau = TRUE,
#'                                    update_sigma_mu = TRUE),
#'                        tgroup = tgroup)
#'
#' # Results
#' out_gcomp$summary_out
#' }
#' @export
BMfits <- function(data,
                   var.type,
                   fixed.regime = NULL,
                   random.regime = NULL,
                   param = NULL,
                   above = FALSE,
                   threshold = NULL,
                   nat_value = FALSE,
                   incremental = FALSE,
                   drop_param = NULL,
                   J = 2000,
                   opts = NULL,
                   Suppress = TRUE,
                   By = 100,
                   tgroup = rep(1,  ncol(data) - 1),
                   Lag1 = FALSE,
                   base_hypers = NULL,
                   ...

) {
    data <- as.matrix(data)

    ## test if the variables in data are continuous or binary etc
    continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1))

    n_Y  <- sum(var.type == "Y")   # Outcome variables
    n_S  <- sum(var.type == "S")   # Survival indicators
    n_Fi <- sum(var.type == "Fi")  # Fixed exposures
    n_Ra <- sum(var.type == "Ra")  # Random exposures

    n_Reg <- if (!is.null(fixed.regime)) length(fixed.regime) else 1

    if (ncol(data) != length(var.type)) stop("Mismatch: data has ", ncol(data), " columns but var.type has length ", length(var.type))

    if (any(is.na(data[, which(var.type == "X0")]))) stop("The baseline covariates includes NAs ")

    if(!is.null(fixed.regime)){
        if (all(sapply(fixed.regime, length) != n_Fi)) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
    }

    if (is.null(opts)) {
      opts <- Opts(
        update_sigma_mu = TRUE, update_s = FALSE,
        update_alpha = FALSE,   # disable longitudinal prior
        update_beta = FALSE, update_gamma = FALSE,
        update_tau = TRUE, update_tau_mean = FALSE,
        update_sigma = TRUE, cache_trees = TRUE
      )
    }

    # opts2 for single time-group (longitudinal prior OFF)
    opts2 <- Opts(
      update_sigma_mu = TRUE, update_s = FALSE,
      update_alpha = FALSE,   # disable longitudinal prior
      update_beta = FALSE, update_gamma = FALSE,
      update_tau = TRUE, update_tau_mean = FALSE,
      update_sigma = TRUE, cache_trees = TRUE
    )

    # If base_hypers is NULL, create defaults
    if (is.null(base_hypers)) {
      base_hypers <- BaseHypers()  # Use your constructor with defaults
    }

    BModels <- vector("list", ncol(data) - 1)

    for (i in 1:(ncol(data) - 1)) {
      id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
      var_i <- var.type[1 + i]

      if (var_i %in% c("Fi", "X0") || (var_i == "Ra" && nat_value == FALSE)) {
        BModels[[i]] <- NA
        next
      }

      X_train <- get_lagged_X(data, var.type, i, id, Lag1)
      if (i == 1) X_train <- as.matrix(cbind(1, X_train))  # Add intercept for first timepoint

      Y_train <- data[id == 1, 1 + i]
      tmp_tg <- tgroup[which(var.type[1:i] != "S" & var.type[1:i] != "D")]
      opts_b <- if (length(unique(tmp_tg)) == 1) opts2 else opts

      # Merge with model-specific arguments
      hypers_args <- c(base_hypers, list(X = X_train, Y = Y_train, tgroup = tmp_tg))

      # Build hypers object
      hypers <- do.call(Hypers, hypers_args)

      #hypers <- make_hypers(X_train, Y_train, tmp_tg, num_tree, alpha, eta, phi, alpha_vec, alpha_shape_1)
      BModels[[i]] <- fit_model(X_train, Y_train, is_binary = !continuous[1 + i], opts_b, hypers, Suppress)

      if (!Suppress) print(i)
    }

    iterations <- seq(from = opts$num_burn + opts$num_thin,
                      to = opts$num_burn + opts$num_thin * opts$n_save,
                      by = opts$num_thin)
    n_save = opts$n_save

    return(list(
      BModels = BModels,
      n_Y = n_Y,
      n_S = n_S,
      n_Reg = n_Reg,
      n_save = n_save,
      continuous = continuous,
      iterations = iterations,
      ...
    ))

  }

#' Performs Bayesian g-computation using SoftBart models to simulate counterfactual outcomes under specified regimes.
#'
#' @param data A matrix or data frame containing confounders, exposures, outcomes, and mortality indicators in temporal order (left to right).
#' @param var.type A character vector specifying the type of each variable in `data`. Options: `"X"` (confounder), `"Fi"` (fixed exposure), `"Ra"` (random exposure), `"Y"` (outcome), `"S"` (survival), `"D"` (dropout).
#' @param fixed.regime Optional list specifying fixed regimes for exposures or confounders.
#' @param random.regime Optional character vector specifying random regime types. Supported: `"uniform"`, `"normal"`, `"triangular"`.
#' @param param List of parameter vectors for each random regime. Length depends on regime type: 1 (binomial), 2 (uniform/normal), 3 (triangular).
#' @param above Logical. If `TRUE`, applies thresholding above a specified value.
#' @param threshold Numeric vector of threshold values for regime application.
#' @param nat_value Logical. If `TRUE`, uses natural value regime.
#' @param incremental Logical. If `TRUE`, applies incremental regime logic.
#' @param drop_param List of parameter vectors for dropout shift modeling using triangular distribution.
#' @param J Integer. Size of pseudo data to generate. Default is 2000.
#' @param opts List of options. Default is NULL.
#' @param Suppress Logical. If `TRUE`, suppresses console output. Default is `TRUE`.
#' @param By Integer. If `Suppress = FALSE`, output is printed every `By` iterations.
#' @param tgroup Integer vector indicating time grouping for each variable. Use `0` for static variables and positive integers for time-varying ones.
#' @param Lag1 Logical. If `TRUE`, applies lag structure to time-varying variables.
#' @param base_hypers  List of hypers. Default is NULL.
#' @param BModels Optional list of pre-fitted SoftBart models.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A named list containing:
#' \describe{
#'   \item{summary_out}{Posterior summaries of predicted outcomes.}
#'   \item{y_hat}{Posterior samples of predicted outcomes.}
#'   \item{summary_surv}{Posterior summaries of survival probabilities.}
#'   \item{s_hat}{Posterior samples of survival probabilities.}
#' }
#'
#' @note Future versions will include an optional argument `return_full = TRUE`
#' to return additional objects such as model options (`opts`), hyperparameters
#' (`base_hypers`), and diagnostics for model fit.
#'
#' @examples
#' \dontrun{
#' n_burn <- 10
#' n_thin <- 1
#' n_save <- 10
#' n_J <- 100
#'
#' set.seed(1234)
#'
#' # Simulate dataset
#' n <- 100
#' times_sim <- rep(1:4, each = 5)
#' lambda <- c(0.25, 0.5, 0.75, 1)
#' a1 <- 0.4; a2 <- 0.2; a3 <- 0.1; sig <- 10
#'
#' Sig <- matrix(c(...), 20, 20, byrow = TRUE) # truncated for brevity
#'
#' vartype_bl <- c(rep("X0", 5), rep("X", 15), "Y")
#' tgroup <- c(rep(1:4, each = 5), 4)
#'
#' training_data <- sim_lfried(n, 20, Sig, sig, lambda)
#'
#' # Fit the model
#' BM <- BMfits(training_data[, 1:21],
#'              var.type = vartype_bl,
#'              opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
#'                          update_s = TRUE, update_alpha = TRUE, update_tau = TRUE,
#'                          update_sigma_mu = TRUE),
#'              tgroup = tgroup)
#'
#' # Perform g-computation
#' out_gcomp <- gcompbart(training_data[, 1:21],
#'                        BModels = BM,
#'                        var.type = vartype_bl,
#'                        J = n_J,
#'                        opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
#'                                    update_s = TRUE, update_alpha = TRUE, update_tau = TRUE,
#'                                    update_sigma_mu = TRUE),
#'                        tgroup = tgroup)
#'
#' # Results
#' out_gcomp$summary_out
#' }
#' @export
gcompbart <- function(
    data,                  # Matrix or data frame with confounders, exposures, outcomes, and mortality indicators in temporal order.
    var.type,              # Vector specifying variable types: X0=baseline confounder, X=confounder, Fi=fixed exposure, Ra=random exposure, Y=outcome, S=survival, D=dropout.
    fixed.regime = NULL,   # List specifying fixed regimes for exposures or confounders.
    random.regime = NULL,  # Character string: "uniform", "normal", or "triangular". Default is NULL.
    param = NULL,          # Parameters for random regime: binomial(a), uniform(a, b), normal(a, b), or triangular(a, b, c).
    above = FALSE,
    threshold = NULL,
    nat_value = FALSE,
    incremental = FALSE,
    drop_param = NULL,
    J = 2000,              # Size of pseudo data. Default is 2000.
    opts = NULL,           # Options passed to SoftBart.
    Suppress = TRUE,       # If FALSE, output is printed every `By` iterations.
    By = 100,              # Iteration interval for printing output (used only if Suppress = FALSE).
    tgroup = rep(1, ncol(data) - 1),  # Time grouping: 0 for static, >0 for time-varying variables.
    Lag1 = FALSE,
    base_hypers = NULL,
    BModels = NULL,
    ...
) {
    ##################################################################
    ## Fit the observed data models using BMfits and save as a list

  # Ensure BModels always holds the full BMfits output
  if (is.null(BModels)) {
    bm_out <- BMfits(...)   # Call BMfits and get full list
  } else {
    bm_out <- BModels       # Rename for clarity
  }

  # Extract components consistently
  BModels     <- bm_out$BModels   # Keep name consistent with rest of code
  n_Y         <- bm_out$n_Y
  n_S         <- bm_out$n_S
  n_Reg       <- bm_out$n_Reg
  continuous  <- bm_out$continuous
  n_save      <- bm_out$n_save
  iterations  <- bm_out$iterations

  ########################################################

  ## Create matrices to store the output. One for the outcome and one for survival probabilities.
  ## Note, the code allows for multiple outcomes/survival, e.g. for different time points.
  nrow_y <- (if (!is.null(fixed.regime) || !is.null(random.regime)) safe_val(n_Reg) else 1) * safe_val(n_Y)
  nrow_s <- (if (!is.null(fixed.regime) || !is.null(random.regime)) safe_val(n_Reg) else 1) * safe_val(n_S)

  y <- matrix(nrow = nrow_y, ncol = n_save)
  s <- matrix(nrow = nrow_s, ncol = n_save)

  b_size <- length(which(tgroup == 0)) + 1

  ## MC-integration over "n_save" iterations
    for(it in 1:n_save) {
        ks <- 1
        ky <- 1
        l <- 1
        n <- 1
        m <- 1
        o <- 1
        d <- 0
        drop <- NULL
        s_hat <- NULL


        # First variable
        if (var.type[1] == "Ra") {
          x <- list(generate_random_matrix(J, random.regime[1], param[[1]]))
          m <- 2

        } else if (var.type[1] == "X0") {
          result <- prepare_X0_data(data, var.type, J)
          x <- result$x
          m <- result$m

        } else if (var.type[1] == "Fi" && var.type[2] == "X0") {
          result <- prepare_Fi_X0_data(data, var.type, fixed.regime, J)
          x <- result$x
          m <- result$m

        } else {
          stop("Unsupported var.type combination at position 1: ", var.type[1])
        }

        for (j in m:ncol(data)) {
          if (var.type[j] == "X") {
            x <- handle_X(x, j, continuous, BModels, iterations, it, Lag1, tgroup, b_size)

          } else if (var.type[j] == "Fi") {
            x <- handle_Fi(x, fixed.regime, l)
            l <- l + 1

          } else if (var.type[j] == "Ra") {
            x <- handle_Ra(x, j, continuous, BModels, iterations, it, random.regime, o, param, above, threshold, nat_value, incremental, Lag1, tgroup, b_size)
            o <- o + 1

          } else if (var.type[j] == "S") {
            surv_hat <- lapply(x, function(x2) {
              pred_surv(BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
            })
            s_tmp <- surv_hat
            s[n:(n + n_Reg - 1), it] <- sapply(s_tmp, mean)
            n <- n + n_Reg

            surv <- lapply(surv_hat, function(a) rbinom(length(a), 1, prob = a))
            x <- Map(function(a, b) a[b == 1, ], x, surv)
            s_hat <- Map(function(a, b) a[b == 1], s_tmp, surv)

            ks <- ks + n_Reg
            if (!is.null(drop)) {
              drop <- Map(function(a, b) a[b == 1], drop, surv)
            }

          } else if (var.type[j] == "D") {
            drop_tmp <- lapply(x, function(x2) {
              pred_drop(BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
            })
            drop <- if (d == 0) drop_tmp else Map(function(e, f) ifelse(e == 0 & f == 0, 0, 1), drop_tmp, drop)
            d <- d + 1
            drop_shift <- if (is.null(drop_param)) {
              lapply(drop, function(f) f * 0)
            } else {
              lapply(drop, function(f) f * EnvStats::rtri(length(f), drop_param[[d]][1], drop_param[[d]][2], drop_param[[d]][3]))
            }

          } else if (var.type[j] == "Y") {
            tmp <- lapply(x, function(x2) {
              pred_y(continuous[j], BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
            })
            if (d > 0) {
              tmp <- Map(`+`, tmp, drop_shift)
            }
            y[ky:(ky + n_Reg - 1), it] <- sapply(tmp, mean)
            ky <- ky + n_Reg

            if (!continuous[j]) {
              tmp <- lapply(tmp, function(tmp2) rbinom(length(tmp2), 1, prob = tmp2))
            } else {
              tmp <- lapply(tmp, function(tmp2) rnorm(length(tmp2), mean = tmp2, sd = mean(BModels[[j - 1]]$sigma)))
            }

            x <- Map(cbind, x, tmp)

            if (Lag1 == TRUE && tgroup[j] > 1) {
              x <- lapply(x, function(d) subset(d, select = -b_size))
            }
          }
        }

        if(it %in% seq(0, n_save, by = By)) {
            print(paste("done ", it, " (out of ", n_save,")", sep=""))
        }
    }


  ## Save posterior samples of the mean predicted outcomes (y_hat) and summaries
  if (n_Y > 0) {
    means_out <- t(apply(y, 1, function(x) {
      c(mean(x), quantile(x, c(0.025, 0.975), na.rm = TRUE))
    }))
    if (is.null(fixed.regime)) {
      rownames(means_out) <- paste0("Y ", 1:n_Y)
      rownames(y) <- paste0("Y ", 1:n_Y)
    } else {
      rownames(means_out) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y ", rep(1:n_Y, each = n_Reg))
      rownames(y) <- paste0("Regime ", rep(1:n_Reg, n_Y), " : Y ", rep(1:n_Y, each = n_Reg))
    }
    colnames(means_out) <- c("mu", "2.5%", "97.5%")
  } else {
    y <- NA
    means_out <- NA
  }

  ## Save posterior samples of the mean survival probabilities and summaries
  if (n_S > 0) {
    means_surv <- t(apply(s, 1, function(x) {
      c(mean(x), quantile(x, c(0.025, 0.975), na.rm = TRUE))
    }))
    if (is.null(fixed.regime)) {
      rownames(means_surv) <- paste0("Survival ", 1:n_S)
      rownames(s) <- paste0("Survival ", 1:n_S)
    } else {
      rownames(means_surv) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival ", rep(1:n_S, each = n_Reg))
      rownames(s) <- paste0("Regime ", rep(1:n_Reg, n_S), " : Survival ", rep(1:n_S, each = n_Reg))
    }
    colnames(means_surv) <- c("mu", "2.5%", "97.5%")
  } else {
    s <- NA
    means_surv <- NA
  }


  out <- list(means_out, y, means_surv, s)
  names(out) <- c("summary_out", "y_hat", "summary_surv", "s_hat")

  return(out)
}

