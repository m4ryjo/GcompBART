#' Helper functions used in BMfits and gcompbart
#' @keywords internal

#' @title Sample from a Dirichlet Distribution
#' @description Generates random samples from a Dirichlet distribution using the Gamma-based transformation. Each sample is obtained by drawing independent Gamma random variables and normalizing them.

#' @param n Integer. Number of samples to draw.
#' @param alpha Numeric vector. Concentration parameters for the Dirichlet distribution.
#'
#' @details
#' The Dirichlet distribution is commonly used for generating random probability vectors.
#' This implementation uses the property that if \eqn{X_i \sim \text{Gamma}(\alpha_i, 1)},
#' then \eqn{X_i / \sum_j X_j \sim \text{Dirichlet}(\alpha)}.
#'
#' @return A numeric matrix of dimension \code{n x length(alpha)}, where each row is a Dirichlet sample.
#'
#' @examples
#' \dontrun{
#' rdirichlet(3, alpha = c(1, 1, 1, 1))
#' }
#'
#' @keywords internal
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  samples <- matrix(0, nrow = n, ncol = k)
  for (i in 1:n) {
    gammas <- rgamma(k, shape = alpha, rate = 1)
    samples[i, ] <- gammas / sum(gammas)
  }
  samples
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
BayesBoot <- function(df, J) {
  # Draw Dirichlet weights
  BBweights <- rdirichlet(1, rep(1, nrow(df)))

  # Multinomial resampling based on Dirichlet weights
  BBdraws <- rmultinom(1, J, BBweights)

  # Repeat rows according to bootstrap draws
  df0 <- as.matrix(data.frame(lapply(df, rep, BBdraws)))
  return(df0)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
expit <- function(x) {
  exp(x)/(1+exp(x))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
safe_val <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(0)
  return(x)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
get_lagged_X <- function(data, var.type, i, id, Lag1) {
  if (Lag1) {
    cols <- which((1:i) == i & var.type[1:i] != "S" & var.type[1:i] != "D")
  } else {
    cols <- which(var.type[1:i] != "S" & var.type[1:i] != "D")
  }
  data[id == 1, cols, drop = FALSE]
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
generate_random_matrix <- function(J, regime, param) {
  switch(regime,
         "uniform" = matrix(runif(J, param[1], param[2])),
         "normal" = matrix(rnorm(J, param[1], param[2])),
         "triangular" = matrix(rtri(J, param[1], param[2], param[3])),
         "binomial" = matrix(rbinom(J, 1, param[1])),
         stop("Unsupported regime type")
  )
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
prepare_X0_data <- function(data, var.type, J) {
  dataC0 <- data.frame(data[, which(var.type == "X0")])
  list_data <- list(BayesBoot(dataC0, J))
  m_val <- ncol(dataC0) + 1
  return(list(x = list_data, m = m_val))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
prepare_Fi_X0_data <- function(data, var.type, fixed.regime, J) {
  data0 <- Map(data.frame, lapply(fixed.regime, function(r) {
    data[which(data[, 1] == r[1]), c(1, which(var.type == "X0"))]
  }))
  list_data <- lapply(data0, function(x) BayesBoot(x, J))
  m_val <- ncol(data0[[1]]) + 1
  return(list(x = list_data, m = m_val))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_X <- function(x, j, continuous, BModels, iterations, it, Lag1, tgroup, b_size) {
  lapply(x, function(x2) {
    pred_comb(continuous[j], BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
  })
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_Fi <- function(x, fixed.regime, l) {
  Map(cbind, x, lapply(fixed.regime, function(r) r[l]))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_Ra <- function(x, j, continuous, BModels, iterations, it, random.regime, o, param, above, threshold, nat_value, incremental, Lag1, tgroup, b_size) {
  lapply(x, function(x2) {
    pred_ra(continuous[j], BModels[[j - 1]], x2, iterations[it], random.regime[o], param[[o]],
            above, threshold[o], nat_value, incremental, Lag1, tgroup[j] - 1, b_size)
  })
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_S <- function(x, j, BModels, iterations, it, Lag1, tgroup, b_size, s_hat, n_Reg, s, ks, drop) {
  surv_hat <- lapply(x, function(x2) {
    pred_surv(BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
  })
  s_tmp <- surv_hat
  s[ks:(ks + n_Reg - 1), it] <- sapply(s_tmp, mean)

  surv <- lapply(surv_hat, function(a) rbinom(length(a), 1, prob = a))
  x <- Map(function(a, b) a[b == 1, ], x, surv)
  s_hat <- Map(function(a, b) a[b == 1], s_tmp, surv)

  ks <- ks + n_Reg
  if (!is.null(drop)) {
    drop <- Map(function(a, b) a[b == 1], drop, surv)
  }

  return(list(x = x, s_hat = s_hat, drop = drop, ks = ks))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_D <- function(x, j, BModels, iterations, it, Lag1, tgroup, b_size, drop, d, drop_param) {
  drop_tmp <- lapply(x, function(x2) {
    pred_drop(BModels[[j - 1]], x2, iterations[it], Lag1, tgroup[j] - 1, b_size)
  })
  drop <- if (d == 0) drop_tmp else Map(function(e, f) ifelse(e == 0 & f == 0, 0, 1), drop_tmp, drop)
  d <- d + 1
  drop_shift <- if (is.null(drop_param)) {
    lapply(drop, function(f) f * 0)
  } else {
    lapply(drop, function(f) f * rtri(length(f), drop_param[[d]][1], drop_param[[d]][2], drop_param[[d]][3]))
  }
  return(list(drop = drop, drop_shift = drop_shift, d = d))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
handle_Y <- function(x, j, continuous, BModels, iterations, it, Lag1, tgroup, b_size, drop_shift, d, y, ky, n_Reg, Lag1_flag) {
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

  if (Lag1_flag && tgroup[j] > 1) {
    x <- lapply(x, function(d) subset(d, select = -b_size))
  }

  return(list(x = x, y = y, ky = ky))
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
fit_model <- function(X, Y, is_binary, opts, hypers, suppress) {
  if (suppress) {
    quiet(if (is_binary) gc_softbart_probit(X, Y, opts = opts, hypers = hypers)
          else gc_softbart_regression(X, Y, opts = opts, hypers = hypers))
  } else {
    if (is_binary) gc_softbart_probit(X, Y, opts = opts, hypers = hypers)
    else gc_softbart_regression(X, Y, opts = opts, hypers = hypers)
  }
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
pred_comb <- function(cont, BM, MCdatatmp, IT, L1, TG, bsize){
    if(ncol(MCdatatmp)==1) {
        MCdata <- as.matrix(cbind(1, MCdatatmp))
    } else {
        MCdata <- MCdatatmp
    }

    X <- MCdata
    if(L1==TRUE & TG>0){
        while(length(BM$ecdfs)<ncol(X)){
            X <- X[,-bsize]
        }
    }

    for(i in 1:ncol(X)) {
        X[,i] <- BM$ecdfs[[i]](X[,i])
    }

    if(cont == FALSE) {
        x_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
        new_MCdata <- cbind(MCdatatmp, rbinom(length(x_hat), 1, prob = x_hat))
    } else if(cont == TRUE) {
        x_hat <- as.numeric(BM$forest$predict_iteration(X, IT)) * BM$sd_Y + BM$mu_Y
        new_MCdata <- cbind(MCdatatmp, rnorm(length(x_hat), mean = x_hat, sd = BM$sigma))
    }

    if(L1==TRUE & TG>1){
        new_MCdata <- new_MCdata[,-bsize]
    } else {
        new_MCdata <- new_MCdata
    }

    return(new_MCdata)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
pred_drop <- function(BM, MCdata, IT, L1, TG, bsize){
  X <- MCdata

  if(L1==TRUE & TG>0){
      while(length(BM$ecdfs)<ncol(X)){
          X <- X[,-bsize]
      }
  }

  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  x_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
  dropout_ind <- rbinom(length(x_hat), 1, prob = x_hat)

  return(dropout_ind)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
pred_ra <- function(cont, BM, MCdatatmp, IT, random.regime, param, above,
                    cutoff, nat_value, incremental, L1, TG, bsize){ #, rank){ # above = NULL
  if(ncol(MCdatatmp)==1) {
    MCdata <- as.matrix(cbind(1, MCdatatmp))
  } else {
    MCdata <- MCdatatmp
  }

  X <- MCdata

  if(L1==TRUE & TG>0){
      while(length(BM$ecdfs)<ncol(X)){
          X <- X[,-bsize]
      }
  }

  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  if(nat_value == TRUE){ # based on the natural value
    if(cont == FALSE) {
      fi_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
      new_fi <- rbinom(length(fi_hat), 1, prob = fi_hat)
    } else if(cont == TRUE) {
      fi_hat <- as.numeric(BM$forest$predict_iteration(X, IT)) * BM$sd_Y + BM$mu_Y
      new_fi <- rnorm(length(fi_hat), mean = fi_hat, sd = BM$sigma)
      if(above == TRUE){
        interv <- which(new_fi > cutoff)
      } else if (above == FALSE) {
        interv <- which(new_fi < cutoff)
      }
      if(random.regime == "uniform") {
        int_dist <- runif(length(interv), param[1], param[2])
      } else if (random.regime == "normal") {
        int_dist <- rnorm(length(interv), param[1], param[2])
      } else if(random.regime == "triangular"){
        int_dist <- rtri(length(interv), param[1], param[2], param[3])
      } else if (random.regime == "binomial") {
        int_dist <- rbinom(length(interv), 1, param[1])
      }
      if (incremental == TRUE){
        new_fi[interv] <- new_fi[interv] + int_dist
        #new_fi[interv] <- ifelse(new_fi[interv]<110, 0, 1)*new_fi[interv] + ifelse(new_fi[interv]<110, 1, 0)*110

      } else {
        new_fi[interv] <- int_dist
      }

    }
  } else if(nat_value == FALSE){ # not based on the natural value
    if(random.regime == "uniform") {
      int_dist <- runif(nrow(MCdata), param[1], param[2])
    } else if (random.regime == "normal") {
      int_dist <- rnorm(nrow(MCdata), param[1], param[2])
    } else if(random.regime == "triangular"){
      int_dist <- rtri(nrow(MCdata), param[1], param[2], param[3])
    } else if (random.regime == "binomial") {
      int_dist <- rbinom(nrow(MCdata), 1, param[1])
    }
    if (incremental == TRUE){
      new_fi <- new_fi + int_dist
    } else {
      new_fi <- int_dist
    }
  }
  new_MCdata <- cbind(MCdata, new_fi)

  if(L1==TRUE & TG>1){
      new_MCdata <- new_MCdata[,-bsize]
  } else {
      new_MCdata <- new_MCdata
  }

  return(new_MCdata)
}


#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
pred_surv <- function(BM, MCdatatmp, IT, L1, TG, bsize){
  if(ncol(MCdatatmp)==1) {
    MCdata <- as.matrix(cbind(1, MCdatatmp))
  } else {
    MCdata <- MCdatatmp
  }

  X <- MCdata

  if(L1==TRUE & TG>0){
      while(length(BM$ecdfs)<ncol(X)){
          X <- X[,-bsize]
      }
  }

  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }

  s_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
  return(s_hat)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
pred_y <- function(cont, BM, MCdata, IT, L1, TG, bsize) {
  X <- MCdata

  if(L1==TRUE & TG>0){
      while(length(BM$ecdfs)<ncol(X)){
          X <- X[,-bsize]
      }
  }

  for(i in 1:ncol(X)) {
    X[,i] <- BM$ecdfs[[i]](X[,i])
  }
  if(cont == FALSE) {
    y_hat <- pnorm(as.numeric(BM$forest$predict_iteration(X, IT)) + BM$offset)
  } else if(cont == TRUE) {
    y_hat <- as.numeric(BM$forest$predict_iteration(X, IT)) * BM$sd_Y + BM$mu_Y
  }
  return(y_hat)
}

#' @title Internal helper function
#' @description This function supports the main gcomputation workflow and is not intended for direct use.
#' @keywords internal
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' @title Internal helper function
#' @description Compute Ratio Adjustment for Survival Probabilities.
#' Returns the maximum of 1 and the ratio \code{a/b}. Used in the computation of the SACE adjustment term.
#'
#' @param a Numeric. Survival probability under intervention.
#' @param b Numeric. Survival probability under natural course.
#'
#' @return A numeric value equal to \code{max(1, a/b)}.
#' @keywords internal
phi <- function(a, b) {
  ratio <- a / b
  ifelse(ratio < 1, 1, ratio)
}

#' @title Draw Samples from a Triangular Distribution
#'
#' @description Generates random samples from a triangular distribution using the inverse transform method.
#'
#' @param n Integer. Number of samples to draw.
#' @param min Numeric. Lower bound of the distribution.
#' @param max Numeric. Upper bound of the distribution.
#' @param mode Numeric. Mode (peak) of the distribution.
#'
#' @return A numeric vector of length \code{n} containing samples from the triangular distribution.
#' @examples
#' rtri(5, min = 0, max = 1, mode = 0.8)
#'
#' @export
rtri <- function(n, min, max, mode) {
  u <- runif(n)
  c <- (mode - min) / (max - min)
  ifelse(u < c,
         min + sqrt(u * (max - min) * (mode - min)),
         max - sqrt((1 - u) * (max - min) * (max - mode)))
}

#' Compute SACE Adjustment Term
#'
#' Computes the Survivor Average Causal Effect (SACE) adjustment term for the SAIE estimand.
#' The adjustment accounts for the "always survivor" stratum under hypothetical interventions.
#'
#' @param Delta_max Numeric. Upper bound for the Delta parameter.
#' @param ps_int Numeric. Survival probability under the intervention regime.
#' @param ps_bl Numeric. Survival probability under the natural course (baseline).
#' @param draw_params Logical. If \code{TRUE}, draws \code{Delta} and \code{lambda} from triangular distributions.
#' If \code{FALSE}, uses fixed values provided via \code{Delta} and \code{lambda}.
#' @param Delta Numeric. Optional fixed value for Delta (required if \code{draw_params = FALSE}).
#' @param lambda Numeric. Optional fixed value for lambda (required if \code{draw_params = FALSE}).
#'
#' @details
#' When \code{draw_params = TRUE}, \code{Delta} is drawn from a triangular distribution with
#' bounds \code{0} and \code{Delta_max}, and \code{lambda} is drawn from a triangular distribution
#' with bounds \code{0.5} and \code{1}.
#'
#' @return A numeric value representing the SACE adjustment term.
#'
#' @examples
#' # Example with random draws
#' sace(Delta_max = 2, ps_int = 0.8, ps_bl = 0.7)
#'
#' # Example with fixed parameters
#' sace(Delta_max = 2, ps_int = 0.8, ps_bl = 0.7, draw_params = FALSE, Delta = 1.5, lambda = 0.8)
#'
#' @export
sace <- function(Delta_max, ps_int, ps_bl, draw_params = TRUE, Delta = NULL, lambda = NULL) {
  U <- phi(ps_int, ps_bl)
  U2 <- 1 / U

  if (draw_params) {
    Delta <- rtri(1, 0, Delta_max, Delta_max - (Delta_max / 100))
    lambda <- rtri(1, 0.5, 1, 1 - (1 / 100))
  } else {
    if (is.null(Delta) || is.null(lambda)) {
      stop("Provide Delta and lambda when draw_params = FALSE")
    }
  }

  Delta * ((ps_int + lambda * (U - ps_int)) * (1 - U2))
}
