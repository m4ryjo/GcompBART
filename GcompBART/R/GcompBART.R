GcompBART <- function(data,
                      var.type,
                      exp.reg,
                      cont.reg,
                      J=2000,
                      Ndraws=200,
                      Suppress = TRUE,
                      By = Ndraws/10,
                      ...
                      ) {

  ## test if continuous or binary
  continous <- apply(data, 2, function(x) !all(x %in% 0:1))
  n_O <- length(which(var.type=="O"))
  n_Fi <- length(which(var.type=="Fi"))

  ## Fit the observed data models and save as a list of models
  BModels <- vector("list", ncol(data) - 1)

  for (i in 1:(ncol(data)-1)) {
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)

    if(var.type[1 + i] == "Fi") {
      BModels[[i]] <- NA
    } else if (continous[1+i] == FALSE & var.type[1 + i] != "Fi" & Suppress == TRUE) {
      quiet(BModels[[i]]  <- pbart(data[id == 1, which(var.type[1:i]!="S")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0))
    } else if (continous[1+i] == TRUE & var.type[1 + i] != "Fi" & Suppress == TRUE) {
      quiet(BModels[[i]]  <- wbart(data[id == 1, which(var.type[1:i] != "S")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0))
    } else if (continous[1+i] == FALSE & var.type[1 + i] != "Fi" & Suppress == FALSE) {
      BModels[[i]]  <- pbart(data[id == 1, which(var.type[1:i] != "S")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0)
    } else if (continous[1+i] == TRUE & var.type[1 + i] != "Fi" & Suppress == FALSE) {
      BModels[[i]]  <- wbart(data[id == 1, which(var.type[1:i] != "S")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0)
    }
  }

  y <- matrix(nrow = 2 * n_O, ncol = Ndraws)
  ## MC-integration over "Ndraws" iterations
  for(it in 1:Ndraws) {
    k <- 1
    l <- 1
    # sample data for the first confounder
    x <- list(data.frame(sample(data[, 1], size = J, replace = TRUE)))
    for (j in 2:ncol(data)) {
      if(var.type[j] == "C") {
        x <- lapply(x, function(x2) {pred_comb(continous[j], BModels[[j - 1]], x2, it)})
      } else if(var.type[j] == "Fi") {
        x <- list(cbind(x[[1]], exp.reg[l]), cbind(x[[length(x)]], cont.reg[l]))
        l <- l + 1
      } else if(var.type[j] == "S") {
        x <- lapply(x, function(x2) {pred_surv(BModels[[j - 1]], x2, it)})
      } else if(var.type[j] == "O") {
        if(all(var.type != "Fi")){
          tmp <- lapply(x, function(x2) {pred_y(continous[j], BModels[[j - 1]], x2, it)})
          y[k, it] <- sapply(tmp, mean)
          if(continous[j] == FALSE) {
            x[[1]] <- cbind(x[[1]], rbinom(length(tmp[[1]]), 1, prob = tmp[[1]]))
          } else if(continous[j] == TRUE) {
            x[[1]] <- cbind(x[[1]], rnorm(length(tmp[[1]]), mean = tmp[[1]], sd = mean(BModels[[j - 1]]$sigma)))
          }
        } else {
          tmp <- lapply(x, function(x2) {pred_y(continous[j], BModels[[j - 1]], x2, it)})
          y[k:(k + 1), it] <- sapply(tmp, mean)
          k <- k + 2
          if(continous[j] == FALSE) {
            x[[1]] <- cbind(x[[1]], rbinom(length(tmp[[1]]), 1, prob = tmp[[1]]))
            x[[2]] <- cbind(x[[2]], rbinom(length(tmp[[2]]), 1, prob = tmp[[2]]))
          } else if(continous[j] == TRUE) {
            x[[1]] <- cbind(x[[1]], rnorm(length(tmp[[1]]), mean = tmp[[1]], sd = mean(BModels[[j-1]]$sigma)))
            x[[2]] <- cbind(x[[2]], rnorm(length(tmp[[2]]), mean = tmp[[1]], sd = mean(BModels[[j-1]]$sigma)))
          }
        }
      }
    }
    if(it %in% seq(0, Ndraws, by = By)) {
      print(paste("done ", it, " (out of ", Ndraws,")", sep=""))
    }
  }
  rownames(y) <- rep(c("Exp", "Contr"), n_O)

  means <- t(apply(y, 1, function(x) {c(mean(x), quantile(x, c(0.025, 0.975)))}))
  rownames(means) <- rep(c("Exp", "Contr"), n_O)
  colnames(means) <- c("mu", "2.5%", "97.5%")


  out <- list(means, y)
  names(out) <- list("summary", "y_hat")
  return(out)
}




