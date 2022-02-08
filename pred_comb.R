pred_comb <- function(cont, BM, MCdata, IT){
  if(cont == FALSE) {
    x_hat <- pnorm(pwbart_it(MCdata, BM$treedraws, mu = BM$binaryOffset, it=IT))
    new_MCdata <- rbind(MCdata, rbinom(length(x_hat), 1, prob = x_hat))
  } else if(cont == TRUE) {
    x_hat <- pwbart_it(MCdata, BM$treedraws, mu = BM$mu, it = IT)
    new_MCdata <- rbind(MCdata, rnorm(length(x_hat), mean = x_hat, sd = mean(BM$sigma)))
  }
  return(new_MCdata)
}
