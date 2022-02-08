pred_y <- function(cont, BM, MCdata, IT) {
  if(cont == FALSE) {
    y_hat <- pnorm(pwbart_it(MCdata, BM$treedraws, mu = BM$binaryOffset, it=IT))
  } else if(cont == TRUE) {
    y_hat <- pwbart_it(MCdata, BM$treedraws, mu = BM$mu, it = IT)
  }
  return(y_hat)
}
