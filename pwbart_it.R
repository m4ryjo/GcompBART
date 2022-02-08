pwbart_it = function(
   x.test,	# x matrix to predict at
   treedraws,	# treedraws from bart
   mu=0,	#mean to add on
   it, #it
   transposed = FALSE
)
{
   if (!transposed) {
      x.test <- t(bartModelMatrix(x.test))
   }
   res = .Call("cpwbart_it",
               treedraws,
               x.test,
               it #thread count
               )
   return(res$yhat.test[it:(it+length(it)-1),]+mu)
}
