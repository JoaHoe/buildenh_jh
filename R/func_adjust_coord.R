##name of script: func_adjust_coord.R
#function for adjustment of corner coordinates
#GNU General Public License (GPL)
adjust_coord <- function(A,b) { 
  At <- t(A)
  AtA <- At %*% A
  Atb <- At %*% b
  x <- rep(0,m)
  AtAinv <- solve(At %*% A)
  x <- AtAinv %*% Atb
  x #adjusted ro-values
  x1 <- x
  
  ## calculation of residuals
  r <- A %*% x1 - b
  r #residuals [pixels]
  r1 <- abs(r)
  max(r1) #maximal residual [pixels]
  rt <- t(r)
  sigma <- sqrt((rt %*% r)/m)
  cat("standard deviation of residuals=",sigma,"[pixel]","\n")
  
  ## test with threshold (thr_r)
  thr_r <- 3*sigma
  if (max(r1) > thr_r) {
    cat("warning: maximum residual exceeds threshold", "\n")
  }
  
  ## calculation of adjusted coordinates
  p <- A%*%x1
  return(p)
} #end of function 'adjust_coord(A,b)'

#end of script 'func_adjust_coord.R'
##########################################################################