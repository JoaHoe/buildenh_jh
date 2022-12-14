## script:'func_w_av.R'
#purpose: weighted average of theta-angle (theta_av)
#used in: intersect_corner_points.R
#GNU General Public License (GPL)

#function
w_av <- function(ang,len) { 
  x <- length(ang)
  y <- 1 : x
  pw <- len/max(len) #weights
  sp <- sum(pw)
  theta <- rep(0,x)
  
  for (i in y) {
    theta[i] <- (ang[i]*pw[i])/sp
  }
  
  theta_av <- sum(theta)
  return(theta_av)
} #end of function 'w_av(ang,len)'

#end of script 'func_w_av.R
