##function for generation of a design matrix
#solution for open polygon
#GNU General Public License (GPL)

design_mat2 <- function(m,phi) { 
  y1 <- 1 : (m-1)
  for (i in y1) {
    k1 <- sin(phi[i]) - cos(phi[i]) * tan(phi[i+1])
    k2 <- cos(phi[i+1]) * tan(phi[i]) - sin(phi[i+1])
    A[2*i-1,i] <- 1/cos(phi[i])-tan(phi[i])/k1
    A[2*i,i] <- 1/k1
    A[2*i-1,i+1] <- tan(phi[i])/k2
    A[2*i,i+1] <- (-1/k2)
  } #end of i-loop
  return(A)
} # end of function 'design_mat2'

