##name of script: func_design_mat.R
cat("version_number= ",v_nr,"\n")
#purpose: solution for closed polygon
#function 'design_mat(m,phi)'
#author: Joachim Höhle
#GNU General Public License (GPL)

design_mat <- function(m,phi) { 
  y <- 1 : (m-1)
  
  for (i in y) {
    k1 <- sin(phi[i]) - cos(phi[i]) * tan(phi[i+1])
    k2 <- cos(phi[i+1]) * tan(phi[i]) - sin(phi[i+1])
    A[2*i-1,i] <<- 1/cos(phi[i])-tan(phi[i])/k1
    A[2*i,i] <<- 1/k1
    A[2*i-1,i+1] <<- tan(phi[i])/k2
    A[2*i,i+1] <<- (-1/k2)
  } #end of loop
  
  i <- m
  k1 <- sin(phi[i]) - cos(phi[i]) * tan(phi[1])
  k2 <- cos(phi[1]) * tan(phi[i]) - sin(phi[1])
  A[2*i-1,i] <<- 1/cos(phi[i]) - tan(phi[i])/k1
  A[2*i,i] <<- 1/k1
  A[2*i-1,1] <<- tan(phi[i])/k2
  A[2*i,1] <<- (-1/k2)
  return(A)
} # end of function 'design_mat(m,phi)'

#end of script 'func_design_mat.R'
