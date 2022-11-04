## name of script: func_dist_v.R
## purpose: calculation of vector with distances 
## function used in: sequence_of_lines.R 
##GNU General Public License (GPL)

dist_v <- function(npx,b_angle_df_seq,refj) {
  #calculation of vector with distances
  x <- rep(5000,n_x)
  j1 <- 1
  
  while (j1 <= npx) {
    dx <- b_angle_df_seq$x_centre[refj] - b_angle_df_seq$x_centre[j1]
    dy <- b_angle_df_seq$y_centre[refj] - b_angle_df_seq$y_centre[j1]
    if (b_angle_df_seq$theta_appr[j1] != b_angle_df_seq$theta_appr[refj]) {
      x[j1]<-as.integer(sqrt(dx^2+dy^2))
    } #end if
    j1 <- j1+1
  } #end loop j1
  
  return(x)
} #end of function

#end of script 'func_dist_v.R'