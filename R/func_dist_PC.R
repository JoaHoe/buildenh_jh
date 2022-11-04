#name of script: func_dist_PC.R
#purpose: derivation of lengths from coordinates of intersected points
#used in: intersect_corner_points.R
#GNU General Public License (GPL)

dist_PC <- function(P_red) {
  x_1 <- (P_red[1,2]) #change to first and last point
  x_last <- (P_red[nrow(P_red),2])
  y_1 <- (P_red[1,3])
  y_last <- (P_red[nrow(P_red),3])
  d_line <- sqrt((x_1-x_last)^2 + (y_1-y_last)^2) #length of line [pixel]
  return(d_line)
} #end function dist_PC

#end of func_dist_PC.R