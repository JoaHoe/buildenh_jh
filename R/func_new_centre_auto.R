#name of script: func_new_centre_auto.R
#purpose: calculation of center of object
#function used in: support_sequence_of_lines.R
#GNU General Public License (GPL)

new_centre_auto <- function() {
  b13_angle_df
  xc_new2 <- mean(b13_angle_df$x_centre)
  yc_new2 <- mean(b13_angle_df$y_centre)
  coo2 <- c(xc_new2, yc_new2)
  cat("new coordinates2:", "xc = ",  coo2[1], "yc = ", coo2[2], "\n")
  points(xc_new2, -yc_new2, asp=1, pch=3, cex=3,col="green") #plot in small scale image 
  return(coo2)
} #end of function new_centre_auto

#end of script 'func_new_centre_auto.R