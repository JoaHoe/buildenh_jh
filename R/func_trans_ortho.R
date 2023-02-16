##name of script: func_trans_ortho.R
cat("version_number= ",v_nr,"\n")
#purpose: calculation of transformation parameter
#         with measurement of check points
#function used in: support_line_detection.R
#author: Joachim Höhle
#GNU General Public License (GPL)


trans_ortho <- function() {
  mar=(wind_x-orig_x)/4
  
  #check point 1
  x1 <- xc - mar 
  y1 <- yc + mar
  points(x1-orig_x, y1-orig_y, pch=16, cex=1.5, asp=1, col="white")
  
  #check point 2
  x2 <- xc + mar
  y2 <- yc - mar
  points(x2-orig_x, y2-orig_y, pch=16, cex=1.5, asp=1, col="white")
  
  #check point 7
  x7=xc
  y7=yc
  points(x7-orig_x, y7-orig_y, pch=16, cex=1.5, asp=1, col="white")
  
  #locator-measurements
  c1 <- locator(1) #left lower control point (1)
  c2 <- locator(1) #right upper control point (2)
  c7 <- locator(1) #check point (center of object )
  #
  #calculation of transformation-parameter (plane)
  dX <- x2 - x1
  dY <- y2 - y1
  dx <- c2$x - c1$x
  dy <- c2$y - c1$y
  N <- dx^2 + dy^2
  #
  a1 <- (dx*dX + dy*dY)/N
  b1 <- (dx*dY - dy*dX)/N
  a0 <- x1 - a1*c1$x + b1*c1$y
  b0 <- y1 - b1*c1$x - a1*c1$y
  #
  x1_ch <- a0 + a1*c1$x - b1*c1$y
  y1_ch <- b0 + b1*c1$x + a1*c1$y
  x2_ch <- a0 + a1*c2$x - b1*c2$y
  y2_ch <- b0 + b1*c2$x + a1*c2$y
  x7_ch <- a0 + a1*c7$x - b1*c7$y
  y7_ch <- b0 + b1*c7$x + a1*c7$y
  points(x7_ch-orig_x, y7_ch-orig_y,pch=3, cex=5,asp=1, col="blue")
  #
  tr_lat <- c(a0,b0)
  D <- matrix(nrow=2, ncol=2)
  D[1,1] <- a1
  D[1,2] <- -b1
  D[2,1] <- b1
  D[2,2] <- a1
  kf2 <- sqrt(a1^2+b1^2)
  L1 <- list(D,tr_lat,kf2)
  return(L1)
} #end of function 'trans_ortho()'

#end of 'func_trans_ortho.R'
