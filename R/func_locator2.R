#name of script: func_locator2.R
#purpose: digitizing and marking of a single pixel in orthoimage-section
#function used in: 'support_line_detection.R'
#instructions: point to pixel and click by mouse 

locator2 <- function() {
  c9 <- locator(1) #standard function
  loc <- c(c9$x,c9$y)
  print(tr_lat)
  print(D)
  pts9 <- tr_lat + D%*%loc #transformation to image-system
  x9 <- pts9[1,1]
  y9 <- pts9[2,1]
  x <<- x9
  y <<- y9
  points(x-orig_x,y-orig_y,pch=3, asp =1, cex=5,asp=1, col="blue")
  cat("x_coordinate= ",x,sep = "","\n")
  cat("y_coordinate= ",y,sep = "","\n")
} #end of function locator2

#end of script 'func_locator2.R'
