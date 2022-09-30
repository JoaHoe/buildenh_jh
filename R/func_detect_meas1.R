#script: func_detect_meas1.R
#purpose: detect lines by measurement of one pixel for rectangular objects
#function used in: support_line_detection.R
#instructions: to be used after locator2()

detect_meas1 <- function() { 
  B5_4_ord <- B5_4[order(B5_4$ro_pixel,decreasing = FALSE),]
  k1 <- nrow(B5_4_ord)
  row.names(B5_4_ord) <- 1 : k1
  pos1 <- list(x=NA,y=NA)
  pos1$x <- x #value of locator2 measurement
  pos1$y <- (-y) #conversion to math-system
  pos1
  
  #lines parallel to ref-line
  theta_ref
  theta_math <- 180 - theta_ref
  theta_math #approximate reference angle
  theta_math_arc <- theta_math/omega
  ro_meas <- cos(theta_math_arc) * pos1$x + sin(theta_math_arc) * pos1$y
  ro_meas <- abs(ro_meas)
  ro_meas <- round(ro_meas)
  cat("ro_meas_ref= ", ro_meas, "\n")
  
  B5_4_ord
  k30 <- nrow(B5_4_ord)
  i=1
  n_line=0
  while (i <= k30) {
    dif_ref <- abs(ro_meas - B5_4_ord$ro_pixel[i])
    if (dif_ref <= 10) { #threshold 10 pixel
      cat("line number=", B5_4_ord$lnr[i],"\n")
      k31 <- i
      print(B5_4_ord[k31,])
      n_line <- n_line+1
    } #end if
    i <- i+1
  } #end while-loop
  cat("n_line_ref= ",n_line,"\n")
  
  #line orthogonal to ref-line
  theta_math <- (-theta_ref) + 180
  theta_orth <- theta_math + 90
  
  if (theta_orth > 180) {
    theta_orth <- theta_orth - 180
  }
  
  theta_orth_arc <- theta_orth / omega
  ro_meas <- cos(theta_orth_arc) * pos1$x + sin(theta_orth_arc) * pos1$y
  ro_meas <- round(abs(ro_meas))
  cat("ro_meas_orth= ", ro_meas, "\n")
  
  n_line <- 0
  i <- 1
  while (i <= k30) {
    dif_orth <- abs(ro_meas - B5_4_ord$ro_pixel[i])
    if (dif_orth <= 10)  {
      cat("line number=", B5_4_ord$lnr[i],"\n")
      k31 <- i
      print(B5_4_ord[k31,])
      n_line <- n_line + 1
    } #end if
    i <- i+1
  } #end while-loop
  #
  cat("n_line_orth= ",n_line,"\n") #select the line where ro_pixel is close to 'ro_meas'
} #end of function 'detect_meas1'

#end of script 'func_detect_meas1.R'