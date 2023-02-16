##name of script: func_detect_thet_len.R
cat("version_number= ",v_nr,"\n")
#purpose: determination of angle theta
#function used in: 'support_line_detection.R'
#author: Joachim Höhle
#GNU General Public License (GPL)

detect_thet_len <- function() { 
  theta_index <- 18 # (=85 degrees)
  w3 <- 15 * k #line-length = 15 pixel
  n_B <- nrow(B0)
 
  i <- 1
  while (i <= n_B) {
    
    if (B0[i,2] == theta_index && B0[i,4] >= w3) { 
      cat("lnr= ",i,"\n")
    }
    
    i <- i + 1
  } #end while

} #end of function 'detect_thet_len()'

#end of script 'func_detect_thet_len.R' 
