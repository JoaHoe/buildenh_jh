##script name: func_rectang_lines.R
#purpose: detecting lines which are parallel/orthogonal 
#         to main orientation
#used in: line_detection.R

rectang_lines <- function(B5_4,theta_ref3,wd3,thr3) { 
  B5_4
  k14 <- length(B5_4$lnr)
  B5_4b <- B5_4
  B5_4b[,1:7] <- 0
  B5_4b[1,] <- B5_4[1,]

  #loop
  i=1
  k1=2
  while (i < k14) {
    i=i+1
    if (B5_4$theta_angle[i] == theta_ref3 && B5_4$n_pixel[i] >= wd3 ||
        B5_4$theta_angle[i] == alph_ref && B5_4$n_pixel[i] >= wd3) {
      B5_4b[k1,] <- B5_4[i,]
      k1 <- k1+1
    } #end if
  } #end loop
  
  B5_4b
  B5_4c <- subset(B5_4b,B5_4b$lnr > 0)
  B5_4c #matrix with lines longer than 'wd3'
  k15 <- length(B5_4c$lnr)
  
  ## Search of lines with theta_ref3
  vec <- 1:k15
  B5_4d <- B5_4c
  B5_4d[,1:7] <- 0
  
  #loop
  j=1
  for (i in vec) {
    if (B5_4c$theta_angle[i] == theta_ref3) {
      B5_4d[j,] <- B5_4c[i,]
      j <- j + 1
    } #end if
  } # end for-loop
  
  B5_4d <- subset(B5_4d,B5_4d$n_pixel >= wd3)
  B5_long_lines <- B5_4d[1:2,]
  B5_4d_ord <- B5_4d[order(B5_4d[,6],decreasing=FALSE),]
  B5_4d_ord2 <- subset(B5_4d_ord,B5_4d_ord$n_pixel >= wd3)
  n_B5_4d_ord2 <- length(B5_4d_ord2$lnr)
  row.names(B5_4d_ord2) <- 1 : n_B5_4d_ord2
  B5_4d_ord2
  
  ##Search of lines with alph_ref
  vec <- 1 : k15
  B5_4dd <- B5_4c
  B5_4dd[,1:7] <- 0
  
  #loop
  j = 1
  for (i in vec) {
    if (B5_4c$theta_angle[i] == alph_ref) {
      B5_4dd[j,] <- B5_4c[i,]
      j = j + 1
    } #end if
  } #end for-loop
  
  B5_4dd
  B5_4dd <- subset(B5_4dd, B5_4dd$n_pixel >= wd3)
  B5_long_lines[3:4,] <- B5_4dd[1:2,] # matrix for weighted mean of angle
  B5_long_lines
  min(B5_long_lines$n_pixel)
  B5_4dd_ord <- B5_4dd[order(B5_4dd[,6],decreasing=FALSE),]
  B5_4dd_ord
  # end of preparation
  
  B5_4d_ord
  B5_4d_ord[,8] <- 0
  names(B5_4d_ord)[8] <- "dif_ro_pixel"
  B5_4d_ord_red <- subset(B5_4d_ord, B5_4d_ord$n_pixel >= wd3)
  k6 <- nrow(B5_4d_ord_red)
  rownames(B5_4d_ord_red) <- 1 : k6
  B5_4d_ord_red
  B5_4d_ord_red[,9:10] <- 0
  names(B5_4d_ord_red)[9] <- "shorter_line"
  names(B5_4d_ord_red)[10] <- "remove_row_nr"
  B5_4d_ord_red
  k5 <- nrow(B5_4d_ord_red)
  vec5 <- 1 : (k5-1)
  
  for (n in vec5) {
    B5_4d_ord_red$dif_ro_pixel[n] <- B5_4d_ord_red$ro_pixel[n+1] - B5_4d_ord_red$ro_pixel[n]
    B5_4d_ord_red$shorter_line[n] <- min(B5_4d_ord_red$n_pixel[n],B5_4d_ord_red$n_pixel[(n+1)])
  } # end for-loop
  B5_4d_ord_red
  
  for (n in vec5) {
    
    if(B5_4d_ord_red$dif_ro_pixel[n] <= thr3 && B5_4d_ord_red$n_pixel[n] == B5_4d_ord_red$shorter_line[n]) {
      B5_4d_ord_red$remove_row_nr[n] <- n 
    } #end if
    
    if(B5_4d_ord_red$dif_ro_pixel[n] <= thr3 && B5_4d_ord_red$n_pixel[n+1] == B5_4d_ord_red$shorter_line[n]) {
      B5_4d_ord_red$remove_row_nr[n+1] <- n+1 
    } #end if
    
  } # end for-loop
  
  B5_4d_ord_red
  B5_4d_ord_red$remove_row_nr
  k7 <- nrow(B5_4d_ord_red)
  vec6 <- 1 : k7
  rem_vec <- rep(0,k7)
  
  for (n in vec6) {
    if(B5_4d_ord_red$remove_row_nr[n] > 0) {
      rem_vec[n] <- B5_4d_ord_red$remove_row_nr[n]
    }
  } # end for-loop
  rem_vec
  #
  rem_vec2 <- NULL
  j = 1
  for (n in vec6) {
    if (rem_vec[n] == 0) {next} 
    else {
      rem_vec2[j] <- rem_vec[n]
      j=j+1
    }
  } #end for-loop
  
  rem_vec2
  le3 <- length(rem_vec2)
  #
  if (le3 == 0) {
    B5_4d_ord_red2 <- B5_4d_ord_red  #solution
  } else {
    B5_4d_ord_red2 <- B5_4d_ord_red[-rem_vec2,] #correction of matrix with rem_vec
    B5_4d_ord_red2
    k8 <- nrow(B5_4d_ord_red2)
    rownames(B5_4d_ord_red2) <- 1:k8
  } #end if-else
  
  B5_4d_ord_red2
  
  if (k8 < 2) {
    stop("error: number of lines is < 2","\n")
  }
  
  if (k8 == 2) { 
    B5_4d_ord_final <- B5_4d_ord_red2 #solution
  } 
  
  if (k8 > 2) {
    B5_4d_ord_red3 <- B5_4d_ord_red2 #new round
    B5_4d_ord_red3[,8:10] <- 0
    vec7 <- 1:(k8-1)
    
    for (n in vec7) {
      B5_4d_ord_red3$dif_ro_pixel[n] <- (B5_4d_ord_red3$ro_pixel[n+1]- B5_4d_ord_red3$ro_pixel[n])
      B5_4d_ord_red3$shorter_line[n] <- min(B5_4d_ord_red3$n_pixel[n],B5_4d_ord_red3$n_pixel[(n+1)])
    } #end for-loop
    
    # select the longer line of two
    
    for (n in vec7) { #condition 1 and 2
      if (B5_4d_ord_red3$dif_ro_pixel[n] <= thr3 && B5_4d_ord_red3$n_pixel[n] == B5_4d_ord_red3$shorter_line[n]) {
        B5_4d_ord_red3$remove_row_nr[n] <- n
      } # end if1
      
      if(B5_4d_ord_red3$dif_ro_pixel[n] <= thr3 && B5_4d_ord_red3$n_pixel[n+1] == B5_4d_ord_red3$shorter_line[n]) {
        B5_4d_ord_red3$remove_row_nr[n+1] <- n+1
      } # end if2
      
    } # end for-loop
  
    k9 <- nrow(B5_4d_ord_red3)
    vec9 <- 1:k9
    rem_vec <- rep(0,k9)
    
    for (n in vec9) {
      if (B5_4d_ord_red3$remove_row_nr[n] > 0) { #if3
        rem_vec[n] <- B5_4d_ord_red3$remove_row_nr[n]
        cat("rem_vec=", rem_vec[n], "\n")
      } # end if3
    } # end for-loop
    
    rem_vec2 <- NULL
    j = 1
    
    for (n in vec9) {
      
      if(rem_vec[n] == 0) {next} else {
        rem_vec2[j] <- rem_vec[n]
        j=j+1
      } #end if-else
      
    } #end for-loop
    
    if (sum(rem_vec) != 0) {
      B5_4d_ord_red4 <- B5_4d_ord_red3[-rem_vec2,]
    } else {
      B5_4d_ord_red4 <- B5_4d_ord_red3
    } #end if-else
    
    k10 <- nrow(B5_4d_ord_red4)
    
    if (k10 < 2) {
      stop("error","\n") 
    }
    
    if (k10 == 2) { 
      B5_4d_ord_final <- B5_4d_ord_red2 #solution
    }
    
    if (k10 >= 2) { #last step
      rownames(B5_4d_ord_red4) <- 1 : k10
      B5_4d_ord_red4[,8:10] <- 0
      vec10 <- 1 : (k10-1)
      
      for (n in vec10) {
        B5_4d_ord_red4$dif_ro_pixel[n] <- B5_4d_ord_red4$ro_pixel[n+1] - B5_4d_ord_red4$ro_pixel[n]
        B5_4d_ord_red4$shorter_line[n] <- min(B5_4d_ord_red4$n_pixel[n], B5_4d_ord_red4$n_pixel[(n+1)])
      } # end for-loop
      
    } #end if last step
    
    B5_4d_ord_red4
    B5_4d_ord_red5 <- subset(B5_4d_ord_red4, n_pixel >= wd3)
    B5_4d_ord_final <- B5_4d_ord_red5
    k11 <- nrow(B5_4d_ord_red5)
  } #end if (k8 >= 2)
  #
  B5_4d_ord_final
  
  ## lines orthogonal to ref
  B5_4dd_ord
  B5_4dd_ord[,8] <- 0
  names(B5_4dd_ord)[8] <- "dif_ro_pixel"
  B5_4dd_ord_red <- B5_4dd_ord
  k6<-nrow(B5_4dd_ord_red)
  rownames(B5_4dd_ord_red) <- 1:k6
  B5_4dd_ord_red[,9:10] <- 0
  names(B5_4dd_ord_red)[9] <- "shorter_line"
  names(B5_4dd_ord_red)[10] <- "remove_row_nr"
  B5_4dd_ord_red
  vec5 <- 1 : (k6-1)
  
  for (n in vec5) {
    B5_4dd_ord_red$dif_ro_pixel[n] <- B5_4dd_ord_red$ro_pixel[n+1] - B5_4dd_ord_red$ro_pixel[n]
    B5_4dd_ord_red$shorter_line[n] <- min(B5_4dd_ord_red$n_pixel[n],B5_4dd_ord_red$n_pixel[(n+1)])
  } #end for-loop
  
  B5_4dd_ord_red
  
  for (n in vec5) {
    
    if(B5_4dd_ord_red$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red$n_pixel[n] == B5_4dd_ord_red$shorter_line[n]) {
      B5_4dd_ord_red$remove_row_nr[n] <- n 
    } #end if 
    
    if(B5_4dd_ord_red$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red$n_pixel[n+1] == B5_4dd_ord_red$shorter_line[n]) {
      B5_4dd_ord_red$remove_row_nr[n+1] <- n+1 
    } #end if 
    
  } #end for-loop
  
  B5_4dd_ord_red
  k7 <- nrow(B5_4dd_ord_red)
  #
  
  if (k7 < 1) {
    stop("error - number of lines is < 1","\n") 
  }
  
  if (k7 == 1) {
    B5_4dd_ord_final <- B5_4dd_ord_red  
  }
  
  if (k7 == 2) {
    B5_4dd_ord_final <- B5_4dd_ord_red
  } else { #else 1
    vec6 = 1:k7
    rem_vec <- rep(0,k7)
    for (n in vec6) {
      
      if(B5_4dd_ord_red$remove_row_nr[n] > 0)  {
        rem_vec[n] <- B5_4dd_ord_red$remove_row_nr[n]
      }
      
    } #end for-loop
    
    rem_vec2 <- subset(rem_vec, rem_vec > 0) 
    
    if (sum(rem_vec) != 0 & nrow(B5_4dd_ord_red) > 2) {
      B5_4dd_ord_red2 <- B5_4dd_ord_red[-rem_vec2,] #some line segments will be removed
    } else { #else2
      B5_4dd_ord_red2 <- B5_4dd_ord_red
    } #end if-else2
  } #end if-else1
  B5_4dd_ord_red2
  #
  k8 <- nrow(B5_4dd_ord_red2)
  rownames(B5_4dd_ord_red2) <- 1 : k8
  B5_4dd_ord_red2[,8:10] <- 0
  
  if (k8 < 2) {
    stop("error","\n") 
  }
  if (k8 == 2) {
    B5_4dd_ord_final <- B5_4dd_ord_red2
  } else {
    vec7 <- 1:(k8-1)
    vec7
    
    for (n in vec7) {
      B5_4dd_ord_red2$dif_ro_pixel[n] <- (B5_4dd_ord_red2$ro_pixel[n+1] - B5_4dd_ord_red2$ro_pixel[n])
      B5_4dd_ord_red2$shorter_line[n] <- min(B5_4dd_ord_red2$n_pixel[n],B5_4dd_ord_red2$n_pixel[(n+1)])
    } #end for-loop
    
  } #end if-else

  k9 <- nrow(B5_4dd_ord_red2)
  
  if (k9 == 2) {
    B5_4dd_ord_final <- B5_4dd_ord_red2
  } else {
    vec7 <- 1 : (k9-1)
    
    for (n in vec7) {
      
      if(B5_4dd_ord_red2$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red2$n_pixel[n] == B5_4dd_ord_red2$shorter_line[n]) {
        B5_4dd_ord_red2$remove_row_nr[n] <- n 
      } #end if 1
      
      if(B5_4dd_ord_red2$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red2$n_pixel[n+1] == B5_4dd_ord_red2$shorter_line[n]) {
        B5_4dd_ord_red2$remove_row_nr[n+1] <- n+1 
      } #end if2
      
    } #end for-loop
  } #end if-else
  k10 <- nrow(B5_4dd_ord_red2)
  vec8 = 1 : k10
  rem_vec <- rep(0,k10)
  
  for (n in vec8) {
    rem_vec[n] <- B5_4dd_ord_red2$remove_row_nr[n]
  } #end for-loop
  
  rem_vec2 <- subset(rem_vec,rem_vec > 0)

  if (sum(rem_vec) != 0) {
    B5_4dd_ord_red3 <- B5_4dd_ord_red2[-rem_vec2,]
  } else {
    B5_4dd_ord_red3 <- B5_4dd_ord_red2
  }
  B5_4dd_ord_red3
  
  B5_4dd_ord_red4 <- subset(B5_4dd_ord_red3, n_pixel >= wd3)
  k11 <- nrow(B5_4dd_ord_red4)
  rownames(B5_4dd_ord_red4) <- 1 : k11
  B5_4dd_ord_red4
  #
  
  #start determination of final line segments
  if (k11 < 1) {
    stop("error - number of lines < 1","\n") 
  }
  
  if (k11 == 1) {
    B5_4dd_ord_final <- B5_4dd_ord_red4
  }
  
  if ( k11 == 2 ) {B5_4dd_ord_final <- B5_4dd_ord_red4} 
  
  if (k11 > 2) {
    B5_4dd_ord_red4[,8:10] <- 0
    rownames(B5_4dd_ord_red4)
    vec9 <- 1 : (k11-1)
    
    for (n in vec9) {
      B5_4dd_ord_red4$shorter_line[n] <- min(B5_4dd_ord_red4$n_pixel[n],B5_4dd_ord_red4$n_pixel[(n+1)])
      B5_4dd_ord_red4$dif_ro_pixel[n] <- (B5_4dd_ord_red4$ro_pixel[n+1] - B5_4dd_ord_red4$ro_pixel[n])
    } #end for-loop
    
    for (n in vec9) {
      if(B5_4dd_ord_red3$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red4$n_pixel[n] == B5_4dd_ord_red4$shorter_line[n]) {
        B5_4dd_ord_red4$remove_row_nr[n] <- n
      } #end if
      
      if(B5_4dd_ord_red4$dif_ro_pixel[n] <= thr3 && B5_4dd_ord_red4$n_pixel[n+1] == B5_4dd_ord_red4$shorter_line[n]) {
        B5_4dd_ord_red4$remove_row_nr[n+1] <- n+1
      } #end if
      
    } #end for-loop
    
    k12 <- nrow(B5_4dd_ord_red4)
    vec10 <- 1 : k12
    rem_vec <- rep(0,k12)
    
    for (n in vec10) {
      rem_vec[n] <- B5_4dd_ord_red4$remove_row_nr[n]
    } #end for-loop
    
    rem_vec
    rem_vec2 <- subset(rem_vec, rem_vec > 0)
    rem_vec2
    
    if (sum(rem_vec) != 0) {
      B5_4dd_ord_red5 <- B5_4dd_ord_red4[-rem_vec2,] #some line segments will be removed
    } else {
      B5_4dd_ord_red5 <- B5_4dd_ord_red4
    } #end if-else 
    
    B5_4dd_ord_red5
    B5_4dd_ord_red6 <- subset(B5_4dd_ord_red5, n_pixel >= wd3)
    k13<-nrow(B5_4dd_ord_red6)
    rownames(B5_4dd_ord_red6) <- 1 : k13
    
    if (k13 < 2) {
      stop("error: number of lines < 2","\n") 
    }
    
    if (k13 >= 2) {
      B5_4dd_ord_final <- B5_4dd_ord_red6
    } #end if
    
  } #end if 
  
  B5_4dd_ord_final
  B5_4d_ord_final
  l1 <- length(B5_4d_ord_final$lnr)
  l2 <- length(B5_4dd_ord_final$lnr)
  
  if (l1 == l2) {
    lnr_det3 <- B5_4d_ord_final$lnr
    lnr_det3[(l1+1) : (l1+l2)] <- B5_4dd_ord_final$lnr
  } else {
    cat("number of line segments is odd - manual interaction is required","\n")
  } #end if-else
  
  # if manual interaction is required
  B5_4d_ord_final
  B5_4dd_ord_final
  l1 <- length(B5_4d_ord_final$lnr)
  l2 <- length(B5_4dd_ord_final$lnr)
  lnr_det3 <- rep(0,(l1+l2))
  lnr_det3 <- B5_4d_ord_final$lnr
  lnr_det3[(l1+1) : (l1+l2)] <- B5_4dd_ord_final$lnr
  lnr_det3 #solution for all line segments of building outline (cas = '100+all')
  l3 <- l1 + 1
  l4 <- l1 + l2
  B5_6 <- matrix(nrow = l4, ncol = 7) 
  B5_6 <- data.frame(B5_6)
  names(B5_6)<-c("lnr","theta_index","ro_index","n","theta_angle","ro_pixel","n_pixel")
  B5_6[,]<-0
  B5_6
  B5_6[1:l1,1] <- B5_4d_ord_final[,1]
  B5_6[l3:l4,1] <- B5_4dd_ord_final[,1]
  B5_6[1:l1,2:7] <- B5_4d_ord_final[,2:7]
  B5_6[l3:l4,2:7] <- B5_4dd_ord_final[,2:7]
  B5_6[,8] <- 0
  names(B5_6)[8] <- "ortho"
  len <- length(B5_6$ortho)
  vec <- 1 : len
  
  for (i in vec) {
    
    if (B5_6$theta_angle[i] == theta_ref3 || B5_6$theta_angle[i] == alph_ref) {  
      B5_6$ortho[i] <- 1
    } else {
      B5_6$ortho[i] <- 0
    } #end if-else
    
  } #end for-loop 
  return(B5_6)
} #end of function 'rectang_lines(B5_4,theta_ref3,wd3,thr3)'

#end of script 'func_rectang_lines.R