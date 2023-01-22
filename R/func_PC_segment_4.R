##script: func_PC_segment_4.R
#purpose: generation of PointClusters (PC) of one line segment
#used in: line_detection.R
#GNU General Public License (GPL)
PC_segment_4 <- function(lnr) { 
  #browser()
  i <- 0
  while (i < n_theta) {
    i <- i + 1
    H[,] <- 0
  }
  #end H null-setting
  
  ratio <- (n_ro - 1)/(ro[n_ro] - ro[1])
  X <- pc2$col
  Y <- pc2$row
  ro2 <- rep(0,n_ro)
  ro_index <- rep(0,n_ro)
  n_X <- length(pc2$col)
  P<-matrix(nrow=n_X, ncol=3)
  P[,] <- 0
  head(P)
  
  ##loop for all points of selectd point cluster(PC)
  k1=0
  k3=1
  while (k1 < n_X) {
    k1 <- k1+1
    i <- 0
    while (i < n_theta) {
      i <- i+1
      ro2[i] <- cos(theta_rad[i])*X[k1] + sin(theta_rad[i])*Y[k1]
      #ro2[i] <- abs(ro2[i]) #change
      ro_index[i] <- round(ratio*(ro2[i]-ro_1)+1) 
      #ro_index[i] <- round(ratio*(ro2[i]+ro_1)+1) 
      if (ro_index[i] >= 1 && ro_index[i] <= n_ro ) {
        k2 <- ro_index[i]
        H[i,k2] <- H[i,k2] + 1
        if (i == B[lnr,1] && k2 == B[lnr,2]) {
          P[k3,1] <- k3
          P[k3,2] <- X[k1]
          P[k3,3] <- Y[k1]
          k3 <- k3+1
        } #end if parameter
      } #end  if ro index
    } #end loop i
  } #end loop k1
  n_P <- k3-1
  P <- P[1 : n_P,]
  setwd(home_dir)
  f<-paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
  write.table(P[1:n_P,],file=f, sep="   ")
  head(P)
  return(list(P,n_P))
} #end function PC_segment_4

#end script 'func_PC_segment_4.R'