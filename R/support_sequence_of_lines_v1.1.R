##name of script: support_sequence_of_lines.R
cat("version_number= ",v_nr,"\n") 
##purpose: supporting scripts to program 'sequence_of_lines.R'
## GNU General Public License (GPL)

##contents:
# 1.digitize and plot new center of object 
# 2.plot of pixel cluster (PC) which represents a line-segment
# 3.plot of all pixel clusters (PC) on orthoimage in small scale
# 4.plot of course of a line (r_dist)
# 5.calculation of ro-value from image coordinates (x, y) and 
#   line orientation (theta_appr) 
# 6.plot of single line pixel cluster (PC) representing a line 
#   segment in small or large scale
# 7.determination of scale
# 8.calculation of center of segment
# 9.correction of line-midpoint-position and calculation of angle
# 10.estimation of object-type (cas)
################################################################################

## 1.digitize and plot new center of object
coo <- new_centre()
#end of script 1
#######################################################################


## 2.plot of pixel cluster (PC) which represents a line segment

plot_PC <- function() { 
  n <- 1 #number of PC (to be changed)
  setwd(home_dir)
  fname=paste("./data/",Img_name,"/b",bnr2,"_",n,".txt", sep="")
  P <- read.table(fname, col.names=c("idx","x","y")) #point cloud
  P_red <- reduce_pointset(P)
  points(P_red[,2] - orig_x, P_red[,3] - orig_y, pch=20, asp=3, cex=0.2, col="blue")
} #end plot_PC
#

plot_PC() #call of function

#end of script 2
################################################################################


## 3.plot of all pixel clusters (PC) on orthoimage in small scale

#function
plot_PC_all <- function() { 

#input of references
  
  if (Img_name == "ISPRS7") {
    OrgImgPathname <- paste(home_dir,"/data",sep = "")
    OrgImgFilename <- "top_mosaic_09cm_area7.tif"
  }
  
  if (Img_name == "ISPRS1") {
    OrgImgPathname <- paste(home_dir,"/data",sep = "")
    OrgImgFilename <- "top_mosaic_09cm_area1.tif"
  }
  
  setwd(OrgImgPathname)
  img_ref <- readImage(OrgImgFilename)
  display(img_ref, method="raster")
  
  for (i in lnr_det3) {
    lnr <- i
    cat("lnr=",lnr,"\n")
    setwd(home_dir)
    fname <- paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
    pc <- read.table(fname, header=TRUE)
    names(pc) <- c("nr","col","row")
    #browser()
    points(pc$col, pc$row, pch=20, asp=1, cex=0.05, col="white")
  } #end for-loop
} #end of function 'plot_PC_all'

plot_PC_all() #call of function

#end of script 3
################################################################################


## 4. plot of course of a line (r_dist)

# histograms of line length (n_pixel)

#function
hist_lin_len <- function() { 
# in ref direction
  B5_4d_ord
  head(B5_4d_ord)
  max(B5_4d_ord$n_pixel)
  min(B5_4d_ord$n_pixel)
  dif_n_pixel <- max(B5_4d_ord$n_pixel) - min(B5_4d_ord$n_pixel)
  hist(B5_4d_ord$n_pixel,nclass=dif_n_pixel)

# in 'orthogonal to ref' direction
  B5_4dd_ord
  head(B5_4dd_ord)
  max(B5_4dd_ord$n_pixel)
  min(B5_4dd_ord$n_pixel)
  dif_n_pixel <- max(B5_4dd_ord$n_pixel) - min(B5_4dd_ord$n_pixel)
  hist(B5_4dd_ord$n_pixel,nclass=dif_n_pixel)
} #end of function 'hist_lin_len'

hist_lin_len() #call of function

#end of script 4
################################################################################

## 5.calculation of ro-value from image coordinates (x,y) 
#    and line-orientation (theta_appr)

ro_from_xy <- function() { 
  theta_appr <- 125 #manual input
  theta_appr_arc <- theta_appr/omega
  X <- 1037.3 #manual input, must be changed
  Y <- 1020.2 #manual input, must be changed
  ro <- cos(theta_appr_arc) * X + sin(theta_appr_arc) * Y
} #end of function 'ro_from_xy'

ro_from_xy #call of function

#end of script 5
################################################################################


## 6.plot of pixel cluster (PC) representing a line segment 
#    in small or large scale

Img_name

#function
plot_PC_2scales <- function() {
  setwd(home_dir)
  lnr <- 1 #manual input of line number
  #browser() #with interaction?
  cat("lnr= ",lnr,"\n")
  H <- array(dim=c(n_theta,n_ro)) # init of Hough matrix
  i <- 0
  while (i < n_theta) {
    i <- i + 1
    H[,] <- 0
  } #end H null-stellung
  
  #generation of H
  ratio <- (n_ro-1)/(ro[n_ro]-ro_1)
  X <- pc3$col
  Y <- pc3$row
  n_X <- length(X)
  P <- matrix(nrow=n_X, ncol=3)
  P[,] <- 0
  
  #loop for all points of selected pixel cluster(PC)
  k1 <- 0
  k3 <- 1
  
  while (k1 < n_X){
    k1 <- k1+1
    i <- 0
    while (i < n_theta) {
      i <- i+1
      PC_segment_4(1) #calculation of Point Cloud of line segment #1
    } #end loop i
  } #end loop k1
  
  n_P <- k3-1
  cat("lnr=",lnr,"n_P=",n_P,"\n")
  setwd(home_dir)
  f <- paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
  write.table(P[1:n_P,],file=f, sep="   ")
  #
  
  #plot point cloud
  #use large scale or small scale presentation
  
  #small scale:
  #points(P[,2],-P[,3], pch=".", asp=1, cex=3.0, col="red") #switch to 'Plots' to see plot (graphics)
  
  #large scale:
  points(P[,2]-orig_x,P[,3]-orig_y,pch=".",asp=1,cex=4.0,col="red") #switch to 'Plots' (enlarged ortho)
  #end plot of line segment
  
} #end of function 'plot_PC_2scales'

plot_PC_2scales() #call of function

#end of script 6
################################################################################


##determination of scale

#function
det_scale <- function() {
  r_max3 <- round(sqrt((2*r_max)^2+(2*r_max)^2))
  co <- locator(2)
  dx <- co$x[1]-co$x[2]
  dy <- co$y[1]-co$y[2]
  s1 <- sqrt(dx^2+dy^2)
  k <- s/r_max3
  1/k
} #end of function 'det_scale'

det_scale() #call of function

#end of script 7
################################################################################

## 8.calculation of center of segment
coo2 <- new_centre_auto() #call of function (contained in: func_loadLib_jh.R)
#end of script 8
################################################################################

## 9.correction of midpoint position and calculation of angle

answ <- readline("Is the position of all midpoints correct? type N if it is not:  ")

if (answ == "N") {
  midpoints
  n_RepPoint <- 8 #row number of midpoint to be corrected (PC_4), must be changed
  n_RepPoint <- as.integer(n_RepPoint)
  b13_angle_df$nr_center <- midpoints[,1]
  b13_angle_df
  all_PC
  r_dist <- dist_v2(n_RepPoint, b13_angle_df, all_PC) #call of function 
  #         function is contained in: func_loadLib_jh.R
  r_dist <- round(r_dist)
  np_r <- max(r_dist)
  
#plot of course of line
  x <- 100
  y <- 100
  plot(x,y,pch = 3,col="red",cex=0.8,xlim=c(0,np_r),ylim=c(np_r,0),xlab="i",
       ylab="r",axes=T,frame.plot=T,main="course of line")
  
  #loop
  i <- 1
  
  while (i < np) {
    x <- i
    y <- r_dist[i]
    points(x, y, pch = 20, col="red", cex = 0.8)
    i <- i+1
  }
  
  i <- 150 #r-value, determined from graph (manual operation)
  x <- i
  y <- r_dist[i]
  points(x,y,pch=20,col="cyan",cex=1.8)
  np <- nrow(all_PC[[n_RepPoint]])
  all_PC2 <- all_PC
  all_PC2[[n_RepPoint]]$dist <- round(r_dist)
  i <- 150 #determined from graph
  r_dist[i]
  x_centre <- all_PC2[[n_RepPoint]]$x[i]
  y_centre <- all_PC2[[n_RepPoint]]$y[i]
  
  #plot
  points(x_centre,-y_centre,pch=20,col="green",cex=1.5)
  b13_angle_df$x_centre[n_RepPoint] <- x_centre
  b13_angle_df$y_centre[n_RepPoint] <- y_centre
  b13_angle_df
}  #end correction of mid_point-position of line segment  
#

#calculation of angle (center of object to new center of segment)

#center of object/building
xc <- plotPar[1]
yc <- plotPar[2]
x_centre <- b13_angle_df2[7,3] #to be transferred to spObj_sequence_of_lines_v1.1.R
y_centre <- b13_angle_df2[7,4] #to be transferred to spObj_sequence_of_lines_v1.1.R

#correction of angle for new midpoint
alpha <- det_of_angle(x_centre, y_centre) #call of function
b13_angle_df2$alpha[7] <- alpha #correction
b13_angle_df2

#end of script 9 (correction of angle for midpoint of segment)
################################################################################

# end of supplementing scripts for program 'sequence of lines'


