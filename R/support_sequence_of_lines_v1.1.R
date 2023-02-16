##name of script: support_sequence_of_lines.R
cat("version_number= ",v_nr,"\n") 
##purpose: supporting scripts to program 'sequence_of_lines.R'
##author: Joachim Höhle
## GNU General Public License (GPL)

##contents:
# 1.digitize and plot center of line 
# 2.plot a pixel cluster (PC) which represents a line-segment
# 3.plot of all pixel clusters (PC) on orthoimage in small scale
# 4.histograms of line-length (n_pixel)
# 5.calculation of ro-value from image coordinates (x, y) and 
#   line orientation (theta_appr) 
# 6.plot of single line pixel cluster (PC) representing a line 
#   segment in small or large scale
# 7.determination of scale
# 8.calculation of new center of object
# 9.correction of line-midpoint-position and calculation of angle
################################################################################

## 1.digitize and plot center of line 
#display enlarged ortho_image and plot of PC of building outline
img_uds <- img_ref[orig_x:wind_x,orig_y:wind_y,1:3]
display(img_uds, method = "raster")
points(xc-orig_x,yc-orig_y,pch=3, asp=1, cex=1.3, col="red")
points(as.integer(pc3$col-orig_x), as.integer(pc3$row-orig_y), 
       pch=20, asp=1, cex=0.3, col="green")

##manual measurement of 3 checkpoints (lower, upper, middle)
L1 <- trans_ortho() 

#derive transformation-matrix
D <- matrix(nrow=2, ncol=2)
D[1,1] <- L1[[1]][1,1]
D[1,2] <- L1[[1]][1,2]
D[2,1] <- L1[[1]][2,1]
D[2,2] <- L1[[1]][2,2]
a0 <- L1[[2]][1]
b0 <- L1[[2]][2]
tr_lat <- c(a0,b0)
kf2 <- L1[[3]]
#

## measurement of one new pixel (center of line)
#results: x,y in image-system of orthoimage)
locator2() #measurement and marking of one pixel's position

#end of script 1.
################################################################################

## 2.plot a pixel cluster (PC) which represents a line segment

plot_PC <- function(PC_number) { 
  n <- PC_number #number of PC (to be changed)
  setwd(home_dir)
  fname=paste("./data/",Img_name,"/b",bnr2,"_",n,".txt", sep="")
  P <- read.table(fname, col.names=c("idx","x","y")) #point cloud
  P_red <- reduce_pointset(P)
  points(P_red[,2] - orig_x, P_red[,3] - orig_y, pch=20, asp=3, cex=0.2, col="red")
} #end plot_PC
#

PC_number
plot_PC(4) #call of function

#end of script 2
################################################################################

## 3.plot all pixel clusters (PC) onto orthoimage in small scale

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
    #plot
    #browser()
    points(pc$col, pc$row, pch=20, asp=1, cex=0.3, col="white")
  } #end for-loop
} #end of function 'plot_PC_all'

#data
Img_name
bnr2
lnr_det3

#plot
plot_PC_all() #call of function

#end of script 3
################################################################################

## 4. histograms of line-length (n_pixel)
#for object types "extr_wd", "4_long", "100_all" only 

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

## 5.calculation of ro-value from image coordinates (X,Y) 
#    and line-orientation (theta_appr)

ro_from_xy <- function() { 
  theta_appr_arc <- theta_appr/omega
  ro <- cos(theta_appr_arc) * X + sin(theta_appr_arc) * Y
} #end of function 'ro_from_xy'

#data
theta_appr <- 125
X <- 1037.3 
Y <- 1020.2  

ro <- ro_from_xy() #call of function
ro

#end of script 5
################################################################################


## 6.plot of pixel cluster (PC) representing a line segment 
#    in small or large scale

#function
plot_PC_2scales <- function(lnr) {
  #input of references
  
  cat("lnr= ",lnr,"\n")
  PC_seg_P_nP <- PC_segment_4(lnr) #call of function
  P <- PC_seg_P_nP[[1]]
  n_P <- PC_seg_P_nP[[2]]
  cat("n_P= ",n_P,"\n")
  #plot point cloud
  #use large scale or small-scale presentation
  setwd(OrgImgPathname)
  
  #small scale:
  img_ref <- readImage(OrgImgFilename)
  display(img_ref, method="raster")
  points(P[,2],P[,3], pch=".", asp=1, cex=1.0, col="green") #switch to 'Plots' to see plot (graphics)
  
  #large scale:
  img_uds <- img_ref[orig_x:wind_x,orig_y:wind_y,1:3]
  display(img_uds, method = "raster")
  points(P[,2]-orig_x,P[,3]-orig_y,pch=".",asp=1,cex=1.0,col="green") #switch to 'Plots' (enlarged ortho)
  #end plot of line segment
  
} #end of function 'plot_PC_2scales'

#data
Img_name
lnr <- 2 

#plot
PC_seg_P_nP <- plot_PC_2scales(lnr) #call of function

#end of script 6
################################################################################

##determination of scale

#function
det_scale <- function() {
  #size of window: 2.4*r_max * 2.4*r_max
  mar <- r_max #specified distance from image center
  #coordinates of check points
  x1 <- xc - mar 
  y1 <- yc + mar
  x2 <- xc + mar
  y2 <- yc - mar
  points((x1-orig_x),(y1-orig_y),pch=16,cex=1.5, col="white") #lower check point
  points((x2-orig_x),(y2-orig_y),pch=16,cex=1.5, col="white") #upper check point
  r_max3 <- sqrt((2*r_max)^2+(2*r_max)^2)
  co <- locator(2)
  dx <- co$x[1]-co$x[2]
  dy <- co$y[1]-co$y[2]
  s1 <- sqrt(dx^2+dy^2)
  k <- r_max3/s1
  cat("scale factor= ",k,"\n") #known distance/measured distance 
} #end of function 'det_scale'

det_scale() #call of function

#end of script 7
################################################################################

## 8.calculation of new center of object

#function: func_new_centre_auto.R 
new_centre_auto <- function() { #function is contained in 'func_loadLib_jh.R'
  b13_angle_df
  xc_new2 <- mean(b13_angle_df$x_centre)
  yc_new2 <- mean(b13_angle_df$y_centre)
  coo2 <- c(xc_new2, yc_new2)
  cat("new coordinates2:", "xc = ",  coo2[1], "yc = ", coo2[2], "\n")
  points(xc_new2, -yc_new2, asp=1, pch=3, cex=3,col="green") #plot in small scale image 
  return(coo2)
} #end of function new_centre_auto

#data
b13_angle_df

#function call
coo2 <- new_centre_auto()

#end of script 8
################################################################################

## 9.correction of line-midpoint position and calculation of angle
#to be used for object-type "100_all+nonortho" and sequence-method "bdr_follow"
#example b11 of ISPRS1

answ <- readline("Is the position of all midpoints correct? ")

if (answ == "N") {
  midpoints
  n_RepPoint <- 8 #row number (index) of midpoint to be corrected (in b13_angle_df), must be changed
  n_RepPoint <- as.integer(n_RepPoint)
  b13_angle_df$nr_center <- midpoints[,1]
  b13_angle_df
  PC_nr #check index of line in 'PC_nr'
  n_RepPoint2 <- 7 #index of midpoint to be corrected (in 'all_PC' and 'PC_nr')
  all_PC[[n_RepPoint2]]
  r_dist <- dist_v2(n_RepPoint2, b13_angle_df, all_PC) #call of function 
  #         function is contained in 'func_loadLib_jh.R'
  r_dist <- round(r_dist) 
  np_r <- max(r_dist) 
  
#plot of course of line
  x <- 100
  y <- 100
  plot(x,y,pch = 3,col="red",cex=0.8,xlim=c(0,np_r),ylim=c(np_r,0),xlab="i",
       ylab="r",axes=T,frame.plot=T,main="course of line")
  
  #loop
  i <- 1
  
  while (i < np2) {
    x <- i
    y <- r_dist[i]
    points(x, y, pch = 20, col="red", cex = 0.8)
    i <- i+1
  }
  
  x <- 100 #x-value, determined from graph (manual operation)
  i <- x
  y <- r_dist[i]
  points(x,y,pch=20,col="blue",cex=1.8)
  np <- nrow(all_PC[[n_RepPoint2]])
  all_PC2 <- all_PC
  all_PC2[[n_RepPoint2]]$dist <- round(r_dist)
  i #determined from graph
  r_dist[i]
  x_centre <- all_PC2[[n_RepPoint2]]$x[i]
  y_centre <- all_PC2[[n_RepPoint2]]$y[i]
  
  #plot
  n_RepPoint
  points(x_centre,-y_centre,pch=20,col="blue",cex=1.5) #change display to graph of object (building)
  b13_angle_df$x_centre[n_RepPoint] <- x_centre
  b13_angle_df$y_centre[n_RepPoint] <- y_centre
  b13_angle_df
}  #end correction of mid_point-position of line segment  
#

#calculation of angle (center of object to new center of segment)

#center of object/building
xc <- plotPar[1]
yc <- plotPar[2]
b13_angle_df3 <- b13_angle_df
# x_centre <- b13_angle_df[6,3] #to be transferred to spObj_sequence_of_lines_v1.1.R
# y_centre <- b13_angle_df[6,4] #to be transferred to spObj_sequence_of_lines_v1.1.R
x_centre <- b13_angle_df3[8,3] #to be transferred to spObj_sequence_of_lines_v1.1.R
y_centre <- b13_angle_df3[8,4] #to be transferred to spObj_sequence_of_lines_v1.1.R

#correction of angle for new midpoint
alpha <- det_of_angle(x_centre, y_centre) #call of function
b13_angle_df3$alpha[6] <- alpha #correction
b13_angle_df3

#end of script 9 
################################################################################

# end of supplementing scripts for program 'sequence of lines'


