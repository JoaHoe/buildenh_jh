## name of script: support_line_detection.R
cat("version_number= ",v_nr,"\n")
##purpose: supporting software for special cases
##ISPRS data (areas #1,#7), DT classification, training by ISPRS orthoimages #26
##instruction: please activate function 'locator' where it is needed
##author: Joachim Höhle
##GNU General Public License (GPL)

## contents:
# 1.improving the ro-range in Hough space
# 2.manual detection of line orientation
# 3.search of lines (lnr) with theta_index and plot of line
# 4.determination of approximate line orientation (theta_appr) from angle (alpha_math) 
# 5.search of lines with 'theta_appr + 90?' and 'theta_appr - 90?' 
# 6.interactive detection of rectangular lines by measurement of one pixel in enlarged image
# 7.estimation of object type (cas)
# 8.determination of line number (lnr) by theta_index and ro_index
# 9.calculation of ro_ind using theta_index and measured point
# 10.plot of detected non-ortho-line onto enlarged orthoimage
################################################################################

## 1.improving the ro-range in Hough-space

# solution after a first calculation with ro_rg = 1
theta_appr_img <- (B[1,1] - 1) * 5 #for theta_step=5 [degrees]
theta1 <- (-theta_appr_img) + 180
theta1_arc <- theta1/omega
theta2 <- theta1 - 90
theta2_arc <- theta2 / omega
d_safety1=150 #ro_min must be negative (manual change)
d_safety2=150
X <- max(pc2$col); Y <- (-max(pc2$row)) 
ro1_max <- cos(theta1_arc) * X + sin(theta1_arc) * Y
ro1_max <- abs(ro1_max)
ro2_max <- cos(theta2_arc) * X + sin(theta2_arc) * Y
ro2_max <- abs(ro2_max)
X <- min(pc2$col); Y <- (-min(pc2$row))
ro1_min <- cos(theta1_arc) * X + sin(theta1_arc) * Y
ro1_min <- abs(ro1_min )
ro2_min <- cos(theta2_arc) * X + sin(theta2_arc) * Y
ro2_min <- abs(ro2_min)
ro_range <- c(ro1_max, ro2_max, ro1_min, ro2_min) #ro1_max, ro1_min are 
#negative when alpha>90 degrees -> check ro-range!

max_ro <- as.integer(max(ro_range))
min_ro <- as.integer(min(ro_range))
#
min_ro2 <- as.integer(min_ro-d_safety1)
max_ro2 <- as.integer(max_ro+d_safety2)
ro <- seq(min_ro2,max_ro2,by=ro_step)
n_ro <- length(ro)
ro_1 <- ro[1]
ro[n_ro]
ro_rg <- 3

#storage of Hough parameters
save(theta_step, ro_step, ro, ro_1, n_theta, n_ro, ro_rg, file=paste("./data/",Img_name,"/H_par",sep=""))
#continue at line 232 in 'line_detection.R'

#end of script 1
######################################################################################################

## 2.manual determination of line-orientation 
#in orthoimage (large scale) 
#determination of theta_ind
#results must be > 0

dir_meas <- locator(2)
x_ang <- (dir_meas$y[1] - dir_meas$y[2]) / (dir_meas$x[1] - dir_meas$x[2])
alpha_meas <- atan(x_ang) * omega
alpha_math <- (-alpha_meas) #change to math-system
theta_math <- alpha_math + 90
theta_img <- (-theta_math) #change to img-system

if(theta_img < 0) {
  theta_img <- theta_img + 180
}

theta_ind <- round(theta_img/theta_step) + 1
theta_ind # theta_index

cat("theta_ind=", theta_ind, "\n")

if (theta_ind < 90/theta_step) {
  alph_ind <- theta_ind + 90/theta_step
} else {
  alph_ind <- theta_ind - 90/theta_step
}

cat("alph_ind=", alph_ind, "\n")

#end of script 2
################################################################################################

## 3.search of lines (lnr) with theta_index and plot of line
theta_ind <- readline("type theta_index= ")
theta_ind <- as.integer(theta_ind)
thr_line_seg = n_pix
ct=0
vec <- 1 : length(B2[,1])
k=1.64 #k=approximate scale factor 

for (i in vec) {
  
  if (B2[i,1] == theta_ind && B2[i,3]/k >= thr_line_seg) { 
    cat("i= ", i,"\n")
    print(B2[i, ])
    ct <- ct + 1
  }
  
}#end search of lines with theta_index

cat("counts= ",ct,"\n")  

##search of lines (lnr) with alph_index
theta_ind

if (theta_ind < 90/theta_step) {
  alph_ind <- theta_ind + 90/theta_step
} else {
  alph_ind <- theta_ind - 90/theta_step
}

cat("alph_ind=", alph_ind, "\n")
vec <- 1 : length(B2[,1])
ct=0
k

for (i in vec) {

  if (B2[i,1] == alph_ind && (B2[i,3]/k) >= thr_line_seg) { #k=approximate scale factor
    cat("i= ", i,"\n")
    print(B2[i,])
    ct <- ct + 1
  }

} #end search of lines with alph_index

cat("counts= ",ct,"\n")
#end search of lines with alph_index

##plot of selected line
lnr <- readline("type line number: ") 
lnr <- as.integer(lnr)

#loop
  L_new  <- PC_segment_4(lnr)  
  P <- L_new[[1]]
  n_P <- L_new[[2]]
  P <- P[1:n_P,]
  P <- as.data.frame(P)
  names(P) <- c("idx","x","y")
  P_red <- reduce_pointset(P) 
  head(P_red)
  x_m <- mean(P_red[,2])
  y_m <- mean(P_red[,3]) 
  points(P[,2]-orig_x,(P[,3]-orig_y), pch=".", asp=1, cex=2.0, col="red") #see 'Plots' (plot))
  points(P_red[,2]-orig_x,(P_red[,3]-orig_y), pch=".", asp=1, cex=2.0, col="black") #see 'Plots' (plot)
  points(x_m-orig_x, y_m-orig_y, pch=16, asp=1, cex=2.0, col="blue")

#plot of selected line onto image extract 
theta_math <- 180 - B4$theta_angle[lnr]

if (theta_math >= 180) {
  theta_math <- theta_math - 180  
}

cat("theta_math= ", theta_math,"\n")
a <- -1/tan(theta_math/omega)
cat("a=",a,"\n")
x <- x_m
y <- y_m 
y_math <- (-y_m) #change to math_system
p2 <- round(x*cos(theta_math/omega) + y_math*sin(theta_math/omega))
b <- round(p2/sin(theta_math/omega))
cat("b= ", b, "\n")
coef = c(b,a)

#calculation of intercept (b2) for image extract
orig_y_math <- (-orig_y) #change to math-system
b_math <- b
y1 <- a * orig_x + b_math
b2_math <- y1 - orig_y_math
cat("b2_math=", b2_math, "\n")

#change to image-system
b2_img <- round(-b2_math)
a_img <- (-a)
coef2 <- c(b2_img,a_img)

if (is.finite(a_img) && is.finite(b2_img)) {
  abline(coef2, col="blue", lty=1, lwd=3, asp=1)
} else {
  ro_l1 <- B4$ro_pixel[lnr]
  ro_l3 <- round(ro_l1 - orig_x)
  lines(c(ro_l3,ro_l3),c(0,wind_y-orig_y),col="blue",lty=1,lwd=3,asp=1)
} #end if-else

#end of script 3
################################################################################################

## 4.determination of approximate line orientation (theta_appr) from angle (alpha_math)

alpha_img <- (-alpha_math)
theta_appr <- alpha_img+90
theta_appr_index <- round(theta_appr/5)+1
theta_appr_index

#end of script 4
###########################################################################################################

##5.search of lines with 'theta_appr + 90' and 'theta_appr - 90'

B_red <- subset(B,B[,3] >= 40) #about 15 pixel
x1 <- nrow(B_red)
vec <- 1 : x1

for (i in vec) {
  
  if (B_red[i,1] == theta_appr_index + 18 
      || B_red[i,1] == theta_appr_index-18 
      || B_red[i,1] == theta_appr_index) {
    cat("i=",i,"B=",B[i,],"\n")
  } #end if
  
} #end for-loop 

#note: there should exist at least two lines with theta_appr_index and theta_appr_index+-18

#end of script 5
###########################################################################################################

## 6.interactive detection of ortho-lines 
# by measurement of one pixel in enlarged orthoimage 

#display enlarged ortho_image and PC of building outline
img_uds <- img_ref[orig_x : wind_x,orig_y:wind_y,1:3]
display(img_uds, method = "raster")
#display(img_uds,method = "browser") #display enables zooming
points(xc-orig_x,yc-orig_y,pch=3, asp=1, cex=1.3, col="red")
points(as.integer(pc3$col-orig_x), as.integer(pc3$row-orig_y), 
       pch=20, asp=1, cex=0.3, col="green")

#determination of transformation parameter by means of orthoimage (large scale)
#measure two control 2 points (left lower, right upper) and one checkpoint (middle)
L1 <- trans_ortho() #

#transformation parameter
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

# measurement of new points (results: x,y)
locator2() #measurement and marking of a pixel's position

#determination of ortholines 
detect_meas1() 

#plot of detected line into enlarged orthoimage
B5_4_ord #use of ref-line (lnr_ref)
i=24 #index in B5_4_ord (value for i has to be changed)
B5_4_ord[i,]
cat("PC_nr=", B5_4_ord$lnr[i], "\n")
y <- (-y) #adapt to math_system #img_system

#coordinates
x
y

#angles
theta_img <- B5_4_ord$theta_angle[i]
theta_math <- 180 - theta_img #change to math-system
cat("theta_math= ", theta_math, "\n")
theta_math_arc <- theta_math/omega
a = (-1/tan(theta_math_arc))
cat("a= ",a,"\n")
B5_4_ord$ro_pixel[i]

p2 <- round(x*cos(theta_math_arc) + y*sin(theta_math_arc)) #sign may be '+' or '-'
b <- p2/sin(theta_math_arc)
orig_y <- (-orig_y) #change to math-system

#calculation by intercept for image extract (math_system)
b2 <- a*orig_x + b - orig_y #original
#b2 <- +a*orig_x + b + orig_y 

#change of parameter to image_system
b2_img <- (-b2)
a_img <- (-a)

coef2 <- c(b2_img,a_img)

# plot
if (is.finite(a)) {
  abline(coef2, col="blue", lty=1, lwd=2, asp=1)
}  else {
  ro_l1 <- B5_4_ord$ro_pixel[i]
  ro_l2 <- ro_l1 + ro_1
  ro_l3 <- round(ro_l2 - orig_x)
  lines(c(ro_l3,ro_l3),c(0, (wind_y - orig_y)),col="red")
} #end if-else

#end of script 6
#########################################################################################################

## 7.estimation of object type (cas)

n_longest_lines <- 10 #number of longest lines after Hough trans (default)
x1 <- B4$theta_index[1 : n_longest_lines]
ce <- matrix(nrow=10, ncol=2)
ce_df <- data.frame(ce)
names(ce_df) <- c("lnr","counts")
ce_df$counts <- 0
#

y1 <- B4$theta_index[1:n_longest_lines]
n1 <- length(y1)
y2 <- rep(0,n1)
y1 <- y1[order(y1,decreasing = F)]
y2[1] <- y1[1]
k <- 2
i <- 2

#loop
while (i <= n1)  {
  
  if (y1[i] != y2[k-1]) {
    y2[k] <- y1[i]
    k <- k + 1
  } #end if
  
  i <- i + 1
} #end while-loop

y2 <- subset(y2,y2>0)
y2
n_y2 <- length(y2)
ce_df <- ce_df[1 : n_y2,]
ce_df$counts
k=1

#loop
for (j in y2) {
  for (i in x1) {
    if (j == i) {
      ce_df$counts[k] <- ce_df$counts[k] + 1
      ce_df$lnr[k] <- j
    } #end if
  } #end loop i
  k <- k+1
} #end loop j
ce_df
ce_df$lnr

#derivation of parameter (ces) for number of orthogonal lines
ces <- 0 #number of orthogonal lines in the first 10 lines at B4 
i <- 1

#loop
while (i <= length(ce_df$lnr)) {
  if (ce_df$lnr[i] == theta_ref_ind || ce_df$lnr[i] == alph_ref_ind) {
    ces <- ces + ce_df$counts[i]
  } #end if
  i <- i+1
} #end while-loop

n_ortholines_1 <- ces
n_nonortholines <- n_longest_lines - n_ortholines_1 #number of non-ortholines at first 10 lines
max_pix <- B4$n_pixel[n_longest_lines] #size of 8th PC

#parameters for estimating of object-type:
cat("max_pix= ", max_pix,"\n") 
cat("n_ortholines_1= ",n_ortholines_1,"\n")
cat("n_nonortholines= ", n_nonortholines,"\n") #number of non-ortholines 
#

##determination of object-type

ty <- fe(max_pix,n_ortholines_1,n_nonortholines) #ty #value is an estimate -> subject of change
cat("suggestion for object type=",ty, "\n")
cat("object types (cas): 1 (extr_wd), 2 (4_long), 3 (100_all), 4 (100_all+nonortho)","\n")

if (ty == 0) {
  cat("suggestion is: 4", "\n") 
}

ty <- readline("type object type= ") #manual input
ty <- as.integer(ty)
cas <- switch(ty,"extr_wd", "4_long", "100_all", "100_all+nonortho")
cat("case= ", cas, "\n")
#end of script 7 
#######################################################################

## 8: determination of line number (lnr) 
#by means of theta_index and ro_index
#by orthoimage (large scale)

img_uds <- img_ref[orig_x : wind_x,orig_y:wind_y,1:3]
display(img_uds, method = "raster")
#display(img_uds,method = "browser") #display enables zooming
points(xc-orig_x,yc-orig_y,pch=3, asp=1, cex=1.3, col="red")
points(as.integer(pc3$col-orig_x), as.integer(pc3$row-orig_y), 
       pch=20, asp=1, cex=0.3, col="green")

#determination of transformation parameter by means of orthoimage (large scale)
#measure two control 2 points (left lower, right upper) and one checkpoint (middle)
L1 <- trans_ortho() #

#transformation parameter
D <- matrix(nrow=2, ncol=2)
D[1,1] <- L1[[1]][1,1]
D[1,2] <- L1[[1]][1,2]
D[2,1] <- L1[[1]][2,1]
D[2,2] <- L1[[1]][2,2]
a0 <- L1[[2]][1]
b0 <- L1[[2]][2]
tr_lat <- c(a0,b0)
kf2 <- L1[[3]]

#measurement of two points to determine theta and ro and line number
c10 <- locator(2) #standard function
c10

pts10 <- c(c10$x[1],c10$y[1],c10$x[2],c10$y[2])
pts10
loc1 <- c(c10$x[1],c10$y[1]) #first point
loc2 <- c(c10$x[2],c10$y[2]) #second point
print(tr_lat)
print(D)
#transformation to image-system (small scale)
pts11 <- tr_lat + D%*%loc1 #first point
pts12 <- tr_lat + D%*%loc2 #second point

x11 <- pts11[1,1] #x of first point
y11 <- pts11[2,1] #y of first point
x <<- x11
y <<- y11
#coordinates of first point, small scale
cat("x_coordinate= ",x,sep = "","\n") 
cat("y_coordinate= ",y,sep = "","\n")

#plot of first point in uds
points(x-orig_x,y-orig_y,pch=3, asp =1, cex=1,asp=1, col="red") #large scale

#second point
x12 <- pts12[1,1]
y12 <- pts12[2,1]
x <<- x12
y <<- y12
cat("x_coordinate= ",x,sep = "","\n") #small scale
cat("y_coordinate= ",y,sep = "","\n") #small scale
#plot of second point in uds
points(x-orig_x,y-orig_y,pch=3, asp =1, cex=1,asp=1, col="red") #large scale
#mean between the two points
x <- (x11+x12)/2
y <- (y11+y12)/2
#plot of mean point in uds 
points(x-orig_x,y-orig_y,pch=3, asp =1, cex=1,asp=1, col="blue") #large scale

#
x_ang <- (y11 - y12)/(x11 - x12)

alpha_meas <- atan(x_ang) * omega
alpha_math <- (-alpha_meas) #change to math-system
theta_math <- alpha_math + 90
theta_img <- 180 - theta_math
theta_ind <- round(theta_img/5) + 1
theta_ind
x
y
theta_img = (theta_ind-1)*5 
theta_math = 180 - theta_img
theta_math_arc=theta_math/omega
y <- -y #change to math_system
ro_math <- round(x*cos(theta_math_arc) + y*sin(theta_math_arc))
ro_ind <- abs(round(ro_math/5)) + 1
cat("ro_index= ", ro_ind, "\n")
b <- ro_math/sin(theta_math_arc)
#
b_img <- (-b) #change to img_system
a = (-1/tan(theta_math_arc))
cat("a= ",a,"\n")
a_img <- (-a) #change to img_system
coef <- c(b_img,a_img)

#calculation by intercept for image extract (math_system)
orig_y_math <- -orig_y
b2 <- a * orig_x + b - orig_y_math 

#change of parameters to image_system
b2_img <- (-b2)
a_img <- (-a)
coef2 <- c(b2_img,a_img)

# plot
#lnr=8 #to be determined in support_line_detection #9

if (is.finite(a_img)) {
  abline(coef2, col="blue", lty=1, lwd=2, asp=1)
}  else {
  ro_l1 <- B2$ro_pixel[lnr]
  ro_l2 <- ro_l1 + ro_1
  ro_l3 <- round(ro_l2 - orig_x)
  lines(c(ro_l3,ro_l3),c(0,(wind_y - orig_y)),col="red")
} #end if-else

#search of line numbers (lnr) with theta_index and ro_ind 
theta_ind <- readline("type theta_index= ")
theta_ind <- as.integer(theta_ind)
ro_ind <- readline("type ro_index= ") #determined evtl. from roof ridge 
ro_ind <- as.integer(ro_ind)
#

head(B2)
vec <- 1 : length(B2[,1])

for (i in vec) {
  if (B2[i,2] == theta_ind && B2[i,3] == ro_ind) {
    print(B2[i, ])
  }
  
} #end search of line numbers with theta_index & ro_index
#

head(B0)
vec <- 1 : length(B0[,1])

for (i in vec) {
  
  if (B0[i,2] == theta_ind && B0[i,3] == ro_ind) {
    print(B0[i, ])
  }
  
} #end search of line numbers with theta_index & ro_index
###############################################################################

## 9:calculation of ro_ind using theta_index and one measured point

#input
theta_ind=15 #change value
x=1261.0 #point (mean, img_system)
y=2124.6 #point (mean, img_system) check!

theta_img = (theta_ind-1)*5 
theta_math = 180 - theta_img
theta_math_arc=theta_math/omega
y <- -y #change to math_system

ro_math <- round(x*cos(theta_math_arc) + y*sin(theta_math_arc))
ro_ind <- abs(round(ro_math/5)) + 1
cat("ro_index= ", ro_ind, "\n")

#end of 9.
###############################################################################

##10.plot of detected non-ortho-line onto enlarged orthoimage
#example ISPRS1, b18, line-index (i) in B: 24

#display enlarged ortho_image and plot of PC of building outline
img_uds <- img_ref[orig_x:wind_x,orig_y:wind_y,1:3]
display(img_uds, method = "raster")
points(xc-orig_x,yc-orig_y,pch=3, asp=1, cex=1.3, col="red")
points(as.integer(pc3$col-orig_x), as.integer(pc3$row-orig_y), 
       pch=20, asp=1, cex=0.3, col="green")

B #all lines after Hought trans
i=24 #index in B (value for i must be adapted)
i=as.integer(i)
B[i,]
cat("PC_nr=", i, "\n")
theta_img <- (B[i,1] - 1) * theta_step

#
lnr <- i #number of point cloud (PC)
PC_seg_P_nP <- PC_segment_4(lnr)
P <- PC_seg_P_nP[[1]] 
n_P <- PC_seg_P_nP[[2]]
setwd(home_dir)
f<-paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
P <- P[1:n_P,]
head(P)
write.table(P,file=f, sep="   ")
#
P <- read.table(f) #point cloud of first line
head(P)
names(P) <- c("idx","x","y")
nrow1 <- nrow(P)
points(P[,2]-orig_x,P[,3]-orig_y, pch=20, asp=1, cex=0.3, col="yellow") #point 
P_red <- reduce_pointset(P)
points(P_red[,2]-orig_x,P_red[,3]-orig_y, pch='.', asp=1, cex=0.3, col="blue")
x_m <- mean(P_red[,2])
y_m <- mean(P_red[,3]) 
y_m_math <- (-y_m) #adapt to math-system
x <- round(x_m)
y <- round(y_m)
points(x-orig_x,y-orig_y, pch=20, asp=1, cex=1.5, col="blue") #point 

#angles
theta_math <- 180 - theta_img #change to math-system
cat("theta_math= ", theta_math, "\n")
theta_math_arc <- theta_math/omega
a = (-1/tan(theta_math_arc))
cat("a= ",a,"\n")
ro_pixel <- (B[i,2]-1)*ro_step + ro_1 

p2 <- round(x*cos(theta_math_arc) + y_m_math*sin(theta_math_arc)) #sign of p2 can be + or -
b <- p2/sin(theta_math_arc)
orig_y_math <- (-orig_y) #change to math-system

#calculation by intercept for image-extract (math_system)
b2 <- a * orig_x + b - orig_y_math 

#change of parameter to image_system
b2_img <- (-b2)
a_img <- (-a)
coef2 <- c(b2_img,a_img)

#plot of line
if (is.finite(a)) {
  abline(coef2, col="red", lty=1, lwd=2, asp=1)
}  else {
  ro_l1 <- B5_4_ord$ro_pixel[i]
  ro_l2 <- ro_l1 + ro_1
  ro_l3 <- round(ro_l2 - orig_x)
  lines(c(ro_l3,ro_l3),c(0, (wind_y-orig_y)),col="red")
} #end if-else

#end 10.
#############################################################################################################


##end of support_line_detection.R
################################################################################