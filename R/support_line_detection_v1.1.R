## name of script: support_line_detection.R
cat("version_number= ",v_nr,"\n")
## purpose: supporting software for special cases
## instruction: please activate function 'locator' where it is needed
## GNU General Public License (GPL)

## contents:
# 1.improving the ro-range in Hough space
# 2.manual detection of line orientation
# 3.search of lines (lnr) with theta_index
# 4.determination of approximate line orientation (theta_appr) from angle (alpha_math) 
# 5.search of lines with 'theta_appr + 90?' and 'theta_appr - 90?' 
# 6.interactive detection of rectangular lines by measurement of one pixel in enlarged image 
# 7.estimation of object type (cas)
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
save(theta_step, ro_step, ro, ro_1, n_theta, n_ro, ro_rg, file="./data/ISPRS7/H_par")

#end of script 1
#continue at line 193 in 'line_detection.R'

######################################################################################################

## 2.manual determination of line orientation in orthoimage (large scale) 
#determination of theta_ind
#result must be > 0

dir_meas <- locator(2)
x_ang <- (dir_meas$y[1] - dir_meas$y[2]) / (dir_meas$x[1] - dir_meas$x[2])
alpha_meas <- atan(x_ang) * omega
alpha_math <- (-alpha_meas) #change to math-system
theta_math <- alpha_math + 90
theta_img <- (-theta_math) #change to img-system

if(theta_img < 0) {
  theta_img <- theta_img + 180
}

theta_ind <- round(theta_img / 5) + 1
theta_ind # theta_ind > 0 ?

if (theta_ind < 0) {
  theta_ind <- theta_ind + 18
}
cat("theta_ind=", theta_ind, "\n")

alph_ind <- theta_ind + 90/5
cat("alph_ind=", alph_ind, "\n")

alph_ind2 <- theta_ind - 90/5 #theta_ind2 > 0 ?
cat("alph_ind2=", alph_ind2, "\n")
#end of script 2
################################################################################################

## 3.search of lines (lnr) with theta_index
theta_ind <- readline("type theta_index= ")
theta_ind <- as.integer(theta_ind)

vec <- 1 : length(B2[,1])
for (i in vec) {
  if (B2[i,2] == theta_ind && B2[i,4]/kf > n_pix) { #kf=scale factor 
    print(B2[i, ])
  }
} #end search of lines with theta_index

##search of lines (lnr) with alph_index
vec <- 1 : length(B2[,1])

for (i in vec) {
  if (B2[i,2] == alph_ind) {
    cat("lnr= ",i,"\n")
  }
} #end search of lines with alph_index

for (i in vec) {
  if (B2[i,2] == alph_ind2) {
    cat("lnr= ",i,"\n")
  }
} #end search of lines with alph_index2

#end of script 3
################################################################################################

## 4.determination of approximate line orientation (theta_appr) from angle (alpha_math)

alpha_img <- (-alpha_math)
theta_appr <- alpha_img+90
theta_appr_index <- round(theta_appr/5)+1
theta_appr_index
#end of script 4
###########################################################################################################

## 5.search of lines with 'theta_appr + 90?' and 'theta_appr - 90?

B_red <- subset(B,B[,3] >= 40) #about 15 pixel
x1 <- nrow(B_red)
vec <- 1 : x1

for (i in vec){
  if (B_red[i,1] == theta_appr_index+18 
      || B_red[i,1] == theta_appr_index-18 
      || B_red[i,1] == theta_appr_index) {
    cat("i=",i,"B=",B[i,],"\n")
  } #end if
} #end for-loop

#note: there should exist at least two lines with theta_appr_index and theta_appr_index+-18

#end of script 5 
###########################################################################################################


## 6.interactive detection of lines by measurement of one pixel in enlarged orthoimage 

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
#

# measurement of new points (results: x,y)
locator2() #measurement and marking of a pixel's position

#determination of lines by measuring one pixel
detect_meas1() 

#plot of detected line into enlarged orthoimage
B5_4_ord #use of ref-line (lnr_ref)
i=35 #index in B5_4_ord (value for i has to be changed)
B5_4_ord[i,]
cat("PC_nr=", B5_4_ord$lnr[i], "\n")
y <- (-y) #adapt to img_system

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
b2 <- a*orig_x + b - orig_y 

#change of parameter to image_system
b2_img <- (-b2)
a_img <- (-a)

coef2 <- c(b2_img,a_img)

# plot
if (is.finite(a)) {
  abline(coef2, col="red", lty=1, lwd=2, asp=1)
}  else {
  ro_l1 <- B5_4_ord$ro_pixel[i]
  ro_l2 <- ro_l1 + ro_1
  ro_l3 <- round(ro_l2 - orig_x)
  lines(c(ro_l3,ro_l3),c(0, (wind_y - orig_y)),col="red")
} #end if-else
#

##plot of detected nonortho-line onto enlarged orthoimage

B2 #all lines after Hought trans
i=136 #index in B2 (value for i has to be changed)
B2[i,]
cat("PC_nr=", B2[i,1], "\n")
theta_img <- (B2[i,2] - 1) * ro_step

#
lnr <- i #longest PC
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
#points(P[,2],-P[,3], pch=20, asp=1, cex=0.5, col="red") #point 
P_red <- reduce_pointset(P)
points(P_red[,2]-orig_x,P_red[,3]-orig_y, pch='.', asp=3, cex=3, col="blue")
x_m <- mean(P_red[,2])
y_m <- mean(P_red[,3]) 
#

y_m <- (-y_m) #adapt to math-system
x <- round(x_m)
y <- round(y_m)

#angles
theta_math <- 180 - theta_img #change to math-system
cat("theta_math= ", theta_math, "\n")
theta_math_arc <- theta_math/omega
a = (-1/tan(theta_math_arc))
cat("a= ",a,"\n")
ro_pixel <- (B2[i,3]-1)*5

p2 <- round(x*cos(theta_math_arc) + y*sin(theta_math_arc)) #sign may be '+' or '-'
b <- p2/sin(theta_math_arc)
orig_y <- (-orig_y) #change to math-system

#calculation by intercept for image extract (math_system)
b2 <- a * orig_x + b - orig_y 

#change of parameter to image_system
b2_img <- (-b2)
a_img <- (-a)
coef2 <- c(b2_img,a_img)

# plot
if (is.finite(a)) {
  abline(coef2, col="red", lty=1, lwd=2, asp=1)
}  else {
  ro_l1 <- B5_4_ord$ro_pixel[i]
  ro_l2 <- ro_l1 + ro_1
  ro_l3 <- round(ro_l2 - orig_x)
  lines(c(ro_l3,ro_l3),c(0, (wind_y - orig_y)),col="red")
} #end if-else

#end of script 6
#############################################################################################################

## 7.estimation of object type (cas)

x1 <- B4$theta_index[1:8]
ce <- matrix(nrow=5, ncol=2)
ce_df <- data.frame(ce)
names(ce_df) <- c("lnr","counts")
ce_df$counts <- 0
ce_df
#

n_longest_lines <- 8 #number of longest lines in Hough trans
y1 <- B4$theta_index[1:n_longest_lines]
n1 <- length(y1)
y2 <- rep(0,n1)
y1 <- y1[order(y1,decreasing = F)]

y2[1] <- y1[1]
k <- 2
i <- 2
while (i <= n1)  {
  if (y1[i] != y2[k-1]) {
    y2[k] <- y1[i]
    k <- k + 1
  } #end if
  i <- i + 1
} #end for-loop

y2 <- subset(y2,y2>0)

ce_df$counts
k=1
for (j in y2) {
  for (i in x1) {
    cat("j=",j,"i=", i,"\n")
    if (j == i) {
      ce_df$counts[k] <- ce_df$counts[k] + 1
      ce_df$lnr[k] <- j
    } #end if
  } #end loop i
  k <- k+1
} #end loop j

ce_df
ce_df$lnr
ces <- 0 #number of ortho-lines in the first 8 lines at B4 
#

i <- 1
while (i <= length(ce_df$lnr)) {
  
  if (ce_df$lnr[i] == theta_ref_ind || ce_df$lnr[i] == alph_ref_ind) {
    ces <- ces + ce_df$counts[i]
  } #end if
  
  i <- i+1
} #end loop
#

n_ortholines_1 <- ces
n_nonortholines <- n_longest_lines - n_ortholines_1
max_pix <- B4$n_pixel[n_longest_lines] #size of 8th PC
#
fe(max_pix,n_ortholines_1,n_nonortholines) #call of function
#

#end of script 7 (estimation of object type (cas))
#######################################################################


##end of support_line_detection.R
################################################################################