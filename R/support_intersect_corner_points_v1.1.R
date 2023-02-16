## name of script: support_intersect_corner_points.R
cat("version_number= ",v_nr,"\n")
## purpose: supporting software for script "intersect_corner_points.R"
##instruction: process first mode 'demo'
##author: Joachim Höhle
##GNU General Public License (GPL)

## contents:
# 1.test of theta average (theta_av)
# 2.interactive generation of theta average
# 3.automated solution for theta average
# 4.manual generation of theta average
# 5.intersection of two lines
# 6.calculation of line segment lengths
################################################################################


## 1.test of theta average
setwd(home_dir)
options(digits=7)
thr_theta_av2 <- 10 #threshold [degrees]
f1 <- paste("./data/",Img_name,"/theta_av_b",bnr2,".txt",sep="")
theta_av <- read.table(f1)
theta_av
#

f2 <- paste("./data/",Img_name,"/th_ref_",bnr2,sep="")
load(f2) #main direction
theta_ref #approximate value

if (abs(theta_av$x[1] - theta_ref) > thr_theta_av2 ) { #threshold theta_ref_appr
  stop("error in theta_av")
}
#end of script 1
################################################################################


## 2.interactive generation of theta average 
# weighted average of theta-angle (theta_av)

#function
w_av <- function(ang,len) { #contained in func_loadLib_jh.R
  x <- length(ang)
  y <- 1 : x
  pw <- len/max(len) #weights
  sp <- sum(pw)
  theta <- rep(0,x)
  
  for (i in y) {
    theta[i] <- (ang[i]*pw[i])/sp
  }
  
  theta_av <- sum(theta)
  return(theta_av)
} #end of function 'w_av(ang,len)'

#data: ortholines of b11 (ISPRS1)
ang_mod <- c(9.4953,5.3238,5.0244,4.8934,4.4975,5.2969)
len <-  c(33,244,179,84,100, 218)

#calculation
theta_average <- w_av(ang_mod,len) #call of function
cat("theta_average= ",theta_average,"degrees","\n")

#end of script 2
################################################################################


## 3.automated solution for theta average
#example b11 (ISPRS1)

setwd(home_dir)
options(digits=7)

##input of table with theta_adj, ro_adj
fname9 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt",sep="")
B6_seq <- read.table(fname9)
B6_seq 

B6_seq[,8] <- 0
names(B6_seq)[8] <- "ortho"
n_B6_seq <- nrow(B6_seq)
z <- 1 : n_B6_seq

#loop 
for (i in z) {
  
  if (B6_seq$theta_ang[i] == theta_ref || B6_seq$theta_ang[i] == alph_ref) {
    B6_seq$ortho[i] <- 1
  } else {
    B6_seq$ortho[i] <- 0
  }

} #end loop i

theta_vec <- rep(0,n_B6_seq)
np_vec <- rep(0,n_B6_seq)

for (i in z) {
  if (B6_seq$ortho[i] == 1) {
    theta_vec[i] <- B6_seq$theta_adj[i]
    np_vec[i] <- B6_seq$np[i]
  }
} #end for-loop

theta_vec_red <- subset(theta_vec, theta_vec != 0)
n_theta_main <- length(theta_vec_red)
vec1 <- 1 : n_theta_main
theta_vec_red2 <- theta_vec_red
# 

for (i in vec1) {
  
  if (theta_vec_red[i] > 90) {
    theta_vec_red2[i] <- theta_vec_red[i] - 90
  }
  
} #end loop i

ang1 <- theta_vec_red2
np_vec_red <- subset(np_vec, np_vec != 0)
len1 <- np_vec_red

#call of function
theta_average <- w_av(ang1,len1) #function contained in 'func_loadLib_jh.R'
cat("theta_average= ",theta_average,"degrees","\n")

##automated solution for theta_average2
B6_seq[,8] <- 0
names(B6_seq)[8] <- "ortho"
n_B6 <- nrow(B6_seq)
z <- 1 : n_B6

#loop
for (i in z) {
  
  if (B6_seq$theta_ang[i] == theta_ref || B6_seq$theta_ang[i] == alph_ref) {
    B6_seq$ortho[i] <- 1
  } else {
    B6_seq$ortho[i] <- 0
  }
  
} #end loop i

##automated solution for a second direction

if (bnr2 == 18) { #line number must be adopted
  theta_vec2 <- rep(0, n_B6_seq)
  np_vec2 <- rep(0,n_B6_seq)
  
  for (i in z) {
    
    if (B6_seq$ortho[i] == 0) {
      theta_vec2[i] <- B6_seq$theta_adj[i]
      np_vec2[i] <- B6_seq$np[i]
    } #end if
    
  } #end loop i
  
  theta_vec_red2 <- subset(theta_vec2, theta_vec2 > 0)
  n_theta_main2 <- length(theta_vec_red2)
  ang2 <- theta_vec_red2
  np_vec2_red <- subset(np_vec2, np_vec2 != 0)
  len2 <- np_vec2_red
  #
  
  #call of function
  theta_average2 <- w_av(ang2,len2) #function contained in 'func_loadLib_jh.R'

} #end if (bnr2 == 18)

cat("theta_average2= ",theta_average2,"degrees","\n")

#end of script 3
#####################################################################


## 4.manual generation of theta average

#input

#example: b36 (ISPRS1)
ang1 <- c(79.9345,81.9270,(171.2564-90)) #theta angles (must be adapted)
len1 <- c(176,56,28) #length of lines
theta_average <- w_av(ang1,len1)

ang1 <- c((161.7411-90),70.0149,72.7405) #theta angles (must be adapted)
len1 <- c(176,56,28) #length of lines
theta_average2 <- w_av(ang1,len1)
#end of script 4 (manual generation of theta average)
###############################################################################


## 5.intersection of two lines
#instruction: change values for other examples

#data
theta_1 <- 89.9684 #input of theta_angle of line 1
theta_1_arc <- theta_1/omega
ro_1 <- 904.94 # input of ro_distance of line 1
theta_2 <- 2.5080 #input of theta_angle of line 2
theta_2_arc <- theta_2/omega
ro_2 <- 673.54 #input of ro_distance of line 2
#

#calculation
N <- (sin(theta_2_arc) - tan(theta_2_arc) * cos(theta_1_arc)) * sin(theta_1_arc)
x <- ((ro_2 * sin(theta_1_arc) - ro_1 * sin(theta_2_arc)) * tan(theta_2_arc))/N
y <- (-1/tan(theta_1_arc)) * x + ro_1 / sin(theta_1_arc)
#

#end of script 5 (intersection of two lines) 
################################################################################


## 6.calculation of line-segment-lengths
#instruction: change paths

#call of function
setwd(home_dir2)
source("func_dist_PC.R") #function contained in 'func_loadLib_jh.R' 
y4 <- B7$lnr
n_pts <- length(y4)
distance <- matrix(nrow=n_pts, ncol=2)
distance[,] <- 0

#plot
plot(xc,-yc, pch=3, cex=2, col="red", asp=1, xlim=c(xc-r_max2,xc+r_max2), ylim=c(-yc-r_max2,
  -yc+r_max2), xlab="col", ylab="row",main=paste("building",bnr2))

#loop
i=0

for (n in y4) {
  setwd(home_dir)
  fname=paste("./data/",Img_name,"/b",bnr2,"_",n,".txt", sep="")#)
  P <- read.table(fname, col.names=c("idx","x","y")) #point cloud
  n_P <- length(P$idx)
  P_red <- reduce_pointset(P) #call of function
  n_P_red <- length(P_red$idx) 
  #browser()
  cat("point cluster nr=", n,"\n")
  points(P_red[,2],-P_red[,3], pch=20, asp=3, cex=1.0, col="red") #with reducing of pixels
  #points(P[,2],-P[,3], pch=20, asp=3, cex=1.0, col="blue") #original PC
  P_red2 <- P_red
  dist <- dist_PC(P_red2) #distance of line
  cat("distance of PC",n," :",dist, "\n")
  i <- i + 1
  distance[i,1:2] <- c(n,dist)
} #end of loop

#result
print(distance)
#

#end script 6 
################################################################################

## end of 'support_intersect_corner_points.R'
