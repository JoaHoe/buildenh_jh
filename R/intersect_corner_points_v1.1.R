#name of program (script): 'intersect_corner_points.R'
cat("version_number= ",v_nr,"\n")
##description of program:
#intersection of successive lines 
#coordinate-system (theta,ro)
#graphics for checking
#derivation of weighted average of main angle (theta)
#using generic design matrix
##author: Joachim Höhle
##GNU General Public License (GPL)
cat("##############################################","\n")

cat("start of program 'intersect_corner_points.R'","\n")
setwd(home_dir)

##input of parameters
options(digits=10)
f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
load(f1) #plot parameter
xc <- plotPar[1]
yc <- plotPar[2]
r_max <- plotPar[3]
#

##input adjusted line parameters (theta,ro)
fname9 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt", sep="")
B6 <- read.table(fname9)
B6
#

n_pts <- length(B6[,1])
PC_nr <- B6[,1]
n_pts <- length(PC_nr) #number of points (vertices)
lnr_seq <- PC_nr

#thresholds
thr_theta_av = 10.0 # [degrees], used for checking of theta-average (theta_av)

#test for odd number of points
if (n_pts %% 2 == 1) { 
  cat("number of points is odd", "\n")
}

## start of program
setwd(home_dir)

#establishing the sequence of line segments in new matrix (B7)
names(B6)[1] <- "lnr"
n_PC <- length(B6$lnr)
B7 <- B6
B7[,] <- 0
B7$lnr <- lnr_seq

#loop
i <- 1
while (i <= n_pts) {
  k=1
  while (k <= n_PC) {
    if(lnr_seq[i] == PC_nr[k]){
    B7[i,] <- B6[k,]
    } #end if
  k <- k + 1
  } #end k-loop
i <- i + 1
} #end i-loop

##plot at large scale 
B7_seq <- B7
y <- 1 : n_pts

plot(xc,-yc, pch=3, cex=2, col="red", asp=1, xlim=c(xc-r_max2, xc+r_max2), 
     ylim=c(-yc-r_max2, -yc+r_max2), xlab="col", ylab="row",
     main=paste("building",bnr2))

y4 <- B7_seq$lnr

#loop
for (n in y4) {
 fname=paste("./data/",Img_name,"/b",bnr2,"_",n,".txt", sep="")
 P_orig <- read.table(fname, col.names=c("idx","x","y")) #point cloud
 P_red <- reduce_pointset(P_orig)
 #browser()
 cat("point cluster nr=", n,"\n")
 #points(P_orig[,2],-P_orig[,3], pch=20, asp=3, cex=0.5, col="blue")
 points(P_red[,2],-P_red[,3], pch=20, asp=3, cex=0.5, col="red") #corrected pixel cluster (PC)
} #end of loop

## intersection of lines
B8 <- B7_seq
B8 #matrix with proper sequence
phi_deg <- (-B8$theta_adj)
phi <- (phi_deg/180)*pi
m <- length(phi_deg) #number of lines in object

#matrix elements for m lines
A <- matrix(nrow=2*m, ncol=m) #matrix
A[,] <- 0
A <- design_mat(m,phi)
x0 <- B8$ro_adj
b0 <- A %*% x0 #calculation of intersections
b0 #approximate coordinates

#correction of coordinates due to tan(90)
z <- 1 : m
for (i9 in z) {
  
  if (B8$theta_adj[i9] == 90) {
    p_pos <- "cor_corner_pts"
    setwd(home_dir2)
    source(paste("spObj_intersect_corner_points_v",v_nr,".R",sep = "")) 
  } 
  
} #end loop i9

b0

#output of approximate coordinates of corners (vertices)
setwd(home_dir)
fname12 <- paste("./data/",Img_name,"/b",bnr2,"_coord_appr.txt",sep="")
write.table(b0,fname12,sep="  ")

##plot of intersected points

#input
b0_1 <- read.table(fname12,col.names=("xy"))
fname <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv",sep="")
idxy <- read.table(fname)
n_r_1 <- length(idxy$x)

#loop
i=1

while(i <= n_r_1) {
  x <- idxy$x[i]
  y <- idxy$y[i]
  points(x, -y, pch=20, cex=0.7, col="blue", asp=1) 
  i <- i+1 
}

#end of plot of all points

##plot of corner coordinates
#separated in x and y

Points_x <- rep(0,n_pts+1)
Points_x[n_pts+1] <- b0_1[1,1]
k <- 1
i <- 1
n2 <- (n_pts*2-1)

while (i <= n2) {
  Points_x[k] <- b0_1[i,1]
  k <- k+1
  i <- i+2
} #end loop Points_x

Points_y <- rep(0,n_pts+1)
Points_y[n_pts+1] <- b0_1[2,1]
k <- 1
i <- 2
n3 <- n_pts*2
while(i <= n3) {
  Points_y[k] <- b0_1[i,1]
  k <- k+1
  i <- i+2
} #end while-loop
#

k1 <- n_pts

#loop
i <- 0
while(i < k1){
  i <- i + 1
  x <- Points_x[i]
  y <- Points_y[i]
  points(x,-y, pch=20, cex=1.5, col="blue", asp=1)
}

b_xy <- cbind(Points_x,Points_y)
dimnames(b_xy)[[2]] <- list("x","y")
cat("approximate coordinates of corners:","\n")
print(b_xy)

#output of approximate corner coordinates
f <- paste("./data/",Img_name,"/b",bnr2,"_xy_appr.txt", sep="")
write.table(b_xy,f)

#plot with line segments
k1 <- n_pts

#loop
i <- 0

while(i < k1) {
  i <- i + 1
  cat("ptnr=",i,"\n")
  #browser()
  x <- Points_x[i]
  y <- Points_y[i]
  points(x, y, pch=20, cex=2.0, col="green", asp=1)
  lines(b_xy,  col="red", asp=1, type="l", lwd=2, lty=1)
}

#end of plot of corner points (derived from intersection with adjusted lines)

##check of orientation angles of the lines (theta_angle)
#difference of angles (theta) at building corners (vertices)
#order of lnr_seq for B6

options(digits=5)
B6_seq <- B8
lnr_seq <- B8$lnr
n_pts <- length(B8$lnr)
n_B6_seq <- n_pts

#loop
x1 <- 1 : n_pts
B6_seq

#difference of angles
angle_dif <- rep(0,n_pts)
y1 <- 1 : (n_pts-1)

for (i in y1) {
  angle_dif[i] <- abs(B6_seq$theta_adj[i] - B6_seq$theta_adj[i+1]) - 90
}

angle_dif[n_pts] <- abs(B6_seq$theta_adj[n_pts] - B6_seq$theta_adj[1]) - 90
angle_dif #differences of angles
max(abs(angle_dif)) #maximal difference
thr_ang <- theta_step/2 #threshold for orthogonality [degree]
b_xy_vertex <- b_xy[1:n_PC,]
b_xy_vertex <- as.data.frame(b_xy_vertex)
b_xy_vertex$ortho <- 0
b_xy_vertex
vec <- 1 : n_PC

#loop
for (i in vec) { 
  
  if (abs(angle_dif[i]) <= thr_ang) { 
    b_xy_vertex$ortho[i] <- 1
  } else {
    b_xy_vertex$ortho[i] <- 0
  } #end if-else

  } #end for-loop  
#

b_xy_vertex
y4 <- B7$lnr
distance <- matrix(nrow=n_pts, ncol=2)
distance[,] <- 0

#loop
i <- 0

for (n in y4) {
  fname <- paste("./data/",Img_name,"/b",bnr2,"_",n,".txt", sep="")
  P <- read.table(fname, col.names=c("idx","x","y")) #point cloud
  P_red <- reduce_pointset(P) #call of function
  #browser()
  points(P_red[,2],-P_red[,3], pch=20, asp=3, cex=1.0, col="red") #reduced PC
  #points(P[,2],-P[,3], pch=20, asp=3, cex=1.0, col="blue") #original PC
  P_red2 <- P_red
  dist <- dist_PC(P_red2) #call of function (distance of line)
  dist <- round(dist)
  i <- i + 1
  distance[i,1:2] <- c(n,dist)
} #end of loop

#distance (length) of line segments to be used for calculation of weights
len <- distance[,2]
np <- len
B6_seq$np <- np
B6_seq2 <- B6_seq
y1 <- 1 : n_pts
B6_seq <- B6_seq2
B6_seq

##output of table with theta_adj, ro_adj
fname9 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt",sep="")
write.table(B6_seq, fname9)
#
#########################################################################

if (cas != "100_all+nonortho") {
  
  #averaging of angles (theta_av)
  theta_vec <- rep(0,n_B6_seq)
  np_vec <- rep(0,n_B6_seq)
  theta_vec2 <- rep(0,n_B6_seq)
  np_vec2 <- rep(0,n_B6_seq)
  z <- 1 : n_B6_seq
  
  #loop
  for (i in z) { 
    
    if (B6_seq$theta_ang[i] == theta_ref || B6_seq$theta_ang[i] == alph_ref) {
      B6_seq$ortho[i] <- 1
    } else {
      B6_seq$ortho[i] <- 0 #other orientation than theta_ref or alph_ref
    } #end if-else
    
  } #end of for-loop
  
  B6_seq
  
  for (i in z) {
    
    if (B6_seq$ortho[i] == 1) {
      theta_vec[i] <- B6_seq$theta_adj[i]
      np_vec[i] <- B6_seq$np[i]
    }
    
  } #end for-loop
  
  np_vec
  theta_vec
  
  for (i in z) {

    if (B6_seq$theta_adj[i] < 0) {
      theta_vec[i] <- B6_seq$theta_adj[i] + 90
      np_vec[i] <- B6_seq$np[i]
    }

  } #end for-loop
  
  theta_vec
  np_vec
  theta_vec_red <- subset(theta_vec, theta_vec >= 0)
  theta_vec_red
  n_theta_main <- length(theta_vec_red)
  n_theta_main
  vec <- 1 : n_theta_main
  ang <- theta_vec_red
  ang_mod <- ang

  for (i in vec) { #reduction to range 0...90
    
    if (ang[i] > 90) {
      ang_mod[i] <- (ang[i] - 90)
    } else {
      ang_mod[i] <- ang[i] #corrected
    }
    
  } #end for-loop
  
  ang_mod
  ang_mod2 <- ang_mod
  len <- np_vec
  ang_mod2
  len
  theta_average <- w_av(ang_mod2,len) #call of function
  cat("theta_average= ",theta_average, "\n")
  theta_av <- theta_average
  
  #output of weighted average of angle
  setwd(home_dir)
  f <- paste("./data/",Img_name,"/theta_av_", bnr2,"_ref.txt",sep="")
  write.table(theta_average,f)
} #end if cas != "100_all+nonortho"

if (cas == "100_all+nonortho") {  
  
  #average of angles
  theta_vec <- rep(0, n_B6_seq)
  np_vec <- rep(0,n_B6_seq)
  theta_vec2 <- rep(0, n_B6_seq)
  np_vec2 <- rep(0,n_B6_seq)
  z <- 1 : n_B6_seq
  
  #loop
  for (i in z) { 
    
    if (B6_seq$theta_ang[i] == theta_ref || B6_seq$theta_ang[i] == alph_ref) {
      B6_seq$ortho[i] <- 1
    } else {
      B6_seq$ortho[i] <- 0 
    } #end if-else
    
  } #end for-loop
  
  B6_seq
  
  for (i in z) {
    
    if (B6_seq$ortho[i] == 1) {
      theta_vec[i] <- B6_seq$theta_adj[i] 
      np_vec[i] <- B6_seq$np[i]
    } 
    
  } #end for-loop
  theta_vec
  
  theta_vec_red <- subset(theta_vec, theta_vec != 0)
  n_theta_main <- length(theta_vec_red)
  vec <- 1 : n_theta_main
  ang <- theta_vec_red
  ang_mod <- ang
  ang_mod
  
  for (i in vec) { #reduction of ang_mod to range 0...90
    
    if (ang_mod[i] > 90) {
      ang_mod[i] <- (ang_mod[i] - 90) 
    }
    
  } #end for-loop

  ang_mod
  np_vec_red <- subset(np_vec, np_vec > 0)
  np_vec_red
  len <- np_vec_red
  #
  theta_average <- w_av(ang_mod,len) #call of function
  cat("theta_average= ",theta_average, "\n")
  
  theta_av <- theta_average
  cat("weighted average of angle=",theta_av," degrees","\n")
  theta_av_mod <- theta_av
  
  #test of theta_av_mod
  theta_av_mod
  theta_ref
  
  if (theta_ref > 90) {
    theta_ref_mod <- theta_ref - 90
  } else {
    theta_ref_mod <- theta_ref
  }
  
  dev_theta <- abs(theta_av_mod - theta_ref_mod)
  
  if (dev_theta > thr_theta_av) {
    cat("theta_av deviates too much from theta_ref:", dev_theta, "degree", "\n")
    p_pos = "cor_theta_av"
    setwd(home_dir2)
    source(paste("spObj_intersect_corner_points_v",v_nr,".R",sep=""))
    theta_av_mod
  } #end if
  
  ##output of weighted average of angle
  
  setwd(home_dir)
  f <- paste("./data/",Img_name,"/theta_av_b", bnr2,".txt",sep="")
  write.table(theta_av_mod,file=f)
  #
  
  ##average for a set of parallel lines of different orientation (theta_ref2)
  
  for (i in z) {
    
    if (B6_seq$ortho[i] == 0) { 
      theta_vec2[i] <- B6_seq$theta_adj[i] 
      np_vec2[i] <- B6_seq$np[i]
    } #end if
    
  } #end for-loop
  
  theta_vec2
  theta_vec_red2 <- subset(theta_vec2, theta_vec2 > 0)
  theta_vec_red2 
  n_theta_main2 <- length(theta_vec_red2)
  
  #
  ang2 <- theta_vec_red2
  np_vec2_red <- subset(np_vec2, np_vec2 != 0)
  len2 <- np_vec2_red
  #
  
  theta_average2 <- w_av(ang2,len2) #call of function
  theta_av2 <- theta_average2
  cat("theta_av2 = ",theta_average2, " degree", sep = "","\n") #second main direction
  f <- paste("./data/",Img_name,"/theta_av2_b", bnr2,".txt",sep="")
  write.table(theta_av2,file=f)
  
  n_ortholines2 <- sum(B6_seq$ortho) 
  cat("n_ortholines2= ",n_ortholines2, "\n")
  
  cat("averaging of non-ortho angles","\n" )
  #two possibilities (2 or 1 line)
  
  #two lines
  if (n_nonortholines == 2) {
    p_pos = "cor_theta_av2"
    setwd(home_dir2)
    source(paste("spObj_intersect_corner_points_v",v_nr,".R",sep=""))
    theta_av2_mod
    
    #output theta_av2
    setwd(home_dir)
    f <- paste("./data/",Img_name,"/theta_av2_b", bnr2,".txt",sep="")
    write.table(theta_av2_mod,file=f)
  } #end if n_nonortholines = 2 
  
  #one line
  if (n_nonortholines == 1) {
      p_pos = "cor_theta_av2"
      
      #loop
      for (i in z) {
        
        if (B6_seq$ortho[i] == 1) {
           i2 <- i
        }
        
      } #end for-loop
      
      setwd(home_dir2)
      source(paste("spObj_intersect_corner_points_v",v_nr,".R",sep=""))
      theta_av2_mod
      
      #output theta_av2
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/theta_av2_b", bnr2,".txt",sep="")
      write.table(theta_av2_mod,file=f)
      
  } #end if n_nonortholines = 1

} #end if cas == "100_all+nonortho" 

## checking of theta_average at objects with orthogonal lines

if (cas != "100_all+nonortho") {
  
  theta_av <- theta_average
  cat("weighted average of angle=",theta_av," degrees","\n")
  
  #test of theta_average by means of theta_ref
  theta_av
  theta_ref
  thr_theta_av #threshold
  
  if (theta_av > 90) {
    theta_av_mod <- theta_av - 90
  } else {
    theta_av_mod <- theta_av
  }
  
  dev_theta <- abs(theta_av_mod - theta_ref)
  cat("deviation: theta_av_mod-theta_ref= ",dev_theta,"\n")
  
  if (dev_theta > thr_theta_av) {
    cat("theta_av deviates too much from theta_ref:", dev_theta, "degree", "\n")
    #stop("interaction required")
  } #end if
  
  #check of average by means of phi-angles
  len_phi_all
  phi_all
  phi_average <- w_av(phi_all,len_phi_all)
  cat("phi_average= ",phi_average,"\n")
  theta_ref
  
  #test of averages derived by different methods
  phi_average
  theta_average
  thr_av=0.5
  
  if (theta_av > 90) {
    theta_av_mod <- theta_av - 90 
  } else {
    theta_av_mod <- theta_av
  }
  
  dif_av <- abs(theta_av_mod - (-phi_average)) #phi into img-system
  cat("differences_in_averages= ",dif_av,"\n")
  
  if (dif_av > thr_av) {
    #stop("error -> check theta_average")
    cat("check theta_average","\n")
  }
  
  ##output of weighted average of angle
  theta_average
  setwd(home_dir)
  f <- paste("./data/",Img_name,"/theta_av_b", bnr2,".txt",sep="")
  write.table(theta_average,file=f)
  
  #test for odd number of points
  if (n_pts %% 2 == 1) {
    cat("warning: number of points is odd", "\n")
  }

} #end if cas != "100_all_nonortho" 

cat("end of program 'intersect_corner_points.R' - continue with program 'adjustment_of_corner_coordinates.R'","\n")
setwd(home_dir2)
source(paste("adjustment_of_corner_coordinates_v",v_nr,".R",sep=""))
##########################################################################


