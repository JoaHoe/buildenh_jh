##name of program (script): adjustment of corner coordinates.R
# description: calculation of final corner coordinates
# weighted average of the main direction (theta_av) is used 
# ro-values are calculated by least-squares adjustment
# average of standard deviation of residuals to be used as quality control
cat("version_number= ",v_nr,"\n")
# author: Joachim Hoehle
# instructions: use supplementing scripts in case of problems
cat("##################################################################","\n")

cat("start of program 'adjustment_of_corner_coordinates.R'","\n")
##inputs

#input of plot-parameter
n_pts # number of corners/lines
xc <- plotPar[1]
yc <- plotPar[2]
r_max <- plotPar[3]

# input of adjusted angle (theta_weighted mean)
setwd(home_dir)
fname10 <- paste("./data/",Img_name,"/theta_av_b",bnr2,".txt", sep="")
theta_av <- as.numeric(read.table(fname10))
cat("theta angle (weighted average) = ", theta_av, "[degrees]", "\n")
theta <- theta_av
#

# input of table with adjusted parameters (theta_adj, ro_adj)
f2 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt",sep="")
B8 <- read.table(f2, header=T)
options(digits = 5)
cat("adjusted parameter of line segments (theta_adj, ro_adj) and length of lines (np):","\n")
print(B8)
n_B8 <- nrow(B8)
z <- 1 : n_B8
B8[,8] <- 0
names(B8)[8] <- "ortho"
# 
for (i in z) {
  if (B8$theta_ang[i] == theta_ref || B8$theta_ang[i] == alph_ref) {
    B8$ortho[i] <- 1
  } else {
    B8$ortho[i] <- 0
  } 
} #end loop i

B8S <- B8 #setup of new data frame with adjusted values 
B8S

#introduction of theta_av
for (i in z) {
  if (B8S$ortho[i] == 1) { #buildings with orthogonal lines
    B8S$theta_adj[i] <- theta_av
  } 
} # end loop
# 

for (i in z) {
  if (B8S$ortho[i] == 1 && B8$theta_ang[i] > 90)
    B8S$theta_adj[i] <- B8S$theta_adj[i] + 90
} # end loop

n_pts <- nrow(B8S)

if (sum(B8S$ortho) < n_pts) { 
  cas <- "100_all+nonortho" #case for objects with non-orthogonal lines
}

B8 <- B8S

if (cas != "100_all+nonortho") {
  n_nonortholines2 <- 0
}

n_nonortholines2 #n_nonortholines2 is the number of lines which belong to the object

##input of theta_av2
if (n_nonortholines2 != 0) { 
  setwd(home_dir)
  fname11 <- paste("./data/",Img_name,"/theta_av2_b", bnr2,".txt",sep="")
  theta_av2 <- as.numeric(read.table(fname11))
  cat("theta angle2 (weighted average) = ", theta_av2, "[degrees]", "\n")
} # end if
#

if (n_nonortholines2 == 0 && soph == 1) {
  for (i in z) {
    if (B8S$ortho[i] == 1) { #buildings with orthogonal lines
      B8S$theta_adj[i] <- theta_av
    }
  } # end loop
  
  for (i in z) {
    if (B8S$ortho[i] == 1 && B8$theta_ang[i] > 90)
      B8S$theta_adj[i] <- B8S$theta_adj[i] + 90
  } # end loop
  n_pts <- nrow(B8S)
  
  if (sum(B8S$ortho) < n_pts) { 
    cas <- "100_all+nonortho" #case for objects with non-orthogonal lines
  }

  B8 <- B8S
} #end if n_nonortholines2 = 0 && soph = 1
B8

if (n_nonortholines2 > 0) { #special object
  p_pos <- "cor_adj_coco"
  setwd(home_dir2)
  source(paste("spObj_adjustment_of_corner_coordinates_v",v_nr,".R",sep = "")) 
} #end if
B8

phi_deg <- B8$theta_adj
phi_deg <- (-phi_deg) #change to math-system
options(digits=6)
cat("sequence of angles [degrees]:","\n")
print(phi_deg)
phi <- (phi_deg*pi/180)
m <- length(phi) # number of lines in object
A <- matrix(nrow=2*m, ncol=m) #matrix for building with m corners
A[,] <- 0
A <- design_mat(m,phi) #call of function
A

# Input of approximate coordinates
setwd(home_dir)
fname12 <- paste("./data/",Img_name,"/b",bnr2,"_coord_appr.txt",sep = "")
b01 <- read.table(fname12,col.names="xy")
b01 # approximate corner coordinates
b <- rep(0,m)
b <- b01$xy
n_b <-length(b01$xy)
b_mat <- matrix(nrow=n_b,ncol=2)
b_mat <- as.data.frame(b_mat)
b_mat[1:n_b,1] <- b
b_mat[1:n_b,2] <- 0
names(b_mat)[1:2] <- c("coo","ortho")

i <- 1
j <- 1
#loop i

while (i <= n_B8) {
    b_mat[j,2] <- b_xy_vertex$ortho[i]
    b_mat[j+1,2] <- b_xy_vertex$ortho[i]
    i <- i + 1
    j <- 2*i - 1
 } #end of loop i

b_mat
b <- b_mat$coo

##adjusted coordinates x,y
p <- adjust_coord(A,b) #function call
p
m1 <- length(p[,1])
m2 <- m1/2


##checks by plotting the results

#plot of adjusted corner coordinates

#separated in x und y
Points_x <- rep(0,(m2+1))
Points_x[m2+1] <- p[1,1] #repeat first point (x-coordinate)
k <- 1
i <- 1

#loop i
while(i <= (2*m2 - 1)) {
  Points_x[k] <- p[i,1]
  k <- k + 1
  i <- i + 2
} #end loop i

Points_y <- rep(0,m2+1)
Points_y[m2+1] <- p[2,1] #repeat first point (y-coordinated)
k <- 1
i <- 2

#loop i
while(i <= 2*m2){
  Points_y[k] <- p[i,1]
  k <- k + 1
  i <- i + 2
} #end loop i

## final coordinates of corners
b_seri_xy2 <- cbind(Points_x,Points_y)
dimnames(b_seri_xy2)[[2]] <- list("x","-y")
row.names(b_seri_xy2) <- 1 : (m2+1)
k1 <- m2

i <- 0

#loop i
while(i < k1) {
  i <- i + 1
  x <- Points_x[i]
  y <- Points_y[i]
  points(x,y, pch=20, cex=1.5, col="green", asp=1)
  lines(b_seri_xy2,  col="red", asp=1, type="l", lwd=2, lty=1)
}

cat("adjusted coordinates:","\n")
print(b_seri_xy2)

## output of final coordinates 

#output of final coordinates for plotting
setwd(home_dir)
fname14 <- paste("./results/",Img_name,"/b",bnr2,"_coord_adj_plot.txt",sep="")
write.table(b_seri_xy2,fname14)

#output of point_numbers (vertices) with intersecting line-pairs
f4 <- paste("./data/",Img_name,"/b",bnr2,"_PC_nr.txt",sep="")
setwd(home_dir)
line_nrs <- read.table(f4, col.names = "lnr")
line_nrs2 <- rep(0,m2)

i=1
while (i < m2) { 
  line_nrs2[i] <- line_nrs$lnr[i+1]
  i <- i + 1
}

line_nrs2[m2] <- line_nrs$lnr[1]
intsec_linepair_vertex_coord <- matrix(nrow=m2, ncol=4)
intsec_linepair_vertex_coord [,1] <- paste(line_nrs$lnr,"_",line_nrs2,sep="")
intsec_linepair_vertex_coord [,2] <- 1 : m2
intsec_linepair_vertex_coord [,3:4] <- b_seri_xy2[1:m2,1:2]
#

f5 <- paste("./results/",Img_name,"/b",bnr2,"_intsec_linepair_vertex_coord.txt",sep="")
write.table (intsec_linepair_vertex_coord, f5)
#

#store bnr2 in a file containing all processed buildings
setwd(home_dir)
fname15 <- paste("./results/",Img_name,"/b_all.txt",sep="")
write.table(bnr2, file= fname15, row.names = F, col.names = F, append=TRUE)

#output of adjusted coordinates
setwd(home_dir)
fname13=paste("./results/",Img_name,"/b",bnr2,"_coord_adj.txt",sep="")
write.table(p,fname13)
#

cat("end of 'adjustment_of_corner_coordinates.R'-continue with 'plot results_on_references.R'","\n")
setwd(home_dir2)
source(paste("plot_results_on_references_v",v_nr,".R",sep=""))
#############################################################################################



