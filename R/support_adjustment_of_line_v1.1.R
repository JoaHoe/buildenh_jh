## name of script: support_adjustment_of_line.R
## purpose: supporting software for special objects
## GNU General Public License (GPL)

cat("version_number= ",v_nr,"\n")

## contents:
# 1.adjustment of all lines stored in matrix 'B6_seq' 
# 2.adjustment of the four longest lines with 'theta_ref','alph_ref' and weighted average (theta_av)
# 3.adjustment of longest line (lnr=1)
################################################################################

## 1.adjustment of all lines stored in matrix 'B6_seq'
# derivation of line parameters (with residuals orthogonal to line)
# calculation in math system

B6 <- B6_seq
B6
x <- length(B6$lnr)
row.names(B6) <- seq(x)
y2 <- 1:x

# loop i
for (i in y2) {
  k4 <- nrow(all_PC[[i]])
  all_PC[[i]][1:k4,]
  
#loop j
  j <- 0
  while(j < k4){
    j <- j + 1
    x <- all_PC[[i]]$x[j]
    y <- all_PC[[i]]$y[j]
    points(x,-y, pch=16, cex=0.4, col="red", asp=1) #plot in graph
  }
  
  x_dat <- all_PC[[i]]$x
  y_dat <- (-all_PC[[i]]$y) #change to math system
  xs <- sum(x_dat)/k4
  ys <- sum(y_dat)/k4
  x_dat_v <- x_dat-xs
  y_dat_v <- y_dat-ys
  N <- (t(x_dat_v)%*%x_dat_v-t(y_dat_v)%*%y_dat_v)
  phi_2 <- (2*t(x_dat_v)%*%y_dat_v)/N
  phi_2
  phi_2_deg <- omega*atan(phi_2) #two solutions (solution 1 or solution 2) are possible for phi_2 (phi_2 or phi_2 +180)
  phi_deg <- 0.5 * (phi_2_deg) #calculation of angle phi (slope angle)
  cat("i=",i, "  phi_deg= ",phi_deg, "xs= ", xs, "ys= ", ys ,"\n")
} #end loop i

#end script 1

################################################################################

## 2.adjustment of four longest lines with theta_ref and alpha_ref + weighted mean
#read in the four long lines with theta_ref and alpha_ref
setwd(home_dir)
points_all <- seq(4)
phi_all <- seq(4)
z=1
y <- B5_long_lines$lnr[1:4]
#

for (i in y) {
  lnr <- i
  cat("line segment =",lnr, "\n")
  f1 <- paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
  PC <- read.table(f1)
  names(PC) <- c("idx","x","y")
  k4 <- nrow(PC)
  
  #adjustment of line
  #loop j
  j <- 0
  while(j < k4){
    j <- j + 1
    x <- PC$x[j]
    y <- PC$y[j]
    points(x,-y, pch=16, cex=0.4, col="red", asp=1)
  }
  x_dat <- PC$x
  y_dat <- (-PC$y) #change to math system
  xs <- sum(x_dat)/k4
  ys <- sum(y_dat)/k4
  x_dat_v <- x_dat-xs
  y_dat_v <- y_dat-ys
  N <- (t(x_dat_v) %*% x_dat_v-t(y_dat_v) %*% y_dat_v)
  phi_2 <- (2*t(x_dat_v) %*% y_dat_v)/N
  phi_2_deg <- omega * atan(phi_2) #two solutions (solution 1 or solution 2) are possible for phi_2 (phi_2 or phi_2 +180)
  phi_deg <- 0.5 * (phi_2_deg) #calculation of angle phi (slope angle)
  cat("line=",lnr, "  phi_deg=",phi_deg, "xs=", xs, "ys=", ys ,"\n")
  phi_all[z] <- phi_deg
  points_all[z] <- B5_long_lines$n_pixel[z]
  z <- z+1
}
phi_all
points_all #length of lines

#weighted mean

max(points_all)
weights <- points_all / max(points_all)
phi_w <- weights * phi_all / sum(weights)
phi_weighted_mean <- sum(phi_w)
cat("weighted mean of 4 longest lines=",phi_weighted_mean,"degrees","\n")
theta_av <- 90 - phi_weighted_mean
fname10 <- paste("./data/",Img_name,"/theta_av_b",bnr2,".txt",sep="")
write.table(theta_av,fname10)
#

#end of script 2
################################################################################

## 3.adjustment of longest line (lnr=1)

#read ref-line
setwd(home_dir)
bnr2 <- 18 #to be changed, for demos: 18(ISPRS7), 11(ISPRS1) 
z <- 1
y <- B6$lnr[6]

for (i in y) {
  lnr <- i
  cat("line segment =",lnr, "\n")
  f1 <- paste("./data/",Img_name,"/b",bnr2,"_",lnr,".txt",sep="")
  PC <- read.table(f1)
  head(PC)
  names(PC) <- c("idx","x","y")
  k4 <- nrow(PC)
  
  #adjustment of line
  #loop j
  j <- 0
  
  while(j < k4){
    j <- j+1
    x <- PC$x[j]
    y <- PC$y[j]
    points(x,-y, pch=16, cex=0.4, col="red", asp=1)
  }
  
  x_dat <- PC$x
  y_dat <- (-PC$y) #change to math system
  xs <- sum(x_dat)/k4
  ys <- sum(y_dat)/k4
  x_dat_v <- x_dat-xs
  y_dat_v <- y_dat-ys
  N <- (t(x_dat_v) %*% x_dat_v-t(y_dat_v) %*% y_dat_v)
  phi_2 <- (2*t(x_dat_v) %*% y_dat_v)/N
  phi_2_deg <- omega*atan(phi_2) #two solutions (solution 1 or 2) 
# are possible for phi_2 (phi_2 or phi_2 + 180)
  phi_deg <- 0.5 * (phi_2_deg) #calculation of angle phi (slope angle)
  cat("line=",lnr, "  phi_deg=",phi_deg, "xs=", xs, "ys=", ys ,"\n")
}

theta_1 <- 90 - phi_deg #solution

#output of result
fname10 <- paste("./data/",Img_name,"/theta_b",bnr2,"_1.txt",sep="")
write.table(theta_1,fname10)
#end of script 3
################################################################################

##end supporting software to script 'adjustment_of_line.R'
