#name of program (script): adjustment_of_line.R
cat("version_number= ",v_nr,"\n")
#description: least-squares adjustment of a straight line 
#parameters are theta (angle) and ro (distance of line from origo)
#residuals are orthogonal to the line
#author: Joachim Hoehle
cat("###########################################################################","\n")

cat("start of program 'adjustment_of_line.R'","\n")
setwd(home_dir)

##preparation
x3 <- length(PC_nr)
all_PC <- list() #generation of a list
for (i in 1:x3){
  all_PC[[i]] <- "PC"
}
all_PC

##plot point clouds of building and origo (small scale)
x <- 0
y <- 0

##Input of plot-parameters (plotPar)
setwd(home_dir)
if (part == "2parts_1" || part == "2parts_2") { #required?
   bnr2 <- bnr2_part
}

##load of plot parameter
f <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
load(f)
xc <- plotPar[1]
yc <- plotPar[2]
r_max <- plotPar[3]
#

if (part != "no_part") {
   bnr2 <- bnr2_part  
}

## plot in large scale
plot(xc,-yc, pch=3, cex=3, col="red", asp=1, xlim=c(xc-r_max2, xc+r_max2), ylim=c(-(yc+r_max2),-(yc-r_max2)),main=paste("b", bnr2, sep="")) #large scale

## input of individual PixelClusters (PC)
x <- as.vector(PC_nr)
x3 <- length(PC_nr)
y2 <- 1 : x3
y2 #sequence of lines (not ordered)

## Input of list PC_all
setwd(home_dir)
i=1
for (n in x) {
  fname8 <- paste("./data/",Img_name,"/all_PC$PC_nr",n,".txt", sep=(""))
  all_PC[[i]] <- read.table(fname8)
  i <- i+1
} #end of input of list PC_all (loop "i")

## plot of separated and corrected point clusters
# loop for each point cluster
palette2 <- c("brown", "red", "gray",  "darkgreen", "blue", "magenta", "black", "cyan")
for (i in y2) {
  cat("i=", i, "\n")
  n <- i
  #browser(text="line-number= ",i)
  points(all_PC[[i]]$x,(-all_PC[[i]]$y), pch='.', asp=1, cex=2.5, col=palette2[n]) #change to math-system
} # end of input PCs

## Input of line parameter (matrix 'B') without proper sequence of lines
fname9 <- paste("./data/",Img_name,"/unsorted_lines_b",bnr2,".txt", sep="")
B6 <- read.table(fname9, header=T)
B6 #line segments of outline
np2 <- nrow(B6)

## generation of proper sequence of lines
B6 <- subset(B6,select=c(lnr,theta_angle,ro_pixel,n_pixel)) 
B6[,5:6] <- 0
names(B6) <- c("PC_nr","theta_ang","ro_pixel","n_pixel", "theta_adj", "ro_adj")

## Storage of B6 for correct sequence
B6_seq <- B6
B6_seq <- B6_seq[1:length(x),]
B6_seq$PC_nr <- x
np <- length(x)
B6_seq[,2:4] <- 0

# loop
i <- 1
while (i <= np) {
  j <- 1
  while (j <= np2) {
    if (B6_seq$PC_nr[i] == B6$PC_nr[j]) {
      B6_seq[i,2:4] <- B6[j,2:4]
    } #end if
    j <- j + 1
  } #end loop j
  i <- i+1
} #end loop i

B6_seq

## derivation of line parameters (with residuals orthogonal to line)
# calculation in math-system

B6 <- B6_seq
x <- length(B6$PC_nr)
y2 <- 1:x
i=1

#loop i
for (i in y2) {
  k4 <- nrow(all_PC[[i]])
  k4
  all_PC[[i]][1 : k4,]
  #loop j
  j <- 0
  while(j < k4) {
    j <- j + 1
    x <- all_PC[[i]]$x[j]
    y <- all_PC[[i]]$y[j]
    points(x, -y, pch=16, cex=0.4, col="red", asp=1)
  }
  x_dat <- all_PC[[i]]$x
  y_dat <- (-all_PC[[i]]$y) #change to math-system
  xs <- sum(x_dat)/k4
  ys <- sum(y_dat)/k4
  x_dat_v <- x_dat-xs
  y_dat_v <- y_dat-ys
  N <- (t(x_dat_v) %*% x_dat_v-t(y_dat_v) %*% y_dat_v)
  phi_2 <- (2*t(x_dat_v) %*% y_dat_v)/N
  phi_2_deg <- omega*atan(phi_2) #two solutions (solution 1 or solution 2) are possible for phi_2   (phi_2 or phi_2 + 180)
  phi_deg <- 0.5*(phi_2_deg) #calculation of angle phi (slope angle)
  
  # solution 1
  th1_img <- 90 - phi_deg #rotation angle in math-system (first axis)
  if ((th1_img < B6[i,2] + 25) && (th1_img > B6[i,2] - 25)) {
    th1 <- phi_deg - 90 #th1 in math system -change
    B6[i,5] <- (-th1) #solution 1
    th1_arc <- (th1/omega) #math system
    ro_test <- cos(th1_arc)*xs + sin(th1_arc)*(ys) #calculation in math system
    B6[i,6] <- (ro_test) #transfer to img_system -change

  #graphical output
    b1 <- (-1/tan(th1_arc))
    a1 <- ro_test/sin(th1_arc)
    abline(a1, b1, col="green",lwd=2) #plot of lines orthogonal to main line
  } else { #solution 2
    th2 <- phi_deg #math_system, rotation of first axis by (phi_degree + 90 degrees)
    th2_arc <- th2/omega
    B6[i,5] <- (-th2) #change of sign due to B6 is in img-system
    ro_test2 <- cos(th2_arc) * xs + sin(th2_arc) * ys #(math-system)
    B6[i,6] <- ro_test2 
    b <- (-1/tan(th2_arc))
    a <- ro_test2/sin(th2_arc)
    abline(a,b, col="blue",lwd=2) #plot of lines parallel to main line
  } #end of if-else
} #end of loop i (generation of line parameters)
cat("line parameters in ordered sequence:","\n")
print(B6)
row.names(B6) <- seq(nrow(B6))

#check of results
i=1
for (i in y2) {
  dif_ang <- abs(B6$theta_ang[i]) - abs(B6$theta_adj[i])
  cat("difference in angle: ", i, "is: ", dif_ang,"\n")
} #end for-loop
print(B6)

cat("manual correction of adjustment-parameter?", "\n")
cat("if proc_mode='demo' -> type N", "\n")
answ4 <- readline("type Y or N: ")

if (answ4 == "Y") {
  p_pos <- "cor_adj_line"
  setwd(home_dir2)
  source(paste("spObj_adjustment_of_line_v",v_nr,".R",sep="")) 
}

print(B6)

## Output of results (B6)
setwd(home_dir)
fname9 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt",sep="") #building parameter
write.table(B6,fname9)
cat("end of 'adjustment_of_line.r'","\n")
#

cat("continue with 'intersect_corner_points_v8.R","\n")
setwd(home_dir2)
source(paste("intersect_corner_points_v",v_nr,".R",sep=""))
################################################################################

