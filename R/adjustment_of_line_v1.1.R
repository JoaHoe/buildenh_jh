#name of program: adjustment_of_line.R
cat("version_number= ",v_nr,"\n")
#description: least-squares adjustment of a straight line 
#unknowns are theta (angle) and ro (distance of line from origin)
#residuals are orthogonal to the line
#author: Joachim Höhle
#GNU General Public License (GPL)
cat("###########################################################################","\n")

cat("start of program 'adjustment_of_line.R'","\n")

setwd(home_dir)

##plot point clouds of building and origin of coordinate system (small scale)
x <- 0
y <- 0

##Input of plot-parameters (plotPar)
setwd(home_dir)

if (part == "2parts_1" || part == "2parts_2") { 
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

##plot at large scale
plot(xc,-yc, pch=3, cex=3, col="red", asp=1, 
     xlim=c(xc-r_max2, xc+r_max2), 
     ylim=c(-(yc+r_max2),-(yc-r_max2)),
     main=paste("b", bnr2, sep="")) 

##input of individual pixel clusters (PC)
x <- PC_nr
x3 <- length(PC_nr)
y2 <- 1 : x3
y2 #sequence of lines (not ordered)

##input of list 'PC_all'
setwd(home_dir)
x
i=1

for (n1 in x) {
  cat("n1= ",n1,"\n")
  fname8 <- paste("./data/",Img_name,"/all_PC$PC_nr",n1,".txt", sep=(""))
  all_PC[[i]] <- read.table(fname8)
  i <- i+1
} #end of input of list 'PC_all'

##plot of separated and corrected point clusters
# loop for each point cluster
palette2 <- c("brown", "red", "gray",  "darkgreen", "blue", "magenta", "black", "cyan")

for (i in y2) {
  cat("i=", i, "\n")
  #browser()
  points(all_PC[[i]]$x,(-all_PC[[i]]$y), pch='.', asp=1, cex=2.5, col="red") #change to math-system
} # end of input PCs

##input of line parameter (matrix 'B') without proper sequence of lines
fname9 <- paste("./data/",Img_name,"/unsorted_lines_b",bnr2,".txt", sep="")
B6 <- read.table(fname9, header=T)
B6 #line segments of outline
np2 <- nrow(B6)

##generation of proper sequence of lines
B6 <- subset(B6,select=c(lnr,theta_angle,ro_pixel,n_pixel)) 
B6[,5:6] <- 0
names(B6) <- c("PC_nr","theta_ang","ro_pixel","n_pixel", "theta_adj", "ro_adj")

## storage of B6 for correct sequence
B6_seq <- B6
B6_seq <- B6_seq[1:length(x),]
B6_seq$PC_nr <- x
np <- length(x)
B6_seq[,2:4] <- 0
B6_seq

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

##derivation of adjusted line parameters 
#with residuals orthogonal to line
#calculation in math-system

B6 <- B6_seq
B6
x <- length(B6$PC_nr)
y2 <- 1 : x
phi_all <- rep(0,x)
a_adj_all <- rep(0,x) 
xs_all <- rep(0,x) 
ys_all <- rep(0,x)
len_phi_all <- rep(0,x)

#loop i
for (i in y2) {
  cat("i= ",i,"\n")
  #browser()
  k4 <- nrow(all_PC[[i]])
  all_PC[[i]][1 : k4,]
  len_phi_all[i] <- k4
  
  #loop j
  j <- 0
  while(j <= k4) { 
    j <- j + 1
    x <- all_PC[[i]]$x[j]
    y <- all_PC[[i]]$y[j]
    points(x, -y, pch=16, cex=0.4, col="red", asp=1)
  }
  
  x_dat <- all_PC[[i]]$x
  y_dat <- (-all_PC[[i]]$y) #change to math-system
  xs <- sum(x_dat)/k4 #math-system
  ys <- sum(y_dat)/k4
  xs_all[i] <- xs
  ys_all[i] <- ys
  x_dat_v <- x_dat - xs
  y_dat_v <- y_dat - ys
  N <- (t(x_dat_v) %*% x_dat_v-t(y_dat_v) %*% y_dat_v)
  phi_2 <- (2*t(x_dat_v) %*% y_dat_v)/N
  phi_2_deg <- omega*atan(phi_2) #two solutions are possible for phi_2  
  phi_deg <- 0.5*(phi_2_deg) #calculation of angle phi (slope angle)
  cat("i= ", i, "phi_deg= ", phi_deg, "\n") #calculation of unknown 1
  phi_all[i] <- phi_deg
  a_adj <- ys - xs*tan(phi_deg/omega)
  cat("a_adj= ", a_adj, "\n")
  a_adj_all[i] <- a_adj
  ro_adj <- (a_adj*cos(phi_deg/omega))
  cat("ro_adj= ", ro_adj, "\n")
  
 
  #two solutions are possible due to tangens-function
  #solution 1 
  #calculation of th1_img
  th1_img <- 90 - phi_deg #th1_img in img-system (adjusted)
  
  if ((th1_img < B6[i,2] + 25) && (th1_img > B6[i,2] - 25)) {

  #update of matrix B6 (img-system)
   
    if(th1_img < 0) {
      th1_img <- th1_img + 180
    }
    
    B6[i,5] <- th1_img #solution 1
    names(B6)[7] <- "solution"
    B6[i,7] <- "1"
    
    #calculation of ro_test (img-system)
    ys_img <- (-ys)
    ro_test_img <- cos(B6[i,5]/omega)*xs + sin(B6[i,5]/omega)*(ys_img)
    #
    cat("i= ",i,"\n")
    cat("ro_test_img= ", ro_test_img,"\n")
    
    #update of table B6
    B6[i,6] <- ro_test_img 
    
    #graphical output 
    # calculation in math-system
    th1 <- phi_deg - 90
    
    if(th1 < 0) {
      th1 <- th1 + 180
    }
    
    #calculation of ro_test (math-system)
    th1_arc <- th1/omega
    ro_test <- cos(th1_arc)*xs + sin(th1_arc)*(ys)
    b1 <- (-1/tan(th1_arc))
    a1 <- ro_test/sin(th1_arc)
    abline(a1, b1, col="green",lwd=2) #plot of lines orthogonal to main line
  
  } else { 
      
    ##solution 2
    #calculation of th2_img
    th2_img <- 180 - phi_deg  
    
    #update of table B6
    if (th2_img < 0) {
      th2_img <- 180 + th2_img 
      B6[i,5] <- th2_img 
    } else {
      B6[i,5] <- th2_img
    }
    
    if (B6[i,5] > 180) {
      B6[i,5] <- B6[i,5] - 180
    }
    
    B6[i,7] <- "2"
    
    #calculation of ro
    ys_img <- (-ys)
    ro_test_img <- cos(B6[i,5]/omega) * xs + sin(B6[i,5]/omega) * ys_img #img-system
    cat("i= ",i,"\n")
    cat("ro_test_img= ", ro_test_img,"\n")
    
    #update of table B6
    B6[i,6] <- ro_test_img
    
    #B6 #agreement with the results of Hough-trans? 
    #see at script-line 442:B4 in 'line_detection.R') 
   
    #plot of adjusted line (graph)
    #calculation in math-system
    th2 <- phi_deg 
    th2_arc <- th2/omega
    b <- (-1/tan(th2_arc))
    ro_test2 <- cos(th2_arc) * xs + sin(th2_arc) * ys 
    a <- ro_test2/sin(th2_arc)
    coef2 <- c(a,b)
    
    if (is.finite(a)) {
      abline(coef2, col="blue",lwd=2) #plot of lines parallel to main line
    } else {
      ro_l1 <- B4$ro_pixel[lnr]
      lines(c(ro_l1,ro_l1),c(0,-2569),col="red",lty=1,lwd=2,asp=1) #ref/th1
      lines(c(ro_l1,ro_l1),c(0,1919),col="red",lty=1,lwd=2,asp=1) #orth/th2
    } #end of if-else
    
  } #end solution2

} #end of loop i (generation of line parameters)

cat("line parameters in ordered sequence:","\n")
print(B6)
row.names(B6) <- seq(nrow(B6))

##checking of results

#standard deviation of residuals
B6
x <- nrow(B6)
y2 <- 1 : x
res <- rep(0,k4)
res2 <- rep(0,k4)
thr_line=3 #threshold for sigma of residuals [pixel] 
thr_res=5*thr_line #threshold for single residual
PC_numb <- B6[,1]
cat("PC_numbers= ", PC_numb,"\n")

#loop i (line_number index)
for (i in y2) {
  #browser()
  PC_number <- B6[i,1]
  cat("i= ",i,"PC_nr= ",PC_number, "\n")
  phi3 <- phi_all[i]
  xs <- xs_all[i]
  ys <- ys_all[i]
  
  if (B6$theta_ang[i] > 90) { 
    phi3 <- phi3 + 90 #solution 2
  } else {
    phi3 <- phi3 #solution 1
  }

  cat("phi3= ",phi3,"\n")
  phi3_arc <- phi3/omega
  a_adj3 <- ys - xs*tan(phi3_arc) #calculation of unknown
  ro_adj <- (a_adj3*cos(phi3_arc))
  cat("ro_adj= ", ro_adj, "\n")
  k4 <- nrow(all_PC[[i]])
  cat("k4= ", k4,"\n")
  res <- rep(0,k4)
  all_PC[[i]][1 : k4,]
  row.names(all_PC[[i]]) <- 1 : k4 
 
  #loop 
  j=1 #pixel number
  
  #calculation of residuals
  while(j <= k4) {
    x <- all_PC[[i]]$x[j]
    y <- -all_PC[[i]]$y[j] #change to math-system
    res[j] <- a_adj3*cos(phi3_arc) + x*sin(phi3_arc) - y*cos(phi3_arc) #math-system
    
    if (res[j] > thr_res) {
      cat("index= ", i ,"pixel= ", j, "res= ", res[j], "to be removed","\n")
    } #end if
    
    j <-j + 1
  } #end loop j

  # output of residuals
  res
  setwd(home_dir)
  fname9 <- paste("./data/",Img_name,"/res_PC_nr_",PC_numb[i],".txt",sep="") 
  write.table(res,fname9,col.names=F)
  
  #calculation of standard deviation
  m <- x 
  res1 <- abs(res)
  res_max <- max(res1) #maximal residual [pixels]
  cat("res_max= ", res_max,"\n")
  rt <- t(res1)
  sigma <- sqrt((rt %*% res1)/(m-2))
  cat("standard deviation of residuals= ",sigma,"[pixel]","\n")
  
  # checking with threshold for sigma (thr_line)
  
  if (sigma > thr_line) {
    cat("warning: sigma exceeds threshold", "\n")
    cat("PC_index= ",i,"\n") #adjustment of this line must be corrected
    cat("PC_nr= ", PC_numb[i],"\n")
  } #end if

} #end loop i

#checking of angles
x <- length(B6$PC_nr)
y2 <- 1 : x

for (i in y2) {
  dif_ang <- abs(B6$theta_ang[i]) - abs(B6$theta_adj[i])
  cat(paste("difference in angle: ", i, "is: ", dif_ang,"degrees","\n"))
} #end for-loop

print(B6)

cat("correction of adjustment-parameter?", "\n")
cat("if proc_mode='demo' -> type N", "\n")
answ4 <- readline("type Y or N: ") #interaction required

if (answ4 == "Y") {
  p_pos <- "cor_adj_line" #correction of adjustment parameter
  setwd(home_dir2)
  source(paste("spObj_adjustment_of_line_v",v_nr,".R",sep=""))
} #end if

print(B6) #final result of line adjustment
B6 <- subset(B6,select= -solution)

# output of results (B6)
setwd(home_dir)
fname9 <- paste("./data/",Img_name,"/param_adj_b",bnr2,".txt",sep="") #building parameter
write.table(B6,fname9)
cat("end of 'adjustment_of_line.r'","\n")
#

cat("continue with 'intersect_corner_points.R'","\n")
setwd(home_dir2)
source(paste("intersect_corner_points_v",v_nr,".R",sep=""))
################################################################################

