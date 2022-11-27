#name of program (script): extract_single_building.R
cat("version_number= ",v_nr,"\n")
#description: extraction of one building  
#from image "buildings of generated land cover map"
#orthoimage: ISPRS data "Vaihingen" of areas: #7, #1
#instruction: use plot of building numbers producible in 'support_extract_single_building' for selecting of object 
#author: Joachim Hoehle
#GNU General Public License (GPL)
cat("###########################################################################","\n")

cat("start of program 'extract_single_building'","\n")
setwd(home_dir)

##mode of processing
cat("select mode of processing? - demo: 1, obj_wise: 2, auto: 3","\n") 
answ <- readline("mode of processing? - type 1, 2, or 3: ") #processing mode

if (answ == "1" && Img_name == "ISPRS7") { #processing of one example
  proc_mode <- "demo" #object b18 (ISPRS7)
  bnr2 <- 18
  part = "no_part"
} 

if (answ == "1" && Img_name == "ISPRS1") { #processing of one example
  proc_mode <- "demo" 
  bnr2 <- 11
  part = "no_part"
}

if (answ == "2") { 
  proc_mode <- "obj_wise" #object-wise processing
}

if (answ == "3") { 
  proc_mode <- "auto" #automatic processing of orthoimage
}

# if (answ == "4") { 
#   stop("terminate automatic processing")
# }

if (proc_mode == "obj_wise") {   
  bnr2 <- readline("type the label of building to be processed: ") #label of building by manual input
  bnr2 <- as.integer(bnr2)
}

if (proc_mode == "auto") { 
  bnr2 <- y_auto[k_y_auto]
}

cat("processing mode= ", proc_mode,"\n")

if (substr(bnr2,3,3) == "1" || substr(bnr2,3,3) == "2") {
    bnr2_part <- bnr2
} else {
    bnr2_orig <- bnr2
    part <- "no_part"
    bnr2_part <- "NA" #setup for partition
} 

#end if_proc_mode

cat("label of building to be extracted= ", bnr2,"\n")

if (part == "2parts_1" || part == "2parts_2") { #repeated? see line 52
  bnr2 <- as.numeric(substr(bnr2,1,2))
  bnr2_orig <- bnr2
} #end if

if (part == "no_part") {
  bnr2_orig <- bnr2
}

cat("label of building to be extracted=", bnr2,"\n") #check if new number is necessary
#

#input of enhanced billede
setwd(home_dir)
LCM_enh_b=readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3_scaled_2.jpg",sep = "")) #classification by method JH, scaled affine
display(LCM_enh_b, method="browser") #use for checking of image
#display(LCM_enh_b, method="raster")

##enhancement of raster image
LCM_enh_b_t <- thresh(LCM_enh_b,2,2,0.01) #thresholding -> white outlines
display(LCM_enh_b_t, method="raster")
LCM_enh_b_t_f <- fillHull(LCM_enh_b_t)
display(LCM_enh_b_t_f,"raster")
LCM_label_A <- bwlabel(LCM_enh_b_t_f) #labeling for area
cat('number of buildings=', max(LCM_label_A),'\n')
display(LCM_label_A, method="raster")
LCM_enh_b_t_f_t2 <- thresh(LCM_enh_b_t_f, 2,2,0.01) #threshholding
display(LCM_enh_b_t_f_t2, method="raster")
LCM_label <- bwlabel(LCM_enh_b_t_f_t2) #labeling for perimeter
display(LCM_label)

##display as negative
LCM_enh_b_t_neg <- (1-LCM_enh_b_t)
#display(LCM_enh_b_t_neg, method="browser")
display(LCM_enh_b_t_neg, method="raster")

##extraction of features (area,radius)
setwd(home_dir)
coor <- computeFeatures.moment(LCM_label) #geometric features (moment)
n9 <- nrow(coor)
cat('number of buildings=',n9,'\n')
shap<-computeFeatures.shape(LCM_label) #geometric features (shape)
shap
n8 <- nrow(shap) #number of buildings to enhance
cat('number of buildings=',n8,'\n')
shap_A<-computeFeatures.shape(LCM_label_A) #geometric features (shape)
shap_A
n7<-nrow(shap_A)
shap1_A<-matrix(nrow=n7,ncol=10)
y<-1:n7
#

for (n in y){
  shap1_A[,2:7]<-shap_A[,1:6]
}
shap1_A[,1] <- y
shap1_A<-data.frame(shap1_A)
shap1_A[,8:9] <- coor[,1:2]
shap1_A[,10] <- coor[,5]
shap2_A <- subset(shap1_A,shap1_A[,2] >= 3086) #removes buildings of area < 5mx5m (3086 pixels)
n8 <- nrow(shap2_A)
cat('number of buildings after area-thresh-holding=',n8,'\n')
rownames(shap2_A) <- 1:n8
names(shap2_A)<-c("bnr","area", "perimeter", "radius.mean", "radius.sd", "radius.min","radius.max","cx","cy","alpha_arc")
shap2_A
y3 <- nrow(shap2_A)
shap2_A_red3<-shap2_A #if no points are removed
shap2_A_red3
#

#output of shape-data
f2=paste("./data/",Img_name,"/shap2_A_red3.csv",sep = "")
write.table(shap2_A_red3,f2)
#

##extraction of one building
cas <- "NA" #case "NA" 
sek <- "NA"
xc <- shap2_A_red3[bnr2,8]
yc <- shap2_A_red3[bnr2,9]
points(xc, yc, col="blue",asp=1, pch=16, cex=1.0)
omega=180/pi
alpha <- shap2_A_red3[bnr2,10]*omega # approximate orientation angle of object
r_max <- shap2_A_red3[bnr2,7]

##write.table plot-parameters of building
plotPar <- c(xc,yc,r_max,alpha)
plotPar_orig <- plotPar

##generation of table (bnr/bnr2)
names(shap2_A_red3) <- c("bnr2","area", "perimeter", "radius.mean", "radius.sd", "radius.min",
                         "radius.max","cx","cy","alpha_arc")
nrow(shap2_A_red3)
y1 <- 1: nrow(shap2_A_red3)
obj_nrs <- matrix(nrow = nrow(shap2_A_red3), ncol=2) #table
colnames(obj_nrs) <- c("bnr","bnr2")
names(shap1_A) <- c("bnr","area", "perimeter", "radius.mean", "radius.sd", "radius.min",
                    "radius.max","cx","cy","alpha_arc")
nrow(shap1_A)
vec <- 1 : nrow(shap1_A)
for (i in y1) {
  for (k in vec) { 
    if (shap1_A$area[k] == shap2_A_red3$area[i] ) {
      obj_nrs[i,1] <- shap1_A$bnr[k]
      obj_nrs[i,2] <- shap2_A_red3$bnr2[i]
    } # end if
  } #end for k
} #end for i

obj_nrs #list with bnr/bnr2 to be used in conversion bnr2<-bnr (and vice versa)
#

##find 'bnr' by 'bnr2'
cat("bnr2= ", bnr2,"\n")

for (i in y1) {
  if (obj_nrs[i,2] == bnr2) {
    bnr <- obj_nrs[[i,1]]
  } #end if
} #end for-loop

cat("bnr= ", bnr,"\n")

##label of building after first labeling (bnr)
is_bnr <- LCM_label@.Data == bnr
display(is_bnr)
display(is_bnr,"raster")
#

#conversion to a vector
coords <- data.frame(x=as.numeric(row(is_bnr)),y=as.numeric(col(is_bnr)), is_bnr=as.numeric(is_bnr))
coords <- coords[coords$is_bnr == 1,] #removal of pixels which do not have the label of the building

##plot of PC and checkpoints
r_max2 <- round(1.1*r_max)
plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
#plot(coords$x, coords$y,pch=16,cex=0.2,col="black",asp=1,xlim=c(1,1887),ylim=c(2557,1), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE) #small scale
points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
#

##output as tiff-image
file1 <-  paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
tiff(file1, width=578, height=578, units="px", bg = "white")
r_max2 <- round(1.1*r_max)
plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
dev.off()
#

##dy-range of window (required for scale factor in program 'sequence_of_lines.R')
ylim <- c(yc+r_max2, yc-r_max2)
ylim[1] - ylim[2]

##calculation of window size
par("usr")
dy_window_plot <- abs(par("usr")[3] - par("usr")[4]) #range of window

##output plot parameter
setwd(home_dir)
plotPar[5] <- dy_window_plot
plorPar_orig <- plotPar
f1 <- paste("./data/",Img_name,"/param_b_",bnr2,".txt",sep="")
write.table(plotPar,file=f1)
f2 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
save(plotPar, file=f2) #parameter xc, yc, r_max, alpha, dy_window_plot
#

#input image
f1 = paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
b_new8 <- readImage(f1)
display(b_new8, method="raster")
#

##partition of object?
cat("partition of object? - no parts: 0, two parts: 1 (first object) or 2 (second object)","\n")

k_part <- readline("select partition-type - if demo - type 0: ") #when object must be parted
if (k_part == "0") { #no partition of object
  part <- "no_part"
}

if (k_part == "1") { #first part
  part <- "2parts_1"
}

if (k_part == "2") { #second part
  part <- "2parts_2"
  #n_y_auto <- n_y_auto + 1
}

if (part == "2parts_1" || part == "2parts_2") {
  bnr2_part <- bnr2
  bnr2 <- as.numeric(substr(bnr2,1,2))
  bnr2_orig <- bnr2
  p_pos <- "cor_sep"
}
#

if(part != "no_part" && p_pos == "cor_sep") {
  setwd(home_dir2)
  source(paste("spObj_extract_single_building_v",v_nr,".R",sep="")) #special object
}

##output of file with PC-coordinates (required in program 'line_detection.R')
coords2 <- subset(coords, select=c(x,y))
head(coords2)
N1 <- length(coords$x)
x_dat <- rep(0,N1)
y_dat <- rep(0,N1)
idx <- rep(0,N1)
idxy <- cbind(idx, x_dat, y_dat)
idxy <- coords2 #use of old object name

#idxy
rownames(idxy) <- 1 : N1
head(idxy)
nrow(idxy)
setwd(home_dir)
fname2 <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv", sep="") #image
write.table(idxy, fname2,  sep= " ", row.names=T) ##output of pixel cluster for one building
#end of output

cat("end of program 'extract_single_building.R' ","\n") 
cat("continue with 'line_detection.R' ","\n")
setwd(home_dir2)
source(paste("line_detection_v",v_nr,".R", sep=""))
###################################################################################


