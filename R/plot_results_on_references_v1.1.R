##plot result on references (orthoimage,ground truth)
#name of program (script): plot_results_on_references.R
cat("version_number= ",v_nr,"\n")
#author: Joachim Hoehle
#examples: ISPRS data: image ISPRS7/LCM1, ISPRS1/LCM2
#instructions: use supplementing scripts for interactions if necessary 
cat("####################################################################","\n")
cat("start of program 'plot_results_on_references.R'","\n")
setwd(home_dir)

##input of table with line-pair, vertex_nr and final coordinates (x,y)
f5 <- paste("./results/",Img_name,"/b",bnr2,"_intsec_linepair_vertex_coord.txt",sep="")
intsec_linepair_vertex_coord2 <- read.table(f5)
names(intsec_linepair_vertex_coord2) <- c("line_pair","vertex_nr","x","y")
cat("table with line-pairs,vertex/corner-number,coordinates(x,y)","\n")
print(intsec_linepair_vertex_coord2)

#plot on orthoimage (small scale)
setwd(OrgImgPathname)
img_ref <- readImage(OrgImgFilename)
display(img_ref, method = "raster")
#

##plot coordinates and connecting lines of object on orthoimage
setwd(home_dir)
fname12 <- paste("./results/",Img_name,"/b",bnr2,"_coord_adj_plot.txt",sep="")
b <- read.table(fname12,header=T)
cat("plot of orthoimage", "\n")
k1 <- nrow(b)
names(b) <- c("Points_x","Points_y")
b$Points_y <- (-b$Points_y) #change to img-system
cat("plot of building outline on top of orthoimage","\n")
#

i <- 0
while(i < k1) {
  i <- i+1
  lines(b, col="white", asp=1, type="l", lwd=2, lty=1)
} #end while

##plot result on orthoimage (large scale)
display(img_uds,method = "raster")
n_x <- length(PC_nr)
vec_y <- 1 : n_x

#loop
for (i in vec_y) {
  cat("i=",i,"\n")
  #browser()
  points(as.integer(pc3$col-orig_x), as.integer(pc3$row-orig_y), pch=20, asp=1, cex=0.3, col="green")
  points(xc-orig_x,yc-orig_y,pch=3, asp=1, cex=1.3, col="red")
  
  #loop
  i <- 0
  
  while(i < k1) {
    i <- i+1
    lines(b$Points_x-orig_x, (b$Points_y - orig_y),  col="blue", asp=1, type="l", lwd=2, lty=1)
  } #end while
  
} #end for-loop
# end of plot large scale

cat("does the result agree with the orthoimage?","\n")

if (proc_mode == "demo") {
  cat("if demo - type Y ","\n")
}

answ <- readline ("type Y or N: ")

if (answ == "N") {
  cat ("start again with this object and select other values for 'cas' and/or 'sek' ","\n")
  setwd(home_dir2)
  source(paste("extract_single_building_v",v_nr,".R",sep = ""))  
}

if (answ == "Y") {
  
#plot on Ground Truth (GT)
  setwd(OrgGtsPathname)
  img_GTS <- readImage(OrgGtsFilename)
  display(img_GTS, method="raster")
  
  #loop
  i <- 0
  while(i < k1) {
    i <- i + 1
    lines(b,col="red", asp=1, type="l", lwd=2, lty=1)
  } #end while
} #end if answ="Y")

cat("does the result agree with the Ground Truth?","\n")

if (proc_mode == "demo") {
  cat("if demo - type Y ","\n")
}

answ <- readline("does the result agree with the Ground Truth? type Y or N: ")

if (answ == "Y") { 

##plot of object(building) as graph
  
  if (Img_name == "ISPRS7") { 
    #plot object onto map (Ground Truth, map_ISPRS7)
    par(mai = c(1.02,0.82,0.82,0.42)) #setup of margins/plot region [inches]
    x=0
    y=0
    plot(x,-y, pch=3, cex=1.5,  cex.axis = 1.2, cex.lab=1.5, col="red", asp=1, xlim=c(1,1887), ylim=c(-2557,-1), 
         axes = TRUE, ann = T, frame.plot = TRUE, main = paste("building #", bnr2," of image '",Img_name,"'",sep = ""))
  } #end image7
  
  if (Img_name == "ISPRS1") { 
    #plot object onto map (Ground Truth, map_ISPRS7)
    par(mai = c(1.02,0.82,0.82,0.42)) #setup of margins/plot region [inches]
    x=0
    y=0
    plot(x,-y, pch=3, cex=1.5,  cex.axis = 1.2, cex.lab=1.5, col="red", asp=1, xlim=c(1,1919), ylim=c(-2569,-1), 
         axes = TRUE, ann = T, frame.plot = TRUE, main = paste("building #", bnr2," of image '",Img_name,"'",sep = ""))
  } #end image1
  
  fname12 <- paste("./results/",Img_name,"/b",bnr2,"_coord_adj_plot.txt",sep="")
  setwd(home_dir)
  b <- read.table(fname12,header=T)
  k1 <- nrow(b)
  names(b) <- c("Points_x","Points_y")
  cat("plot of building-outline","\n")

  #loop
  i <- 0
  while(i < k1) {
    i <- i+1
    lines(b, col="black", asp=1, type="l", lwd=1, lty=1)
  } #end while
    
} else { #poor agreement
  cat ("start again with this object and select other values for 'cas' and/or 'sek' ","\n")
  setwd(home_dir2)
  source(paste("extract_single_building_v",v_nr,".R",sep = ""))  
} #end if-else
  
##processing of other objects (buildings)

answ2 <- readline("other buildings to process? type Y or N: ")

if (answ2 == "Y" && proc_mode == "auto") {
  k_y_auto <- k_y_auto + 1 #next building
  
  if (k_y_auto < n_y_auto) {
    setwd(home_dir2)
    source(paste("extract_single_building_v",v_nr,".R",sep = ""))
  } #end if (k_y_auto < n_y_auto)
    
  if (k_y_auto >= n_y_auto) {
    cat(paste("end of processing object ",bnr2, sep = ""),"\n")
    cat("end of processing mode 'auto' ","\n")
    proc_mode <- "NA"
  } #end if (k_y_auto > n_y_auto)
    
} #end if answ2 = "Y" && proc_mode = "auto"  

if (answ2 == "Y" && proc_mode != "auto") { 
  setwd(home_dir2)
  source(paste("extract_single_building_v",v_nr,".R",sep = ""))
}  #end (answ2 = "Y" && proc_mode = "obj_wise")

if (answ2 == "N") {
    cat("end of program 'plot_results_on_references.R'","\n")
    stop("end of program package 'buildenh' ","\n")
} #end if answ2="N"

#end of program 'plot_results_on_references.R'

#end of program package 'buildenh'
################################################################################

