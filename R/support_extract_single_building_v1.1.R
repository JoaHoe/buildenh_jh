## name of script: support_extract_single_building.R
## purpose: scripts for special cases 
cat("version_number= ",v_nr,"\n")

##contents:

# 1. plot numbers of all buildings (before and after applying a threshold for minimum size of area)
# 2. plot number of one building (before and after applying a threshold for minimum size of area)
# 3. find number 'bnr2' by number 'bnr' and vice versa
# 4. generation of table with bnr/bnr2 columns
# 5. find bnr2 by bnr and vice versa by using the generated table 
###########################################################################################

## 1.plot numbers (label) of buildings (before and after 
#    applying a threshold for minimum size of area)

#numbering without area thresholding: object label is 'bnr'
#use proper image under 'Plots'  (apply arrows: <- or -> )

#input of LCM with all outlines of objects
setwd(home_dir)
LCM_enh_b=readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3_scaled_2.jpg",sep="")) #classification by JH, scaled affine

#display of enhanced billede
#display(LCM_enh_b, method="browser") #for checking of image
display(LCM_enh_b, method="raster")

y1 <- 1 : nrow(coor)
for (i in y1) {
  x <- coor[i,1]
  y <- coor[i,2]
  text(x,y,i,cex=1.2,col="red")
} #end for loop

#end of script plot numbers (label) of buildings
#############################################################################################################

## 2.plot numbers of buildings after applying the threshold 
#    for minimum size (3086 pixels)

#input of enhanced image
LCM_enh_b <- readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3_scaled_2.jpg",sep="")) #classification by JH, scaled affine
#display(LCM_enh_b, method="browser") #for checking of image
display(LCM_enh_b, method="raster")
#
y1 <- 1:nrow(shap2_A_red3)

for (i in y1) {
  nr <- i 
  xc <- shap2_A_red3[i,8]
  yc <- shap2_A_red3[i,9]
  text(xc,yc,nr,cex=1.2,col="red")
} #end for loop

#end of script plot numbers of buildings after applying the threshold for minimum size
#############################################################################################################

## 3.plot the number of one building (before and after applying 
#    a threshold for minimum size of area)

display(LCM_enh_b, method="raster")
i <- 18 #select the number of building (bnr2), here b18
x <- coor[i,1]
y <- coor[i,2]
text(x,y,i,cex=1.2,col="red")

#plot number of one building (bnr2)
#numbering with thresholding using area >= 3086 pixels: bnr2
#use proper image under 'Plots' (use arrows: <- or -> )

display(LCM_enh_b, method="raster")
names(shap2_A_red3) <- c("bnr2","area","perimeter","radius.mean","radius.sd","radius.min","radius.max","cx","cy","alpha_arc")
y1 <- 1 : nrow(shap2_A_red3)
i <- 18 #type number (bnr2)
xc <- shap2_A_red3[i,8]
yc <- shap2_A_red3[i,9]
text(xc,yc,i,cex=1.2,col="red")
#
#end of script 3
###############################################################################################################

## 4.find number 'bnr2' by number 'bnr' and vice versa

bnr=18 #type building number, here b18
cat("bnr=", bnr, "\n")

y1 <- 1 : nrow(shap2_A_red3)

for (n in y1) {
  if (shap2_A_red3$bnr[n] ==  bnr) {
    bnr2 <- n
    cat("bnr2=",bnr2,"\n")
  } #end if
} #end for-loop

#end of script
##############################################################################################################

## 5.generation of table with bnr/bnr2 columns

names(shap2_A_red3) <- c("bnr2","area","perimeter","radius.mean","radius.sd","radius.min",
                       "radius.max","cx","cy","alpha_arc")
nrow(shap2_A_red3)
y1 <- 1 : nrow(shap2_A_red3)
obj_nrs <- matrix(nrow=nrow(shap2_A_red3), ncol=2)
colnames(obj_nrs) <- c("bnr","bnr2")
names(shap1_A) <- c("bnr","area", "perimeter", "radius.mean", "radius.sd", "radius.min",
                    "radius.max","cx","cy","alpha_arc")
vec <- 1 : nrow(shap1_A)
for (i in y1) {
  for (k in vec) { 
    if (shap1_A$area[k] == shap2_A_red3$area[i] ) {
      obj_nrs[i,1] <- shap1_A$bnr[k]
      obj_nrs[i,2] <- shap2_A_red3$bnr2[i]
    } #end if
  } #end for k
} #end for i
obj_nrs
#

#end of script 5
######################################################################################################

## 6.find bnr2 by bnr and vice versa using a table

answ <- readline("type bnr: ")
bnr <- as.integer(answ)
cat("bnr= ",bnr,"\n")
y1 <- 1 : nrow(shap2_A_red3)

for (i in y1) {
  if (obj_nrs[i,1] == bnr) {
    bnr2 <- obj_nrs[i,2]
  } #end if
} #end for-loop

cat("bnr2= ", bnr2,"\n")

#find 'bnr' by 'bnr2'
answ <- readline("type bnr2: ")
bnr2 <- as.integer(answ)
cat("bnr2= ",bnr2,"\n")

for (i in y1) {
  if (obj_nrs[i,2] == bnr2) {
    bnr <- obj_nrs[i,1]
  } #end if
} #end for-loop

cat("bnr= ", bnr,"\n")
#

#end of script 6
#########################################################################################################

##end of 'support_extract_single_building.R'

