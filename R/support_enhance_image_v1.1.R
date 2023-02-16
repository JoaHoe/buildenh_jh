## name of script: support_enhance_image_v1.1.R
cat("version_number= ",v_nr,"\n")
## supporting software for program "enhance_image.R"
## author: Joachim Höhle
## instruction: run all programs in "demo" mode before using 'support_enhance_image.R'
## GNU General Public License (GPL)

##contents:
# 1.checking of scale
# 2.scaling of image (simple and affine)
# 3.cutting out of net-image
# 4.affine scaling with manual measurement in browser
# 5.affine transformation without measuring

###################################################################################

## 1.checking of scale

#orthoimages of ISPRS-data

#ISPRS1
setwd(home_dir)
LCM_b_1 <-  readImage(paste("./data/", Img_name,"/images/LCM_cart_enh_b3.jpg",sep = "")) #land cover map (ISPRS1)
display(LCM_b_1, method = "raster")
colorMode(LCM_b_1)
LCM_b_2 <- imageData(LCM_b_1)
co_xy <- locator(4,"p") #digitize of sides: left (1), right (2), upper (3), lower (4)
co_xy #check
dcol <- (co_xy$x[2]-co_xy$x[1])
drow <- (co_xy$y[4]-co_xy$y[3])

if (abs(round(dcol - img_x_max)) < 3 && abs(round(drow - img_y_max)) < 3) { #orthoimage #1
  cat("scaling is OK")
} else {
  cat("x_deviation=",abs(round(dcol-img_x_max)),"\n")
  cat("y_deviation=",abs(round(dcol-img_y_max)),"\n")
  stop("decide if scaling must be repeated")
}

#end 1.checking of scale
#############################################################################################################

## 2.scaling of image

#orthoimages of ISPRS-data

#ISPRS7

#measurement in raster image
LCM_b=readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3.jpg",sep=""))  #unscaled image
#display(LCM_b, method="browser") #find coordinates of corners
display(LCM_b, method="raster") #find coordinates of corners
co_xy <- locator(4) #digitize of corners
co_xy
x_left_side <- round(co_xy$x[1])
x_right_side <- round(co_xy$x[2])
y_upper_side <- round(co_xy$y[3])
y_lower_side <- round(co_xy$y[4])
x_left_side; x_right_side
y_upper_side; y_lower_side
LCM_b_netto <- LCM_b[x_left_side:x_right_side,y_upper_side:y_lower_side] #cut out net image 
#

#simple scaling
LCM_b_mod=resize(LCM_b_netto,output.dim=c(1887,2557),output.origin=c(0,0),w=1887,h=2557) #resize orthoimage #7
#end of simple scaling

#affine scaling
LCM_b_1887_2557 <- readImage("./data/ISPRS7/images/LCM_cart_enh_b3_1887_2557_new.jpg")
display(LCM_b_1887_2557, method = "raster")
#measurement of sides
co_xy<-locator(4,"p") #digitize of corners

#determination of scale factors
mx<-1887/(co_xy$x[2]-co_xy$x[1]) #orthoimage ISPRS7
my<-2557/(co_xy$y[4]-co_xy$y[3])
m1=matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
m1

display(affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557))) #image ISPRS7
LCM_b_mod2 <- affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557)) #accurate solution
LCM_b_2 <- affine(LCM_b_1, m1, output.dim=c(1887,2557)) #accurate solution
LCM_b_2 <- affine(LCM_b, m1, output.dim=c(1887,2557)) #accurate solution
display(LCM_b_2)
LCM_b_mod2 <- affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557)) #accurate solution
writeImage(LCM_b_2,"./data/ISPRS7/LCM_cart_enh_b3_affin_new.jpg")

#display
LCM_b_3 <- readImage("./data/ISPRS7/LCM_cart_enh_b3_affin_new.jpg")
display(LCM_b_3)
display(LCM_b_3,"raster")
###############################################

##ISPRS1
co_xy$x[1] <- 130 #to be checked and modified
co_xy$x[2] <- 1819
co_xy$y[3] <- 147
co_xy$y[4] <- 2408

#determination of scaling factors
mx <- 1919/(co_xy$x[2]-co_xy$x[1]) 
my <- 2569/(co_xy$y[4]-co_xy$y[3]) 
m1 <- matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
m1
#end scale factors for ISPRS1

########################################################################################################

## 3.cutting out of net-image

#ISPRS7
LCM_b=readImage("./data/ISPRS7/LCM_cart_enh_b3.jpg")
#display(LCM_b, method="browser") #find coordinates of corners
display(LCM_b, method="raster") #find coordinates of corners
co_xy <- locator(4) #digitizing of side-lines in raster image
co_xy
x_left_side <- round(co_xy$x[1]) #left side
x_right_side <- round(co_xy$x[2]) #right side
y_upper_side <- round(co_xy$y[3]) #upper side
y_lower_side <- round(co_xy$y[4]) #lower side
#
x_left_side
x_right_side
y_upper_side
y_lower_side
LCM_b_netto <- LCM_b[x_left_side:x_right_side,y_upper_side:y_lower_side] #cutout net-image #7
#

##ISPRS1
setwd(home_dir)
LCM_b=readImage("./data/ISPRS1/images/LCM_cart_enh_b3.jpg")
display(LCM_b, method="browser") #find coordinates of sides (Zoom:200%)
#

co_xy$x[1] <- 130 #to be checked
co_xy$x[2] <- 1819
co_xy$y[3] <- 147
co_xy$y[4] <- 2408

 x_left_side <- co_xy$x[1] #left side
x_right_side <- co_xy$x[2] #right side
y_upper_side <- co_xy$y[3] #upper side
y_lower_side <- co_xy$y[4] #lower side

LCM_b_netto <- LCM_b[x_left_side:x_right_side,y_upper_side:y_lower_side] #cutout net-image #1
display(LCM_b_netto)

#end 3.cutting-out of net-image

###############################################################################################

## 4.affine scaling with manual measurement in browser

#ISPRS7

#results of script 'enhance_image'
setwd(home_dir)
LCM_b <- readImage(paste("./data/", Img_name,"/images/LCM_cart_enh_b3.jpg",sep = "")) #land cover map theme 'building'
display(LCM_b, method = "browser")

#measurement of sides
co_xy <- matrix(rep(0,8),nrow = 4)
co_xy <- data.frame(co_xy)
names(co_xy) <- c("x","y")
co_xy

#ISPRS7
co_xy$x[1] <- 128 #to be checked
co_xy$x[2] <- 1789
co_xy$y[3] <- 146
co_xy$y[4] <- 2396

x_left_side <- co_xy$x[1] #left side
x_right_side <- co_xy$x[2] #right side
y_upper_side <- co_xy$y[3] #upper side
y_lower_side <- co_xy$y[4] #lower side

LCM_b_netto <- LCM_b[x_left_side:x_right_side,y_upper_side:y_lower_side] #cutout net-image #7
display(LCM_b_netto)

#determination of scaling factors
mx <- 1887/(co_xy$x[2]-co_xy$x[1]) 
my <- 2557/(co_xy$y[4]-co_xy$y[3]) 
m1 <- matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
m1

#affine transformation
LCM_b_2 <- EBImage::affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557)) #orthoimage7
display(LCM_b_2) #display on browser
##

##ISPRS1
setwd(home_dir)
LCM_b=readImage("./data/ISPRS1/images/LCM_cart_enh_b3.jpg")
display(LCM_b, method="browser") #find coordinates of sides (Zoom:200%)
#

#measurement of sides
co_xy <- matrix(rep(0,8),nrow = 4)
co_xy <- data.frame(co_xy)
names(co_xy) <- c("x","y")
co_xy

co_xy$x[1] <- 130 #to be checked
co_xy$x[2] <- 1819
co_xy$y[3] <- 147
co_xy$y[4] <- 2408

x_left_side <- co_xy$x[1] #left side
x_right_side <- co_xy$x[2] #right side
y_upper_side <- co_xy$y[3] #upper side
y_lower_side <- co_xy$y[4] #lower side

LCM_b_netto <- LCM_b[x_left_side:x_right_side,y_upper_side:y_lower_side] #cutout net-image #1
display(LCM_b_netto)

#determination of scaling factors
mx <- 1919/(co_xy$x[2]-co_xy$x[1]) 
my <- 2569/(co_xy$y[4]-co_xy$y[3]) 
m1 <- matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
m1
#end scale-factors for ISPRS1

#affine transformation
LCM_b_2 <- EBImage::affine(LCM_b_netto,m1,filter='bilinear',output.dim=c(1919,2569)) #ISPRS1
display(LCM_b_2) #display on browser
#

#end 4.affine scaling with manual measurement in browser
##########################################################################################################

## 5.affine transformation without measuring

#ISPRS7
setwd(home_dir)
LCM_b_1 <- readImage("./data/ISPRS7/images/LCM_cart_enh_b3.jpg")
display(LCM_b_1,"raster")
par("mai") #margins in inches
par("usr") 
LCM_b_netto <- LCM_b_1[128:1789,146:2396]  #cutting-out net image 
display(LCM_b_netto)
LCM_b_netto2 <- imageData(LCM_b_netto)
mx <- 1.136063 #determined by manual measurements 
my <- 1.136444 #determined by manual measurements
m1 = matrix(c(mx, 0, 0, 0, my, 0), nrow=3)
m1
LCM_b_2 <- EBImage::affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557))
display(LCM_b_2,"raster")
setwd(home_dir)
writeImage(LCM_b_2,"./data/ISPRS7/images/LCM_cart_enh_b3_scaled_2.jpg") #scaled affine
#

##ISPRS1
setwd(home_dir)
LCM_b_1 <- readImage("./data/ISPRS1/images/LCM_cart_enh_b3.jpg")
display(LCM_b_1,"raster")
par("mai") #margins in inches
par("usr") 
LCM_b_netto <- LCM_b_1[130:1818,147:2407]  #cutting-out net image 
display(LCM_b_netto)
LCM_b_netto2 <- imageData(LCM_b_netto)
mx <- 1.70739 #determined by manual measurements 
my <- 1.708185 #determined by manual measurements
m1 = matrix(c(mx, 0, 0, 0, my, 0), nrow=3)
m1
LCM_b_2 <- EBImage::affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1919,2569))
display(LCM_b_2,"raster")
setwd(home_dir)
writeImage(LCM_b_2,"./data/ISPRS1/images/LCM_cart_enh_b3_scaled_2.jpg") #scaled affine
#end ISPRS1

#end of 5.affine transformation without measuring
################################################################################

##end of 'support_enhance_image.R'


