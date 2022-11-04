## name of script: support_enhance_image_v1.1.R
## supporting software for program "enhance_image.R"
## GNU General Public License (GPL)
cat("version_number= ",v_nr,"\n")
## instruction: run all programs in "demo" mode before using 'support_enhance_image.R'

##contents:
# -checking of scale
# -scaling of image (simple and affine)
# -cutting out of net-image
# -affine scaling with manual measurement in browser
# -affine transformation without measuring

###################################################################################

##checking of scale

#orthoimage ISPRS7
setwd(home_dir)
#LCM_b_1 <-  readImage("./data/ISPRS7/images/LCM_cart_enh_b3.jpg") #land cover map 1
LCM_b_1 <-  readImage(paste("./data/", Img_name,"/images/LCM_cart_enh_b3.jpg",sep = "")) #land cover map
display(LCM_b_1, method = "raster")
colorMode(LCM_b_1)
LCM_b_2 <- imageData(LCM_b_1)
co_xy <- locator(4,"p") #digitize of sides: left (1), right (2), upper (3), lower (4)
co_xy #check
dcol <- (co_xy$x[2]-co_xy$x[1])
drow <- (co_xy$y[4]-co_xy$y[3])

# if (abs(round(dcol - 1887)) < 3 && abs(round(drow - 2557)) < 3) { #orthoimage #7
#   cat("scaling is OK")
# } else {
#   cat("x_deviation=",abs(round(dcol-1887)),"\n")
#   cat("y_deviation=",abs(round(dcol-2557)),"\n")
#   stop("decide if scaling must be repeated")
# }

if (abs(round(dcol - 1919)) < 3 && abs(round(drow - 2569)) < 3) { #orthoimage #1
  cat("scaling is OK")
} else {
  cat("x_deviation=",abs(round(dcol-1919)),"\n")
  cat("y_deviation=",abs(round(dcol-2569)),"\n")
  stop("decide if scaling must be repeated")
}

#############################################################################################################

##scaling of image

#orthoimage ISPRS7

#measurement in raster image
LCM_b=readImage("./data/ISPRS7/LCM_cart_enh_b3.jpg")  #not-to-scale image
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
########################################################################################################

## cutting out of net-image

LCM_b=readImage("./data/ISPRS7/LCM_cart_enh_b3.jpg")
#display(LCM_b, method="browser") #find coordinates of corners
display(LCM_b, method="raster") #find coordinates of corners
co_xy <- locator(4) #digitize of side-lines in raster image
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
###############################################################################################

##affine scaling with manual measurement in browser
setwd(home_dir)
#LCM_b=readImage("./data/ISPRS7/ISPRS_#7_b.tiff")
#LCM_b=readImage("./data/ISPRS1/images/ISPRS_#1_b.tiff")
LCM_b <- readImage(paste("./data/", Img_name,"/images/LCM_cart_enh_b3.jpg",sep = "")) #land cover map
display(LCM_b, method = "browser")

#measurement of sides
co_xy <- matrix(rep(0,8),nrow = 4)
co_xy <- data.frame(co_xy)
names(co_xy) <- c("x","y")
co_xy

#ortho image ISPRS7
# co_xy$x[1] <- 128 #to be checked and modified
# co_xy$x[2] <- 1789
# co_xy$y[3] <- 146
# co_xy$y[4] <- 2396

#determination of scaling factors
# mx <- 1887/(co_xy$x[2]-co_xy$x[1]) #image 7
# my <- 2557/(co_xy$y[4]-co_xy$y[3]) #image 7
# m1 <- matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
# m1

#affine transformation
#LCM_b_2 <- affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557)) #orthoimage7
#display(LCM_b_2) #display on browser

################################

##image ISPRS1
co_xy$x[1] <- 130 #to be checked and modified
co_xy$x[2] <- 1819
co_xy$y[3] <- 147
co_xy$y[4] <- 2408

#determination of scaling factors
mx <- 1919/(co_xy$x[2]-co_xy$x[1]) #orthoimage ISPRS1
my <- 2569/(co_xy$y[4]-co_xy$y[3]) #orthoimage ISPRS1
m1 <- matrix(c(mx, 0, 0, 0,my, 0), nrow=3)
m1
##############################

##output
# setwd(home_dir)
# writeImage(LCM_b_2,paste("./data/",Img_name,"/images/LCM_cart_enh_b3_affin_new.jpg",sep = ""))
# 
# 
##display
# LCM_b_3 <- readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3_affin_new.jpg",sep = ""))
# display(LCM_b_3)
#display(LCM_b_3,"raster")
##########################################################################################################

##affine transformation without measuring
setwd(home_dir)
LCM_b_1 <- readImage("./data/ISPRS7/images/LCM_cart_enh_b3.jpg")
display(LCM_b_1,"raster")
par("mai") #margins in inches
par("usr") 
LCM_b_netto <- LCM_b_1[128:1789,146:2396]  #cut out net image ('ISPRS#7)
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
#end of affine transformation without measuring

##end of 'support_enhance_image.R'


