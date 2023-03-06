## name of program: spObj_enhance_image.R
cat("version_number= ",v_nr,"\n")
# purpose: scaling of image 
# ISPRS data "Vaihingen", area#7, area#1
# author: Joachim Höhle
# GNU General Public License (GPL)

## ISPRS_#7

if (Img_name == "ISPRS7") { 
  #affine transformation without measuring (applied)
  setwd(home_dir)
  LCM_b_1 <- readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3.jpg",sep=""))
  #par("mai") #margins in inches
  #par("usr") #users plotting region (default: x1=0, x2=1, y1=0, y2=1)
  LCM_b_netto <- LCM_b_1[128:1789,146:2396]  #cut out of net-image 
  display(LCM_b_netto)
  LCM_b_netto2 <- imageData(LCM_b_netto)
  mx <- 1.136063 #determined by manual measurements in 'support_enhance_image.R' 
  my <- 1.136444 #determined by manual measurements in 'support_enhance_image.R' 
  m1 = matrix(c(mx, 0, 0, 0, my, 0), nrow=3)
  LCM_b_2 <- EBImage::affine(LCM_b_netto, m1, filter='bilinear', output.dim=c(1887,2557))
  display(LCM_b_2,"raster")
  setwd(home_dir)
  writeImage(LCM_b_2,paste("./data/",Img_name,"/images/LCM_cart_enh_b3_scaled_2.jpg",sep = "")) #scaled affine
} #end of scaling orthoimage ISPRS7
#


## ISPRS_#1 

if (Img_name == "ISPRS1") { 
  #affine transformation without measuring 
  setwd(home_dir)
  LCM_b_1 <- readImage(paste("./data/",Img_name,"/images/LCM_cart_enh_b3.jpg",sep = ""))
  #par("mai") #margins in inches
  #par("usr") #users plotting region (default: x1=0, x2=1, y1=0, y2=1)
  display(LCM_b_1)
  LCM_b_netto <- LCM_b_1[130:1819,147:2408]  #cut out of net-image
  display(LCM_b_netto)
  mx <- 1.136175 #determined by manual measurements 
  my <- 1.136223 #determined by manual measurements
  m1 = matrix(c(mx, 0, 0, 0, my, 0), nrow=3)
  LCM_b_2 <- EBImage::affine(LCM_b_netto, m1, filter=c("bilinear","none"), output.dim=c(1919,2569))
  display(LCM_b_2,"raster")
  #display(LCM_b_2)
  setwd(home_dir)
  writeImage(LCM_b_2,paste("./data/",Img_name,"/images/LCM_cart_enh_b3_scaled_2.jpg",sep = "")) #scaled affine
} #end of scaling orthoimage ISPRS_#1

#end of script 'spObj_enhance_image.R'

