##name of program: startup_buildenh.R
#description: program(script) starts the package 'buildenh'
v_nr = "1.1" #version number of the program
cat("version_number= ",v_nr,"\n")
#examples: extracted buildings from land cover maps derived by classification programs
#data: ISPRS test "Vaihingen": orthoimages of area #7, #1
#author: Joachim Hoehle
#instructions: change directories for input;
#instructions: input project title and image name
#instructions: save your home directory
#instructions: type 'Ctrl+A'(select all) and 'Source'
#instructions: new users may start by examples (processing mode = demo)
#depends: R-4.2.1
#Copyright(C) 2022 Joachim Hoehle
#GNU General Public License (GPL)
###################################################################################
cat("start of software package 'buildenh_jh' ","\n")

cat("first program/script 'startup_buildenh.R' ","\n")

#save your home directory
old_dir <- setwd("./")
getwd()
#
#home_dir <- "C:/Users/Joachim/OneDrive/Documents/GitHub/buildenh"
#home_dir2 <- "C:/Users/Joachim/OneDrive/Documents/GitHub/buildenh/R"
home_dir <- "C:/Users/Joachim/R_programs/buildenh_jh/clone1"
home_dir2 <- "C:/Users/Joachim/R_programs/buildenh_jh/clone1/R"

###################################################################################

## title of project (activate manually)

setwd(home_dir)
#prj_title <- "ISPRS7_LCM1" #orthoimage ISPRS#7, 
#classification method: DT/LCM1 by 17 attributes
#training by orthoimage #26
prj_title <- "ISPRS1_LCM2" #orthoimage ISPRS#1, classification method: DT/LCM2 by 5 attributes
cat("project title is = ", prj_title,"\n")

#select orthoimage (activate manually)
#Img_name <- "ISPRS7" #name of orthoimage to be processed, change for other image
Img_name <- "ISPRS1" #name of orthoimage to be processed, change for other image

if (Img_name == "ISPRS7") {
##setting of path- & file-name for original data:
  OrgClassResFilename <- "ISPRS_#7_b.tiff" #extracted buildings
  OrgClassResPathname <- paste(home_dir,"/data",sep = "")
  OrgImgPathname <- paste(home_dir,"/data",sep = "")
  OrgImgFilename <- "top_mosaic_09cm_area7.tif"  #GSD=0.09m
  OrgGtsPathname <- paste(home_dir,"/data",sep = "")
  OrgGtsFilename <- "gts_top_mosaic_09cm_area7.tif" #GSD=0.09m
  #GSD=Ground Sampling Distance
} #end of image7

if (Img_name == "ISPRS1") {
  ##setting of path- & file-name for original data:
  setwd(home_dir)
  OrgClassResPathname <- paste(home_dir,"/data",sep = "")
  OrgClassResFilename <- "ISPRS_#1_b.tiff" #extracted buildings
  OrgImgPathname <- paste(home_dir,"/data",sep = "")
  OrgImgFilename <- "top_mosaic_09cm_area1.tif" #GSD=0.09m
  OrgGtsPathname <- paste(home_dir,"/data",sep = "") 
  OrgGtsFilename <- "gts_top_mosaic_09cm_area1.tif" #GSD=0.09m
  #GSD=Ground Sampling Distance
} #end of image1

proc_mode <- "NA" #mode of processing

## install packages
# install.packages('EBImage')
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("EBImage")
#
# install.packages('roxygen2')
# install.packages('devtools')
# install.packages('testthat')
# install.packages('knitr')
# install.packages('spatstat')
# install.packages('tiff')
# install.packages('stringr')
# install.packages('usethis')
#

##loading of libraries 
setwd(home_dir2)
source("func_loadLib_op.R") #load of other R-packages, v1.1
source("func_loadLib_jh.R") #load of functions for the R-package 'buildenh', v1.1
#
loadLib_op() #call of function
loadLib_jh() #call of function

#other functions
display = function(...) if (interactive()) EBImage::display(...)

#setup for processing mode "auto"

if (Img_name == "ISPRS1") {
  y_auto <- c(4,5,7,18) #file for automatic processing (orthoimage #7)  
}

if (Img_name == "ISPRS7") {
  y_auto <- c(20,22,23,5) #file for automatic processing (orthoimage #7)  
}

n_y_auto <- length(y_auto)
k_y_auto <- 1

#setup of parameter/variables
par(mai = c(1.02,0.82,0.82,0.42)) #setup of margins/plot region [inches]
bnr2_part <- "NA" #partition of object

cat("end of 'startup_buildenh.R' - continue with 'extract_all_buildings.R' ", "\n")
#end of 'startup_buildenh.R'

##start the next program ("enhance_image.R")
setwd(home_dir2)
source(paste("enhance_image_v",v_nr,".R",sep=""))
#

