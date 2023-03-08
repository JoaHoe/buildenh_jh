## name of script: func_loadLib_op.R
cat("version_number= ",v_nr,"\n")
## purpose: required open source R-packages 
#author: Joachim Höhle
## GNU General Public License (GPL)

loadLib_op <- function() { #required R-packages
  library("EBImage")
  library("spatstat")
  library("tiff")
  library("rpart")
  library("nlme")
} #end of function 'loadLib_op()'

#end of script 'func_loadLib_op.R