## name of script: func_loadLib_op.R
## purpose: required open source R-packages 

loadLib_op <- function() { #required R-packages
  library("devtools")
  library("EBImage")
  library("knitr")
  library("roxygen2")
  library("spatstat")
  library("stringr")
  library("testthat")
  library("tiff")
  library("usethis")
} #end of function 'loadLib_op()'

#end of script 'func_loadLib_op.R