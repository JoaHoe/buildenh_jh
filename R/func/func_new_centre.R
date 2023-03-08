##name of script: func_new_centre.R
cat("version_number= ",v_nr,"\n")
#purpose: digitize new center in plot of single object
#function used in: script 'support_sequence_of_lines.R'
#author: Joachim Höhle
#GNU General Public License (GPL)

new_centre <- function() {
  cat("digitize new center in plot of single object","\n")
  pos <- locator(n=1, type="n")
  xc_new <- pos$x
  yc_new <- pos$y
  coo <<- c(xc_new, yc_new)
  points(xc_new, yc_new, asp=1, pch=3, cex=3,col="green")
  cat("new coordinates:", "xc = ",  coo[1], "yc = ", coo[2], "\n")
  return(coo)
} #end of function new_centre

#end of script 'func_new_centre.R'


