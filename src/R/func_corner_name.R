##name of script: func_corner_name.R
cat("version_number= ",v_nr,"\n")
#purpose: gives name/label to corners (vertex) of object
#function used in script 'sequence of lines'
#author: Joachim Höhle
#GNU General Public License (GPL)

corner_name <- function(nPC1,nPC2){
  name_of_corner <- paste("P",nPC1,"_","P",nPC2,sep="")
  return(name_of_corner)
} #end of function 'corner_name'

#end of script 'func_corner_name.R'
