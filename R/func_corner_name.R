##name of script: func_corner_name.R
#purpose: gives name/label to corners (vertex) of object
#function used in script 'sequence of lines'

corner_name <- function(nPC1,nPC2){
  name_of_corner <- paste("P",nPC1,"_","P",nPC2,sep="")
  return(name_of_corner)
} #end of function 'corner_name'

#end of script 'func_corner_name.R'
