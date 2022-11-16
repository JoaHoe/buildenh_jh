#name of script: func_reduce_pointset.R
#purpose: eliminating of pixels which do not belong to PC
#function used in: script 'sequence_of_lines.R'
##GNU General Public License (GPL)

reduce_pointset <- function(P) {
  xgaps <- test_gaps_onedirection(hist(P$x,plot=FALSE)) #modified
  ygaps <- test_gaps_onedirection(hist(P$y,plot=FALSE)) #modified
  
  if (!is.null(xgaps) | !is.null(ygaps)) {

    if (!is.null(xgaps) & !is.null(ygaps)) {
      if (xgaps$nGaps > ygaps$nGaps) { ygaps <- NULL} else { xgaps <- NULL}
    }
    
    if (is.null(xgaps)) {
      P <- if (ygaps$direction == "lower") subset(P, P$y <= ygaps$cut) else subset(P, P$y > ygaps$cut)
    }
    
    if (is.null(ygaps)) {
      P <- if (xgaps$direction == "lower") subset(P, P$x <= xgaps$cut) else subset(P, P$x > xgaps$cut)
    }
  }
  return(P)
} #end of function "reduce_pointset(P)"

#end of script 'func_reduce_pointset.R'
