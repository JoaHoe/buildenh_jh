## name of script: func_det_of_angle.R
# purpose: determination of angle at corner points
# function used in: spObj_sequence_of_lines.R & support_sequence_of_lines.R

det_of_angle <- function(corner_x,corner_y) { 
  alph <- (atan2(-(corner_y-yc),(corner_x-xc)))*omega # -dy/dx because alpha 
    #in math-coor-system
  
  if(alph < 0) {
    alph3 <- alph + 360
  } else {
    alph3 <- alph
  } #end if-else
  
  angle_center_corner <- alph3
  return(angle_center_corner)
} #end of function 'det_of_angle'

#end of script "func_det_of_angle.R"