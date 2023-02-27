##name of script: 'spObj_intersect_corner_coordinates.R'
cat("version_number= ",v_nr,"\n")
#purpose: average angle for objects with non-orthogonal lines for special objects
#author: Joachim Höhle
#GNU General Public License (GPL)

## ISPRS_#7

if (Img_name == "ISPRS7") { 
  #no corrections
} #end ISPRS7

################################################################################

if (Img_name == "ISPRS1") { 
  
  #b36
  if (bnr2==36 && p_pos == "cor_theta_av2") { 
    theta_av2_mod <- 71.4769 #manual calculation by script 'support_intersect_corner_points', script #4
  } #end if
  
  #b372
  
  if (bnr2==372 && p_pos == "cor_theta_av") {
    theta_av_mod <- theta_ref_adj
    theta_av_mod
  } #end if
  
  if (bnr2==372 && p_pos == "cor_theta_av2") {
    cat("i2 =", i2,"\n") 
    theta_av2_mod <- B6_seq$theta_adj[i2]
  } #end if
  
  #end b372
  
  #b39
  
  if (bnr2 == 39 && p_pos == "cor_corner_pts") { #object type:"extr_wd"
    b0[5,1] <- 1.71790e+03 #correction due to tan(90)
  } #end if
  
  #end b39
  
  
  #b45
  
  if (bnr2==45 && p_pos == "cor_theta_av") {
    theta_av_mod <- theta_av_mod
  } #end if
  
  #end b45
  
} #end ISPRS1

################################################################################
## end of spObj_intersect_corner_points.R
