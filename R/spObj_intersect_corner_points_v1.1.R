##name of script 'spObj_intersect_corner_coordinates.R'
cat("version_number= ",v_nr,"\n")
#GNU General Public License (GPL)
#purpose: average angle for objects with non-orthogonal lines 

## ISPRS_#7
if (Img_name == "ISPRS7") { 
  #no corrections
} #end ISPRS7
#

if (Img_name == "ISPRS1") { 
  
  #b36
  if (bnr2==36 && p_pos == "cor_theta_av2") { 
    theta_av2_mod <- 71.4769 #manual calculation by support_intersect_corner_points #4
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
  
  #b45
  if (bnr2==45 && p_pos == "cor_theta_av") {
    theta_av_mod <- theta_av_mod
  } #end if
  
  

} #end ISPRS1
#######################################################
## end of spObj_intersect_corner_points.R
