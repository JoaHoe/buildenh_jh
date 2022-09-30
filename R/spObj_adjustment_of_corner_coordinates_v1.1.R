##name of program: spObj_adjustment_of_corner_coordinates.R
cat("version_number= ",v_nr,"\n")


##buildings of orthoimage ISPRS_#7

if (Img_name == "ISPRS7") {  

  #b5
  if (bnr2 == 5 && p_pos == "cor_adj_coco") { #b5 (with nonortholines)
    #loop
    B8S <- B8
    
    for (i in z) {
      
      if (B8$ortho[i] == 1) { #buildings with orthogonal lines
        B8S$theta_adj[i] <- theta_av
      } else {
        B8S$theta_adj[i] <- theta_av2
      }
      
    } # end loop
  } #end bnr2=5

  #b18
  if (bnr2 == 18 && p_pos == "cor_adj_coco") { #b18 (with nonortholines)
    #loop
    B8S <- B8
    
    for (i in z) {
      if (B8$ortho[i] == 1) { #buildings with orthogonal lines
        B8S$theta_adj[i] <- theta_av
      } else {
        B8S$theta_adj[i] <- theta_av2
      }
      
    } # end loop
    
    for (i in z) {
      
      if (B8S$ortho[i] == 1 && B8$theta_ang[i] > 90) {
        B8S$theta_adj[i] <- B8S$theta_adj[i] + 90 
      }
        
    } # end loop
    
    B8S
    n_pts <- nrow(B8S)
  
    if (sum(B8S$ortho) < n_pts) { 
      cas <- "100_all+nonortho" #case for objects with non-orthogonal lines
    }
    B8 <- B8S
    
  } #end b18
  
} #end of ISPRS7
##############################################################

##buildings of orthoimage ISPRS_#1

if (Img_name == "ISPRS1") { 
  #no corrections
} #end of ISPRS1

##end of script 'spObj_adjustment_of_corner_coordinates.R' 