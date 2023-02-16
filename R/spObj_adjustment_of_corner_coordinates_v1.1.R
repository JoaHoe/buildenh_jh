##name of script: spObj_adjustment_of_corner_coordinates.R
cat("version_number= ",v_nr,"\n")
#author: Joachim Höhle
#GNU General Public License (GPL)

##objects (buildings) of orthoimage ISPRS_#7

if (Img_name == "ISPRS7") {  

  #b5
  if (bnr2 == 5 && p_pos == "cor_adj_coco") { #b5 (with non-ortholines)
    #loop
    B8S <- B8
    
    for (i in z) {
      
      if (B8$ortho[i] == 1) { #label of orthogonal lines
        B8S$theta_adj[i] <- theta_av
      } else {
        B8S$theta_adj[i] <- heta_av2
      } #end if
      
    } # end loop
    
  } #end bnr2=5

  #b18
  if (bnr2 == 18 && p_pos == "cor_adj_coco") { #b18 (with nonortholines)
    
    B8S <- B8
    #loop
    
    for (i in z) {
      
      if (B8$ortho[i] == 1) { #label of orthogonal lines
        B8S$theta_adj[i] <- theta_av
      } else {
        B8S$theta_adj[i] <- theta_av2
      }
      
    } # end loop
    
    for (i in z) {
      
      if (B8S$ortho[i] == 1 && B8$theta_ang[i] > 90) {
        B8S$theta_adj[i] <- B8S$theta_adj[i] + 90 
      }
        
    } # end for-loop
    
    B8S
    n_pts <- nrow(B8S)
  
    if (sum(B8S$ortho) < n_pts) { 
      cas <- "100_all+nonortho" #object with non-orthogonal lines
    }
    
    B8 <- B8S
    
  } #end b18
  
} #end of ISPRS7
################################################################################


##buildings of orthoimage ISPRS1

if (Img_name == "ISPRS1") { 
  
  #b36
  #correction of non-ortholines
  
  if (bnr2 == 36 && p_pos == "cor_adj_coco") { 
    
    B8
    B8S <- B8
    n_B8S <- nrow(B8S)
    z <- 1 : n_B8S
    
    #loop
    for (i in z) {
      
      if (B8S$ortho[i] == 0) { #buildings with non-orthogonal lines
        B8S$theta_adj[i] <- theta_av2
      }

    } # end of for-loop
    
    B8S
    for (i in z) {
      
      if (B8S$ortho[i] == 0 && B8S$theta_ang[i] > 90) {
        B8S$theta_adj[i] <- B8S$theta_adj[i] + 90 
      }
      
    } # end for-loop
    
  } #end b36
  
  
  #b372
  if (bnr2 == 372 && p_pos == "cor_adj_coco") { #b372 (with non-ortholines)
    
    z <- 1 : nrow(B8S)
    
    for (n2 in z) {
      
      if (B8S$ortho[n2] == 1) { #B8S$ortho[n2]==0 (if other line than theta_ref) 
        cat("n2= ", n2, "\n")
        B8S$theta_adj[n2] <- theta_av # theta_av=theta_ref_adj
      } #end if
      
      B8S  
    } #end for-loop
    
  } #end b372
 
   
  #b45
  
  if (bnr2 == 45 && p_pos == "cor_adj_coco") { #b45 (object with non-ortholines)
    B8 <- B8 #dummy line (no correction)
  } #end b45
  
} #end of orthoimage ISPRS1

##end of script' spObj_adjustment_of_corner_coordinates.R' 