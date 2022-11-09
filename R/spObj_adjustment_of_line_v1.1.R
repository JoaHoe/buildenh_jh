##name of script: spObj_adjustment_of_line.R
#GNU General Public License (GPL)
cat("version_number= ",v_nr,"\n")


##buildings of orthoimage ISPRS_#7

if (Img_name == "ISPRS7") {  
  
  #b4
  
  if (bnr2 == 4 && p_pos == "cor_adj_line") {
    B6$theta_adj[3] <- B6$theta_ang[3] #correction of theta
    B6$ro_adj[3] <- B6$ro_pixel[3] #correction of ro
  }

  #b5
  
  if (bnr2 == 5 && p_pos == "cor_adj_line") {   
    B6$ro_adj[3] <- (-B6$ro_adj[3]) #correction of ro
    B6$theta_adj[3] <- 180 + B6$theta_adj[3] #angle must be positive
  } 

  #b24
  if (bnr2 == 24 && p_pos == "cor_adj_line") {   
    B6$ro_adj[5] <- B6$ro_pixel[5] #correction of ro (line315)
    B6$theta_adj[5] <- B6$theta_ang[5] #correction of theta (line315)
  } 

} #end of Img_name = "ISPRS7"


##buildings of orthoimage ISPRS_#1

if (Img_name == "ISPRS1") {  
 #no corrections 
} #end of Img_name = "ISPRS1"

##end of script 'spObj_adjustment_of_line.R' 


