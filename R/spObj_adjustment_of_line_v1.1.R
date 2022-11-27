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
  
  #b21
  if (bnr2 == 21 && p_pos == "cor_adj_line") {   
    B6$ro_adj[2] <- B6$ro_pixel[2] #correction of ro (line 524)
    B6$theta_adj[2] <- B6$theta_ang[2] #correction of theta (line 524)
  } 
  
  #b271
  if (bnr2 == 271 && p_pos == "cor_adj_line") { 
    #B6$theta_adj[4] <- B6$theta_ang[4] #correction of theta (line 982)
    #B6$ro_adj[4] <- B6$ro_pixel[4] #correction of ro (line 982)
    #B6$theta_adj[10] <- B6$theta_ang[10] #correction of theta (line 1053)
    #B6$ro_adj[10] <- B6$ro_pixel[10] #correction of ro (line 1053)
    #B6$theta_adj[11] <- B6$theta_ang[11] #correction of theta (line 283)
    #B6$ro_adj[11] <- B6$ro_pixel[11] #correction of ro (line 283)
  } #end b271
  
} #end of Img_name = "ISPRS1"

##end of script 'spObj_adjustment_of_line.R' 


