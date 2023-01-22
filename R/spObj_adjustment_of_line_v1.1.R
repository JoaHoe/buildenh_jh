##name of script: spObj_adjustment_of_line.R
#GNU General Public License (GPL)
cat("start of spObj_adjustment of line ","\n")
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
    B6 <- B6 #dummy line
  } #end b271
  
  #b36
  if (bnr2==36 && p_pos == "cor_adj_line") {
    i # index of PC, which must be corrected
    cat("PC_nr= ",PC_numb[i],"\n") #PC_nr, which must be corrected
    res1
    residuals <- matrix(nrow=length(res1), ncol=2)
    residuals[,1] <- 1 : length(res1)
    residuals[,2] <- round(t(res1),digits=1)
    residuals <- as.data.frame(residuals)
    #
    cat("PC_number= ",PC_numb[i],"\n")
    nrow(all_PC[[i]])
    vec <- 1 : nrow(residuals)
    thr_res #3xthr_line
    
    #remove gross errors
    nrs_errors <- rep(0,nrow(residuals))
    if (nrow(residuals) == nrow(all_PC[[i]])) {
      n2 <- 1
      for (n1 in vec) { 
        #cat("n1= ",n1,"\n")
        if (residuals$V2[n1] > thr_res) {
          cat("index of wrong pixel= ",n1,"\n")
          nrs_errors[n2] <- residuals$V1[n1]
          n2 <- n2 + 1
        } #end if 
      } #end for
    } #end if
    nrs_errors
    nrs_errors <- nrs_errors[nrs_errors > 0]
    nrs_errors #numbers of errors, to be transfered into correction vector (not done yet!)
    all_PC[[i]] <- all_PC[[i]][-nrs_errors,]
    all_PC[[i]]
    
    #store corrected PC
    setwd(home_dir)
    fname8 <- paste("./data/",Img_name,"/all_PC$PC_nr",PC_numb[i],".txt", sep=(""))
    write.table(all_PC[[i]],fname8)
    
    #start again adjustment of line
    setwd(home_dir2)
    source(paste("adjustment_of_line_v",v_nr,".R",sep="")) 
   } #end of bnr2=36
  
  #b372
  if (bnr2==372 && p_pos == "cor_adj_line") {
   B6$theta_adj[2] <- B6$theta_ang[2]
   B6$ro_adj[2] <- B6$ro_pixel[2]
   B6$theta_adj[3] <- B6$theta_ang[3]
   B6$ro_adj[3] <- B6$ro_pixel[3]
   B6 
  } #end of b372
  
  #b38
  if (bnr2 == 38 && p_pos == "cor_adj_line") {
    B6 <- B6 #dummy line 
  } #end of b38

  #b46
  if (bnr2 == 46 && p_pos == "cor_adj_line") {
    B6 <- B6 #(dummy line) 
  } #end of b46
  
} #end of Img_name = "ISPRS1"

##end of script 'spObj_adjustment_of_line.R' 


