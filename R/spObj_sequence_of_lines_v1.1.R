##name of script: spObj_sequence_of_lines.R
cat("version_number= ",v_nr,"\n")
#GNU General Public License (GPL)
#purpose: dealing with special objects 
#correction of position of line-centers, angle of line and sequence of lines 
#instruction: if corrections in the positions are required -> use interactive detection of lines in 'support_line_detection.R' 
#orthoimage ISPRS_#7, #1


##orthoimage ISPRS_#7

if (Img_name == "ISPRS7") { 

  #b4
  if (bnr2 == 4 && p_pos == "cor_pos") {  
      b13_angle_df2 <- b13_angle_df
      b13_angle_df2[8,3:4] <- c(1360.9,181.8)
      b13_angle_df2[11,3:4] <- c(1355.1,211.2)
      b13_angle_df2[12,3:4] <- c(1509.9,132.6)
      b13_angle_df2
  } #end of b4
  
  #b5
  if (bnr2 == 5 && p_pos == "cor_sek") {  
      sequence_seg <- c(7,23,98,264,70,49,34) # when selected "manually"
  } #end of b5
  
  #b8
  if (bnr2 == 8 && p_pos == "cor_pos") { 
      b13_angle_df2
      b13_angle_df2[10,3:4] <- c(1157,532)
    } #end of b8
  
  #b16
  if (bnr2 == 16 && p_pos == "cor_pos") { 
    
    #correction of position of line-center
    b13_angle_df2
    b13_angle_df2[8,3:4] <- c(368,956) #measured manually in 'browser'
    b13_angle_df2
    #
    
    #correction of angle
    
    #center of object/building
    xc <- plotPar[1]
    yc <- plotPar[2]
    x_centre <- b13_angle_df2[8,3] 
    y_centre <- b13_angle_df2[8,4]
    
    #correction of angle for new midpoint
    alpha <- det_of_angle(x_centre, y_centre) #call of function
    b13_angle_df2$alpha[8] <- alpha #manual correction
    b13_angle_df2
  } #end of b16
  
  #b17
  if (bnr2 == 17 && p_pos == "cor_pos") {  
    b13_angle_df2 <- b13_angle_df
    b13_angle_df2[7,3:4] <- c(839.2,976.3) #measured manually in 'browser'
    b13_angle_df2[6,3:4] <- c(914.2,903.7)
    b13_angle_df2[5,3:4] <- c(877.1,980.3)
    b13_angle_df2
  } #end of b17
  
  #b18
  if (bnr2 == 18 && p_pos == "cor_pos") {  
    b13_angle_df[6,3:4] <- c(953.2,921.5)
    b13_angle_df[8,3:4] <- c(1017,1084)
    b13_angle_df[9,3:4] <- c(1131.6,1044.1)
    b13_angle_df[10,3:4] <- c(789.8,1144.4)
    b13_angle_df2 <- b13_angle_df
  } #end of b18
  
  #b23
  if (bnr2 == 23 && p_pos == "cor_det") { 
    lnr_det3[6] <- c(36) 
    lnr_det5 <- lnr_det3
    n_PCR <- length(lnr_det5)
    B5_6R <- B5_6
    B5_6R <- B5_6R[1:n_PCR,]
    B5_6R$lnr[1:n_PCR] <- lnr_det5 
    row.names(B5_6R) <- 1:n_PCR
    B5_6R
    B5_4
    n_B5_4 <- length(B5_4$lnr)
    
    # loop
    i <- 1
    
    while (i <= n_B5_4) {
      if (B5_4$lnr[i] == lnr_det3[6]) { 
        B5_6R[6,2:7] <- B5_4[i,2:7]
        B5_6R[6,8] <- 1
      } #end if 
      i <- i + 1 
    } #end of loop
    
    B5_6R4 <- B5_6R
  } #end of b23
  
  #b24
  
  if (bnr2 == 24 && p_pos == "cor_pos") { 
    b13_angle_df
    b13_angle_df[1,3:4] <- c(1346.5,1293.0)
    b13_angle_df[5,3:4] <- c(1399.9,1458.2)
    b13_angle_df[6,3:4] <- c(1291.2,1422.2)
    b13_angle_df2 <- b13_angle_df
    b13_angle_df2
  } #end of b24
  
  ##b27
  #first part
  
  if (bnr2 == 27 && p_pos == "cor_sep") { 
    #removal of points of the other part by separating the line
    B5_6
    theta_ref <- B5_6$theta_angle[13]
    alph_ref <- theta_ref - 90
    B5_6R2 <- B5_6
    B5_6R2[,] <- 0
    n_PCR2 <- length(B5_6R2$lnr)
    vec_1 <- 1 : n_PCR2
    k2 = 1
    #loop
    for (n in vec_1) {
      if (B5_6$theta_angle[n] == alph_ref && B5_6$ro_pixel[n] <= 1430) {
        B5_6R2[k2,] <- B5_6[n,]
        k2 <- k2+1
      } #end if
    } #end for-loop
    
    B5_6R2 <- B5_6R2[1:(k2-1),]
    B5_6R3 <- B5_6R2
    n_B6R3 <- n_PCR2-(k2-1)
    B5_6R3 <- B5_6R3[1:(n_B6R3),] 
    B5_6R3[,] <- 0
    B5_6R3
    #
    k=1
    for (n in vec_1) {
      if (B5_6$theta_angle[n] == theta_ref) {
        B5_6R3[k,] <- B5_6[n,]
        k <- k+1
      } #end if
    } #end for-loop
    
    B5_6R3[15:19,] <- B5_6R2[1:5,] 
    B5_6R3 <- B5_6R3[1:19,]
    row.names(B5_6R3) <- 1:19
    B5_6R3
    B5_6R4 <- B5_6R3[-c(1,3:5,7:8,10:12,14:15),]
    B5_6R <- B5_6R4
    B5_6R
  } #end b27, first part
  
  #second part
  
  if (bnr2 == 27 && p_pos == "cor_sep") { 
    #separation of parts
    B5_6
    theta_ref
    B5_6R2 <- B5_6
    B5_6R2[,] <- 0
    n_PCR2 <- length(B5_6R2$lnr)
    vec_1 <- 1 : n_PCR2
    k2 = 1
    for (n in vec_1) {
      if (B5_6$theta_angle[n] == theta_ref && B5_6$ro_pixel[n] >= 1325) {
        B5_6R2[k2,] <- B5_6[n,]
        k2 <- k2+1
      } #end if
    } #end for-loop
    B5_6R2
    k3 <- k2-1
    B5_6R2 <- B5_6R2[1:k3,]
    B5_6R3 <- B5_6R2
    n_B6R3 <- n_PCR2 - k3
    B5_6R3 <- B5_6R3[1 : n_B6R3,] 
    row.names(B5_6R3) <-1: n_B6R3
    #
    B5_6R3[,] <- 0
    k=1
    for (n in vec_1) {
      if (B5_6$theta_angle[n] == alph_ref) {
        B5_6R3[k,] <- B5_6[n,]
        k <- k+1
      } #end if
    } #end for-loop
    
    B5_6R3[14:21,] <- B5_6R2[1:8,]
    B5_6R4 <- B5_6R3[-c(2,4:10,13,15,19),]
    row.names(B5_6R4) <- 1:10
    lnr_det3 <- B5_6R4$lnr
    
    #supplement of double points
    B5_6R4[11,] <- B5_6R4[2,]
    B5_6R4[12,] <- B5_6R4[3,]
    lnr_det3 <- B5_6R4$lnr
    
    #solution by manual operation 
    lnr_det5 <- lnr_det3
    lnr_det5
  } #end if (b27,p_pos="cor_sep)
  
  ##corrections for postion of line-midpoints
  
  if (bnr2 == 271 && p_pos == "cor_pos") { #b27, first part 
    b13_angle_df[8,3:4] <- c(588.2,1818.5)
    b13_angle_df2 <- b13_angle_df
  }
  
  if (bnr2 == 272 && p_pos == "cor_pos") { #b27, second part
    #correction of position
    b13_angle_df
    b13_angle_df[2,3:4] <- c(739.7,1636.5)
    b13_angle_df[5,3:4] <- c(824.6,1691.0)
    b13_angle_df[8,3:4] <- c(732.1,1704.0)
    b13_angle_df[9,3:4] <- c(612.2,1911.5)
    b13_angle_df[11,3:4] <- c(787.1,1728.8)
    b13_angle_df[12,3:4] <- c(703.3,1948.8)
    b13_angle_df2 <- b13_angle_df
  } #end b272 / "cor_pos"
  
  
  #b30
  if (bnr2 == 30 && p_pos == "cor_sek") {
    sequence_seg2 <- sequence_seg
    sequence_seg2 <- c(1,5,16,3,98)
    sequence_seg2
  } #end b30
  
  #b32
  if (bnr2 == 321 && p_pos == "cor_pos") {
    b13_angle_df2[8,3:4] <- c(293,2042) #new position to be measured in system b32/viewer
  } #end b321
  
  #b34
  if (bnr2 == 34 && p_pos == "cor_pos") {  
    b13_angle_df[3,3:4] <- c(1796,2288)
    b13_angle_df[5,3:4] <- c(1692,2071)
  } #end b34

} #end of ISPRS7
##########################################################

##orthoimage ISPRS_#1

if (Img_name == "ISPRS1") {
  
  #b22
  if (bnr2 == 22 && p_pos == "cor_pos") {  
    #determine new position by 'support_line_detection.R', #6
    #or with angle by 'support_sequence_of_lines.R', #9
    b13_angle_df[6,3:4] <- c(1487, 1423) #manually derived
    b13_angle_df2 <- b13_angle_df
  } #end b22
  
  #b271
  if (bnr2 == 271 && p_pos == "cor_pos") {  
    #determine new position by 'support_line_detection.R', #6
    #or with angle by 'support_sequence_of_lines.R', #9
    b13_angle_df
    b13_angle_df[1,3:4] <- c(582,1672)
    b13_angle_df[4,3:4] <- c(571,1721) #manually derived
    b13_angle_df[5,3:4] <- c(456,2199) #manually derived
    b13_angle_df[6,3:4] <- c(546,1808) #manually derived
    b13_angle_df[10,3:4] <- c(651,2006) #manually derived
    b13_angle_df[12,3:4] <- c(626,1840) #manually derived
    b13_angle_df[13,3:4] <- c(633,1749) #manually derived
    b13_angle_df2 <- b13_angle_df
    b13_angle_df2
  } #end b271
  
} #end of ISPRS1

##end of script 'spObj_sequence_of_lines.R'
################################################################################
