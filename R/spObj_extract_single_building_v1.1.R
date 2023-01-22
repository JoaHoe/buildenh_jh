##name of script: spObj_extract_single_building.R
#GNU General Public License (GPL)
cat("version_number= ",v_nr,"\n")
#purpose: partition of special objects 
#instruction: measure 2 pixels which will separate object (zoom:200%)
#             use 'display(is_bnr)'
#####################################################################

if (Img_name == "ISPRS7") {

  #b27
  if (bnr2 == 27 && p_pos == "cor_sep") {
      cat("start of spObj_extract_single_building.R -first part",sep = "")
      display(is_bnr)
      is_bnr2 <- is_bnr
      imageData(is_bnr2)[726:730, 1639:1643] = FALSE
      imageData(is_bnr2)[715:719, 1673:1677] = FALSE
      display(is_bnr2)
      is_bnr2_label <- bwlabel(is_bnr2)
      cat("number of objects= ",max(is_bnr2_label),"\n")
      coor_part <- computeFeatures.moment(is_bnr2_label) #geometric features (moment)
      shap_part <- computeFeatures.shape(is_bnr2_label) #geometric features (shape)
      is_label_1 <- is_bnr2_label@.Data == 1 #left part
      imageData(is_label_1)[720:735, 1630:1648]
      display(is_label_1)
      #
      
      is_label_2 <- is_bnr2_label@.Data == 2 #second part
      display(is_label_2)
      #save original values
      xc_orig <- xc
      yc_orig <- yc
      #
      if (proc_mode == "obj_wise" && part == "2parts_1") {
          bnr2 <- 271
      }
      
      if (proc_mode == "obj_wise" && part == "2parts_2") {
        bnr2 <- 272
      }
  
      if (bnr2 == 271) { 
        is_label_1
        display(is_label_1)
        #storage
        setwd(home_dir)
        f1<-paste("./data/",Img_name,"/param_b_",bnr2,sep="")
        save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
        #
        
        coords <- data.frame(x=as.numeric(row(is_label_1)),
                             y=as.numeric(col(is_label_1)), 
                             is_label_1=as.numeric(is_label_1))
        coords <- coords[coords$is_label_1 == 1,] #removal of pixels 
                    #which do not have the label of the building
        head(coords)
        
        #calculation of new center of object from connected components
        xc <- coor_part[1,1]
        yc <- coor_part[1,2]
        alpha <- coor_part[1,5]*omega
        r_max <- shap_part[1,6]
      
        #storage of plot parameters
        setwd(home_dir)
        plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
        f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
        save(plotPar, file=f1) #parameters xc,yc,r_max,alpha
        
        #plot of PCs and checkpoints
        dev.set(2)
        r_max2 <- round(1.1*r_max)
        plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
             ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
             " - left",sep = ""), axes=TRUE)
        points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #center of PC
        #
        
        #generation of image
        setwd(home_dir)
        f=paste("./data/",Img_name,"/b_nr",sep = "")
        save(bnr2,file=f)
        file1 = paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
        tiff(file1, width=578, height=578, units="px", bg = "white")
        r_max2 <- round(1.1*r_max)
        plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
        points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
        points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        dev.off()
        #
        
        #output of file with PC-coordinates (required in program 'line_detection.R')
        coords2 <- subset(coords, select=c(x,y))
        head(coords2)
        N1 <- length(coords$x)
        x_dat <- rep(0,N1)
        y_dat <- rep(0,N1)
        idx <- rep(0,N1)
        idxy <- cbind(idx, x_dat, y_dat)
        idxy <- coords2 #use of old object name
        rownames(idxy) <- 1 : N1
        fname2 <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv", sep="") 
        write.table(idxy, fname2,  sep= " ", row.names=T) #output of PCs for one building
        #end of output
        
        setwd(home_dir2)
        source(paste("line_detection_v",v_nr,".R", sep=""))
      } #end bnr2 = 271
  
      if (bnr2 == 272) { #second part
        cat("start of spObj_extract_single_building.R -second part",sep = "", "\n")
        bnr2_part <- bnr2 
        is_label_2 <- is_bnr2_label@.Data == 2 
        display(is_label_2)
        is_label_2
        coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
        coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
        xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=272
        yc <- coor_part[2,2]
        alpha <- coor_part[2,5]*omega
        r_max <- shap_part[2,6]
        dy_window_plot
        
        #output plot parameters
        plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
        setwd(home_dir)
        f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
        save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
        f <- paste("./data/",Img_name,"/b_nr",sep = "")
        save(bnr2,file=f)
        
        #plot of PC and checkpoints
        r_max2 <- round(1.1 * r_max)
        plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
             asp=1,xlim=c(xc-r_max2,xc+r_max2),
             ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, 
             main=paste("b", bnr2," - second part",sep = ""), axes=TRUE)
        points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
        
        #generation of image
        file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
        tiff(file1, width=578, height=578, units="px", bg = "white")
             r_max2 <- round(1.1*r_max)
        plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
             asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), 
             xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2," - right",
             sep = ""), axes=TRUE)
        #plot(coords$x, coords$y,pch=16,cex=0.2,col="black",asp=1,xlim=c(1,1887),ylim=c(2557,1), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE) #small scale
        points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
        points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #center of PC
        dev.off()
        #windows() #to establish new window if necessary
        
        #output of file with PC-coordinates (required in program 'line_detection.R')
        coords2 <- subset(coords, select=c(x,y))
        head(coords2)
        N1 <- length(coords$x)
        x_dat <- rep(0,N1)
        y_dat <- rep(0,N1)
        idx <- rep(0,N1)
        idxy <- cbind(idx, x_dat, y_dat)
        idxy <- coords2 #use of old object name
        rownames(idxy) <- 1 : N1
        fname2 <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv", sep="")
        write.table(idxy,fname2,sep= " ",row.names=T) ##output of pixel cluster for one building
        #end of output
        
        cat("end of program 'extract_single_building.R' - 
            continue with 'line_detection.R' ","\n")
        setwd(home_dir2)
        source(paste("line_detection_v",v_nr,".R", sep=""))
      } #end of bnr2=272
  } # end bnr2 = 27
  ############################################
  
  #b32
  if (bnr2 == 32 && p_pos == "cor_sep") {
    cat("start of spObj_extract_single_building.R - 
        first part",sep = "")
    display(is_bnr)
    is_bnr2 <- is_bnr
    imageData(is_bnr2)[462:466, 2284:2288] = FALSE
    imageData(is_bnr2)[494:498, 2230:2234] = FALSE
    display(is_bnr2)
    is_bnr2_label <- bwlabel(is_bnr2)
    cat("Number of objects= ",max(is_bnr2_label),"\n")
    coor_part<-computeFeatures.moment(is_bnr2_label) #geometric features (moment)
    shap_part<-computeFeatures.shape(is_bnr2_label) #geometric features (shape)
    
    is_label_1 <- is_bnr2_label@.Data == 1 #left part
    imageData(is_label_1)[460:493, 2287:2230]
    display(is_label_1)
    is_label_2 <- is_bnr2_label@.Data == 2 #second part
    display(is_label_2)
    
    #save original values
    xc_orig <- xc
    yc_orig <- yc
    
    if (proc_mode == "obj_wise" && part == "2parts_1") {
      bnr2 <- 321
      bnr2_part <- bnr2 
    }
    
    if (proc_mode == "obj_wise" && part == "2parts_2") {
      bnr2 <- 322 
    }
    
    if (bnr2 == 321) { 
      is_label_1
      display(is_label_1)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      coords <- data.frame(x=as.numeric(row(is_label_1)),y=as.numeric(col(is_label_1)), is_label_1=as.numeric(is_label_1))
      coords <- coords[coords$is_label_1 == 1,] #removal of pixels which do not have the label of the building
      
      #calculation of new centre of object from connected components
      xc <- coor_part[1,1]
      yc <- coor_part[1,2]
      alpha <- coor_part[1,5] * omega
      r_max <- shap_part[1,6]
      
      #storage of plot parameters of separated object:
      # xc, yc, r_max, alpha, dy_window_plot
      
      setwd(home_dir)
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar,file=f1) 
      
      #plot of PC and checkpoints
      dev.set(2)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
           " - left", sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      
      #generation of image
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
           r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
      points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      dev.off()

      #output of file with PC-coordinates (required in program 'line_detection.R')
      coords2 <- subset(coords, select=c(x,y))
      head(coords2)
      N1 <- length(coords$x)
      x_dat <- rep(0,N1)
      y_dat <- rep(0,N1)
      idx <- rep(0,N1)
      idxy <- cbind(idx, x_dat, y_dat)
      idxy <- coords2 #use of old object name
      
      #idxy
      rownames(idxy) <- 1 : N1
      nrow(idxy)
      fname2 <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv", sep="")
      write.table(idxy, fname2,  sep= " ", row.names=T) ##output of pixel cluster for one building
      #end of output
      
      setwd(home_dir2)
      source(paste("line_detection_v",v_nr,".R", sep=""))
    } #end of bnr2=321
    
    #second part 
    if (bnr2 == 322) { 
      cat("start of spObj_extract_single_building.R -second part",sep = "", "\n")
      part <-  "2parts_2"
      bnr2_part <- bnr2 
      is_label_2 <- is_bnr2_label@.Data == 2 
      display(is_label_2)
      is_label_2
      coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
      coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
      xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=272
      yc <- coor_part[2,2]
      alpha <- coor_part[2,5]*omega
      r_max <- shap_part[2,6]
      dy_window_plot
      
      #output plot parameter
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameters xc,yc,r_max,alpha
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      
      #plot of PC and checkpoints
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
        ylim=c(yc+r_max2,yc-r_max2),xlab = NULL,ylab=NULL,ann=T,main=paste("b", bnr2," - second",sep = ""),axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      ##generation of image
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
           " - second part",sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      dev.off()
      #windows() #if necessary
      cat("end of spObj_extract_single_building.R - first part", sep = "","\n")
      
      #output of file with PC-coordinates (required in program 'line_detection.R')
      coords2 <- subset(coords, select=c(x,y))
      head(coords2)
      N1 <- length(coords$x)
      x_dat <- rep(0,N1)
      y_dat <- rep(0,N1)
      idx <- rep(0,N1)
      idxy <- cbind(idx, x_dat, y_dat)
      idxy <- coords2 #use of old object name
      
      #idxy
      rownames(idxy) <- 1 : N1
      nrow(idxy)
      fname2 <- paste("./data/",Img_name,"/idxy_LCM_b",bnr2,".csv", sep="")
      write.table(idxy, fname2,  sep= " ", row.names=T) ##output of pixel cluster for one building
      #end of output
      
      cat("end of program 'spObj_extract_single_building.R' - continue with 'line_detection.R' ","\n")
      
      setwd(home_dir2)
      source(paste("line_detection_v",v_nr,".R", sep=""))
    } #end of bnr2=322
  } # end bnr2 = 32

} #end ISPRS7
#################################################################################


if (Img_name == "ISPRS1") {
  
  #b6
  if (bnr2 == 6 && p_pos == "cor_sep") {
    cat("start of spObj_extract_single_building.R -first part",sep = "")
    display(is_bnr)
    is_bnr2 <- is_bnr
    imageData(is_bnr2)[92:96, 262:266] = FALSE
    imageData(is_bnr2)[184:188, 434:438] = FALSE
    display(is_bnr2)
    is_bnr2_label <- bwlabel(is_bnr2)
    cat("Number of objects= ",max(is_bnr2_label),"\n")
    coor_part<-computeFeatures.moment(is_bnr2_label) #geometric features (moment)
    shap_part<-computeFeatures.shape(is_bnr2_label) #geometric features (shape)
    
    #first part
    is_label_1 <- is_bnr2_label@.Data == 1 
    imageData(is_label_1)[92:188, 262:438]
    display(is_label_1)
    
    #second part
    is_label_2 <- is_bnr2_label@.Data == 2 
    display(is_label_2)
    
    #save original values
    xc_orig <- xc
    yc_orig <- yc
    #
    
    if (proc_mode == "obj_wise" && part == "2parts_1") {
      bnr2 <- 61
      bnr2_part <- bnr2 
    }
    
    if (proc_mode == "obj_wise" && part == "2parts_2") {
      bnr2 <- 62 
    }
    
    #first part
    if (bnr2 == 61) { 
      is_label_1
      display(is_label_1)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      coords <- data.frame(x=as.numeric(row(is_label_1)),y=as.numeric(col(is_label_1)), is_label_1=as.numeric(is_label_1))
      coords <- coords[coords$is_label_1 == 1,] #removal of pixels which do not have the label of the building
      
      #calculation of new centre of object from connected components
      xc <- coor_part[1,1]
      yc <- coor_part[1,2]
      alpha <- coor_part[1,5]*omega
      r_max <- shap_part[1,6]
      
      #storage of plot parameter of separated object
      setwd(home_dir)
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar,file=f1) #parameter xc, yc, r_max, alpha
      
      #plot of PC and checkpoints (large scale)
      dev.set(2)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   " - left",sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      
      #generation of image
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
      points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      dev.off()
      cat("end of spObj_extract_single_building.R -first part",sep = "")
    } #end bnr2 = 61 (first part)
    
    #second part 
    if (bnr2 == 62) { 
      cat("start of spObj_extract_single_building.R - second part",sep = "", "\n")
      part <-  "2parts_2"
      bnr2_part <- bnr2 
      is_label_2 <- is_bnr2_label@.Data == 2 #right part
      display(is_label_2)
      is_label_2
      coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
      coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
      xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=342
      yc <- coor_part[2,2]
      alpha <- coor_part[2,5]*omega
      r_max <- shap_part[2,6]
      dy_window_plot
      
      #output plot parameter
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      
      #plot of PC and checkpoints (large scale)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
           asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, 
                                                    yc-r_max2), xlab = NULL, ylab=NULL, ann=T, 
           main=paste("b", bnr2,sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      
      #generation of image
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      dev.off()
      #windows() #if necessary
      cat("end of spObj_extract_single_building.R - second part", sep = "","\n")
    } #end of bnr2=62
  } # end bnr2 = 6
  
  
  
  ##################################################
  #27
  if (bnr2 == 27 && p_pos == "cor_sep") {
    cat("start of spObj_extract_single_building.R -first part",sep = "")
    display(is_bnr)
    is_bnr2 <- is_bnr
    imageData(is_bnr2)[649:654, 1877:1881] = FALSE
    imageData(is_bnr2)[663:667, 1886:1892] = FALSE
    display(is_bnr2)
    is_bnr2_label <- bwlabel(is_bnr2)
    cat("Number of objects= ",max(is_bnr2_label),"\n")
    coor_part<-computeFeatures.moment(is_bnr2_label) #geometric features (moment)
    shap_part<-computeFeatures.shape(is_bnr2_label) #geometric features (shape)
    
    #first part
    is_label_1 <- is_bnr2_label@.Data == 1 
    imageData(is_label_1)[649:663, 1877:1886]
    display(is_label_1)
    
    #second part
    is_label_2 <- is_bnr2_label@.Data == 2 
    display(is_label_2)
    
    #save original values
    xc_orig <- xc
    yc_orig <- yc
    #
    
    if (proc_mode == "obj_wise" && part == "2parts_1") {
      bnr2 <- 271
      bnr2_part <- bnr2 
    }
    
    if (proc_mode == "obj_wise" && part == "2parts_2") {
      bnr2 <- 272 
    }
    
    #first part
    if (bnr2 == 271) { 
      is_label_1
      display(is_label_1)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      coords <- data.frame(x=as.numeric(row(is_label_1)),y=as.numeric(col(is_label_1)), is_label_1=as.numeric(is_label_1))
      coords <- coords[coords$is_label_1 == 1,] #removal of pixels which do not have the label of the building
      
      #calculation of new centre of object from connected components
      xc <- coor_part[1,1]
      yc <- coor_part[1,2]
      alpha <- coor_part[1,5]*omega
      r_max <- shap_part[1,6]
      
      #storage of plot parameter of separated object
      setwd(home_dir)
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar,file=f1) #parameter xc, yc, r_max, alpha
      
      #plot of PC and checkpoints (large scale)
      dev.set(2)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   " - left",sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      
      #generation of image
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
      points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      dev.off()
      cat("end of spObj_extract_single_building.R -first part",sep = "")
    } #end bnr2 = 321 (first part)
    
    #second part 
    if (bnr2 == 272) { 
      cat("start of spObj_extract_single_building.R - second part",sep = "", "\n")
      part <-  "2parts_2"
      bnr2_part <- bnr2 
      is_label_2 <- is_bnr2_label@.Data == 2 #right part
      display(is_label_2)
      is_label_2
      coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
      coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
      xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=342
      yc <- coor_part[2,2]
      alpha <- coor_part[2,5]*omega
      r_max <- shap_part[2,6]
      dy_window_plot
      
      #output plot parameter
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      
      #plot of PC and checkpoints (large scale)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
           asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, 
                                                    yc-r_max2), xlab = NULL, ylab=NULL, ann=T, 
           main=paste("b", bnr2,sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      
      #generation of image
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      dev.off()
      #windows() #if necessary
      cat("end of spObj_extract_single_building.R - second part", sep = "","\n")
    } #end of bnr2=272
  } # end bnr2 = 27
  
  ##########################################
  
  #b34
  if (bnr2 == 34 && p_pos == "cor_sep") {
    cat("start of spObj_extract_single_building.R -first part",sep = "")
    display(is_bnr)
    is_bnr2 <- is_bnr
    imageData(is_bnr2)[1524:1528, 2055:2059] = FALSE
    imageData(is_bnr2)[1544:1548, 2095:2099] = FALSE
    display(is_bnr2)
    is_bnr2_label <- bwlabel(is_bnr2)
    cat("Number of objects= ",max(is_bnr2_label),"\n")
    coor_part<-computeFeatures.moment(is_bnr2_label) #geometric features (moment)
    shap_part<-computeFeatures.shape(is_bnr2_label) #geometric features (shape)
    
    #first part
    is_label_1 <- is_bnr2_label@.Data == 1 
    imageData(is_label_1)[1524:1548, 2055:2099]
    display(is_label_1)
    
    #second part
    is_label_2 <- is_bnr2_label@.Data == 2 
    display(is_label_2)
    
    #save original values
    xc_orig <- xc
    yc_orig <- yc
    #
    
    if (proc_mode == "obj_wise" && part == "2parts_1") {
      bnr2 <- 341
      bnr2_part <- bnr2 
    }
    
    if (proc_mode == "obj_wise" && part == "2parts_2") {
      bnr2 <- 342 
    }
    
    #first part
    if (bnr2 == 341) { 
      is_label_1
      display(is_label_1)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      coords <- data.frame(x=as.numeric(row(is_label_1)),y=as.numeric(col(is_label_1)), is_label_1=as.numeric(is_label_1))
      coords <- coords[coords$is_label_1 == 1,] #removal of pixels which do not have the label of the building
      
      #calculation of new centre of object from connected components
      xc <- coor_part[1,1]
      yc <- coor_part[1,2]
      alpha <- coor_part[1,5]*omega
      r_max <- shap_part[1,6]
      
      #storage of plot parameter of separated object
      setwd(home_dir)
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar,file=f1) #parameter xc, yc, r_max, alpha
      
      #plot of PC and checkpoints (large scale)
      dev.set(2)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
           " - left",sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      
      #generation of image
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
      points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      dev.off()
      cat("end of spObj_extract_single_building.R -first part",sep = "")
    } #end bnr2 = 321 (first part)
    
    #second part 
    if (bnr2 == 342) { 
      cat("start of spObj_extract_single_building.R - second part",sep = "", "\n")
      part <-  "2parts_2"
      bnr2_part <- bnr2 
      is_label_2 <- is_bnr2_label@.Data == 2 #right part
      display(is_label_2)
      is_label_2
      coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
      coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
      xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=342
      yc <- coor_part[2,2]
      alpha <- coor_part[2,5]*omega
      r_max <- shap_part[2,6]
      dy_window_plot
      
      #output plot parameter
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      
      #plot of PC and checkpoints (large scale)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
       asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, 
       yc-r_max2), xlab = NULL, ylab=NULL, ann=T, 
       main=paste("b", bnr2,sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      
      #generation of image
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
           sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      dev.off()
      #windows() #if necessary
      cat("end of spObj_extract_single_building.R - second part", sep = "","\n")
    } #end of bnr2=342
  } # end bnr2 = 34
  
  #b37
  
  if (bnr2 == 37 && p_pos == "cor_sep") {
    cat("start of spObj_extract_single_building.R -first part",sep = "")
    display(is_bnr)
    is_bnr2 <- is_bnr
    imageData(is_bnr2)[1278:1284, 2258:2262] = FALSE #first cut
    imageData(is_bnr2)[1332:1336, 2215:2219] = FALSE #second cut
    display(is_bnr2)
    is_bnr2_label <- bwlabel(is_bnr2)
    cat("Number of objects= ",max(is_bnr2_label),"\n")
    coor_part<-computeFeatures.moment(is_bnr2_label) #geometric features (moment)
    shap_part<-computeFeatures.shape(is_bnr2_label) #geometric features (shape)
    
    #first part
    is_label_1 <- is_bnr2_label@.Data == 1 
    imageData(is_label_1)[1282:1334, 2260:2217]
    display(is_label_1)
    
    #second part
    is_label_2 <- is_bnr2_label@.Data == 2 
    display(is_label_2)
    
    #save original values
    xc_orig <- xc
    yc_orig <- yc
    #
    
    if (proc_mode == "obj_wise" && part == "2parts_1") {
      bnr2 <- 371
      bnr2_part <- bnr2 
    }
    
    if (proc_mode == "obj_wise" && part == "2parts_2") {
      bnr2 <- 372 
    }
    
    #first part
    if (bnr2 == 371) { 
      is_label_1
      display(is_label_1)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      coords <- data.frame(x=as.numeric(row(is_label_1)),y=as.numeric(col(is_label_1)), is_label_1=as.numeric(is_label_1))
      coords <- coords[coords$is_label_1 == 1,] #removal of pixels which do not have the label of the building
      
      #calculation of new centre of object from connected components
      xc <- coor_part[1,1]
      yc <- coor_part[1,2]
      alpha <- coor_part[1,5]*omega
      r_max <- shap_part[1,6]
      
      #storage of plot parameter of separated object
      setwd(home_dir)
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar,file=f1) #parameter xc, yc, r_max, alpha
      
      #plot of PC and checkpoints (large scale)
      dev.set(2)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   " - left",sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      #
      
      #generation of image
      setwd(home_dir)
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann= FALSE, main=paste("b", bnr2), axes=TRUE)
      points(xc, yc, pch = 16, cex=1.5, col = "black", asp=1) #centre of PC
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      dev.off()
      cat("end of spObj_extract_single_building.R -first part",sep = "")
    } #end bnr2 = 371 (first part)
    
    #second part 
    if (bnr2 == 372) { 
      cat("start of spObj_extract_single_building.R - second part",sep = "", "\n")
      part <-  "2parts_2"
      bnr2_part <- bnr2 
      is_label_2 <- is_bnr2_label@.Data == 2 #right part
      display(is_label_2)
      is_label_2
      coords <- data.frame(x=as.numeric(row(is_label_2)),y=as.numeric(col(is_label_2)), is_label_2=as.numeric(is_label_2))
      coords <- coords[coords$is_label_2 == 1,] #removal of pixels which do not have the label of the building
      xc <- coor_part[2,1] #generate xc,yc,r_max,alpha and store with bnr2=342
      yc <- coor_part[2,2]
      alpha <- coor_part[2,5]*omega
      r_max <- shap_part[2,6]
      dy_window_plot
      
      #output plot parameter
      plotPar <- c(xc,yc,r_max,alpha,dy_window_plot)
      setwd(home_dir)
      f1 <- paste("./data/",Img_name,"/param_b_",bnr2,sep="")
      save(plotPar, file=f1) #parameter xc, yc, r_max, alpha
      f <- paste("./data/",Img_name,"/b_nr",sep = "")
      save(bnr2,file=f)
      
      #plot of PC and checkpoints (large scale)
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",
           asp=1,xlim=c(xc-r_max2,xc+r_max2),ylim=c(yc+r_max2, 
                                                    yc-r_max2), xlab = NULL, ylab=NULL, ann=T, 
           main=paste("b", bnr2,sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      
      #generation of image
      file1 <- paste('./data/',Img_name,'/images/b',bnr2,'_new8.tif',sep = "")
      tiff(file1, width=578, height=578, units="px", bg = "white")
      r_max2 <- round(1.1*r_max)
      plot(coords$x, coords$y, pch=16, cex=0.2,col="black",asp=1,xlim=c(xc-r_max2,xc+r_max2),
           ylim=c(yc+r_max2, yc-r_max2), xlab = NULL, ylab=NULL, ann=T, main=paste("b", bnr2,
                                                                                   sep = ""), axes=TRUE)
      points(xc+r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc+r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc-r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc+r_max, yc-r_max, pch=16, cex=1.5, col="black", asp=1) #point for scaling
      points(xc, yc, pch = 3, cex=1.5, col = "red", asp=1) #centre of PC
      dev.off()
      #windows() #if necessary
      cat("end of spObj_extract_single_building.R - second part", sep = "","\n")
    } #end of bnr2=372
  } # end bnr2 = 37
  
} #end ISPRS1
###################################################################
cat("end of program 'spObj_extract_single_building.R' ",sep = "","\n")

###################################################################

