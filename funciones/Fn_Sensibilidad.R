


Sensibilidad<-function(dir.0,dir.1,Carpeta,admb){
  
  dir.6<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  exe_admb<-paste(admb,".exe",sep="")
  
  unlink(dir.6,recursive=T) 
  dir.create(file.path(dir.0,Carpeta))
  
  setwd(dir.1);file.copy(c(dat_admb,exe_admb),dir.6) 
  setwd(dir.6);system(paste(admb,"-nox",sep=" "))
  
  data        <- lisread(paste("admb",".dat",sep="")) 
  names(data) <- str_trim(names(data), side="right")
  dat         <- data
  
  #==========================================================================
  #######################     CREA Y CORRE ESCENARIOS #####################
  #==========================================================================
  setwd(dir.8)
  #--------------
  # escenario 0: Caso 0	Igual al caso base de septiembre 2019 (MAE0919)
  #--------------
  # caso base
  dat<- data
  dat$opRec <- 1
  writeData(paste("admb","s",1,".dat",sep=""), dat, append=FALSE)
  ###########################################################################
  # **Actualizaci?n 2019**
  ###########################################################################
  #--------------
  # escenario 1:  Caso 0 + Captura 2018 + descarte 2018(2%)
  #--------------
  dat<- data
  dat$opRec <- 1
  dat$Ind[29,8] <- 62187
  dat_mar$NMs<-c(10, 17, 6, 0)
  writeData(paste("admb","s",2,".dat",sep=""), dat, append=FALSE)
  #--------------
  # escenario 2**: Caso 1 + Captura 2019 + descarte 2019(1%)
  #--------------
  dat<- data
  dat$opRec <- 1
  dat$Ind[30,8] <- 160968
  dat_mar$NMs<-c(10, 17, 6, 0)
  writeData(paste("admb","s",3,".dat",sep=""), dat, append=FALSE)
  #--------------
  # escenario 3**: Caso 2 + Composici?n de edad de la flota 2019
  #--------------
  dat<- data
  dat$opRec <- 1
  dat$Ind[29,8] <- 62187
  dat$Ind[30,8] <- 160968
  dat$Captura_edad[30,] <- c(686, 1820, 1412, 726, 160)
  dat_mar$NMs<-c(10, 17, 6, 0)
  writeData(paste("admb","s",4,".dat",sep=""), dat, append=FALSE)
  
  
}