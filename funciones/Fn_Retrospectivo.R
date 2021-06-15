Retrospectivobase<-function(admb,dir.0,dir.1,Carpeta,system){

  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "Retospectivobase"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta "Retospectivobase"
  setwd(dir) 
  
  

  data        <- lisread(paste(dir,dat_admb, sep='/'))
  names(data) <- str_trim(names(data), side="right")
  data.1      <- data
  retros      <- c(0:4)
  for(i in 1:length(retros)){
    data.1$nanos               <- data$nanos-retros[i]
    data.1$Ind                 <- data$Ind[1:(data$nanos-retros[i]),]
    data.1$Frecuencia_flota    <- data$Frecuencia_flota[1:(data$nanos-retros[i]),]
    data.1$Frecuencia_cruceros <- data$Frecuencia_cruceros[1:(data$nanos-retros[i]),]
    
    writeData(paste(admb,"s",i,".dat",sep=""), data.1, append=F)
  
  setwd(dir.1)
  file.copy(c(paste(admb,".tpl",sep="")),dir)
  
  setwd(dir)
  file.rename(paste(admb,".tpl",sep=""),paste(admb,"s",i,".tpl",sep="")) 
  
  if(system=="mac"){
    system(paste("~/admb-12.2/admb ",admb,"s",i,sep=""))
    system(paste("./",admb,"s",i,sep="")) }
  
  if(system=="windows"){
    system(paste("admb ",admb,"s",i,sep=""))
    system(paste(admb,"s",i,sep="")) }
  
  file.remove(paste(admb,"s",i,".htp", sep=""),
              paste(admb,"s",i,".cpp", sep=""),
              paste(admb,"s",i,".obj", sep=""),
              paste(admb,"s",i,".p01", sep=""),
              paste(admb,"s",i,".b01", sep=""),
              paste(admb,"s",i,".r01", sep=""),
              paste(admb,"s",i,".p02", sep=""),
              paste(admb,"s",i,".b02", sep=""),
              paste(admb,"s",i,".r02", sep=""),
              paste(admb,"s",i,".p03", sep=""),
              paste(admb,"s",i,".b03", sep=""),
              paste(admb,"s",i,".r03", sep=""),
              paste(admb,"s",i,".p04", sep=""),
              paste(admb,"s",i,".b04", sep=""),
              paste(admb,"s",i,".r04", sep=""),
              paste(admb,"s",i,".p05", sep=""),
              paste(admb,"s",i,".b05", sep=""),
              paste(admb,"s",i,".r05", sep=""),
              paste(admb,"s",i,".p06", sep=""),
              paste(admb,"s",i,".b06", sep=""),
              paste(admb,"s",i,".r06", sep=""),
              paste(admb,"s",i,".par", sep=""),
              paste(admb,"s",i,".bar", sep=""),
              paste(admb,"s",i,".eva", sep=""),
              paste(admb,"s",i,".cor", sep=""),
              paste(admb,"s",i,".std", sep=""),
              paste(admb,"s",i,".log", sep=""),
              paste(admb,"s",i,".tpl", sep=""),
              paste(admb,"s",i,".exe", sep=""))
  
  }
}




Retrospectivoalternativo<-function(admb,dir.0,dir.1,Carpeta,system,Fase_dev_Rt){
  
  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "Retospectivobase"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta "Retospectivobase"
  setwd(dir) 
  
 
  
  data        <- lisread(paste(dir,dat_admb, sep='/'))
  names(data) <- str_trim(names(data), side="right")
  data.1      <- data
  retros      <- c(0:5)
  for(i in 1:length(retros)){
    data.1$nanos               <- data$nanos-retros[i]
    data.1$Ind                 <- data$Ind[1:(data$nanos-retros[i]),]
    data.1$Frecuencia_flota    <- data$Frecuencia_flota[1:(data$nanos-retros[i]),]
    data.1$Frecuencia_cruceros <- data$Frecuencia_cruceros[1:(data$nanos-retros[i]),]
    data.1$Fase_dev_Rt_devNo[1] <- Fase_dev_Rt[i]
    
    writeData(paste(admb,"s",i,".dat",sep=""), data.1, append=F)
  
    setwd(dir.1)
    file.copy(c(paste(admb,".tpl",sep="")),dir)
    
    setwd(dir)
    file.rename(paste(admb,".tpl",sep=""),paste(admb,"s",i,".tpl",sep="")) 
    
    if(system=="mac"){
      system(paste("~/admb-12.2/admb ",admb,"s",i,sep=""))
      system(paste("./",admb,"s",i,sep="")) }
    
    if(system=="windows"){
      system(paste("admb ",admb,"s",i,sep=""))
      system(paste(admb,"s",i,sep="")) }
    
    file.remove(paste(admb,"s",i,".htp", sep=""),
                paste(admb,"s",i,".cpp", sep=""),
                paste(admb,"s",i,".obj", sep=""),
                paste(admb,"s",i,".p01", sep=""),
                paste(admb,"s",i,".b01", sep=""),
                paste(admb,"s",i,".r01", sep=""),
                paste(admb,"s",i,".p02", sep=""),
                paste(admb,"s",i,".b02", sep=""),
                paste(admb,"s",i,".r02", sep=""),
                paste(admb,"s",i,".p03", sep=""),
                paste(admb,"s",i,".b03", sep=""),
                paste(admb,"s",i,".r03", sep=""),
                paste(admb,"s",i,".p04", sep=""),
                paste(admb,"s",i,".b04", sep=""),
                paste(admb,"s",i,".r04", sep=""),
                paste(admb,"s",i,".p05", sep=""),
                paste(admb,"s",i,".b05", sep=""),
                paste(admb,"s",i,".r05", sep=""),
                paste(admb,"s",i,".p06", sep=""),
                paste(admb,"s",i,".b06", sep=""),
                paste(admb,"s",i,".r06", sep=""),
                paste(admb,"s",i,".par", sep=""),
                paste(admb,"s",i,".bar", sep=""),
                paste(admb,"s",i,".eva", sep=""),
                paste(admb,"s",i,".cor", sep=""),
                paste(admb,"s",i,".std", sep=""),
                paste(admb,"s",i,".log", sep=""),
                paste(admb,"s",i,".tpl", sep=""),
                paste(admb,"s",i,".exe", sep=""))
    
  }
}

Figura1_Retrospectivo<-function(dir.2,admb,ylabdif){
  setwd(dir.2)
  rep.0       <- reptoRlist(paste(admb,"s1.rep",sep=""))
  retros      <- c(0:4)
  
  # =====================================================================================#
  # AN?LISIS DE PATRONES RETROSPECTIVOS "MOHN RHO"
  # =====================================================================================#		
  # Reclutamientos
  retroR      <- matrix(0,nrow=length(rep.0$Years),ncol=length(retros)+1)
  retroR[,1]  <- rep.0$Years;		
  retroBD     <- matrix(0,nrow=length(rep.0$Years),ncol=length(retros)+1)
  retroBD[,1] <- rep.0$Years;		
  retroF      <- matrix(0,nrow=length(rep.0$Years),ncol=length(retros)+1)
  retroF[,1]  <- rep.0$Years;	
  
  
  for(i in 1:length(retros)){
    rep<- reptoRlist(paste(admb,"s",i,".rep",sep=""))
    retroR[,i+1] <- c(rep$Reclutamiento,rep(0,i-1))
    retroBD[,i+1] <- c(rep$Biomasa_desovante,rep(0,i-1))
    retroF[,i+1]  <- c(rep$F,rep(0,i-1))
  }
  
  #Código de Francisco y Fernando 
  # Mohn rho
  years      <- retroBD[,1]
  nyears     <- length(years)
  nyrs.retro <- length(retros)
  
  mohn.r     <- rep(NA, nyrs.retro) 
  rel.diff.r <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  mohn.ssb     <- rep(NA, nyrs.retro) 
  rel.diff.ssb <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  mohn.f       <- rep(NA, nyrs.retro)
  rel.diff.f   <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  
  for (j in 1:nyrs.retro) { 
    rel.diff.r[,j] <- (retroR[,(j+1)]-retroR[,2])/retroR[,2]
    mohn.r[j]      <- rel.diff.r[(nyears-j),j]
    rel.diff.ssb[,j] <- (retroBD[,(j+1)]-retroBD[,2])/retroBD[,2]
    mohn.ssb[j]      <- rel.diff.ssb[(nyears-j),j]
    rel.diff.f[,j]   <- (retroF[,(j+1)]-retroF[,2])/retroF[,2]
    mohn.f[j]        <- rel.diff.f[(nyears-j),j]}
  
  #reclutamientos
  ave.mohn.r  <- mean(mohn.r)
  a0r            <- which(rel.diff.r==-1,arr.ind=TRUE)
  junkr          <- rel.diff.r
  a1r            <- which(rel.diff.r==-1,arr.ind=TRUE);
  ffr            <- dim(a1r)
  for(l in 1:ffr[1]){rel.diff.r[a1r[l,1],a1r[l,2]]<-NA}
  
  #biomasa desovante
  ave.mohn.ssb  <- mean(mohn.ssb)
  a0b            <- which(rel.diff.ssb==-1,arr.ind=TRUE)
  junkb          <- rel.diff.ssb
  a1b            <- which(rel.diff.ssb==-1,arr.ind=TRUE);
  ffb            <- dim(a1b)
  for(l in 1:ffb[1]){rel.diff.ssb[a1b[l,1],a1b[l,2]]<-NA}
  
  #mortalidad por pesca
  ave.mohn.f    <- mean(mohn.f)
  a0f           <- which(rel.diff.f==-1,arr.ind=TRUE)
  junkf         <- rel.diff.f
  a1f           <- which(rel.diff.f==-1,arr.ind=TRUE)
  fff           <- dim(a1f)
  for(l in 1:fff[1]){rel.diff.f[a1f[l,1],a1f[l,2]]<-NA}
  
  ##########################
  par(mfcol=c(3,2),mar=c(4,4,2,1)+0.05)
  #reclutamientos
  plot(rep$Years,rep$Reclutamiento/10^3,las=1,type="n",xlim=c(2000,2020),ylim=c(0,25),
       xaxp=c(2000,2020,20),ylab="Reclutamientos*10^3",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$Years,rep$Reclutamiento/10^3,type="l",col=i,lwd=1)}
  #Biomasa desovante
  plot(rep$Years,rep$Biomasa_desovante/10^3,las=1,type="n",xlim=c(2000,2020),ylim=c(0,100),
       xaxp=c(2000,2020,20),ylab="Biomasa desovante*10^3",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$Years,rep$Biomasa_desovante/10^3,type="l",col=i,lwd=1)}
  legend (2002,45,rev(seq(2015,2019,1)),lty=1,col=seq(1,8,1),bty="n",cex=1,lwd=1)
  #Mortalidad por pesca
  plot(rep$Years,rep$F,las=1,type="n",xlim=c(2000,2020),ylim=c(0,1.1),
       xaxp=c(2000,2020,20),ylab="Mortalidad por Pesca",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$Years,rep$F,type="l",col=i,lwd=1)}
  #reclutamientos
  plot(years,rel.diff.r[,1],type="l",lty=2,ylim=c(-2,2),xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Reclutamientos",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.r[,j],col=j,lwd=1) }
  text(2003,-0.5,paste("Rho=",round(ave.mohn.r,2),sep=""),cex=1.2)
  #biomasa desovante
  plot(years,rel.diff.ssb[,1],type="l",lty=2,ylim=c(-2,2),xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Biomasa desovante",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.ssb[,j],col=j,lwd=1)}
  text(2003,-0.5,paste("Rho=",round(ave.mohn.ssb,2),sep=""),cex=1.2)
  legend (2010,-0.1,rev(seq(2015,2018,1)),lty=1,col=seq(2,8,1),bty="n",cex=1,lwd=1)
  #mortalidad por pesca 
  plot(years,rel.diff.f[,1],type="l",lty=2,ylim=c(-2,2),xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Mortalidad por pesca",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.f[,j],col=j,lwd=1)}
  text(2003,0.5,paste("Rho=",round(ave.mohn.f,2),sep=""),cex=1.2)
 }

Figura2_Retrospectivo<-function(dir.3,admb,ylabdif){
  setwd(dir.3)
  rep.0       <- reptoRlist(paste(admb,"s1.rep",sep=""))
  retros      <- c(0:4)
  
  # =====================================================================================#
  # AN?LISIS DE PATRONES RETROSPECTIVOS "MOHN RHO"
  # =====================================================================================#		
  # Reclutamientos
  retroR      <- matrix(0,nrow=length(rep.0$YRS),ncol=length(retros)+1)
  retroR[,1]  <- rep.0$YRS;		
  retroBD     <- matrix(0,nrow=length(rep.0$YRS),ncol=length(retros)+1)
  retroBD[,1] <- rep.0$YRS;		
  retroF      <- matrix(0,nrow=length(rep.0$YRS),ncol=length(retros)+1)
  retroF[,1]  <- rep.0$YRS;	
  
  
  for(i in 1:length(retros)){
    rep<- reptoRlist(paste(admb,"s",i,".rep",sep=""))
    retroR[,i+1] <- c(rep$Reclutas,rep(0,i-1))
    retroBD[,i+1] <- c(rep$BD,rep(0,i-1))
    retroF[,i+1]  <- c(rep$F,rep(0,i-1))
  }
  
  #Código de Francisco y Fernando 
  # Mohn rho
  years      <- retroBD[,1]
  nyears     <- length(years)
  nyrs.retro <- length(retros)
  
  mohn.r     <- rep(NA, nyrs.retro) 
  rel.diff.r <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  mohn.ssb     <- rep(NA, nyrs.retro) 
  rel.diff.ssb <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  mohn.f       <- rep(NA, nyrs.retro)
  rel.diff.f   <- matrix(NA, nrow=nyears, ncol=(nyrs.retro ))
  
  for (j in 1:nyrs.retro) { 
    rel.diff.r[,j] <- (retroR[,(j+1)]-retroR[,2])/retroR[,2]
    mohn.r[j]      <- rel.diff.r[(nyears-j),j]
    rel.diff.ssb[,j] <- (retroBD[,(j+1)]-retroBD[,2])/retroBD[,2]
    mohn.ssb[j]      <- rel.diff.ssb[(nyears-j),j]
    rel.diff.f[,j]   <- (retroF[,(j+1)]-retroF[,2])/retroF[,2]
    mohn.f[j]        <- rel.diff.f[(nyears-j),j]}
  
  #reclutamientos
  ave.mohn.r  <- mean(mohn.r)
  a0r            <- which(rel.diff.r==-1,arr.ind=TRUE)
  junkr          <- rel.diff.r
  a1r            <- which(rel.diff.r==-1,arr.ind=TRUE);
  ffr            <- dim(a1r)
  for(l in 1:ffr[1]){rel.diff.r[a1r[l,1],a1r[l,2]]<-NA}
  
  #biomasa desovante
  ave.mohn.ssb  <- mean(mohn.ssb)
  a0b            <- which(rel.diff.ssb==-1,arr.ind=TRUE)
  junkb          <- rel.diff.ssb
  a1b            <- which(rel.diff.ssb==-1,arr.ind=TRUE);
  ffb            <- dim(a1b)
  for(l in 1:ffb[1]){rel.diff.ssb[a1b[l,1],a1b[l,2]]<-NA}
  
  #mortalidad por pesca
  ave.mohn.f    <- mean(mohn.f)
  a0f           <- which(rel.diff.f==-1,arr.ind=TRUE)
  junkf         <- rel.diff.f
  a1f           <- which(rel.diff.f==-1,arr.ind=TRUE)
  fff           <- dim(a1f)
  for(l in 1:fff[1]){rel.diff.f[a1f[l,1],a1f[l,2]]<-NA}
  
  ##########################
  par(mfcol=c(3,2),mar=c(4,4,2,1)+0.05)
  #reclutamientos
  plot(rep$YRS,rep$Reclutas/10^3,las=1,type="n",xlim=c(2000,2020),ylim=c(0,25),
       xaxp=c(2000,2020,20),ylab="Reclutamientos*10^3",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$YRS,rep$Reclutas/10^3,type="l",col=i,lwd=1)}
  #Biomasa desovante
  plot(rep$YRS,rep$BD/10^3,las=1,type="n",xlim=c(2000,2020),ylim=c(0,100),
       xaxp=c(2000,2020,20),ylab="Biomasa desovante*10^3",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$YRS,rep$BD/10^3,type="l",col=i,lwd=1)}
  legend (2004,45,rev(seq(2015,2019,1)),lty=1,col=seq(1,8,1),bty="n",cex=1,lwd=1)
  #Mortalidad por pesca
  plot(rep$YRS,rep$F,las=1,type="n",xlim=c(2000,2020),ylim=c(0,1.5),
       xaxp=c(2000,2020,20),ylab="Mortalidad por Pesca",xlab="",cex.lab=1,cex.axis=1)
  for(i in 1:length(retros)){rep<-reptoRlist(paste(admb,"s",i,".rep",sep=""))
  lines(rep$YRS,rep$F,type="l",col=i,lwd=1)}
  #reclutamientos
  plot(years,rel.diff.r[,1],type="l",lty=2,ylim=c(-2.5,2.5),xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Reclutamientos",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.r[,j],col=j,lwd=1) }
  text(2003,-0.5,paste("Rho=",round(ave.mohn.r,2),sep=""),cex=1.2)
  #biomasa desovante
  plot(years,rel.diff.ssb[,1],type="l",lty=2,ylim=c(-1.5,1.5),xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Biomasa desovante",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.ssb[,j],col=j,lwd=1)}
  text(2003,-0.5,paste("Rho=",round(ave.mohn.ssb,2),sep=""),cex=1.2)
  legend (2010,-0.1,rev(seq(2015,2018,1)),lty=1,col=seq(2,8,1),bty="n",cex=1,lwd=1)
  #mortalidad por pesca 
  plot(years,rel.diff.f[,1],type="l",lty=2,xlim=c(2000,2020),xaxp=c(2000,2020,20),main="Mortalidad por pesca",ylab=ylabdif,xlab="")
  for (j in 1:nyrs.retro){lines(years,rel.diff.f[,j],col=j,lwd=1)}
  text(2003,0.5,paste("Rho=",round(ave.mohn.f,2),sep=""),cex=1.2)
}



