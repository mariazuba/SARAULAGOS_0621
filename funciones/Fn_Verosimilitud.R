
Verosimilitudbase<-function(admb,dir.0,dir.1,Carpeta,l_log_Ro,l_opt_Ro,l_Fase_desRt,l_Fase_devNo,system){
  
  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "Retospectivobase"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta "Retospectivobase"
  setwd(dir) 
  
  #if(system=="mac"){
   # system(paste("~/admb-12.2/admb",admb,sep=" "))
   # system(paste("./",admb,sep=""))}
  
  #if(system=="windows"){
  #  system(paste("admb",admb,sep=" "))
  #  system(admb)}
  
  #===========================================================================================================
  # PRIMER PASO: LEE ARCHIVO DE DATOS .DAT (CASOBASE)
  #===========================================================================================================
  data.0 <- readLines(dat_admb,encoding="UTF-8")
  data.1 <- data.0
  
  #===========================================================================================================
  # SEGUNDO PASO: CREA VECTOR PARA R0 (CASOS O ESCENARIOS A PROBAR)
  #===========================================================================================================
  log_Ro    <- as.numeric(data.1[l_log_Ro])
  prior_Ro  <- seq(log_Ro*0.90,log_Ro*1.10,0.05)
  casos     <- length(prior_Ro);casos
  #===========================================================================================================
  # TERCER PASO: CREA LOS ARCHIVOS .DAT DE CADA CASO
  #===========================================================================================================
  for(i in 1:casos){
    data.1[l_log_Ro]  <- prior_Ro[i]
    data.1[l_opt_Ro]  <- -1
    data.1[l_Fase_desRt]  <- 1
    data.1[l_Fase_devNo]  <- 2
    
    cat(data.1,file=(can<-file(paste(admb,"s",i,".dat", sep=''),"wb",encoding="UTF-8")),sep="\n")
    close(can)
    
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



Verosimilitudalternativo<-function(admb,dir.0,dir.1,Carpeta,l_log_Ro,l_opt_Ro,l_Fase_desRt_devNo,system){
  
  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "Retospectivobase"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta "Retospectivobase"
  setwd(dir) 
  

  
  #===========================================================================================================
  # PRIMER PASO: LEE ARCHIVO DE DATOS .DAT (CASOBASE)
  #===========================================================================================================
  data.0 <- readLines(dat_admb,encoding="UTF-8")
  data.1 <- data.0
  
  #===========================================================================================================
  # SEGUNDO PASO: CREA VECTOR PARA R0 (CASOS O ESCENARIOS A PROBAR)
  #===========================================================================================================
  log_Ro    <- as.numeric(data.1[l_log_Ro])
  prior_Ro  <- seq(log_Ro*0.90,log_Ro*1.10,0.05)
  casos     <- length(prior_Ro);casos
  #===========================================================================================================
  # TERCER PASO: CREA LOS ARCHIVOS .DAT DE CADA CASO
  #===========================================================================================================
  for(i in 1:casos){
    data.1[l_log_Ro]  <- prior_Ro[i]
    data.1[l_opt_Ro]  <- -1
    data.1[l_Fase_desRt_devNo]  <- paste(1,2,sep="\t") #Fase_dev_Rt_devNo

    
    cat(data.1,file=(can<-file(paste(admb,"s",i,".dat", sep=''),"wb",encoding="UTF-8")),sep="\n")
    close(can)
    
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









Verosimilitudseptiembre<-function(admb,dir.0,dir.1,dir.5,Carpeta){
  
  dat_admb<-paste(admb,".dat",sep="")
  exe_admb<-paste(admb,".exe",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  unlink(dir.5,recursive=T) #borra la carpeta "CBA2016"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
  setwd(dir.1);file.copy(c(dat_admb,exe_admb),dir.5) #copia los archivos de la carpeta MAE0316
  setwd(dir.5) 
  system(paste(admb,"-nox",sep=" "))
  #===========================================================================================================
  # PRIMER PASO: LEE ARCHIVO DE DATOS .DAT (CASOBASE)
  #===========================================================================================================
  data.0 <- readLines(dat_admb,encoding="UTF-8")
  data.1 <- data.0
  
  #===========================================================================================================
  # SEGUNDO PASO: CREA VECTOR PARA R0 (CASOS O ESCENARIOS A PROBAR)
  #===========================================================================================================
  log_Ro    <- as.numeric(data.1[155])
  prior_Ro  <- seq(log_Ro*0.90,log_Ro*1.10,0.05)
  casos     <- length(prior_Ro);casos
  #===========================================================================================================
  # TERCER PASO: CREA LOS ARCHIVOS .DAT DE CADA CASO
  #===========================================================================================================
  for(i in 1:casos){
    data.1[155]  <- prior_Ro[i]
    data.1[152]  <- -1
    data.1[161]  <- 1
    data.1[164]  <- 2
    
    cat(data.1,file=(can<-file(paste(admb,"s",i,".dat", sep=''),"wb",encoding="UTF-8")),sep="\n")
    close(can)}
  #===========================================================================================================
  # CUARTO PASO: CORRE TODOS LOS CASOS CREADOS EN .DAT Y GUARDA LOS ARCHIVOS .REP, .STD Y .PAR DE CADA CASO
  #===========================================================================================================
  # Nota: es necesario tener en el directorio los archivos por lotes "correall" y "run", modificado seg?n corresponda
  dat_run<-paste(admb, "-ind %1.dat",sep=" ")
  rep_run<-paste("copy", rep_admb,"%1.rep",sep=" ")
  std_run<-paste("copy", std_admb,"%1.std",sep=" ")
  
  run<-rbind(dat_run,rep_run,std_run)
  cat(run,file=(can<-file("run.bat","wb",encoding="UTF-8")),sep="\n")
  close(can)
  
  n<-casos 
  corre <- rep(NA,n)
  s     <- seq(1,n,1)
  for(i in 1:n){corre[i]<-paste("call run ",admb,"s",s[i],sep="")}
  cat(corre,file=(can<-file("corre.bat","wb",encoding="UTF-8")),sep="\n");close(can)
  
  system("corre") 
  
}

Verosimilitudjunio<-function(admb,dir.0,dir.1,dir.5,Carpeta){
  
  dat_admb<-paste(admb,".dat",sep="")
  exe_admb<-paste(admb,".exe",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  unlink(dir.5,recursive=T) #borra la carpeta "CBA2016"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
  setwd(dir.1);file.copy(c(dat_admb,exe_admb),dir.5) #copia los archivos de la carpeta MAE0316
  setwd(dir.5) 
  system(paste(admb,"-nox",sep=" "))
  #===========================================================================================================
  # PRIMER PASO: LEE ARCHIVO DE DATOS .DAT (CASOBASE)
  #===========================================================================================================
  data.0 <- readLines(dat_admb,encoding="UTF-8")
  data.1 <- data.0
  
  #===========================================================================================================
  # SEGUNDO PASO: CREA VECTOR PARA R0 (CASOS O ESCENARIOS A PROBAR)
  #===========================================================================================================
  log_Ro    <- as.numeric(data.1[158])
  prior_Ro  <- seq(log_Ro*0.90,log_Ro*1.10,0.05)
  casos     <- length(prior_Ro);casos
  #===========================================================================================================
  # TERCER PASO: CREA LOS ARCHIVOS .DAT DE CADA CASO
  #===========================================================================================================
  for(i in 1:casos){
    data.1[158]  <- prior_Ro[i]
    data.1[155]  <- -1
    data.1[164]  <- 1
    data.1[167]  <- 2
    
    cat(data.1,file=(can<-file(paste(admb,"s",i,".dat", sep=''),"wb",encoding="UTF-8")),sep="\n")
    close(can)}
  #===========================================================================================================
  # CUARTO PASO: CORRE TODOS LOS CASOS CREADOS EN .DAT Y GUARDA LOS ARCHIVOS .REP, .STD Y .PAR DE CADA CASO
  #===========================================================================================================
  # Nota: es necesario tener en el directorio los archivos por lotes "correall" y "run", modificado seg?n corresponda
  dat_run<-paste(admb, "-ind %1.dat",sep=" ")
  rep_run<-paste("copy", rep_admb,"%1.rep",sep=" ")
  std_run<-paste("copy", std_admb,"%1.std",sep=" ")
  
  run<-rbind(dat_run,rep_run,std_run)
  cat(run,file=(can<-file("run.bat","wb",encoding="UTF-8")),sep="\n")
  close(can)
  
  n<-casos 
  corre <- rep(NA,n)
  s     <- seq(1,n,1)
  for(i in 1:n){corre[i]<-paste("call run ",admb,"s",s[i],sep="")}
  cat(corre,file=(can<-file("corre.bat","wb",encoding="UTF-8")),sep="\n");close(can)
  
  system("corre") 
  
}
#===========================================================================================================

Figura_verosimilitud<-function(dir.4,admb){
setwd(dir.4)
casos <-35
logRo    <- rep(0,casos)
likeval  <- matrix(ncol=9,nrow=casos)
slikeval <- matrix(ncol=10,nrow=casos)

for(i in 1:casos){
  rep         <- reptoRlist(paste(admb,"s",i,".rep",sep=""))
  data        <- readLines(paste(admb,"s",i,".dat", sep=''),encoding="UTF-8")
  logRo[i]    <- as.numeric(data[139])
  likeval[i,] <- rep$likeval}

#===========================================================================================================
# SEXTO PASO: ESTANDARIZAR VEROSIMILITUD
#===========================================================================================================
like    <- data.frame(round(likeval,3),Total=apply(likeval,1,sum))
minLik  <- apply(like,2,min)                         # busca el m?nimo
for(i in 1:10){slikeval[,i]<-like[,i]-minLik[i]}    # Estandarizaci?n
#===========================================================================================================
# ULTIMO PASO: GUARDAR TABLAS Y FIGURA
#===========================================================================================================
names<-c("Ro","cpue",	"cru","mph", 	"desemb",	"proflo",	"procru",	"desvRo",	"desNo",	"Lo",
         "Total")
#--------------------------------------------------------------------------------------------
# Tabla verosimilitud
TLk1 <- data.frame(exp(logRo),like);
colnames(TLk1)<-names
# Tabla estandarizada
TLk2 <- data.frame(exp(logRo),slikeval);
colnames(TLk2)<-names

par(mar=c(4,4,1,1)+0.5)
plot(TLk2$Ro,TLk2$Total,type="l",lwd=3,ylim=c(0,4),xlim=c(800,16000), xaxs= "i",	ylab="L-min(L)",xlab="Ro",
     las=1,main="Modelo alternativo")
#lines(c(0,TLk2$Ro),rep(2,casos+1),lty=2,lwd=2)
for(i in 2:7){lines(TLk2$Ro,TLk2[,i],col=i,lty=2,lwd=2)}
#for(i in 8:10){lines(TLk2$Ro,TLk2[,i],col=i,lty=3,lwd=2)}
legend(12000,4,names[c(11,2:7)],col=1:8,lty=c(1,rep(2,7)),lwd=2,bty="n",cex=0.8)
}

#===========================================================================================================

Figura_verosimilitudbase<-function(dir.5,admb,titulo){
  setwd(dir.5)
  casos <-35
  logRo    <- rep(0,casos)
  likeval  <- matrix(ncol=9,nrow=casos)
  slikeval <- matrix(ncol=10,nrow=casos)
  
  for(i in 1:casos){
    rep         <- reptoRlist(paste(admb,"s",i,".rep",sep=""))
    data        <- readLines(paste(admb,"s",i,".dat", sep=''),encoding="UTF-8")
    logRo[i]    <- as.numeric(data[155])
    likeval[i,] <- rep$Likeval}
  
  #===========================================================================================================
  # SEXTO PASO: ESTANDARIZAR VEROSIMILITUD
  #===========================================================================================================
  like    <- data.frame(round(likeval,3),Total=apply(likeval,1,sum))
  minLik  <- apply(like,2,min)                         # busca el minimo
  for(i in 1:10){slikeval[,i]<-like[,i]-minLik[i]}    # Estandarizaci?n
  #===========================================================================================================
  # ULTIMO PASO: GUARDAR TABLAS Y FIGURA
  #===========================================================================================================
  #CPUE   BCru   Bmph   Desemb  pfFlota  pfCru  dev_R  devNo LR
  names<-c("Ro","cpue",	"cru","mph", 	"desemb",	"proflo",	"procru",	"desvRo",	"desNo",	"Lo",
           "Total")
  #--------------------------------------------------------------------------------------------
  # Tabla verosimilitud
  TLk1 <- data.frame(exp(logRo),like);
  colnames(TLk1)<-names
  # Tabla estandarizada
  TLk2 <- data.frame(exp(logRo),slikeval);
  colnames(TLk2)<-names
  
  par(mar=c(4,4,1,1)+0.5)
  plot(TLk2$Ro,TLk2$Total,type="l",lwd=3,ylim=c(0,4),xlim=c(800,16000), xaxs= "i",	ylab="L-min(L)",xlab="Ro",
       las=1,main=titulo)
   for(i in 2:7){lines(TLk2$Ro,TLk2[,i],col=i,lty=2,lwd=2)}
   legend(12000,4,names[c(11,2:7)],col=1:8,lty=c(1,rep(2,7)),lwd=2,bty="n",cex=0.8)
}


Figura_verosimilitudjunio<-function(dir.5,admb,titulo){
  setwd(dir.5)
  casos <-35
  logRo    <- rep(0,casos)
  likeval  <- matrix(ncol=9,nrow=casos)
  slikeval <- matrix(ncol=10,nrow=casos)
  
  for(i in 1:casos){
    rep         <- reptoRlist(paste(admb,"s",i,".rep",sep=""))
    data        <- readLines(paste(admb,"s",i,".dat", sep=''),encoding="UTF-8")
    logRo[i]    <- as.numeric(data[158])
    likeval[i,] <- rep$Likeval}
  
  #===========================================================================================================
  # SEXTO PASO: ESTANDARIZAR VEROSIMILITUD
  #===========================================================================================================
  like    <- data.frame(round(likeval,3),Total=apply(likeval,1,sum))
  minLik  <- apply(like,2,min)                         # busca el minimo
  for(i in 1:10){slikeval[,i]<-like[,i]-minLik[i]}    # Estandarizaci?n
  #===========================================================================================================
  # ULTIMO PASO: GUARDAR TABLAS Y FIGURA
  #===========================================================================================================
  #CPUE   BCru   Bmph   Desemb  pfFlota  pfCru  dev_R  devNo LR
  names<-c("Ro","cpue",	"cru","mph", 	"desemb",	"proflo",	"procru",	"desvRo",	"desNo",	"Lo",
           "Total")
  #--------------------------------------------------------------------------------------------
  # Tabla verosimilitud
  TLk1 <- data.frame(exp(logRo),like);
  colnames(TLk1)<-names
  # Tabla estandarizada
  TLk2 <- data.frame(exp(logRo),slikeval);
  colnames(TLk2)<-names
  
  par(mar=c(4,4,1,1)+0.5)
  plot(TLk2$Ro,TLk2$Total,type="l",lwd=3,ylim=c(0,4),xlim=c(800,16000), xaxs= "i",	ylab="L-min(L)",xlab="Ro",
       las=1,main=titulo)
  for(i in 2:7){lines(TLk2$Ro,TLk2[,i],col=i,lty=2,lwd=2)}
  legend(12000,4,names[c(11,2:7)],col=1:8,lty=c(1,rep(2,7)),lwd=2,bty="n",cex=0.8)
}



