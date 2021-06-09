######################################################################
######################################################################
#                                                                    #
#                      TAMA?O DE MUESTRA                             #
#                                                                    #
######################################################################
######################################################################
rm(list=ls())   
library(stringr)
library(knitr) # para generar reporte Rmarkdown
library(stringr)
library(reshape) 
library(dplyr) 
library(ggplot2)
library(ggthemes) # para ggplot
library(patchwork) # para unir gráficos de ggplot
library(strucchange) # libreria utilizada para análisis de quiebres

dir.Fig     <-"Figuras/" # carpeta de las figuras utilizadas y generadas en este estudio
fig         <-c("pdf") # formato de figuras generadas por este código
dir.0       <-getwd() # directorio de trabajo 
dir.1       <-paste(dir.0,"/codigos_admb",sep="") # carpeta de códigos ADMB 
dir.2       <-paste(dir.0,"/Retrospectivobase",sep="") # carpeta de códigos ADMB 
dir.3       <-paste(dir.0,"/Retrospectivoalternativo",sep="") # carpeta de códigos ADMB 
dir.4       <-paste(dir.0,"/Verosimilitudalternativo",sep="") # carpeta de códigos ADMB 
dir.5       <-paste(dir.0,"/Verosimilitudbase",sep="") # carpeta de códigos ADMB 

dir.fun     <-paste(dir.0,"/funciones/",sep="") # carpeta de funciones utilizadas en este informe
source(paste(dir.fun,"functions.R",sep="")) # funciones para leer .dat y .rep
source(paste(dir.fun,"Fn_PBRs.R",sep="")) # funciones para leer .dat y .rep


setwd(dir.1)
system("~/admb-12.2/admb MTT0621")
system("./MTT0621")
#------------------------------------------------
data.0  <- lisread(paste(dir.1,"MTT0621.dat", sep='/'));
names(data.0)<-str_trim(names(data.0), side="right")
data.1        <- data.0
rep     <- reptoRlist("MTT0621.rep")                                               
std     <- read.table("MTT0621.std",header=T,sep="",na="NA",fill=T) 

#============================================================#
# II. COMPOSICI?N EDAD DE LAS CAPTURAS                       #
#============================================================#
years  <- data.1$Ind[,1]                                                              
nyears <- data.1$nanos                                              
age  <-seq(5.5,20,0.5)                                            
nage<-length(age)                                            
#Proporci?n observada                                        
pobsF<-rep$pf_obs                                              
pobsR<-rep$pobs_RECLAN                                       
                                
#Proporci?n predicha                                         
ppredF<-rep$pf_pred                                          
ppredR<-rep$ppred_RECLAN                                   
                                   
#=================================================================#
# M?TODO de Francis
#=================================================================#
Nf1 <-60
Nr1 <-34
#-------------------------------------------#
#FLOTA
fanos<-years
fobs <-pobsF
fpre <-ppredF 
#RECLAS
ranos<-years
robs <-pobsR[rowSums(pobsR)>0,]
rpre <-ppredR[rowSums(pobsR)>0,]
#composicion de edad Flota
Of  <- rep(0,length(fanos))
Ef  <- rep(0,length(fanos))
vf  <- rep(0,length(fanos))
vNf <- rep(0,length(fanos))
#composicion de edad crucero de verano reclas
Or  <- rep(0,length(robs[,1]))
Er  <- rep(0,length(robs[,1]))
vr  <- rep(0,length(robs[,1]))
vNr <- rep(0,length(robs[,1]))
#-------------------------------------------#
#composicion de edad Flota
for(i in 1:length(fanos)){
	Of[i]  <- sum(fobs[i,]*age)
	Ef[i]  <- sum(fpre[i,]*age)
	vf[i]  <- sum(fpre[i,]*age^2)-Ef[i]^2
	vNf[i] <- vf[i]/Nf1}
#composicion de edad crucero de verano reclas
for(i in 1:length(robs[,1])){
	Or[i]  <- sum(robs[i,]*age)
	Er[i]  <- sum(rpre[i,]*age)
	vr[i]  <- sum(rpre[i,]*age^2)-Er[i]^2
	vNr[i] <- vr[i]/Nr1}
#--------------------------------------------#
wf  <- 1/var((Of-Ef)/sqrt(vNf)) #Flota
wr  <- 1/var((Or-Er)/sqrt(vNr)) #Reclas
Nf2 <- Nf1*wf                   # NM FLOTA
Nr2 <- Nr1*wr                   # NM RECLAS

#-----------------------------------------------------------------#
NM_Fran <- data.frame(nmF=c(Nf1,Nf2),nmR=c(Nr1,Nr2));NM_Fran
#-----------------------------------------------------------------#



#=================================================================#
# M?todo de Ianelli 2002
#=================================================================#
#Composici?n de edad de la FLOTA
Ofl <-ppredF[rowSums(pobsF)>0,]*(1-ppredF[rowSums(pobsF)>0,])
Efl <-(pobsF[rowSums(pobsF)>0,]-ppredF[rowSums(pobsF)>0,])^2
wfl <-rep(0,length(Ofl[,1]))
for(i in 1:length(Ofl[,1])){
	wfl[i] <-sum(Ofl[i,])/sum(Efl[i,])}

nmf_ari <-mean(wfl)                      # MEDIA ARITMETICA
nmf_geo <-exp(sum(log(wfl))/length(wfl)) # MEDIA GEOM?TRICA
nmf_arm <-1/mean(1/wfl)                  # MEDIA ARM?NICA

#------------------------------------------------------------
#Composici?n de edad Crucero de verano RECLAS
Ore <-ppredR[rowSums(pobsR)>0,]*(1-ppredR[rowSums(pobsR)>0,])
Ere <-(pobsR[rowSums(pobsR)>0,]-ppredR[rowSums(pobsR)>0,])^2
wre <-rep(0,length(Ore[,1]))
for(i in 1:length(Ore[,1])){	
	wre[i] <-sum(Ore[i,])/sum(Ere[i,])}
nmr_ari <-mean(wre)                      # MEDIA ARITMETICA
nmr_geo <-exp(sum(log(wre))/length(wre)) # MEDIA GEOM?TRICA
nmr_arm <-1/mean(1/wre)                  # MEDIA ARM?NICA
#------------------------------------------------------------
#------------------------------------------------------------
NM_Ian <- data.frame(nmF=c(nmf_ari,nmf_geo,nmf_arm),nmR=c(nmr_ari,nmr_geo,nmr_arm));NM_Ian
#------------------------------------------------------------
