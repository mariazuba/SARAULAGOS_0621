library(MASS)
library(tweedie)
library(statmod)
library(Hmisc)
library(alr4)
library(car)
library(nlme)
library(nortest) # Para realizar test de nomalidad

rm(list=ls())           # Limpieza de variables del ambiente, ?til cada vez que empezamos un nuevo an?lisis
setwd("C:/Users/macristina.perez/Documents/Recursos/Sardina Austral/CPUE_Los_Lagos/estandarizacion/elson")
#=========================================================================================================
# CARGAR ARCHIVO DE DATOS (.csv)
#=========================================================================================================
matriz<-read.csv("cpue019.csv",header=T,sep=",")
head(matriz)  # veo un poco los datos (t?tulos + 6 primeras filas)
str(matriz) 

# VARIABLES DEL ARCHIVO DE DATOS (Nombres de las columnas)
#---------------------------------------------------------
# cb        = capacidad de bodega (variables continua)
# rangocb   = categorias de rangos de la capacidad de bodega (1-6) 
# rangocb2   = categorias de rangos de la capacidad de bodega (1-4) 
# dfp       = d?as fuera de puerto (variable continua)
# ano       = a?os de datos ()
# trim      = trimestres (1=ene-marz, 2=abr-jun, 3=jul-sept, 4=oct-dic)
# cexito    = presencia/ausencia de captura (1=presencia, 0=ausencia)
matriz$cexito    <-ifelse(matriz$CPUE1>0,1,0)  # si CPUE1 es mayor a "0" escriba "1" de lo contrario escriba "0"

# zonapesca = 3 zonas en sentido norte-sur
matriz$zonapesca[matriz$AREA%in%c("A","B","C")] <- "1"    # Zona 1 (norte)
matriz$zonapesca[matriz$AREA%in%c("D","E","F")] <- "2"    # Zona 2 (centro)
matriz$zonapesca[matriz$AREA%in%c("G","H","I")] <- "3"    # Zona 3 (sur)

# zonapesca = 3 zonas en sentido Oeste-ESte (solo para explorar)
matriz$zonapesca2[matriz$AREA%in%c("A","D","G","I")] <- "1"    # Zona 1 (Isla)
matriz$zonapesca2[matriz$AREA%in%c("B","E","H")]     <- "2"    # Zona 2 (centro)
matriz$zonapesca2[matriz$AREA%in%c("C","F")]         <- "3"    # Zona 3 (continente)

matriz$cb2 <- cut(matriz$cb,c(0,25,50,100),labels=c("1","2","3"))

head(matriz)  #


tb1 <- ftable(matriz$ano,matriz$zonapesca);tb1               # n? REGISTROS X ZONA
tb2 <- ftable(matriz$ano,matriz$cexito);tb2    

####################################################################################
#===================================================================================
#
# 1. ELECCI?N DE LA ESTRUCTURA DE ERRORES Y FUNCI?N DE VINCULO
#===================================================================================
####################################################################################

--------------------------------------------------------------
  #Grafico para visualizar el tipo de distribuci?n de los datos
  #--------------------------------------------------------------

x11(width=7,height=5)
par(mfrow=c(1,2))
hist(matriz$CPUE1,xlab="CPUE (ton/viaje)",ylab="Frecuencia", main='',ylim=c(0,2500))
title("Histograma CPUE",cex.main=1)
box()
hist(matriz$lncpue2,nclass=25,xlab=" Log CPUE ",ylab="Frecuencia",main='')
title("Histograma log (CPUE)", cex.main=1 )
box()



#PARA ANALISIS DESCARTANDO LOS VIAJES SIN CAPTURA
CPUEpos   <- subset(matriz,CPUE1>1);names(CPUEpos)

x11(width=7,height=5)
par(mfrow=c(1,2))
hist(CPUEpos$CPUE1,xlab="CPUE (ton/vcp)",ylab="Frecuencia", main='',ylim=c(0,2500))
title("Histograma CPUE",cex.main=1)
box()
hist(CPUEpos$lncpue2,nclass=25,xlab=" Log CPUE ",ylab="Frecuencia",main='')
title("Histograma log (CPUE)", cex.main=1 )
box()

# Re-definamos la naturaleza de las variables seg?n se requiera (factores)
matriz$cb   <-factor(matriz$rangocb)
#matriz$cb2   <-factor(matriz$cb2)
matriz$mes <-factor(matriz$mes)
matriz$trim <-factor(matriz$trim)
matriz$zona <-factor(matriz$zonapesca)
#matriz$zona2 <-factor(matriz$zonapesca2)
matriz$Año <-factor(matriz$ano)
matriz$exipes <- !is.na(matriz$lncpue1)
matriz$cpue <- matriz$cpt
# Convertimos a factores en el caso de considerar solo viajes positivos (descartar sin captura)
CPUEpos$CB   <-factor(CPUEpos$rangocb)
CPUEpos$TRIM <-factor(CPUEpos$trim)
CPUEpos$ZONA <-factor(CPUEpos$zonapesca)
CPUEpos$ZONA2 <-factor(CPUEpos$zonapesca2)
CPUEpos$AÑO <-factor(CPUEpos$ano)
CPUEpos$exipes <- !is.na(CPUEpos$lncpue1)

#=========================================================================
#MODELOS
#=========================================================================

#---modelo DELTA-LOGNORMAL #----
options(contrasts=c(factor="contr.treatment","contr.sum"))
normal1 <- formula('lncpue2 ~ Año + mes  + cb + zona')

#Modelo solop viajes positivos
modelo1  <- glm(normal1, na.action=na.exclude, data=CPUEpos, family=gaussian(link = "identity"))
gamma2 <- formula('CPUE1 ~ Año + mes + cb + zona')
mod.tweedie <- formula('cpue ~ Año + mes + cb + zona')
mod.tweedie2 <- formula('cpue ~ Año + mes + cb + zona + zona:mes + zona:Año')


#Modelo con todos los viajes (incluye sin pesca)
modelo1  <- glm(normal3, na.action=na.exclude, data=matriz, subset= CPUE1>0, family=gaussian(link = "identity"))

# ---modelo DELTA-GAMMA #----
modelo2  <- glm(gamma2, family=Gamma(link=log),na.action=na.exclude, data=matriz,subset=cexito==1 & CPUE1>0)
# modelo2  <- glm(gamma1, family=Gamma(link=log),na.action=na.exclude, data=matriz, subset=cexito==1 & lncpue1>0 & !is.na(cb) & !is.na(zona))
# ---modelo BINOMIAL #----
modelo3  <- glm(binomial, family=binomial(link="logit"),na.action=na.exclude, data=matriz)
# ---modelo TWEEDIE #----

library(tweedie)
library(statmod)
modelo4  <- glm(mod.tweedie, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=matriz)
modelo5  <- glm(mod.tweedie2, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=matriz)

prom <-tapply(matriz$cpue,list(matriz$Año),mean,na.rm=TRUE);prom
var  <-tapply(matriz$cpue,list(matriz$Año),var,na.rm=TRUE);var

n <-tapply(matriz$cpue,list(matriz$Año),length)
x11()
plot(log(prom),log(var))
media=as.vector(prom)
varianza=as.vector(var)
model<-lm(log(varianza)~log(media));model
anova(model,test="F")


gamma2 <- formula('CPUE1 ~ Año + mes + cb + zona')
mod.tweedie <- formula('cpue ~ Año + mes + cb + zona')

PseudoR <- function (objeto)
{
  salida <-  1-(objeto$deviance/objeto$null.deviance)
  return(salida)
  end
}

modelo1$resumen  <- summary(modelo1); modelo1$resumen              # entrega los coeficientes del GLM
modelo1$anova    <- anova(modelo1,test="Chisq"); modelo1$anova     # An?lisis de devianza podemos elegir el mejor modelo
PseudoR(modelo1)

modelo2$resumen  <- summary(modelo2); modelo2$resumen              # entrega los coeficientes del GLM
modelo2$anova    <- anova(modelo2,test="Chisq"); modelo2$anova     # An?lisis de devianza podemos elegir el mejor modelo
PseudoR(modelo2)

modelo3$resumen  <- summary(modelo3); modelo3$resumen              # entrega los coeficientes del GLM
modelo3$anova    <- anova(modelo3,test="Chisq"); modelo3$anova     # An?lisis de devianza podemos elegir el mejor modelo
PseudoR(modelo3)

modelo4$resumen  <- summary(modelo4); modelo4$resumen              # entrega los coeficientes del GLM
modelo4$anova    <- anova(modelo4,test="Chisq"); modelo4$anova     # An?lisis de devianza podemos elegir el mejor modelo
PseudoR(modelo4)

modelo5$resumen  <- summary(modelo5); modelo5$resumen              # entrega los coeficientes del GLM
modelo5$anova    <- anova(modelo5,test="Chisq"); modelo5$anova   
PseudoR(modelo5)

#Despliegue de coefficientes
summary(modelo4)
names(modelo4)
win.graph()
par(mfrow=c(2,2))
termplot(modelo4,se=T,ylim="free")
x11()
par(mfcol=c(2,2))
plot(modelo4)