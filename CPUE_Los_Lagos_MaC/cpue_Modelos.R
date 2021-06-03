
library(MASS)
library(tweedie)
library(statmod)
library(Hmisc)
library(alr3)
library(car)
library(nlme)
library(nortest) # Para realizar test de nomalidad

rm(list=ls())           # Limpieza de variables del ambiente, ?til cada vez que empezamos un nuevo an?lisis
setwd("D:/RESPALDOS/Proyectos/AUSTRAL SARDINA/cpue/2019") # le indico el directorio de trabajo 
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




# Test de normalidad para lncpue

#pearson.test(matriz$cpue)
#ad.test(matriz$lncpue1)
#cvm.test(matriz$lncpue1)


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
#CPUEpos$CB   <-factor(CPUEpos$rangocb)
#CPUEpos$TRIM <-factor(CPUEpos$trim)
#CPUEpos$ZONA <-factor(CPUEpos$zonapesca)
#CPUEpos$ZONA2 <-factor(CPUEpos$zonapesca2)
#CPUEpos$A?O <-factor(CPUEpos$ano)
#CPUEpos$exipes <- !is.na(CPUEpos$lncpue1)


# Formulas
# ---------
#Modelo Solo viajes positivos
#normal1 <- formula('lncpue1 ~ A?O + TRIM  + CB + ZONA')

#Modelo todos los viajes

binomial <- formula('cexito ~  Año + mes + cb + zona')

normal <- formula('lncpue2 ~ Año + mes  + cb + zona')
normal2 <- formula('lncpue2 ~ Año + trim  + cb + zona + cb:zona')               #incluye una interacci?n
normal3 <- formula('lncpue2 ~ Año + trim  + cb + zona + zona:mes + zona:Año')   #incluye dos interacci?nes


gamma <- formula('CPUE1 ~ mes:Año + Año + cb + zona')
gamma2 <- formula('CPUE1 ~ Año + mes + cb + zona')
gamma3 <- formula('CPUE1 ~ Año + trim  + cb + zona + zona:mes + zona:Año')
gamma4 <- formula('CPUE1 ~ Año + mes + cb + zona + Año:zona + cb:zona')


mod.tweedie <- formula('cpue ~ Año + mes + cb + zona')
mod.tweedie2 <- formula('cpue ~ Año + mes + cb + zona + zona:mes + zona:Año')
mod.tweedie3 <- formula('cpue ~ Año + mes + cb + zona +  cb:zona + trim:zona')

# Codigo para encontrar el valor optimo (var.power) en el modelo tweedie
#metodo interpolacion
out <- tweedie.profile( matriz$CPUE1 ~ 1, p.vec=seq(1.2, 1.965, length=10),
do.plot=TRUE, method="interpolation", do.smooth=TRUE, do.ci=TRUE)


#=========================================================================
#MODELOS
#=========================================================================

#---modelo DELTA-LOGNORMAL #----
options(contrasts=c(factor="contr.treatment","contr.sum"))

#Modelo solop viajes positivos
# modelo1  <- glm(normal1, na.action=na.exclude, data=CPUEpos, family=gaussian(link = "identity"))

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

# Se debe prestar atenci?n a:
# a) test de normalidad para la variable respuesta
# b) los test de significaci?n para los estimadores del modelo
# c) La cantidad de varianza explicada por el modelo (devianza D2)
# d) analisis de residuos (normalidad)

# Test para el tipo de distribuci?n de la varaible respuesta (Kolmogorv & Smirnov)

ks.test(matriz$lncpue2,matriz$Año,"pgamma")

# (McCullangh & Nelde,1989, Gunnar 1996) para gamma pendiente cercana a 2 entre log mean y log var
prom <-tapply(matriz$cpue,list(matriz$Año),mean,na.rm=TRUE);prom
var  <-tapply(matriz$cpue,list(matriz$Año),var,na.rm=TRUE);var

n <-tapply(matriz$cpue,list(matriz$Año),length)
x11()
plot(log(prom),log(var))
#help(lm)
media=as.vector(prom)
varianza=as.vector(var)
model<-lm(log(varianza)~log(media));model
anova(model,test="F")

##########################################################################################################################
# DEVIANZA
#-----------
# La devianza nos da una idea de la variabilidad de los datos
# para obtener una medida de variabilidad explicada por el modelo, se compara la "Devianza Nulo" (Dnull)
# con a devianza resiudal (Dfull)
# Esto es una medida de cu?nto de la variabilidad de la variable respuesta (CPUE) no es explicada por el modelo.
# % devianza explicada= (Dnull ? Dfull)/(Dnull)
# =(202932-164584)/202932 (ejemplo jurel)
# =18,9%
##########################################################################################################################
# Analisis de devianza: 
# esta basada la seleccion del mejor conjuntos de predictores
# -------------------------------------------------------------
source('anadev.r')
anadev(model.obj=modelo2, model.form=gamma, stat='F', type=c(3, 'LR'))

anova(modelo2,test="F")



#Despliegue de coefficientesmodelo
################################################################################################
# Almacenamos todo el an?lisis estad?stico en el mismo objeto pero dentro de 
# otras variables (resumen y anova)
################################################################################################
modelo1$resumen  <- summary(modelo1); modelo1$resumen              # entrega los coeficientes del GLM
modelo1$anova    <- anova(modelo1,test="Chisq"); modelo1$anova     # An?lisis de devianza podemos elegir el mejor modelo

modelo2$resumen  <- summary(modelo2); modelo2$resumen              # entrega los coeficientes del GLM
modelo2$anova    <- anova(modelo2,test="Chisq"); modelo2$anova     # An?lisis de devianza podemos elegir el mejor modelo

modelo3$resumen  <- summary(modelo3); modelo3$resumen              # entrega los coeficientes del GLM
modelo3$anova    <- anova(modelo3,test="Chisq"); modelo3$anova     # An?lisis de devianza podemos elegir el mejor modelo

modelo4$resumen  <- summary(modelo4); modelo4$resumen              # entrega los coeficientes del GLM
modelo4$anova    <- anova(modelo4,test="Chisq"); modelo4$anova     # An?lisis de devianza podemos elegir el mejor modelo

modelo5$resumen  <- summary(modelo5); modelo5$resumen              # entrega los coeficientes del GLM
modelo5$anova    <- anova(modelo5,test="Chisq"); modelo5$anova   

#Despliegue de coefficientes
summary(modelo4)
names(modelo4)
win.graph()
par(mfrow=c(2,2))
termplot(modelo4,se=T,ylim="free")
x11()
par(mfcol=c(2,2))
plot(modelo4)


##########################################################################################################################################
#=========================================================================================================================================
#
# 3. AN?LISIS DE LOS RESIDUOS
#=========================================================================================================================================
##########################################################################################################################################

# Los residuos son las diferencias entre los valores estimados por el modelo y los valores observados.
# Sin embargo, muchas veces se utilizan los residuos estandarizados, que tienen que seguir una distribuci?n normal
# Conviene analizar los siguientes gr?ficos:
#------------------------------------------
# a) Histograma de residuos
# b) Gr?fico de residuos frente a valores estimados.  Estos gr?ficos pueden indicar falta de linealidad,
#    heterocedasticidad (varianza no constante) y valores at?picos)
# c) El gr?fico probabil?stico de normalidad (qqplot), que permite constrastar la normalidad (simetr?a) de la distribuci?n de residuos.
#   OPCIONALMENTE:
# a1) Gr?ficos de residuos frente a variables explicativas.  Pueden ayudar a identificar 
#     si la falta de linealidad o la heterocedasticidad es debida a alguna variable explicativa
# b1) Gr?fico de los residuos frente al tiempo (u orden de medida).  Permite detectar cambios sistem?ticos en el muestreo.
# c1) Gr?fico de valores at?picos.  Existen tests que permiten detectar valores at?picos. 
#     Los ?ndices m?s comunes son el ?ndice de Cook y el apalancamiento o leverage
#------------------------------------------------------------------------------------------------------------------------------------------------

# Todos estos gr?ficos (y opcionalmente algunos tests estad?sticos complementarios) nos pueden ayudar en la evaluaci?n del modelo utilizado.
# En caso necesario, ser? preciso volver a plantear el modelo, tal vez utilizando una estructura de errores m?s adecuada, otra funci?n de v?nculo
# o incluso eliminando ciertos datos que pueden estar "sobre-influenciando" nuestro an?lisis.
#-------------------------------------------------------------------------------------------------------------------------------------------------

# Gr?ficas residuos
# ------------------
modelo4$qres <- rstudent(modelo5)
res <- residuals(modelo5,type="deviance")


x11(width=7,height=5)
par(mfrow=c(1,2))
hist(modelo4$qres, freq=FALSE, xlab='Residuales Estdar', ylab='Densidad', main=' ', xlim=c(-5,5),ylim=c(0,0.6))
curve(dnorm(x),col="gray25", lty=1, lwd=1.5, add=TRUE)
box()
qqnorm(res,ylab="Residuales",xlab="cuantiles Normal Estandar", main=' ')
qqline(modelo4$qres)



# Grafico Serie cpue total Log-Normal
# ------------------------------
modelo <- modelo4
predict.indice <- data.frame(A?o=levels(matriz$A?o),Normal=NA) 
year <- c(0,summary(modelo4)$coefficients[2:dim(predict.indice)[[1]]])
predict.indice$Normal <- exp(summary(modelo1)$coefficients[1] + year + (sigmaHat(modelo4)/2))
predict.indice$Metodo <- "Tradicional"
rm(year)
graphics.off()
plot(as.numeric(as.character(predict.indice$A?o)), predict.indice$Normal, type="b", xlab='A?os', ylab='CPUE')



#Grafico Serie cpue total
# -----------------
modelo <- modelo1
tspan <- seq(min(as.numeric(as.character(matriz$A?o))), max(as.numeric(as.character(matriz$A?o)))) #Intervalo de tiempo
year <- tspan; year[1] <- 0
year[2:length(tspan)] <- summary(modelo)$coefficients[2:length(tspan)]
years <- exp(summary(modelo)$coefficients[1] + year + (sigmaHat(modelo)/2))
years <- exp(summary(modelo)$coefficients[1] + year)
predict.indice$Gamma <- years
length(years)
graphics.off()
plot(tspan, years, type="b", xlab='A?os', ylab='CPUE',ylim=c(0,35))
#plot(tspan[2:length(tspan)], years[2:length(years)], type="b", xlab='A?os', ylab='CPUE')
lines(tspan,predict.indice$Normal,col="blue",type="b")
legend("topright",fill=c("black","blue"),legend=c("Gamma","Normal"),bty="n")
