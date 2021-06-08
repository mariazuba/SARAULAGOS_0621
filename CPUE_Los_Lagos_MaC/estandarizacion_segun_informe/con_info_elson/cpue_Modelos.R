
library(MASS)
library(tweedie)
library(statmod)
library(Hmisc)
library(alr4)
library(car)
library(nlme)
library(nortest) # Para realizar test de nomalidad

rm(list=ls())           # Limpieza de variables del ambiente, ?til cada vez que empezamos un nuevo an?lisis
setwd("C:/Users/macristina.perez/Documents/GitHub/SARAULAGOS_0621/CPUE_Los_Lagos_MaC/estandarizacion_segun_informe/con_info_elson")
#=========================================================================================================
# CARGAR ARCHIVO DE DATOS (.csv)
#=========================================================================================================
matriz<-read.csv("cpue019.csv",header=T,sep=",")
names(matriz)
head(matriz)  # veo un poco los datos (t?tulos + 6 primeras filas)
str(matriz) 

#matriz$logcpue <- log(matriz$CPUE1)
#matriz$logcpue2 <- log(matriz$CPUE2)

#write.csv(matriz$logcpue2, file = "matriz$logcpue2.csv", row.names = T)

# TABLAS

# ------------ agrupamos los dias fuera de puerto por ano ----------------
vcp <- matriz %>% group_by(zonapesca, ano) %>% #filter(COD_PESQUERIA %in% c(5,92)) %>%
  summarise (VCPn = n())
vcp
#write.csv(vcp, file = "vcp_zona.csv", row.names = T)

ncatch <- matriz %>% group_by(zonapesca, ano) %>% 
  summarise (catch = sum(cpt, na.rm = TRUE),
             catchsd = sd(cpt, na.rm = TRUE),
             catchn = n(),
             catche = qnorm(0.975)*catchsd/sqrt(catchn),
             catchl = catch - catche,
             catchr = catch + catche)
ncatch
#write.csv(ncatch, file = "catch_prom_zona.csv", row.names = T)

ncpue <- matriz %>% group_by(ano) %>% 
  summarise (f1 = sum(cpt, na.rm = TRUE),
             f2 = n(),
             cpue = f1/f2,
             cpuen = n())
ncpue
#write.csv(ncpue, file = "ncpue.csv", row.names = T)
cpue<-ncpue$cpue
vcp<-ncpue$f2

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
#write.csv(matriz$cexito, file = "matriz$cexito.csv", row.names = T)

# zonapesca = 3 zonas en sentido norte-sur
matriz$zonapesca[matriz$AREA%in%c("A","B","C")] <- "1"    # Zona 1 (norte)
matriz$zonapesca[matriz$AREA%in%c("D","E","F")] <- "2"    # Zona 2 (centro)
matriz$zonapesca[matriz$AREA%in%c("G","H","I")] <- "3"    # Zona 3 (sur)

# zonapesca = 3 zonas en sentido Oeste-ESte (solo para explorar)
matriz$zonapesca2[matriz$AREA%in%c("A","D","G","I")] <- "1"    # Zona 1 (Isla)
matriz$zonapesca2[matriz$AREA%in%c("B","E","H")]     <- "2"    # Zona 2 (centro)
matriz$zonapesca2[matriz$AREA%in%c("C","F")]         <- "3"    # Zona 3 (continente)

matriz$cb2 <- cut(matriz$cb,c(0,25,50,100),labels=c("1","2","3"))
#write.csv(matriz$cb2, file = "matriz$cexito.csv", row.names = T)
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

#x11(width=7,height=5)
#par(mfrow=c(1,2))

ppi<-300
png("figuras/hist_todos.png", width=10*ppi, height=4.5*ppi, res=ppi)
par(mfrow=c(1,2))
nbb = 30
hist(matriz$CPUE1,xlab="CPUE (ton/viaje)",ylab="Frecuencia", main='',ylim=c(0,2500))
title("Histograma CPUE",cex.main=1)
box()
hist(matriz$lncpue2,nclass=25,xlab=" Log CPUE ",ylab="Frecuencia",main='')
title("Histograma log (CPUE)", cex.main=1 )
box()
dev.off()

#plot_2 <- ggarrange(p1, p2, p4, ncol = 1, nrow = 3, align = "v", common.legend = F, legend = "right")
#ggsave(plot_2, filename = "figuras/figura2_b.png", width=7, height=8, dpi=300)

#PARA ANALISIS DESCARTANDO LOS VIAJES SIN CAPTURA
CPUEpos   <- subset(matriz,CPUE1>1)
summary(CPUEpos)

#x11(width=7,height=5)
ppi<-300
png("figuras/hist_positivos.png", width=10*ppi, height=4.5*ppi, res=ppi)
par(mfrow=c(1,2))
nbb = 30
hist(CPUEpos$CPUE1,xlab="CPUE (ton/vcp)",ylab="Frecuencia", main='',ylim=c(0,2500))
title("Histograma CPUE",cex.main=1)
box()
hist(CPUEpos$lncpue2,nclass=25,xlab=" Log CPUE ",ylab="Frecuencia",main='')
title("Histograma log (CPUE)", cex.main=1 )
box()
dev.off()

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
CPUEpos$cb   <-factor(CPUEpos$rangocb)
CPUEpos$trim <-factor(CPUEpos$trim)
CPUEpos$zona <-factor(CPUEpos$zonapesca)
#CPUEpos$ZONA2 <-factor(CPUEpos$zonapesca2)
CPUEpos$Año <-factor(CPUEpos$ano)
CPUEpos$exipes <- !is.na(CPUEpos$lncpue1)


# Formulas
# ---------
#Modelo Solo viajes positivos
#normal1 <- formula('lncpue1 ~ A?O + TRIM  + CB + ZONA')

#Modelo todos los viajes
#solo datos positivos
#binomial <- formula('cexito ~  Año + mes + cb + zona') # %0.09

normal1 <- formula('lncpue1 ~ Año + mes  + cb + zona') #0.19%
normal2 <- formula('lncpue1 ~ Año + trim  + cb + zona + cb:zona') #0.17%             #incluye una interacci?n
normal3 <- formula('lncpue1 ~ Año + trim  + cb + zona + zona:mes + zona:Año') #0.25%  #incluye dos interacci?nes


gamma1 <- formula('CPUE1 ~ mes:Año + Año + cb + zona') #0.10
gamma2 <- formula('CPUE1 ~ Año + mes + cb + zona') #0.12
gamma3 <- formula('CPUE1 ~ Año + trim  + cb + zona + zona:mes + zona:Año') #0.14
gamma4 <- formula('CPUE1 ~ Año + mes + cb + zona + Año:zona + cb:zona') #0.12

#datos totales positivos y ceros
mod.tweedie1 <- formula('cpue ~ Año + mes + cb + zona') # 0.16%
mod.tweedie2 <- formula('cpue ~ Año + mes + cb + zona + zona:mes + zona:Año') # 0.18%
#mod.tweedie3 <- formula('cpue ~ Año + mes + cb + zona +  cb:zona + trim:zona')
#mod.tweedie4 <- formula('cpue ~ Año + trim + cb + zona') # 0.16%

# Codigo para encontrar el valor optimo (var.power) en el modelo tweedie
#metodo interpolacion
out <- tweedie.profile( matriz$CPUE1 ~ 1, p.vec=seq(1.2, 1.965, length=10),
do.plot=TRUE, method="interpolation", do.smooth=TRUE, do.ci=TRUE)


#=========================================================================
#MODELOS
#=========================================================================

#---modelo DELTA-LOGNORMAL #----
options(contrasts=c(factor="contr.treatment","contr.sum"))

#Modelo solo viajes positivos
#modelo1  <- glm(normal3, na.action=na.exclude, data=CPUEpos, family=gaussian(link = "identity")) # 0.19%

#Modelo con todos los viajes (incluye sin pesca)
modelo1  <- glm(normal3, na.action=na.exclude, data=matriz, family=gaussian(link = "identity"))

 # ---modelo DELTA-GAMMA #----
#Modelo solo viajes positivos
modelo2  <- glm(gamma4, family=Gamma(link=log),na.action=na.exclude, data=matriz, subset=cexito==1 & CPUE1>0)
# modelo2  <- glm(gamma1, family=Gamma(link=log),na.action=na.exclude, data=matriz, subset=cexito==1 & lncpue1>0 & !is.na(cb) & !is.na(zona))

# ---modelo BINOMIAL #----
#modelo3  <- glm(binomial, family=binomial(link="logit"),na.action=na.exclude, data=matriz)

library(tweedie)
library(statmod)

# ---modelo TWEEDIE #----
modelo4  <- glm(mod.tweedie1, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=matriz) #0.16
modelo5  <- glm(mod.tweedie2, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=matriz) #0.18
#modelo6  <- glm(mod.tweedie3, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=matriz) #0.17
#modelo7  <- glm(mod.tweedie2, family=tweedie(var.power=1.635903, link.power=0), na.action=na.omit, data=CPUEpos) #0.20

years<-seq(2007,2020,1)

cpue_mod2   <- exp(c(modelo2$resumen$coefficients[1,1],  
                     modelo2$resumen$coefficients[1,1] + modelo2$resumen$coefficients[2:length(years),1]))


cpue_mod4   <- exp(c(modelo4$resumen$coefficients[1,1],  
                     modelo4$resumen$coefficients[1,1] + modelo4$resumen$coefficients[2:length(years),1]))

cpue_mod5   <- exp(c(modelo5$resumen$coefficients[1,1],  
                     modelo5$resumen$coefficients[1,1] + modelo5$resumen$coefficients[2:length(years),1]))

#MODELO 2, 4 Y 5 


# Se debe prestar atenci?n a:
# a) test de normalidad para la variable respuesta
# b) los test de significaci?n para los estimadores del modelo
# c) La cantidad de varianza explicada por el modelo (devianza D2)
# d) analisis de residuos (normalidad)

# Test para el tipo de distribuci?n de la varaible respuesta (Kolmogorv & Smirnov)

#ks.test(matriz$lncpue2,matriz$Año,"pgamma")

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

PseudoR <- function (objeto)
{
  salida <-  1-(objeto$deviance/objeto$null.deviance)
  return(salida)
  end
}


#Despliegue de coefficientesmodelo
################################################################################################
# Almacenamos todo el an?lisis estad?stico en el mismo objeto pero dentro de 
# otras variables (resumen y anova)
################################################################################################
modelo1$resumen  <- summary(modelo1); modelo1$resumen              # entrega los coeficientes del GLM
modelo1$anova    <- anova(modelo1,test="Chisq"); modelo1$anova     # An?lisis de devianza podemos elegir el mejor modelo
PseudoR(modelo1)
AIC(modelo1)

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

modelo6$resumen  <- summary(modelo6); modelo6$resumen              # entrega los coeficientes del GLM
modelo6$anova    <- anova(modelo6,test="Chisq"); modelo6$anova   
PseudoR(modelo6)

modelo7$resumen  <- summary(modelo7); modelo7$resumen              # entrega los coeficientes del GLM
modelo7$anova    <- anova(modelo7,test="Chisq"); modelo7$anova   
PseudoR(modelo7)


#Despliegue de coefficientes
summary(modelo1)
names(modelo1)
win.graph()

ppi<-300
png("figuras/mod5_con dataposi.png", width=10*ppi, height=6*ppi, res=ppi)
par(mfrow=c(2,2))
nbb = 30
par(mfrow=c(2,2))
termplot(modelo7,se=T,ylim="free")
dev.off()

x11()
par(mfcol=c(2,2))
plot(modelo1)


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
modelo7$qres <- rstudent(modelo7)
res <- residuals(modelo7,type="deviance")


x11(width=7,height=5)
par(mfrow=c(1,2))

ppi<-300
png("figuras/mod5_res_condataposi.png", width=10*ppi, height=4.5*ppi, res=ppi)
par(mfrow=c(1,2))
nbb = 30
hist(modelo7$qres, freq=FALSE, xlab='Residuales Estdar', ylab='Densidad', main=' ', xlim=c(-5,5),ylim=c(0,0.6))
curve(dnorm(x),col="gray25", lty=1, lwd=1.5, add=TRUE)
box()
qqnorm(res,ylab="Residuales",xlab="cuantiles Normal Estandar", main=' ')
qqline(modelo7$qres)
dev.off()


# Grafico Serie cpue total Log-Normal
# ------------------------------
modelo <- modelo4
predict.indice <- data.frame(Año=levels(matriz$Año),Normal=NA) 
year <- c(0,summary(modelo4)$coefficients[2:dim(predict.indice)[[1]]])
predict.indice$Normal <- exp(summary(modelo1)$coefficients[1] + year + (sigmaHat(modelo4)/2))
predict.indice$Metodo <- "Tradicional"
rm(year)
graphics.off()
plot(as.numeric(as.character(predict.indice$Año)), predict.indice$Normal, type="b", xlab='Años', ylab='CPUE')



#Grafico Serie cpue total
# -----------------
modelo <- modelo4
tspan <- seq(min(as.numeric(as.character(matriz$Año))), max(as.numeric(as.character(matriz$Año)))) #Intervalo de tiempo
year <- tspan; year[1] <- 0
year[2:length(tspan)] <- summary(modelo)$coefficients[2:length(tspan)]
years <- exp(summary(modelo)$coefficients[1] + year + (sigmaHat(modelo)/2))
years <- exp(summary(modelo)$coefficients[1] + year)
predict.indice$Gamma <- years
length(years)
graphics.off()
plot(tspan, years, type="b", xlab='Años', ylab='CPUE',ylim=c(0,4))
#plot(tspan[2:length(tspan)], years[2:length(years)], type="b", xlab='A?os', ylab='CPUE')
lines(tspan,predict.indice$Normal,col="blue",type="b")
legend("topright",fill=c("black","blue"),legend=c("Gamma","Normal"),bty="n")

