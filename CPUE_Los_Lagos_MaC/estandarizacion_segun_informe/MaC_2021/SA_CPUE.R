# DATOS SARDINA AUSTRAL  ####
library(dplyr)
library(ggplot2)
devtools::source_url("https://github.com/ale-yanez/RFunctions/blob/master/multiplot.R?raw=TRUE")

rm(list=ls())

# Directorio, Librerias, funciones, data inicial ####
setwd("C:/Users/macristina.perez/Documents/Recursos/Sardina Austral/CPUE_Los_Lagos_MaC/estandarizacion")
getwd()

# Carga data sin filtrar
data <- read.csv('Bitacora_comercial_2020.csv', header=T, sep=',')
names(data)
dim(data)
#7570   29
#7752   30
temp <- data %>% tbl_df() %>% 
  select(NOMBRE_BARCO,CAPACIDAD_BODEGA,TIPO_EMB,FECHA_HORA_RECALADA,FECHA_HORA_ZARPE,FECHA_LANCE, LATITUD, LONGITUD, PESO_ton, AÑO, MES, AREA,TRIMESTRE)# %>% 
  #filter(COD_PESQUERIA==50)#filter(cpue>0 & cpue<10000) # !(year==2015 &
dim(temp)
#  7570  13 
#  7752  13 
temp$AREA <- as.factor(temp$AREA)
summary(temp$AREA)
unique(temp$AREA)


temp$Zonapesca[temp$AREA%in%c("A","B","C")] <- "1"    # Zona 1 (norte)      
temp$Zonapesca[temp$AREA%in%c("D","E","F")] <- "2"    # Zona 2 (centro)
temp$Zonapesca[temp$AREA%in%c("G","H","I")] <- "3"    # Zona 3 (sur)

temp$Zonapesca<- as.factor(temp$Zonapesca)    # Define la zona como factor

class(temp$Zonapesca)

# TABLAS
# ------------ agrupamos los dias fuera de puerto por ano ----------------
vcp <- temp %>% group_by(Zonapesca, AÑO) %>% #filter(COD_PESQUERIA %in% c(5,92)) %>%
  summarise (VCPn = n())
vcp
#write.csv(vcp, file = "vcp_zona.csv", row.names = T)

ncatch <- temp %>% group_by(Zonapesca, AÑO) %>% 
  summarise (catch = mean(PESO_ton, na.rm = TRUE),
             catchsd = sd(PESO_ton, na.rm = TRUE),
             catchn = n(),
             catche = qnorm(0.975)*catchsd/sqrt(catchn),
             catchl = catch - catche,
             catchr = catch + catche)
ncatch
#write.csv(ncatch, file = "catch_prom_zona.csv", row.names = T)

ncpue <- temp %>% group_by(AÑO) %>% 
  summarise (f1 = sum(PESO_ton, na.rm = TRUE),
             f2 = n(),
             cpue = f1/f2,
             cpuen = n())
ncpue
#write.csv(ncpue, file = "ncpue.csv", row.names = T)
cpue<-ncpue$cpue
vcp<-ncpue$f2

# FACTOR Embarcacion ####
temp$embarcacion[temp$CAPACIDAD_BODEGA <= 30] <- 1
temp$embarcacion[temp$CAPACIDAD_BODEGA > 30 & temp$CAPACIDAD_BODEGA == 40] <- 2
temp$embarcacion[temp$CAPACIDAD_BODEGA > 40 & temp$CAPACIDAD_BODEGA == 50] <- 3
temp$embarcacion[temp$CAPACIDAD_BODEGA > 50 & temp$CAPACIDAD_BODEGA == 60] <- 4
temp$embarcacion[temp$CAPACIDAD_BODEGA > 60 & temp$CAPACIDAD_BODEGA == 70] <- 5
temp$embarcacion[temp$CAPACIDAD_BODEGA >= 70] <-6

temp$embarcacion<- as.factor(temp$embarcacion)
class(temp$embarcacion)
summary(temp)

# AS.FACTOR ####
temp$embarcacion <- as.factor(temp$embarcacion)
temp$MES <- as.factor(temp$MES)
temp$TRIMESTRE <- as.factor(temp$TRIMESTRE)
temp$Zonapesca <- as.factor(temp$Zonapesca)
temp$AÑO <- as.factor(temp$AÑO)
temp$cpue <- as.factor(temp$PESO_ton)

############# ESTANDARIZACIÓN #############

temp2 <- temp %>% tbl_df() %>% 
  select(AÑO, TRIMESTRE, cpue, Zonapesca, embarcacion, MES, PESO_ton) %>% 
  filter(!is.na(TRIMESTRE)) %>% 
  filter(!is.na(Zonapesca)) %>% 
  filter(!is.na(embarcacion))

dim(temp2)
#7570 
#7570
#7570
#3257
#3416  

############### modelacion cpue ###################

PseudoR <- function (objeto)
{
  salida <-  1-(objeto$deviance/objeto$null.deviance)
  return(salida)
  end
}

mod1 <- glm(log(cpue) ~ as.factor(year) + as.factor(Trim) + 
             as.factor(Zona) + as.factor(embarcacion), data=temp2, family=gaussian)

mod1$resumen <- summary(mod1)
mod1$anova   <- anova(mod1,test="F")
mod1$resumen
mod1$anova
PseudoR(mod1)

x11()
par(mfrow=c(2,3))
termplot(mod1, se=T, ylim="free")

drop1(mod1, test="F")
AIC(mod1)


mod2<-glm( cpue ~ as.factor(year) + as.factor(Trim)+as.factor(Zona)+
             as.factor(embarcacion), data=temp2, family=Gamma(link = "log"))

mod2$resumen <- summary(mod2)
mod2$anova   <- anova(mod2,test="F")
mod2$resumen
mod2$anova
PseudoR(mod2)

x11()
par(mfrow=c(2,2))
termplot(mod2, se=T, ylim="free")

drop1(mod2, test="F")
AIC(mod2)

mod3<-glm( cpue ~ as.factor(year) + as.factor(Trim)+as.factor(Zona)+
             as.factor(year)*as.factor(Zona)  , data=temp2, family=Gamma(link = "log"))

mod3$resumen <- summary(mod3)
mod3$anova   <- anova(mod3,test="F")
mod3$resumen
mod3$anova
PseudoR(mod3)

x11()
par(mfrow=c(2,2))
termplot(mod3, se=T, ylim="free")

drop1(mod3, test="F")
AIC(mod3)



x11()
par(mfcol=c(2,2))
Residuals <- resid(mod3)
xfit<-seq(min(Residuals),max(Residuals),length=40) 
yfit<-dnorm(xfit) 
hist(Residuals, freq=FALSE, xlab="Std. residuals", ylab="Frequency")
lines(xfit, yfit)
qqnorm(Residuals,xlab="Theoric quantiles", ylab="Std. residuals")
lines(xfit, xfit)
plot(fitted(mod3),Residuals, xlab="Predicted values", ylab="Std. residuals")

x11()
par(mfcol=c(2,2))
plot(mod3)

years<-seq(1998,2019,1)
cpue_prom<- ncpue2$cpue
cpue_mod1   <- exp(c(mod1$resumen$coefficients[1,1],  
                     mod1$resumen$coefficients[1,1] + mod1$resumen$coefficients[2:length(years),1]))
cpue_mod2 <- exp(c(mod2$resumen$coefficients[1,1],  
                     mod2$resumen$coefficients[1,1] + mod2$resumen$coefficients[2:length(years),1]))
cpue_mod3 <- exp(c(mod3$resumen$coefficients[1,1],  
                     mod3$resumen$coefficients[1,1] + mod3$resumen$coefficients[2:length(years),1]))

X11()
plot(years,cpue_prom,type="l",ylab="CPUE",xlab="Año",ylim=c(0,3000),lwd=2,lty=2)
lines(years,cpue_mod1,type="l",col=5,lwd=2,lty=2)
lines(years,cpue_mod2,type="l",col=10,lwd=2,lty=2)
lines(years,(cpue_mod3/3),type="l",col=7,lwd=2,lty=2)
points(year,cpue_prom,col=1,lwd=2)
legend (2015,3000,c("CPUE promedio","CPUE Mod1","CPUE Mod2","CPUE Mod3"),
        lwd=c(1,1,1,1),lty=c(2,2,2,2),col=c(1,5,10,7))

