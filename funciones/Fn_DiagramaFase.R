#********************************************************
# DIAGRAMA DE FASE  EQUILIBRIO 
#********************************************************
DiagramaFase<-function(name,FRMS,BRMS,SpB,SpBSE,ln_Fyr,ln_FSE,Year){


Fval     <- exp(ln_Fyr)/FRMS
Bval     <- SpB/BRMS
lastB    <- tail(SpB,1)
lastF    <- tail(Fval,1)

Fvalue   <- exp(ln_Fyr+(ln_FSE*ln_FSE)/2)
# Calculate confidence intervals
Qmult    <- -qnorm((1-(80/100))/2.0)
sbSE     <- tail(SpBSE,1)
sb95     <- c(lastB-Qmult*sbSE,lastB+Qmult*sbSE)
B95      <- sb95/BRMS

FvSE     <- tail(ln_FSE,1)
F95      <- c(lastF*exp(-Qmult*FvSE),lastF*exp(Qmult*FvSE))

###############################################################################
# Define the regions of the Phase diagram
BlimC    <- 0.5         #l?mite de colapso o agotamiento
BTargC   <- 1.0                #Brms
FTarg    <- 1.0                #Frms
#-------------------------------------------------------------------------------------------------
maxaxes  <- length(Year)
ifelse(length(maxaxes)>1,xupdown<-c(0,maxaxes[1]),xupdown<-c(0,3))
ifelse(length(maxaxes)>1,yupdown<-c(0,maxaxes[3]),yupdown<-c(0,6))
#-------------------------------------------------------------------------------------------------
greenx   <- c(BTargC,xupdown[2],xupdown[2],BTargC,BTargC)   # l?mites de sub-explotaci?n
greeny   <- c(0,0,FTarg,FTarg,0)                            # l?mites de sub-explotaci?n
yellowx  <- c(BlimC,BTargC,BTargC,BlimC,BlimC)              
yellowy  <- c(0,0,yupdown[2],yupdown[2],0)
yellowx2 <- c(BTargC,xupdown[2],xupdown[2],BTargC,BTargC)
yellowy2 <- c(BTargC,BTargC,yupdown[2],yupdown[2],BTargC)
redx     <- c(0,BlimC,BlimC,0,0)                             # l?mites de agotamiento y/o colapso
redy     <- c(0,0,yupdown[2],yupdown[2],0)                   # l?mites de agotamiento y/o colapso
#-----------------------------------------------------------------
leftFE   <- 0.90    # l?mite izquierdo de plenaExplotaci?n
rightFE  <- 1.35    # l?mite derecho de plenaExplotaci?n 
upFE     <- 1.10    # l?mite superior de plenaExplotaci?n
downFE   <- 0.45    # l?mite inferior de plenaExplotaci?n
#-----------------------------------------------------------------
orangex     <- c(BlimC,leftFE,leftFE,BlimC,BlimC)  # l?mite de sobre-explotaci?n
orangey     <- c(0,0,yupdown[2],yupdown[2],0)      # l?mite de sobre-explotaci?n
orangex2    <- c(BlimC,leftFE,leftFE,BlimC,BlimC)  # l?mite de sobre-explotaci?n y sobrepesca
orangey2    <- c(0,0,upFE,upFE,0)                  # l?mite de sobre-explotaci?n y sobrepesca
#-----------------------------------------------------------------
if (downFE < 0.1*FTarg) downFE <- 0.1*FTarg                                #l?mite inferior de plenaExplotaci?n
fullyEx  <- c(leftFE,leftFE,xupdown[2],xupdown[2],rightFE,rightFE,leftFE)  #l?mites de plenaExplotaci?n
fullyEy  <- c(0,upFE,upFE,downFE,downFE,0,0)                               #l?mites de plenaExplotaci?n
#-----------------------------------------------------------------
lastB    <- tail(Bval,1)
lastF    <- tail(Fval,1)
###############################################################################
col1   <-"gray90"; col2<-"gray80"; col3<-"#C0C0C0"; col4<-"gray47";col5<-"gray70"
usefont <- 2
#par(mfrow = c(1,1))
ifelse((nchar(name) > 1),upspace <- 1.5,upspace <- 0)
par(mai=c(0.45,0.45,0.05,0.05),oma=c(0.0,0,upspace,0.0))
par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=usefont,font=usefont)

plot(Bval,Fval,type="n",pch=16,cex=0.8,lwd=1,xlim=c(0,3),xaxs="i",
       ylim=c(0,4),yaxs="i",xlab="",ylab="")
#----------------------------------
polygon(redx,redy,col=col4)
polygon(yellowx2,yellowy2,col=col2)
polygon(yellowx,yellowy,col=col2)
polygon(greenx,greeny,col=col1)
polygon(fullyEx,fullyEy,col="white")  
polygon(orangex,orangey,col=col5)
polygon(orangex2,orangey2,col=col5)
#----------------------------------
abline(v=BTargC,lty=1,col="white")  
abline(v=BTargC,lty=2,col=1)
abline(h=FTarg,lty=1,col="white")
abline(h=FTarg,lty=2,col=1)
abline(v=BlimC,lty=1,col="white")  
abline(v=BlimC,lty=2,col=1)
#----------------------------------
title(xlab=list(expression("BD/BD"[RMS]), cex=1., font=usefont),
      ylab=list(expression("F/F"[RMS]), cex=1., font=usefont))
if (upspace > 0) mtext(name,side=3,outer=T,cex=1,font=usefont)

  lines(Bval,Fval,lwd=1,col=1)
  points(Bval,Fval,pch=16,cex=1.0)
  arrows(x0=B95[1],y0=lastF,x1=B95[2],y1=lastF,length=0.05,angle=90,col=4,lwd=1,code=3,lty=2)
  arrows(x0=lastB,y0=F95[1],x1=lastB,y1=F95[2],length=0.05,angle=90,col=4,lwd=1,code=3,lty=2)
  points(c(tail(Bval,1),tail(Bval[1])),c(tail(Fval,1),tail(Fval[1])),pch=19,cex=1,col=c(4,3))
  text(c(Bval[1],Bval[maxaxes]),c(Fval[1],Fval[maxaxes])+0.1,c(Year[1],Year[length(Year)]),cex=1)
  
  box()
  
}




