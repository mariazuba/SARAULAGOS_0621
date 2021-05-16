SPRFmort<-function(R0,Fmort,talla,edad,Dat)
{
ntalla   <- length(talla)
nedad    <- length(edad)
nspr     <- length(Fmort) 
npro     <- matrix(nrow=ntalla,ncol=nedad)
npro[,1] <- Dat$Pre*R0
npr      <- matrix(nrow=ntalla,ncol=nedad)
npr[,1]  <- Dat$Pre*R0
spro 	 	 <- matrix(nrow=ntalla,ncol=nedad)
spr 	 	 <- matrix(nrow=ntalla,ncol=nedad)
ypr 	 	 <- matrix(nrow=ntalla,ncol=nedad)
SPR      <- rep(0,nspr)
YPR      <- rep(0,nspr)

for(j in 1:nspr)
	{

  for(i in 2:nedad)
   { 
  	npro[,i] 	<- Ta%*%(npro[,i-1]*exp(-Dat$M)) 
  	npr[,i] 	<- Ta%*%(npr[,i-1]*exp(-((Fmort[j]*Dat$Sel)+Dat$M))) 
   }
	   for(i in 1:nedad)
	  		{
	   	  spro[,i] 	<- npro[,i]*Dat$Mad*Dat$Wmed*exp(-Dat$M*Dat$Tspw)
	 		  spr[,i] 	<- npr[,i]*Dat$Mad*Dat$Wmed*exp(-((Fmort[j]*Dat$Sel)+Dat$M)*Dat$Tspw)
   	    ypr[,i] 	<- npr[,i]*(Fmort[j]*Dat$Sel)*Dat$Wmed*(1-exp(-((Fmort[j]*Dat$Sel)+Dat$M)))/((Fmort[j]*Dat$Sel)+Dat$M)
   	    }
	      		SPRo<-sum(colSums(spro))	
        		SPR[j]<-sum(colSums(spr))	
	      		YPR[j]<-sum(colSums(ypr))	
  }
  Pspr<-SPR/SPRo
	out 		<- cbind(Fmort,SPR,YPR,Pspr=round(Pspr,2))
	assign("out",out,pos=1)
}

#SPRFpbr <- function (Fpbrs)

 # F85 <- Fpbrs[1]
#  F80 <- Fpbrs[2]
 # F60 <- Fpbrs[3]
  #F55 <- Fpbrs[4]
  #F45 <- Fpbrs[5]
  #F30 <- Fpbrs[6]

  #SBF0 <- 0.
  #SBF85 <- 0.
  #SBF80 <- 0.
  #SBF60 <- 0.
  #SBF55 <- 0.
  #SBF45 <- 0.
  #SBF30 <- 0.
  
  #nsbf <- 6

 # Nspr[,1] <- matrix (0, nsbf, Dat$Pre*R0)





