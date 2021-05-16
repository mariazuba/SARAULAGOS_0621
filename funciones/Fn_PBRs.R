
#============================================================================
# EQUILIBRIUM SSB PER RECRUIT  "%SPR(F)"
#============================================================================
SPRFmort<-function(R0,Fmort,Amax,Dat)
{
	nspr 	 	<- length(Fmort)
	npr 	 	<- rep(0,Amax)
	spr 	 	<- rep(0,nspr)
	ypr 	 	<- rep(0,nspr)
	npr[1] 	<- R0
	
	for(j in 1:nspr)
	{
		for(i in 2:Amax)
		{
		npr[i] 	<- npr[i-1]*exp(-(Dat$M+Fmort[j]*Dat$Sel[i-1]))
		npr[Amax] 	<- npr[Amax]/(1-exp(-(Dat$M+Fmort[j]*Dat$Sel[Amax])))
		}
	  		for(i in 1:Amax)
	  		{
	 		  spr[j] 	<- spr[j]+npr[i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+Fmort[j]*Dat$Sel[i])*Dat$Tspw)
	  		ypr[j] 	<- ypr[j]+Fmort[j]*Dat$Sel[i]*npr[i]*Dat$Wmed[i]*(1-exp(-(Dat$M+Fmort[j]*Dat$Sel[i])))/(Dat$M+Fmort[j]*Dat$Sel[i]) 
	  		}
	}
	Pspr 		<- spr/spr[1]
	out 		<- cbind(Fmort,spr,ypr,Pspr=round(Pspr,3))
	assign("out",out,pos=1)

}
#============================================================================================
SPR1<-function(R0,Fmort,Amax,Dat)
{
  nf      <- length(Fmort)
  npr 	 	<- rep(0,Amax)
  spr 	 	<- rep(0,nf)
  ypr     <- rep(0,nf)
  npr[1] 	<- R0
  
  for(j in 1:nf)
  {
      for(i in 2:Amax)
    {
      npr[i] 	<- npr[i-1]*exp(-(Dat$M+Fmort[j]*Dat$Sel[i-1]))
      npr[Amax] 	<- npr[Amax]/(1-exp(-(Dat$M+Fmort[j]*Dat$Sel[Amax])))
    }
    for(i in 1:Amax)
    {
      spr[j] 	<- spr[j]+npr[i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+Fmort[j]*Dat$Sel[i])*Dat$Tspw)
      ypr[j]   <- ypr[j]+Fmort[j]*Dat$Sel[i]*npr[i]*Dat$Wmed[i]*(1-exp(-(Dat$M+Fmort[j]*Dat$Sel[i])))/(Dat$M+Fmort[j]*Dat$Sel[i]) 
    }
  }
  sal<-NULL
  sal$R0   <- R0
  sal$Fobj <- Fmort
  sal$spr <-  spr
  sal$ypr <-  ypr
  sal$pspr <- spr/spr[1]
  return(npr)
}

####################################################################################
#    ESTIMACI?N DE Fpbr
####################################################################################
SPRFpbr <- function(Fpbr)
{
	F85 <- Fpbr[1]
	F80 <- Fpbr[2]
	F60 <- Fpbr[3]
	F55 <- Fpbr[4]
	F52 <- Fpbr[5]
	F50 <- Fpbr[6]
  F45 <- Fpbr[7]
	F40 <- Fpbr[8]
	F30 <- Fpbr[9]
	F325<-Fpbr[10]
	F425<-Fpbr[11]
	
	SBF0 <- 0.
	SBF85 <- 0.
	SBF80 <- 0.
	SBF60 <- 0.
	SBF55 <- 0.
	SBF52 <- 0.
	SBF50 <- 0.
  SBF45 <- 0.
	SBF40 <- 0.
	SBF30 <- 0.
	SBF325 <- 0.
	SBF425 <- 0.
	
	
	nsbf	<- 12
	
	Nspr  <- matrix(0,nsbf,Amax)

	for(j in 1:nsbf)
	{
	 Nspr[j,1] <- R0
	}
  
	for(i in 2:Amax)
	{
  	 Nspr[1,i] 	  <- Nspr[1,i-1]*exp(-Dat$M)
	 Nspr[1,Amax] <- Nspr[1,Amax]/(1-exp(-Dat$M))

	 Nspr[2,i] 	  <- Nspr[2,i-1]*exp(-(Dat$M+F85*Dat$Sel[i-1]))
	 Nspr[2,Amax] <- Nspr[2,Amax]/(1-exp(-(Dat$M+F85*Dat$Sel[Amax])))
		
	 Nspr[3,i] 	  <- Nspr[3,i-1]*exp(-(Dat$M+F80*Dat$Sel[i-1]))
	 Nspr[3,Amax] <- Nspr[3,Amax]/(1-exp(-(Dat$M+F80*Dat$Sel[Amax])))
		
	 Nspr[4,i] 	  <- Nspr[4,i-1]*exp(-(Dat$M+F60*Dat$Sel[i-1]))
	 Nspr[4,Amax] <- Nspr[4,Amax]/(1-exp(-(Dat$M+F60*Dat$Sel[Amax])))

	 Nspr[5,i] 	  <- Nspr[5,i-1]*exp(-(Dat$M+F55*Dat$Sel[i-1]))
 	 Nspr[5,Amax] <- Nspr[5,Amax]/(1-exp(-(Dat$M+F55*Dat$Sel[Amax])))
	 
	 Nspr[6,i] 	  <- Nspr[6,i-1]*exp(-(Dat$M+F52*Dat$Sel[i-1]))
	 Nspr[6,Amax] <- Nspr[6,Amax]/(1-exp(-(Dat$M+F52*Dat$Sel[Amax])))
		
	 Nspr[7,i] 	  <- Nspr[7,i-1]*exp(-(Dat$M+F50*Dat$Sel[i-1]))
	 Nspr[7,Amax] <- Nspr[7,Amax]/(1-exp(-(Dat$M+F50*Dat$Sel[Amax])))
		
   Nspr[8,i]     <- Nspr[8,i-1]*exp(-(Dat$M+F45*Dat$Sel[i-1]))
	 Nspr[8,Amax] <- Nspr[8,Amax]/(1-exp(-(Dat$M+F45*Dat$Sel[Amax])))
	 
	 Nspr[9,i] 	  <- Nspr[9,i-1]*exp(-(Dat$M+F40*Dat$Sel[i-1]))
	 Nspr[9,Amax] <- Nspr[9,Amax]/(1-exp(-(Dat$M+F40*Dat$Sel[Amax]))) 
		
	 Nspr[10,i] 	  <- Nspr[10,i-1]*exp(-(Dat$M+F30*Dat$Sel[i-1]))
 	 Nspr[10,Amax] <- Nspr[10,Amax]/(1-exp(-(Dat$M+F30*Dat$Sel[Amax])))
	
	 Nspr[11,i] 	  <- Nspr[11,i-1]*exp(-(Dat$M+F325*Dat$Sel[i-1]))
	 Nspr[11,Amax] <- Nspr[11,Amax]/(1-exp(-(Dat$M+F325*Dat$Sel[Amax])))
 
	 Nspr[12,i] 	  <- Nspr[12,i-1]*exp(-(Dat$M+F425*Dat$Sel[i-1]))
	 Nspr[12,Amax] <- Nspr[12,Amax]/(1-exp(-(Dat$M+F425*Dat$Sel[Amax])))
			 		
	 }
	for(i in 1:Amax)
	{
      SBF0  <- SBF0+Nspr[1,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-Dat$M*Dat$Tspw)
		  SBF85 <- SBF85+Nspr[2,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F85*Dat$Sel[i])*Dat$Tspw)
		  SBF80 <- SBF80+Nspr[3,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F80*Dat$Sel[i])*Dat$Tspw)
    	SBF60 <- SBF60+Nspr[4,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F60*Dat$Sel[i])*Dat$Tspw)
    	SBF55 <- SBF55+Nspr[5,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F55*Dat$Sel[i])*Dat$Tspw)
	   	SBF52 <- SBF52+Nspr[6,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F52*Dat$Sel[i])*Dat$Tspw)
	  	SBF50 <- SBF50+Nspr[7,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F50*Dat$Sel[i])*Dat$Tspw)
      SBF45 <- SBF45+Nspr[8,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F45*Dat$Sel[i])*Dat$Tspw)
		  SBF40 <- SBF40+Nspr[9,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F40*Dat$Sel[i])*Dat$Tspw)
      SBF30 <- SBF30+Nspr[10,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F30*Dat$Sel[i])*Dat$Tspw)
	  	SBF325 <- SBF325+Nspr[11,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F325*Dat$Sel[i])*Dat$Tspw)
		  SBF425 <- SBF425+Nspr[12,i]*Dat$Mad[i]*Dat$Wmed[i]*exp(-(Dat$M+F425*Dat$Sel[i])*Dat$Tspw)
		      	}
    sprpen <- (SBF85/SBF0-0.85)^2
	  sprpen <- sprpen+(SBF80/SBF0-0.80)^2
	  sprpen <- sprpen+(SBF60/SBF0-0.60)^2
    sprpen <- sprpen+(SBF55/SBF0-0.55)^2
	  sprpen <- sprpen+(SBF52/SBF0-0.52)^2
	  sprpen <- sprpen+(SBF50/SBF0-0.50)^2
  	sprpen <- sprpen+(SBF45/SBF0-0.45)^2
	  sprpen <- sprpen+(SBF40/SBF0-0.40)^2
    sprpen <- sprpen+(SBF30/SBF0-0.30)^2
	  sprpen <- sprpen+(SBF325/SBF0-0.325)^2
	  sprpen <- sprpen+(SBF425/SBF0-0.425)^2
	 }


