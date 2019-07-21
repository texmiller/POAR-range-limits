## Fit size-dependent vital rates and build MPM
library(R2jags)
library(mcmcplots)
library(popbio)

## load the Bayesian model output and other raw data miscellany
load("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/analysis/poar.demography.BAYES.RData")

longitude<- -100:-90

#### identify the JAGS object to draw final params
final.params<-poar.out$BUGSoutput$mean
### Paramater vectors for females and males
F.params<-c()
##SURVIVAL#########################
F.params[1]<-final.params$alpha.surv  #surv intercpt
F.params[2]<-final.params$beta.surv.size  #slope log size
F.params[3]<-final.params$beta.surv.long  #slope long
##GROWTH#########################
F.params[11]<-final.params$alpha.grow  #grow int
F.params[12]<-final.params$beta.grow.size  #slope log size
F.params[13]<-final.params$beta.grow.long  #slope long
F.params[14]<-final.params$alpha.gamma.grow  #NegBin size param
##GROWTH#########################
F.params[21]<-final.params$alpha.flow  #flow int
F.params[22]<-final.params$beta.flow.size  #slope log size
F.params[23]<-final.params$beta.flow.long  #slope long
##PANIC#########################
F.params[31]<-final.params$alpha.panic  #panic int
F.params[32]<-final.params$beta.panic.size  #slope log size
F.params[33]<-final.params$beta.panic.long  #slope long
##GERMINATION#########################
F.params[41]<-final.params$germ.mu  #F-dominant seed germination (twosex==F)
F.params[42]<-final.params$germ.beta1 #This is slope of SR dependence for 'maybe' tetrazolium viability (twosex==T)
#F.params[41]<-0.4253165 
#F.params[42]<- 9.6583863
#F.params[43]<- -12.9141399
##MISC#############################
F.params[51]<-200  #seeds per panicle ## hard coded from Fig 2 of resp surf paper, low density
F.params[52]<-0.5  #seed sex ratio  # assumed 1:1 re: Jason Goldman
F.params[61]<-quantile(poar$tillerN_t1,na.rm=T,probs=0.95)  #max size (# tillers), based on real data but not abs max

M.params<-c()
##SURVIVAL#########################
M.params[1]<-final.params$alpha.surv+final.params$beta.surv.sex  #surv intercpt
M.params[2]<-final.params$beta.surv.size  #slope log size
M.params[3]<-final.params$beta.surv.long+final.params$beta.surv.sex.long  #slope long
##GROWTH#########################
M.params[11]<-final.params$alpha.grow+final.params$beta.grow.sex  #grow int
M.params[12]<-final.params$beta.grow.size  #slope log size
M.params[13]<-final.params$beta.grow.long+final.params$beta.grow.sex.long  #slope long
M.params[14]<-final.params$alpha.gamma.grow  #NegBin size paraM
##GROWTH#########################
M.params[21]<-final.params$alpha.flow+final.params$beta.flow.sex  #flow int
M.params[22]<-final.params$beta.flow.size  #slope log size
M.params[23]<-final.params$beta.flow.long+final.params$beta.flow.sex.long  #slope long
##PANIC#########################
M.params[31]<-final.params$alpha.panic+final.params$beta.panic.sex  #panic int
M.params[32]<-final.params$beta.panic.size  #slope log size
M.params[33]<-final.params$beta.panic.long+final.params$beta.panic.sex.long  #slope long
##GERMINATION#########################
M.params[41]<-F.params[41]  #F-dominant seed germination (twosex==F)
#Hard-coded for now. Value from Fig 1 of response surf ms (max assuming no mate limitation)
M.params[42]<-F.params[42]  #intercept for sex ratio dependence (twosex==T)
#M.params[43]<-F.params[43] 
##MISC#############################
M.params[51]<-F.params[51]  #seeds per panicle ## hard coded from Fig 2 of resp surf paper, low density
M.params[52]<-F.params[52]  #seed sex ratio  # assumed 1:1 re: Jason Goldman
M.params[61]<-F.params[61]  #max size (# tillers), based on real data but not abs max

plot(seq(0,1,.01),invlogit(F.params[41]+F.params[42]*seq(0,1,.01)),ylim=c(0,1))
lines(seq(0,1,.01),invlogit(10+-20*seq(0,1,.01)))
#########################################################################################
## Define vital rate functions
#SURVIVAL AT SIZE X.
sx<-function(x,params,long){
  surv.mean<-params[1] + params[2]*log(x) + params[3]*long
  return(invlogit(surv.mean))
}
#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,long){
  grow.mean<-params[11] + params[12]*log(x) + params[13]*long
  grow<-dnbinom(x=y,mu=exp(grow.mean),size=params[14],log=F)
  truncLower<-dnbinom(x=0,mu=exp(grow.mean),size=params[14],log=F)
  truncUpper<-sum(dnbinom(x=params[61]:1000,mu=exp(grow.mean),size=params[14],log=F))
  return(grow/(1-(truncLower+truncUpper)))
}
#SURVIVAL*GROWTH
pxy<-function(x,y,params,long){
  sx(x,params,long) * gxy(x,y,params,long)
}
# PROBABILITY OF FLOWERING
Pfx<-function(x,params,long){
  flow.mean<-params[21] + params[22]*log(x) + params[23]*long
  return(invlogit(flow.mean))
}
#NUMBER OF PANICLES
Nfx<-function(x,params,long){
  panic.mean<-params[31] + params[32]*log(x) + params[33]*long
  return(exp(panic.mean))
}
#SEED GERMINATION
germ<-function(params,twosex,SR){
  if(twosex==F){return(invlogit(params[41]))}
  #if(twosex==F){return(invlogit(params[41]+params[42]*0.5+params[43]*(0.5^2)))}
  if(twosex==T){return(invlogit(params[41]+params[42]*SR))}
  #if(twosex==T){return(invlogit(params[41]+params[42]*SR+params[43]*(SR^2)))}
}
#FERTILITY--returns number of recruits
##Female offspring
Fertx.F<-function(x,params,long,twosex,SR=NULL){
  seedlings<-Pfx(x,params,long)*Nfx(x,params,long)*params[51]*germ(params,twosex,SR)*params[52]
  return(seedlings)
}
##Male offspring
Fertx.M<-function(x,params,long,twosex,SR=NULL){
  seedlings<-Pfx(x,params,long)*Nfx(x,params,long)*params[51]*germ(params,twosex,SR)*(1-params[52])
  return(seedlings)
}
################################################################################
## Put it all together
megamatrix<-function(F.params,M.params,long,twosex,SR=NULL){   
  
  matdim<-F.params[61]         ## bigmatrix dimension is max size
  y<-1:F.params[61]
  
  ## F-to-F growth/survival transition
  F.Tmat<-matrix(0,matdim,matdim)
  F.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=F.params,long=long))
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim,matdim)
  F.Fmat[1,1:matdim]<-Fertx.F(x=y,params=F.params,long=long,twosex,SR=SR) 
  
  ## M-to-M growth/survival transition
  M.Tmat<-matrix(0,matdim,matdim)
  M.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=M.params,long=long))
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim,matdim)
  M.Fmat[1,1:matdim]<-Fertx.M(x=y,params=F.params,long=long,twosex,SR=SR) 
  
  #M-to-F
  zero.mat<-matrix(0,matdim,matdim)

  # Put it all together
  MEGAmat<-cbind(rbind(F.Tmat+F.Fmat,  ##Female growth/survival + recruitment[1,1]
                       M.Fmat), ##Male recruitment [2,1]
                 rbind(zero.mat,   ##Females from males [1,2]
                       M.Tmat))   ##Male growth/survival

  return(list(MEGAmat=MEGAmat))
}
################################################################################
## Simulation for 2-sex model
#####################################################################
lambdaSim=function(F.params,M.params,long,max.yrs){

  #simulate and store lambdas
  matdim<-F.params[61]
  y<-1:matdim

  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/(matdim*2),(matdim*2))

  for(t in 1:max.yrs){ #Start loop
    ##Estimate panicle SR
    ## scalar multiplication to weight females by flowering prob
    flowering.females<-n0[1:matdim]*Pfx(x=y,param=F.params,long=long)
    ##Vector operation to sum female panicles
    F.panicles<-flowering.females%*%Nfx(x=y,param=F.params,long=long)
    ##Same for males
    flowering.males<-n0[(matdim+1):(matdim*2)]*Pfx(x=y,param=M.params,long=long)
    M.panicles<-flowering.males%*%Nfx(x=y,param=M.params,long=long)
    ##Panicle sex ratio (proportion female)
    OSRtracker[t]<-F.panicles/(F.panicles+M.panicles)
    SRtracker[t]<-sum(n0[1:matdim])
    
    #Store matrix
    MEGAmat<-megamatrix(F.params=F.params,M.params=M.params,long=long,twosex=T,SR=SRtracker[t])$MEGAmat

    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  
  #discard initial values (to get rid of transient)
  #burnin    <- round(max.yrs*0.1)
  #rtracker  <- rtracker[-c(1:burnin)]
  
  #Finish and return
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker))
  
}

## LTRE function
## store LTRE output in a matrix of columns for melt date and rows for parameters
poar.LTRE<-function(F.params,M.params,long,max.yrs,perturbation,
                    LTRE.params,dparam.dlong){
  
  #This is the output matrix
  F.dlam.dparam<-matrix(NA,ncol=length(long),nrow=length(LTRE.params))
  M.dlam.dparam<-matrix(NA,ncol=length(long),nrow=length(LTRE.params))
  F.LTRE.out<-matrix(NA,ncol=length(long),nrow=length(LTRE.params))
  M.LTRE.out<-matrix(NA,ncol=length(long),nrow=length(LTRE.params))
  lambda<-c()
  #Loop over longitudes
  for(i in 1:length(long)){
    ##estimate lambda for this longitude
    lambda[i]<-lambdaSim(F.params=F.params,M.params=M.params,long=long[i],max.yrs=max.yrs)$lambdatracker[max.yrs]
    
    ## loop over long-dependent parameters
    for(j in 1:length(LTRE.params)){
      ##estimate sensitivity of lambda to this parameter
      F.params.perturb<-F.params
      F.params.perturb[LTRE.params[j]]<-F.params[LTRE.params[j]]+perturbation
      F.lambda.perturb<-lambdaSim(F.params=F.params.perturb,M.params=M.params,long=long[i],max.yrs=max.yrs)$lambdatracker[max.yrs]
      F.dlam.dparam[j,i]<-(F.lambda.perturb-lambda[i])/perturbation
      F.LTRE.out[j,i]<-F.dlam.dparam[j,i]*F.params[dparam.dlong[j]]

      M.params.perturb<-M.params
      M.params.perturb[LTRE.params[j]]<-M.params[LTRE.params[j]]+perturbation
      M.lambda.perturb<-lambdaSim(F.params=F.params,M.params=M.params.perturb,long=long[i],max.yrs=max.yrs)$lambdatracker[max.yrs]
      M.dlam.dparam[j,i]<-(M.lambda.perturb-lambda[i])/perturbation
      M.LTRE.out[j,i]<-M.dlam.dparam[j,i]*M.params[dparam.dlong[j]]
    }
  }
  return(list(lambda=lambda,F.LTRE.out=F.LTRE.out,M.LTRE.out=M.LTRE.out,
              F.dlam.dparam=F.dlam.dparam,M.dlam.dparam=M.dlam.dparam))
}


test<-lambdaSim(F.params=F.params,M.params=M.params,long=0,max.yrs=100)
plot(1:100,test$lambdatracker)
plot(1:100,test$SRtracker,ylim=c(0,1))
lines(1:100,test$OSRtracker)

## Asymptotic results for the F-dominant model. This one is pretty easy, no simulation.

long.seq<-seq(min(poar$long.center),max(poar$long.center)*3,0.5)

lambda.long.Fdom<-c()
lambda.long.2sex<-c()
SR.long.2sex<-c()
OSR.long.2sex<-c()
max.yrs<-50
for(i in 1:length(long.seq)){
  lambda.long.Fdom[i]<-lambda(megamatrix(F.params=F.params,M.params=M.params,long=long.seq[i],twosex=F)$MEGAmat)
  lambda.long.2sex[i]<-lambdaSim(F.params=F.params,M.params=M.params,long=long.seq[i],max.yrs=max.yrs)$lambdatracker[max.yrs]
  SR.long.2sex[i]<-lambdaSim(F.params=F.params,M.params=M.params,long=long.seq[i],max.yrs=max.yrs)$SRtracker[max.yrs]
  OSR.long.2sex[i]<-lambdaSim(F.params=F.params,M.params=M.params,long=long.seq[i],max.yrs=max.yrs)$OSRtracker[max.yrs]
}

plot(long.seq+mean(poar$Longitude),(lambda.long.Fdom),type="l",ylim=c(0.8,7))
lines(long.seq+mean(poar$Longitude),(lambda.long.2sex),lty=2)
abline(1,0)

plot(long.seq+mean(poar$Longitude),SR.long.2sex,type="l",ylim=c(0,1))
lines(long.seq+mean(poar$Longitude),OSR.long.2sex,lty=2)

### LTRE analysis
LTRE.out<-poar.LTRE(F.params=F.params,
          M.params=M.params,
          long=long.seq,
          max.yrs=50,
          perturbation=.1,
          LTRE.params=c(1,11,21,31), ## intercepts of long-dependent vital rates
          dparam.dlong=c(3,13,23,33))  ## slopes of long-dependent vital rates


plot(long.seq+mean(poar$Longitude),LTRE.out$lambda,ylim=c(0,6))
abline(h=1)
abline(v=max(poar$long.center)+mean(poar$Longitude))

long.seq[which.min(abs(LTRE.out$lambda-1))]+mean(poar$Longitude)

dlam.dlong<-(LTRE.out$lambda[2:length(long.seq)]-LTRE.out$lambda[1:(length(long.seq)-1)])/0.5
plot(long.seq[-1],dlam.dlong)

tot.LTRE<-rbind(LTRE.out$F.LTRE.out,LTRE.out$M.LTRE.out)
lines(long.seq,colSums(tot.LTRE))

LTRE.out$M.dlam.dparam

real.long<-long.seq+mean(poar$Longitude)
plot(real.long,LTRE.out$F.LTRE.out[1,],type="l",ylim=c(-1,0))
abline(h=0,col="gray")
lines(real.long,LTRE.out$F.LTRE.out[2,],lty=2)
lines(real.long,LTRE.out$F.LTRE.out[3,],lty=3)
lines(real.long,LTRE.out$F.LTRE.out[4,],lty=4)
lines(real.long,colSums(LTRE.out$F.LTRE.out),lwd=3)
lines(real.long,LTRE.out$M.LTRE.out[1,],col="blue")
lines(real.long,LTRE.out$M.LTRE.out[2,],col="blue",lty=2)
lines(real.long,LTRE.out$M.LTRE.out[3,],col="blue",lty=3)
lines(real.long,LTRE.out$M.LTRE.out[4,],col="blue",lty=4)
lines(real.long,colSums(LTRE.out$M.LTRE.out),col="blue",lwd=3)

par(mfrow=c(2,1))
plot(real.long,LTRE.out$F.LTRE.out[1,],type="l",ylim=c(-.7,0),
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff ,"Long.")))
abline(h=0,col="gray")
lines(real.long,LTRE.out$F.LTRE.out[2,],lty=2)
lines(real.long,LTRE.out$F.LTRE.out[3,],lty=3)
lines(real.long,LTRE.out$F.LTRE.out[4,],lty=4)
lines(real.long,colSums(LTRE.out$F.LTRE.out),lwd=3)

plot(real.long,LTRE.out$M.LTRE.out[1,],col="blue",type="l",ylim=c(-.7,0),
xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff ,"Long.")))
abline(h=0,col="gray")
lines(real.long,LTRE.out$M.LTRE.out[2,],col="blue",lty=2)
lines(real.long,LTRE.out$M.LTRE.out[3,],col="blue",lty=3)
lines(real.long,LTRE.out$M.LTRE.out[4,],col="blue",lty=4)
lines(real.long,colSums(LTRE.out$M.LTRE.out),col="blue",lwd=3)

###########################################################
###########################################################
## Posterior sampling
##main goal is to get posteriors for one-sex and two-sex lambdas
poar.posterior<-poar.out$BUGSoutput$sims.list
long.seq<-seq(min(poar$long.center),max(poar$long.center)*3,0.5)
N.draws<-200
draw<-sample.int(poar.out$BUGSoutput$n.sims,size=N.draws)
lambda.onesex.posterior<-matrix(NA,nrow=N.draws,ncol=length(long.seq))
lambda.twosex.posterior<-matrix(NA,nrow=N.draws,ncol=length(long.seq))

for(i in 1:N.draws){
  F.params<-c()
  ##SURVIVAL#########################
  F.params[1]<-poar.posterior$alpha.surv[draw[i]]  #surv intercpt
  F.params[2]<-poar.posterior$beta.surv.size[draw[i]]  #slope log size
  F.params[3]<-poar.posterior$beta.surv.long[draw[i]]  #slope long
  ##GROWTH#########################
  F.params[11]<-poar.posterior$alpha.grow[draw[i]]  #grow int
  F.params[12]<-poar.posterior$beta.grow.size[draw[i]]  #slope log size
  F.params[13]<-poar.posterior$beta.grow.long[draw[i]]  #slope long
  F.params[14]<-poar.posterior$alpha.gamma.grow[draw[i]]  #NegBin size param
  ##GROWTH#########################
  F.params[21]<-poar.posterior$alpha.flow[draw[i]]  #flow int
  F.params[22]<-poar.posterior$beta.flow.size[draw[i]]  #slope log size
  F.params[23]<-poar.posterior$beta.flow.long[draw[i]]  #slope long
  ##PANIC#########################
  F.params[31]<-poar.posterior$alpha.panic[draw[i]]  #panic int
  F.params[32]<-poar.posterior$beta.panic.size[draw[i]]  #slope log size
  F.params[33]<-poar.posterior$beta.panic.long[draw[i]]  #slope long
  ##GERMINATION#########################
  F.params[41]<-poar.posterior$germ.mu[draw[i]]  #F-dominant seed germination (twosex==F)
  F.params[42]<-poar.posterior$germ.beta1[draw[i]] #This is slope of SR dependence for 'maybe' tetrazolium viability (twosex==T)
  #F.params[41]<-0.4253165 
  #F.params[42]<- 9.6583863
  #F.params[43]<- -12.9141399
  ##MISC#############################
  F.params[51]<-200  #seeds per panicle ## hard coded from Fig 2 of resp surf paper, low density
  F.params[52]<-0.5  #seed sex ratio  # assumed 1:1 re: Jason Goldman
  F.params[61]<-quantile(poar$tillerN_t1,na.rm=T,probs=0.95)  #max size (# tillers), based on real data but not abs max
  
  M.params<-c()
  ##SURVIVAL#########################
  M.params[1]<-poar.posterior$alpha.surv[draw[i]]+poar.posterior$beta.surv.sex[draw[i]]  #surv intercpt
  M.params[2]<-poar.posterior$beta.surv.size[draw[i]]  #slope log size
  M.params[3]<-poar.posterior$beta.surv.long[draw[i]]+poar.posterior$beta.surv.sex.long[draw[i]]  #slope long
  ##GROWTH#########################
  M.params[11]<-poar.posterior$alpha.grow[draw[i]]+poar.posterior$beta.grow.sex[draw[i]]  #grow int
  M.params[12]<-poar.posterior$beta.grow.size[draw[i]]  #slope log size
  M.params[13]<-poar.posterior$beta.grow.long[draw[i]]+poar.posterior$beta.grow.sex.long[draw[i]]  #slope long
  M.params[14]<-poar.posterior$alpha.gamma.grow[draw[i]]  #NegBin size paraM
  ##GROWTH#########################
  M.params[21]<-poar.posterior$alpha.flow[draw[i]]+poar.posterior$beta.flow.sex[draw[i]]  #flow int
  M.params[22]<-poar.posterior$beta.flow.size[draw[i]]  #slope log size
  M.params[23]<-poar.posterior$beta.flow.long[draw[i]]+poar.posterior$beta.flow.sex.long[draw[i]]  #slope long
  ##PANIC#########################
  M.params[31]<-poar.posterior$alpha.panic[draw[i]]+poar.posterior$beta.panic.sex[draw[i]]  #panic int
  M.params[32]<-poar.posterior$beta.panic.size[draw[i]]  #slope log size
  M.params[33]<-poar.posterior$beta.panic.long[draw[i]]+poar.posterior$beta.panic.sex.long[draw[i]]  #slope long
  ##GERMINATION#########################
  M.params[41]<-F.params[41]  #F-dominant seed germination (twosex==F)
  #Hard-coded for now. Value from Fig 1 of response surf ms (max assuming no mate limitation)
  M.params[42]<-F.params[42]  #intercept for sex ratio dependence (twosex==T)
  #M.params[43]<-F.params[43] 
  ##MISC#############################
  M.params[51]<-F.params[51]  #seeds per panicle ## hard coded from Fig 2 of resp surf paper, low density
  M.params[52]<-F.params[52]  #seed sex ratio  # assumed 1:1 re: Jason Goldman
  M.params[61]<-F.params[61]  #max size (# tillers), based on real data but not abs max
  
  for(j in 1:length(long.seq)){
    lambda.onesex.posterior[i,j]<-lambda(megamatrix(F.params=F.params,M.params=M.params,long=long.seq[j],twosex=F)$MEGAmat)
    lambda.twosex.posterior[i,j]<-lambdaSim(F.params=F.params,M.params=M.params,long=long.seq[j],max.yrs=max.yrs)$lambdatracker[max.yrs]
  }
  
}


##grab confidence limits
lambda.onesex.posterior.CI<-lambda.twosex.posterior.CI<-matrix(0,2,length(long.seq))
for(j in 1:length(long.seq)){
  lambda.onesex.posterior.CI[,j]<-quantile(lambda.onesex.posterior[,j],probs=c(0.025,0.975))
  lambda.twosex.posterior.CI[,j]<-quantile(lambda.twosex.posterior[,j],probs=c(0.025,0.975))
}

par(mfrow=c(1,1))
plot(long.seq+mean(poar$Longitude),
     colMeans(lambda.onesex.posterior),type="l",ylim=c(0,10),
     xlab="Longitude",ylab=expression(lambda),cex.lab=1.5)
abline(1,0,lty=3,col="black")
polygon(x=c(long.seq+mean(poar$Longitude),rev(long.seq+mean(poar$Longitude))),
        y=c(lambda.onesex.posterior.CI[1,],rev(lambda.onesex.posterior.CI[2,])),
        col=adjustcolor("lightgray",alpha.f=0.8),border=NA,lwd=2
)
lines(long.seq+mean(poar$Longitude),colMeans(lambda.onesex.posterior),lwd=2,col="black")
polygon(x=c(long.seq+mean(poar$Longitude),rev(long.seq+mean(poar$Longitude))),
        y=c(lambda.twosex.posterior.CI[1,],rev(lambda.twosex.posterior.CI[2,])),
        col=adjustcolor("tomato",alpha.f=0.3),border=NA,lwd=2
)
lines(long.seq+mean(poar$Longitude),colMeans(lambda.twosex.posterior),lwd=2,col="red")

abline(v=max(poar$long.center)+mean(poar$Longitude))

abline(v=long.seq[which.min(abs(colMeans(lambda.onesex.posterior)-1))]+mean(poar$Longitude))
abline(v=long.seq[which.min(abs(colMeans(lambda.twosex.posterior)-1))]+mean(poar$Longitude),lty=2)

##NICE FIGURE
plot(long.seq+mean(poar$Longitude),colMeans(lambda.twosex.posterior),type="n",ylim=c(0,8),
     xlab="Longitude",ylab=expression(lambda),cex.lab=1.5)
polygon(x=c(long.seq+mean(poar$Longitude),rev(long.seq+mean(poar$Longitude))),
        y=c(lambda.twosex.posterior.CI[1,],rev(lambda.twosex.posterior.CI[2,])),
        col=adjustcolor("tomato",alpha.f=0.3),border=NA,lwd=2
)
lines(long.seq+mean(poar$Longitude),colMeans(lambda.twosex.posterior),lwd=3,col="red")
abline(1,0,lty=3,col="black")
abline(v=max(poar$long.center)+mean(poar$Longitude))


par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(real.long,LTRE.out$F.LTRE.out[1,],type="l",ylim=c(-.7,0),col="darkblue",cex.lab=1.4,
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff ,"Long.")))
abline(h=0,col="gray")
lines(real.long,LTRE.out$F.LTRE.out[2,],lty=2,col="darkblue")
lines(real.long,LTRE.out$F.LTRE.out[3,],lty=3,col="darkblue")
lines(real.long,LTRE.out$F.LTRE.out[4,],lty=4,col="darkblue")
lines(real.long,colSums(LTRE.out$F.LTRE.out),lwd=4,col="darkblue")
legend("bottomright",legend=c("Survival","Growth","Flowering","Fertility","Total"),
       lty=c(1,2,3,4,1),lwd=c(1,1,1,1,4),col="darkblue",bty="n",cex=.8)

plot(real.long,LTRE.out$M.LTRE.out[1,],col="darkorange3",type="l",ylim=c(-.7,0),cex.lab=1.4,
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff ,"Long.")))
abline(h=0,col="gray")
lines(real.long,LTRE.out$M.LTRE.out[2,],col="darkorange3",lty=2)
lines(real.long,LTRE.out$M.LTRE.out[3,],col="darkorange3",lty=3)
lines(real.long,LTRE.out$M.LTRE.out[4,],col="darkorange3",lty=4)
lines(real.long,colSums(LTRE.out$M.LTRE.out),col="darkorange3",lwd=4)
legend("bottomright",legend=c("Survival","Growth","Flowering","Fertility","Total"),
       lty=c(1,2,3,4,1),lwd=c(1,1,1,1,4),col="darkorange3",bty="n",cex=.8)