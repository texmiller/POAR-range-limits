
# vital rate and megamatrix functions ---------------------------------------------------

#SURVIVAL AT SIZE X.
sx<-function(x,params,long,rfx){
  surv_mean<-params$surv_mu + 
    params$surv_size*log(x) + 
    params$surv_long*long + 
    params$surv_size_long*log(x)*long +
    params$surv_long2*(long^2) + 
    params$surv_size_long2*log(x)*(long^2) + 
    rfx["surv","site"] + rfx["surv","block"] + rfx["surv","source"]
  return(invlogit(surv_mean))
}

#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,long,rfx){
  grow_mean<-params$grow_mu + 
    params$grow_size*log(x) + 
    params$grow_long*long + 
    params$grow_size_long*log(x)*long +
    params$grow_long2*(long^2) + 
    params$grow_size_long2*log(x)*(long^2) + 
    rfx["grow","site"] + rfx["grow","block"] + rfx["grow","source"]
  grow<-dnbinom(x=y,mu=exp(grow_mean),size=params$phi_g,log=F)
  truncLower<-dnbinom(x=0,mu=exp(grow_mean),size=params$phi_g,log=F)
  truncUpper<-sum(dnbinom(x=params$max_size:10000,mu=exp(grow_mean),size=params$phi_g,log=F))
  return(grow/(1-(truncLower+truncUpper)))
}

#SURVIVAL*GROWTH
pxy<-function(x,y,params,long,rfx){
  sx(x,params,long,rfx) * gxy(x,y,params,long,rfx)
}

# PROBABILITY OF FLOWERING
pfx<-function(x,params,long,rfx){
  flow_mean<-params$flow_mu + 
    params$flow_size*log(x) + 
    params$flow_long*long + 
    params$flow_size_long*log(x)*long +
    params$flow_long2*(long^2) + 
    params$flow_size_long2*log(x)*(long^2) + 
    rfx["flow","site"] + rfx["flow","block"] + rfx["flow","source"]
  return(invlogit(flow_mean))
}

#NUMBER OF PANICLES
nfx<-function(x,params,long,rfx){
  panic_mean<-params$panic_mu + 
    params$panic_size*log(x) + 
    params$panic_long*long + 
    params$panic_size_long*log(x)*long +
    params$panic_long2*(long^2) + 
    params$panic_size_long2*log(x)*(long^2) + 
    rfx["panic","site"] + rfx["panic","block"] + rfx["panic","source"]
  return(exp(panic_mean))
}

#SEED VIABILITY
viab<-function(params,twosex,OSR=NULL){
  if(twosex==F){return(params$v0)}
  if(twosex==T){return(params$v0 * (1 - OSR ^ params$a_v))}
}

#FERTILITY--returns number of recruits
##Female offspring
fertx_F<-function(x,params,rfx,long,twosex,OSR=NULL){
  seedlings<-pfx(x,params,long,rfx)*nfx(x,params,long,rfx)*params$ov_per_inf*viab(params,twosex,OSR)*params$germ*params$PSR
  return(seedlings)
}

##Male offspring
fertx_M<-function(x,params,rfx,long,twosex,OSR=NULL){
  seedlings<-pfx(x,params,long,rfx)*nfx(x,params,long,rfx)*params$ov_per_inf*viab(params,twosex,OSR)*params$germ*(1-params$PSR)
  return(seedlings)
}

## assmble vital rate functions into a projection matrix
megamatrix<-function(F_params,M_params,long,twosex,OSR=NULL,rfx){  
  matdim<-F_params$max_size         
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim,matdim)
  F.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=F_params,long=long,rfx=rfx))
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim,matdim)
  F.Fmat[1,1:matdim]<-fertx_F(x=y,params=F_params,rfx=rfx,long=long,twosex=twosex,OSR=OSR) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim,matdim)
  M.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=M_params,long=long,rfx=rfx))
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim,matdim)
  M.Fmat[1,1:matdim]<-fertx_M(x=y,params=F_params,rfx=rfx,long=long,twosex=twosex,OSR=OSR) 
  
  #M-to-F
  zero.mat<-matrix(0,matdim,matdim)
  
  # Put it all together as a megamatrix
  MEGAmat<-cbind(rbind(F.Tmat+F.Fmat,  ##Female growth/survival + recruitment[1,1]
                       M.Fmat), ##Male recruitment [2,1]
                 rbind(zero.mat,   ##Females from males [1,2]
                       M.Tmat))   ##Male growth/survival
  
  return(list(MEGAmat=MEGAmat,y=y))
}

# Analysis of 2sex model --------------------------------------------------
# this needs to be done by simulation

lambdaSim<-function(F_params,M_params,long,rfx,max.yrs){
  matdim<-F.params$max_size         
  y<-1:F.params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/(matdim*2),(matdim*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[1:matdim]*pfx(x=y,param=F_params,long=long,rfx=rfx) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,long=long,rfx=rfx,rfx=rfx) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+1):(matdim*2)]*pfx(x=y,param=M_params,long=long)
    M_panicles<-flowering_males%*%nfx(x=y,param=M_params,long=long,rfx=rfx)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:matdim])
    #assmble matrix
    MEGAmat<-megamatrix(F.params=F_params,M_params=M.params,long=long,twosex=T,OSR=OSRtracker[t],rfx=rfx)$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker))
}


# LTRE function -----------------------------------------------------------
## store LTRE output in a matrix of columns for longitude and rows for parameters
POAR_LTRE<-function(F_params,M_params,long,max.yrs,perturbation,
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

