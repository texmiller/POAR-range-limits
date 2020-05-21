
# vital rate and megamatrix functions ---------------------------------------------------

#SURVIVAL AT SIZE X.
sx<-function(x,params,long,rfx,surv_perturb=0){
  surv_mean<-params$surv_mu + 
    params$surv_size*log(x) + 
    params$surv_long*long + 
    params$surv_size_long*log(x)*long +
    params$surv_long2*(long^2) + 
    params$surv_size_long2*log(x)*(long^2) + 
    rfx["surv","site"] + rfx["surv","block"] + rfx["surv","source"]
  return(invlogit(surv_mean) + surv_perturb)
}

#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,long,rfx,grow_perturb=0){
  grow_mean<-exp(params$grow_mu + 
    params$grow_size*log(x) + 
    params$grow_long*long + 
    params$grow_size_long*log(x)*long +
    params$grow_long2*(long^2) + 
    params$grow_size_long2*log(x)*(long^2) + 
    rfx["grow","site"] + rfx["grow","block"] + rfx["grow","source"]) + grow_perturb
  
  grow<-dpoisinvgauss(x=y,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  truncLower<-dpoisinvgauss(x=0,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-dpoisinvgauss(x=params$max_size:10000,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
  return(grow/(1-(truncLower+truncUpper)))
}

#SURVIVAL*GROWTH
pxy<-function(x,y,params,long,rfx,grow_perturb=0,surv_perturb=0){
  sx(x,params,long,rfx,surv_perturb) * gxy(x,y,params,long,rfx,grow_perturb)
}

# PROBABILITY OF FLOWERING
pfx<-function(x,params,long,rfx,flow_perturb=0){
  flow_mean<-params$flow_mu + 
    params$flow_size*log(x) + 
    params$flow_long*long + 
    params$flow_size_long*log(x)*long +
    params$flow_long2*(long^2) + 
    params$flow_size_long2*log(x)*(long^2) + 
    rfx["flow","site"] + rfx["flow","block"] + rfx["flow","source"]
  return(invlogit(flow_mean) + flow_perturb)
}

#NUMBER OF PANICLES
nfx<-function(x,params,long,rfx,fert_perturb=0){
  panic_mean<-params$panic_mu + 
    params$panic_size*log(x) + 
    params$panic_long*long + 
    params$panic_size_long*log(x)*long +
    params$panic_long2*(long^2) + 
    params$panic_size_long2*log(x)*(long^2) + 
    rfx["panic","site"] + rfx["panic","block"] + rfx["panic","source"]
  return(exp(panic_mean) + fert_perturb)
}

#SEED VIABILITY
viab<-function(params,twosex,OSR=NULL,viab_perturb=0){
  if(twosex==F){return(params$v0 + viab_perturb)}
  if(twosex==T){return((params$v0 * (1 - OSR ^ params$a_v)) + viab_perturb)}
}

#FERTILITY--returns number of recruits
##Female offspring
fertx_F<-function(x,params,rfx,long,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,long,rfx,flow_perturb)*nfx(x,params,long,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*params$PSR
  return(seedlings)
}

##Male offspring
fertx_M<-function(x,params,rfx,long,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,long,rfx,flow_perturb)*nfx(x,params,long,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*(1-params$PSR)
  return(seedlings)
}

megamatrix_delay<-function(F_params,M_params,long,twosex,OSR=NULL,rfx,
                           grow_perturb=0,surv_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){  
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim+1,matdim+1)
  F.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=F_params,long=long,rfx=rfx,grow_perturb=grow_perturb,surv_perturb=surv_perturb))
  # surviving seedlings emerge in continuous population
  F.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=F_params,long=long,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim+1,matdim+1)
  # seedlings in top row
  F.Fmat[1,2:(matdim+1)]<-fertx_F(x=y,params=F_params,rfx=rfx,long=long,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim+1,matdim+1)
  M.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=M_params,long=long,rfx=rfx,grow_perturb=grow_perturb,surv_perturb=surv_perturb))
  M.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=M_params,long=long,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim+1,matdim+1)
  M.Fmat[1,2:(matdim+1)]<-fertx_M(x=y,params=F_params,rfx=rfx,long=long,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
  #M-to-F
  zero.mat<-matrix(0,matdim+1,matdim+1)
  
  # Put it all together as a megamatrix
  MEGAmat<-cbind(rbind(F.Tmat+F.Fmat,  ##Female growth/survival + recruitment[1,1]
                       M.Fmat), ##Male recruitment [2,1]
                 rbind(zero.mat,   ##Females from males [1,2]
                       M.Tmat))   ##Male growth/survival
  
  return(list(MEGAmat=MEGAmat,F.Tmat=F.Tmat,M.Tmat=M.Tmat,y=y))
}

# Analysis of 2sex model --------------------------------------------------
# this needs to be done by simulation

lambdaSim_delay<-function(F_params,M_params,long,rfx,max.yrs,
                          grow_perturb=0,surv_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/((matdim+1)*2),((matdim+1)*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[2:(matdim+1)]*pfx(x=y,param=F_params,long=long,rfx=rfx,flow_perturb=flow_perturb) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,long=long,rfx=rfx,fert_perturb=fert_perturb) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+3):((matdim+1)*2)]*pfx(x=y,param=M_params,long=long,rfx=rfx,flow_perturb=flow_perturb)
    M_panicles<-flowering_males%*%nfx(x=y,param=M_params,long=long,rfx=rfx,fert_perturb=fert_perturb)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:(matdim+1)])
    #assmble matrix
    MEGAmat<-megamatrix_delay(F_params=F_params,M_params=M_params,long=long,twosex=T,OSR=OSRtracker[t],rfx=rfx,
                              grow_perturb=grow_perturb,surv_perturb=surv_perturb,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb)$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker,n0=n0))
}


# RFX fun -----------------------------------------------------------------
rfx_fun <- function(site_tau_s=0,block_tau_s=0,source_tau_s=0,
                    site_tau_g=0,block_tau_g=0,source_tau_g=0,
                    site_tau_f=0,block_tau_f=0,source_tau_f=0,
                    site_tau_p=0,block_tau_p=0,source_tau_p=0){
  rfx <- data.frame(site = rnorm(4,0,c(site_tau_s,site_tau_g,site_tau_f,site_tau_p)),
                    block = rnorm(4,0,c(block_tau_s,block_tau_g,block_tau_f,block_tau_p)),
                    source = rnorm(4,0,c(source_tau_s,source_tau_g,source_tau_f,source_tau_p)))
  rownames(rfx) <- c("surv","grow","flow","panic")
  return(rfx)
}

# basement garbage --------------------------------------------------------


## assmble vital rate functions into a projection matrix
## This was the original megamatrix function, was not updated to include vital rate perturbations (because we are not using it)
megamatrix<-function(F_params,M_params,long,twosex,OSR=NULL,rfx){  
  matdim<-F_params$max_size         
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim,matdim)
  F.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=F_params,long=long,rfx=rfx,grow_perturb=grow_perturb,surv_perturb=surv_perturb))
  
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
  
  return(list(MEGAmat=MEGAmat,F.Tmat=F.Tmat,M.Tmat=M.Tmat,y=y))
}


lambdaSim<-function(F_params,M_params,long,rfx,max.yrs){
  matdim<-F_params$max_size         
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/(matdim*2),(matdim*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[1:matdim]*pfx(x=y,param=F_params,long=long,rfx=rfx) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,long=long,rfx=rfx) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+1):(matdim*2)]*pfx(x=y,param=M_params,long=long,rfx=rfx)
    M_panicles<-flowering_males%*%nfx(x=y,param=M_params,long=long,rfx=rfx)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:matdim])
    #assmble matrix
    MEGAmat<-megamatrix(F_params=F_params,M_params=M_params,long=long,twosex=T,OSR=OSRtracker[t],rfx=rfx)$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker,n0=n0))
}


