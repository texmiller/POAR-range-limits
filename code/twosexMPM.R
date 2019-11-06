library(tidyverse)

# vital rate and megamatrix functions ---------------------------------------------------

#SURVIVAL AT SIZE X.
sx<-function(x,params,long,rfx){
  surv_mean<-params$surv_mu + 
    params$surv_size*log(x) + 
    params$surv_long*long + 
    params$surv_size_long*log(x)*long +
    params$surv_long2*(long^2) + 
    params$surv_size_long2*log(x)*(long^2) + 
    rfx$site["surv"] + rfx$block["surv"] + rfx$source["surv"]
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
    rfx$site["grow"] + rfx$block["grow"] + rfx$source["grow"]
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
    rfx$site["flow"] + rfx$block["flow"] + rfx$source["flow"]
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
    rfx$site["panic"] + rfx$block["panic"] + rfx$source["panic"]
  return(exp(panic_mean))
}

#SEED VIABILITY
viab<-function(params,twosex,OSR){
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

## Put it all together
megamatrix<-function(F.params,M.params,long,twosex,SR=NULL,proc_err=F,rand_seed=NULL){  
  
  ## draw random effects for site, block, source
  rfx <- tibble(site = rep(0,4),
                block = rep(0,4),
                source = rep(0,4))
  rownames(rfx) <- c("surv","grow","flow","panic")
  if(proc_err==T){
    set.seed(rand_seed)
    rfx["surv",] <- rnorm(3,0,c(F.params$site_sd_s,F.params$block_sd_s,F.params$source_sd_s))
    rfx["grow",] <- rnorm(3,0,c(F.params$site_sd_g,F.params$block_sd_g,F.params$source_sd_g))
    rfx["flow",] <- rnorm(3,0,c(F.params$site_sd_f,F.params$block_sd_f,F.params$source_sd_f))
    rfx["panic",] <- rnorm(3,0,c(F.params$site_sd_p,F.params$block_sd_p,F.params$source_sd_p))
  }

  matdim<-F.params$max_size         
  y<-1:F.params$max_size
  
  ## stopped here 11/6/2019
  
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
# demographic parameters --------------------------------------------------
dir <- "C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits"
dir <- "C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits"
demog_par <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds"))

