library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
library(countreg)
## useful functions
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
invlogit<-function(x){exp(x)/(1+exp(x))}

dir <- "C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits"
dir <- "C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits"


# First, check posteriors of all parameters in full model -------------------------------------------------------------------------


## Load the Stan output for full vital rate model
#fit_dropsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_full.rds"))
fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds"))

##survival coefs
mcmc_dens_overlay(fit_allsites_full,par=quote_bare(b0_s,bsize_s,bsex_s,blong_s,bsizesex_s,bsizelong_s,blongsex_s,bsizelongsex_s,
                                                   blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s))
#growth coefs
mcmc_dens_overlay(fit_allsites_full,par=quote_bare(b0_g,bsize_g,bsex_g,blong_g,bsizesex_g,bsizelong_g,blongsex_g,bsizelongsex_g,
                                                   blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,phi_g))
#flowering coefs
mcmc_dens_overlay(fit_allsites_full,par=quote_bare(b0_f,bsize_f,bsex_f,blong_f,bsizesex_f,bsizelong_f,blongsex_f,bsizelongsex_f,
                                                   blong2_f,bsizelong2_f,blong2sex_f,bsizelong2sex_f))
#panicle coefs
mcmc_dens_overlay(fit_allsites_full,par=quote_bare(b0_p,bsize_p,bsex_p,blong_p,bsizesex_p,bsizelong_p,blongsex_p,bsizelongsex_p,
                                                   blong2_p,bsizelong2_p,blong2sex_p,bsizelong2sex_p,phi_p))
## these generally look good. I don' think there is anything wrong with the fit.


# Second, check the lambda-long relationship at mean params -------------------------------------------------------------------------

# pull out stan coefficients
mean_coef <- lapply(rstan::extract(fit_full, pars = quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                                        bsizesex_s, bsizelong_s,blongsex_s,bsizelongsex_s,
                                                        blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s,
                                                        site_tau_s,block_tau_s,source_tau_s,
                                                        
                                                        b0_g,bsize_g,bsex_g,blong_g,
                                                        bsizesex_g, bsizelong_g,blongsex_g,bsizelongsex_g,
                                                        blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,
                                                        site_tau_g,block_tau_g,source_tau_g,
                                                        phi_g,
                                                        
                                                        b0_f,bsize_f,bsex_f,blong_f,
                                                        bsizesex_f, bsizelong_f,blongsex_f,bsizelongsex_f,
                                                        blong2_f,bsizelong2_f,blong2sex_f,bsizelong2sex_f,
                                                        site_tau_f,block_tau_f,source_tau_f,
                                                        
                                                        b0_p,bsize_p,bsex_p,blong_p,
                                                        bsizesex_p, bsizelong_p,blongsex_p,bsizelongsex_p,
                                                        blong2_p,bsizelong2_p,blong2sex_p,bsizelong2sex_p,
                                                        site_tau_p,block_tau_p,source_tau_p,
                                                        
                                                        v0,a_v,m,lambda_d))
                    ,mean)
                    #,median)

# set up param vectors-------------------------------------------------------------------------
F_params <- M_params <- list()
## survival
F_params$surv_mu <- mean_coef$b0_s
F_params$surv_size <- mean_coef$bsize_s
F_params$surv_long <- mean_coef$blong_s 
F_params$surv_size_long <- mean_coef$bsizelong_s 
F_params$surv_long2 <- mean_coef$blong2_s 
F_params$surv_size_long2 <- mean_coef$bsizelong2_s 
M_params$surv_mu <- mean_coef$b0_s + mean_coef$bsex_s  
M_params$surv_size <- mean_coef$bsize_s + mean_coef$bsizesex_s 
M_params$surv_long <- mean_coef$blong_s + mean_coef$blongsex_s 
M_params$surv_size_long <- mean_coef$bsizelong_s + mean_coef$bsizelongsex_s 
M_params$surv_long2 <- mean_coef$blong2_s + mean_coef$blong2sex_s 
M_params$surv_size_long2 <- mean_coef$bsizelong2_s + mean_coef$bsizelong2sex_s 
## growth
F_params$grow_mu <- mean_coef$b0_g 
F_params$grow_size <- mean_coef$bsize_g 
F_params$grow_long <- mean_coef$blong_g 
F_params$grow_size_long <- mean_coef$bsizelong_g 
F_params$grow_long2 <- mean_coef$blong2_g 
F_params$grow_size_long2 <- mean_coef$bsizelong2_g 
F_params$phi_g <- mean_coef$phi_g 
M_params$grow_mu <- mean_coef$b0_g + mean_coef$bsex_g 
M_params$grow_size <- mean_coef$bsize_g + mean_coef$bsizesex_g 
M_params$grow_long <- mean_coef$blong_g + mean_coef$blongsex_g 
M_params$grow_size_long <- mean_coef$bsizelong_g + mean_coef$bsizelongsex_g 
M_params$grow_long2 <- mean_coef$blong2_g + mean_coef$blong2sex_g 
M_params$grow_size_long2 <- mean_coef$bsizelong2_g + mean_coef$bsizelong2sex_g 
M_params$phi_g <- mean_coef$phi_g 
## flowering
F_params$flow_mu <- mean_coef$b0_f 
F_params$flow_size <- mean_coef$bsize_f 
F_params$flow_long <- mean_coef$blong_f 
F_params$flow_size_long <- mean_coef$bsizelong_f 
F_params$flow_long2 <- mean_coef$blong2_f 
F_params$flow_size_long2 <- mean_coef$bsizelong2_f 
M_params$flow_mu <- mean_coef$b0_f + mean_coef$bsex_f 
M_params$flow_size <- mean_coef$bsize_f + mean_coef$bsizesex_f 
M_params$flow_long <- mean_coef$blong_f + mean_coef$blongsex_f 
M_params$flow_size_long <- mean_coef$bsizelong_f + mean_coef$bsizelongsex_f 
M_params$flow_long2 <- mean_coef$blong2_f + mean_coef$blong2sex_f 
M_params$flow_size_long2 <- mean_coef$bsizelong2_f + mean_coef$bsizelong2sex_f 
## panicles
F_params$panic_mu <- mean_coef$b0_p 
F_params$panic_size <- mean_coef$bsize_p 
F_params$panic_long <- mean_coef$blong_p 
F_params$panic_size_long <- mean_coef$bsizelong_p 
F_params$panic_long2 <- mean_coef$blong2_p 
F_params$panic_size_long2 <- mean_coef$bsizelong2_p 
M_params$panic_mu <- mean_coef$b0_p + mean_coef$bsex_p 
M_params$panic_size <- mean_coef$bsize_p + mean_coef$bsizesex_p 
M_params$panic_long <- mean_coef$blong_p + mean_coef$blongsex_p 
M_params$panic_size_long <- mean_coef$bsizelong_p + mean_coef$bsizelongsex_p 
M_params$panic_long2 <- mean_coef$blong2_p + mean_coef$blong2sex_p 
M_params$panic_size_long2 <- mean_coef$bsizelong2_p + mean_coef$bsizelong2sex_p 
## seed viability and misc fertility params
F_params$v0 <- mean_coef$v0 
F_params$a_v <- mean_coef$a_v 
F_params$ov_per_inf <- mean_coef$lambda_d 
F_params$germ <- mean_coef$m 
F_params$PSR <- 0.5
## set max size equal between the sexes
F_params$max_size <- quantile(na.omit(poar$tillerN_t0),probs=0.99) #max(na.omit(poar$tillerN_t0)); 
M_params$max_size <- F_params$max_size

long_seq <- seq(min(poar_surv_binned$long),max(poar_surv_binned$long),0.2)
lambda_long_mean<-SR_long_mean<-OSR_long_mean<-c()
for(l in 1:length(long_seq)){
  lambda_run <- lambdaSim(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)
  lambda_long_mean[l] <- lambda_run$lambdatracker[max_yrs]
  SR_long_mean[l] <- lambda_run$SRtracker[max_yrs]
  OSR_long_mean[l] <- lambda_run$OSRtracker[max_yrs]
}

plot(long_seq,lambda_long_mean,type="b")

## the problem is panicle production. See where the fitted function flips in concavity
sizes <- 1:F_params$max_size

plot(long_seq,nfx(x=F_params$max_size,params=F_params,long=long_seq,rfx=rfx_fun()),type="n",ylim=c(0,5))
for(i in 1:length(sizes)){
  lines(long_seq,nfx(x=sizes[i],params=F_params,long=long_seq,rfx=rfx_fun()))
}

par(mfrow=c(1,3))
with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    plot(long[size_bin==i],mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()   

      points(long[sex==1 & size_bin==i],mean_grow[sex==1 & size_bin==i],
             bg="black",pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      points(long[sex==2 & size_bin==i],mean_grow[sex==2 & size_bin==i],
             bg="white",pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)


  }
})
