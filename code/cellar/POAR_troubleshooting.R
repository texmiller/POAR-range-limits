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
#fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds"))
fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full_noLong2intx.rds"))

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
fit_full <- fit_allsites_full  #fit_full <- fit_dropsites_full #

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


# Re-fit panicles alone ---------------------------------------------------


# data for model
data_panicles <- list( n_sites    = poar_allsites.panic$site %>% n_distinct,
                           n_sources  = poar_allsites.panic$source %>% n_distinct(),
                           
                          # panicle data
                           n_blocks_p = poar_allsites.panic$block %>% n_distinct,
                           site_p   = poar_allsites.panic$site,
                           source_p =  poar_allsites.panic$source,
                           block_p  = poar_allsites.panic$block,
                           site_block_p = data.frame( site_i  = poar_allsites.panic$site,
                                                      block_i = poar_allsites.panic$block ) %>% 
                             unique %>% .$site_i,
                           male_p   = poar_allsites.panic$sex-1,
                           long_p   = poar_allsites.panic$long.center,
                           size_p   = poar_allsites.panic$log_size_t1,
                           y_p      = poar_allsites.panic$panic_t1,
                           n_p      = nrow(poar_allsites.panic)
                          )
# simulation parameters
sim_pars <- list(
  warmup = 2000, 
  iter = 10000, 
  thin = 3, 
  chains = 3
)
# fit the "big" model 
fit_panicles_noZT <- stan(
  file = 'code/stan/panicle_noZT.stan',
  data = data_panicles,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

predP <- rstan::extract(fit_panicles_noZT, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_panicles_noZT, pars = c("phi_p"))$phi_p
n_post_draws <- 100
post_draws <- sample.int(dim(predP)[1], n_post_draws)

y_p_sim1 <- y_p_sim2 <- matrix(NA,n_post_draws,length(data_panicles$y_p))

for(i in 1:n_post_draws){
  mu_i <- exp(predP[i,])
  ## sample panicle data (zero-truncated NB)
  ##y_p_sim[i,] <- rztnbinom(n=length(data_allsites_all$y_p), mu = exp(predP[i,]), size=phi_P[i])
  y_p_sim1[i,] <- rnbinom(n=length(data_panicles$y_p), mu = mu_i, size=phi_P[i])
  ## this is adjusting the mean according to: https://data.princeton.edu/wws509/notes/countmoments
  y_p_sim2[i,] <- rnbinom(n=length(data_panicles$y_p), mu = mu_i / (1 - (1+mu_i*phi_P[i])^(-1/phi_P[i])), size=phi_P[i])
}

ppc_dens_overlay(data_panicles$y_p, y_p_sim1)+xlim(0, 50)
ppc_dens_overlay(data_panicles$y_p, y_p_sim2)+xlim(0, 50)

mean_panicle_coef <- lapply(rstan::extract(fit_panicles_noZT, pars = quote_bare(b0_p,bsize_p,bsex_p,blong_p,
                                                               bsizesex_p, bsizelong_p,blongsex_p,bsizelongsex_p,
                                                               blong2_p,bsizelong2_p,blong2sex_p,bsizelong2sex_p))
                    ,mean)
test_params<-c()
test_params$panic_mu <- mean_panicle_coef$b0_p 
test_params$panic_size <- mean_panicle_coef$bsize_p 
test_params$panic_long <- mean_panicle_coef$blong_p 
test_params$panic_size_long <- mean_panicle_coef$bsizelong_p 
test_params$panic_long2 <- mean_panicle_coef$blong2_p 
test_params$panic_size_long2 <- mean_panicle_coef$bsizelong2_p 

plot(long_seq,nfx(x=F_params$max_size,params=test_params,long=long_seq,rfx=rfx_fun()),type="n",ylim=c(0,5))
for(i in 1:length(sizes)){
  lines(long_seq,nfx(x=sizes[i],params=test_params,long=long_seq,rfx=rfx_fun()))
}


# No long2 interactions ---------------------------------------------------

fit_panicles_noLong2intx <- stan(
  file = 'code/stan/panicle_noLong2interactions.stan',
  data = data_panicles,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

predP <- rstan::extract(fit_panicles_noLong2intx, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_panicles_noLong2intx, pars = c("phi_p"))$phi_p
n_post_draws <- 100
post_draws <- sample.int(dim(predP)[1], n_post_draws)

y_p_sim <- matrix(NA,n_post_draws,length(data_panicles$y_p))

for(i in 1:n_post_draws){
  mu_i <- exp(predP[i,])
  ## sample panicle data (zero-truncated NB)
  ##y_p_sim[i,] <- rztnbinom(n=length(data_allsites_all$y_p), mu = exp(predP[i,]), size=phi_P[i])
  y_p_sim[i,] <- rnbinom(n=length(data_panicles$y_p), mu = mu_i, size=phi_P[i])
}

ppc_dens_overlay(data_panicles$y_p, y_p_sim)+xlim(0, 50)

mean_panicle_coef <- lapply(rstan::extract(fit_panicles_noLong2intx, pars = quote_bare(b0_p,bsize_p,bsex_p,blong_p,
                                                                                bsizesex_p, bsizelong_p,blongsex_p,bsizelongsex_p,
                                                                                blong2_p))
                            ,mean)
test_params<-c()
test_params$panic_mu <- mean_panicle_coef$b0_p 
test_params$panic_size <- mean_panicle_coef$bsize_p 
test_params$panic_long <- mean_panicle_coef$blong_p 
test_params$panic_size_long <- mean_panicle_coef$bsizelong_p 
test_params$panic_long2 <- mean_panicle_coef$blong2_p 
test_params$panic_size_long2 <- 0 

plot(long_seq,nfx(x=F_params$max_size,params=test_params,long=long_seq,rfx=rfx_fun()),type="n",ylim=c(0,5))
for(i in 1:length(sizes)){
  lines(long_seq,nfx(x=sizes[i],params=test_params,long=long_seq,rfx=rfx_fun()))
}

## see what happens when long2 interactions are dropped from full model

## Load the Stan output for full vital rate model
#fit_dropsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_full.rds"))
fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full_noLong2intx.rds"))


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
F_params$surv_size_long2 <- mean_coef$bsizelong2_s  #0#
M_params$surv_mu <- mean_coef$b0_s + mean_coef$bsex_s  
M_params$surv_size <- mean_coef$bsize_s + mean_coef$bsizesex_s 
M_params$surv_long <- mean_coef$blong_s + mean_coef$blongsex_s 
M_params$surv_size_long <- mean_coef$bsizelong_s + mean_coef$bsizelongsex_s 
M_params$surv_long2 <- mean_coef$blong2_s + mean_coef$blong2sex_s #0#
M_params$surv_size_long2 <- mean_coef$bsizelong2_s + mean_coef$bsizelong2sex_s  #0#
## growth
F_params$grow_mu <- mean_coef$b0_g #+4
F_params$grow_size <- mean_coef$bsize_g 
F_params$grow_long <- mean_coef$blong_g 
F_params$grow_size_long <- mean_coef$bsizelong_g 
F_params$grow_long2 <- mean_coef$blong2_g 
F_params$grow_size_long2 <- mean_coef$bsizelong2_g #0#
F_params$phi_g <- mean_coef$phi_g 
M_params$grow_mu <- mean_coef$b0_g + mean_coef$bsex_g #+4
M_params$grow_size <- mean_coef$bsize_g + mean_coef$bsizesex_g 
M_params$grow_long <- mean_coef$blong_g + mean_coef$blongsex_g 
M_params$grow_size_long <- mean_coef$bsizelong_g + mean_coef$bsizelongsex_g 
M_params$grow_long2 <- mean_coef$blong2_g + mean_coef$blong2sex_g #0#
M_params$grow_size_long2 <- mean_coef$bsizelong2_g + mean_coef$bsizelong2sex_g #0#
M_params$phi_g <- mean_coef$phi_g 
## flowering
F_params$flow_mu <- mean_coef$b0_f 
F_params$flow_size <- mean_coef$bsize_f 
F_params$flow_long <- mean_coef$blong_f 
F_params$flow_size_long <- mean_coef$bsizelong_f 
F_params$flow_long2 <- mean_coef$blong2_f 
F_params$flow_size_long2 <- mean_coef$bsizelong2_f  #0#
M_params$flow_mu <- mean_coef$b0_f + mean_coef$bsex_f 
M_params$flow_size <- mean_coef$bsize_f + mean_coef$bsizesex_f 
M_params$flow_long <- mean_coef$blong_f + mean_coef$blongsex_f 
M_params$flow_size_long <- mean_coef$bsizelong_f + mean_coef$bsizelongsex_f 
M_params$flow_long2 <- mean_coef$blong2_f + mean_coef$blong2sex_f #0#
M_params$flow_size_long2 <- mean_coef$bsizelong2_f + mean_coef$bsizelong2sex_f #0#
## panicles
F_params$panic_mu <- mean_coef$b0_p 
F_params$panic_size <- mean_coef$bsize_p 
F_params$panic_long <- mean_coef$blong_p 
F_params$panic_size_long <- mean_coef$bsizelong_p 
F_params$panic_long2 <- mean_coef$blong2_p 
F_params$panic_size_long2 <- mean_coef$bsizelong2_p #0#
M_params$panic_mu <- mean_coef$b0_p + mean_coef$bsex_p 
M_params$panic_size <- mean_coef$bsize_p + mean_coef$bsizesex_p 
M_params$panic_long <- mean_coef$blong_p + mean_coef$blongsex_p 
M_params$panic_size_long <- mean_coef$bsizelong_p + mean_coef$bsizelongsex_p 
M_params$panic_long2 <- mean_coef$blong2_p + mean_coef$blong2sex_p #0#
M_params$panic_size_long2 <- mean_coef$bsizelong2_p + mean_coef$bsizelong2sex_p #0#
## seed viability and misc fertility params
F_params$v0 <- mean_coef$v0 
F_params$a_v <- mean_coef$a_v 
F_params$ov_per_inf <- mean_coef$lambda_d 
F_params$germ <- mean_coef$m 
F_params$PSR <- 0.5
## set max size equal between the sexes
F_params$max_size <- quantile(na.omit(poar$tillerN_t1),probs=0.95) #max(na.omit(poar$tillerN_t0)); 
M_params$max_size <- F_params$max_size
## survival of 1yo's
F_params$surv_1yo <- M_params$surv_1yo <- 0.1 ## can I get a number from Poa autumnalis?

long_seq <- seq(min(poar_surv_binned$long),max(poar_surv_binned$long),0.3)
lambda_long_mean<-SR_long_mean<-OSR_long_mean<-c()
n0<-matrix(NA,(F_params$max_size+1)*2,length(long_seq))
for(l in 1:length(long_seq)){
  #lambda_run <- lambdaSim_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)
  #lambda_long_mean[l] <- lambda_run$lambdatracker[max_yrs]
  #SR_long_mean[l] <- lambda_run$SRtracker[max_yrs]
  #OSR_long_mean[l] <- lambda_run$OSRtracker[max_yrs]
  #n0[,l] <- lambda_run$n0
  mat <- megamatrix_delay(F_params,M_params,twosex=F,long=long_seq[l],rfx=rfx_fun())$MEGAmat
  lambda_long_mean[l] <- lambda(mat)
  ssd <- stable.stage(mat)
  SR_long_mean[l] <- sum(ssd[1:(F_params$max_size+1)])
  n0[,l] <- ssd
  fem_panic <- ssd[2:(F_params$max_size+1)] %*% 
    (pfx(x=1:F_params$max_size,param=F_params,long=long_seq[l],rfx=rfx_fun()) *
    nfx(x=1:F_params$max_size,param=F_params,long=long_seq[l],rfx=rfx_fun()))
  male_panic <- ssd[(M_params$max_size+3):((M_params$max_size+1)*2)] %*% 
    (pfx(x=1:M_params$max_size,param=M_params,long=long_seq[l],rfx=rfx_fun()) *
       nfx(x=1:M_params$max_size,param=M_params,long=long_seq[l],rfx=rfx_fun()))
  OSR_long_mean[l] <- fem_panic / (fem_panic + male_panic)
}

par(mfrow=c(3,1))
plot(long_seq,lambda_long_mean,type="b",main=F_params$max_size)
plot(long_seq,SR_long_mean,type="b",ylim=c(0,1));abline(h=0.5)
plot(long_seq,OSR_long_mean,type="b",ylim=c(0,1));abline(h=0.5)

barplot(n0[2:(F_params$max_size+1),1])

plot(1:F_params$max_size,
     pfx(x=1:F_params$max_size,param=F_params,long=min(long_seq),rfx=rfx_fun()),type="l",col="hotpink")
lines(1:F_params$max_size,
      pfx(x=1:F_params$max_size,param=F_params,long=max(long_seq),rfx=rfx_fun()),lty=2,col="hotpink")
lines(1:M_params$max_size,
     pfx(x=1:M_params$max_size,param=M_params,long=min(long_seq),rfx=rfx_fun()),col="dodgerblue")
lines(1:M_params$max_size,
      pfx(x=1:M_params$max_size,param=M_params,long=max(long_seq),rfx=rfx_fun()),lty=2,col="dodgerblue")

plot(1:F_params$max_size,
     nfx(x=1:F_params$max_size,param=F_params,long=min(long_seq),rfx=rfx_fun()),type="l",col="hotpink")
lines(1:F_params$max_size,
      nfx(x=1:F_params$max_size,param=F_params,long=max(long_seq),rfx=rfx_fun()),lty=2,col="hotpink")
lines(1:M_params$max_size,
      nfx(x=1:M_params$max_size,param=M_params,long=min(long_seq),rfx=rfx_fun()),col="dodgerblue")
lines(1:M_params$max_size,
      nfx(x=1:M_params$max_size,param=M_params,long=max(long_seq),rfx=rfx_fun()),lty=2,col="dodgerblue")


# check size transitions --------------------------------------------------
library(plot.matrix)

poar.grow <- poar %>% 
  subset( tillerN_t0 > 0 & tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          tillerN_t0, tillerN_t1 ) %>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          unique.block = unique.block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>% 
  rename( sex   = Sex,
          block = unique.block ) %>% 
  mutate( log_size_t0   = log(tillerN_t0),
          log_size_t0_z = log(tillerN_t0) %>% scale %>% .[,1] )

size_trans <- as.matrix(table(poar.grow$tillerN_t0,poar.grow$tillerN_t1))
image(size_trans[1:F_params$max_size,1:F_params$max_size])

plot(density(poar.grow$tillerN_t0))

image(
t(outer(1:F_params$max_size,
        1:F_params$max_size,
        gxy,params=F_params,long=0,rfx=rfx_fun()))
)
