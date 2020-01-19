library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
library(countreg)

dir <- "C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits"
dir <- "C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits"
dir <- "C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits"

## read data
poar_dropsites <- read.csv(paste0(dir,"/Experiment/Demography/POAR-range-limits/data/demography.csv"), stringsAsFactors = F)
poar_allsites <- read.csv(paste0(dir,"/Experiment/Demography/POAR-range-limits/data/demography_allsites.csv"), stringsAsFactors = F)
viabVr <- read.csv(paste0(dir,"/Experiment/Demography/POAR-range-limits/data/viability.csv"))

## Load the Stan output for full vital rate model
#this model drops the three bad sites
#fit_dropsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_full.rds"))

#this model includes the three bad sites
fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds"))

#this model includes the three bad sites but drops the long^2 interactions terms (but keeps long^2)
#fit_allsites_full <- readRDS(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/fit_allsites_full_noLong2intx.rds"))


## read in the site latlong file to un-scale longitude
latlong <- read.csv(paste0(dir,"/Experiment/Demography/data/SiteLatLong.csv"))

## useful functions
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
invlogit<-function(x){exp(x)/(1+exp(x))}

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
predS <- rstan::extract(fit_allsites_full, pars = c("predS"))$predS
predG <- rstan::extract(fit_allsites_full, pars = c("predG"))$predG
phi_G <- rstan::extract(fit_allsites_full, pars = c("phi_g"))$phi_g
predF <- rstan::extract(fit_allsites_full, pars = c("predF"))$predF
predP <- rstan::extract(fit_allsites_full, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_allsites_full, pars = c("phi_p"))$phi_p
predV <- rstan::extract(fit_allsites_full, pars = c("predV"))$predV
predM <- rstan::extract(fit_allsites_full, pars = c("predM"))$predM

n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws)

y_s_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_s))
y_g_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_g))
y_f_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_f))
y_p_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_p))
y_v_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_v))
y_m_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_m))

for(i in 1:n_post_draws){
  ## sample survival data (bernoulli)
  y_s_sim[i,] <- rbinom(n=length(data_allsites_all$y_s), size=1, prob = invlogit(predS[i,]))
  ## sample growth data (zero-truncated NB)
  #y_g_sim[i,] <- rnbinom(n=length(data_allsites_all$y_g), mu = exp(predG[i,]), size=phi_G[i])
  ## sample flowering data (bernoulli)
  y_f_sim[i,] <- rbinom(n=length(data_allsites_all$y_f), size=1, prob = invlogit(predF[i,]))
  ## sample panicle data (zero-truncated NB)
  #y_p_sim[i,] <- rnbinom(n=length(data_allsites_all$y_p), mu = exp(predP[i,]), size=phi_P[i])
  ## sample viability data (binomial)
  y_v_sim[i,] <- rbinom(n=length(data_allsites_all$y_v), size=data_allsites_all$tot_seeds_v, prob = predV[i,])
  ## sample germination data (binomial)
  y_m_sim[i,] <- rbinom(n=length(data_allsites_all$y_m), size=data_allsites_all$tot_seeds_m, prob = predM[i,])
  
  for(j in 1:length(data_allsites_all$y_g)){
    y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predG[i,j]), size=phi_G[i]) / (1 - dnbinom(0, mu = exp(predG[i,j]), size=phi_G[i])))
  }
  for(j in 1:length(data_allsites_all$y_p)){
    y_p_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predP[i,j]), size=phi_P[i]) / (1 - dnbinom(0, mu = exp(predP[i,j]), size=phi_P[i])))
  }
}

ppc_dens_overlay(data_allsites_all$y_s, y_s_sim)
ppc_dens_overlay(data_allsites_all$y_g, y_g_sim)+xlim(0, 100)
ppc_dens_overlay(data_allsites_all$y_f, y_f_sim)
ppc_dens_overlay(data_allsites_all$y_p, y_p_sim)+xlim(0, 50)
ppc_dens_overlay(data_allsites_all$y_v, y_v_sim) ## maybe need beta-binomial?
ppc_dens_overlay(data_allsites_all$y_m, y_m_sim) ## maybe need beta-binomial?

mcmc_dens_overlay(fit_allsites_full,par=quote_bare(blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,phi_g))
mcmc_trace(fit_full,par=quote_bare(blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g))
mcmc_dens_overlay(fit_full,par=quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                          bsizesex_s, bsizelong_s,blongsex_s,bsizelongsex_s,
                                          blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s))

# Core vital rates --------------------------------------------------------
# estimate sex- and longitude-specific vital rates for small, medium, and large plants

##set which dataset to use
poar <- poar_allsites  #poar <- poar_dropsites #
fit_full <- fit_allsites_full  #fit_full <- fit_dropsites_full #

# first bin size groups
size_bin_num <- 3


## Survival
poar.surv <- poar %>% 
  subset(tillerN_t0>0 ) %>%
  select(year, Code, site, unique.block, Sex, 
         long.center, long.scaled, 
         tillerN_t0, surv_t1) %>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          unique.block = unique.block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric) %>% 
  rename( sex   = Sex,
          block = unique.block ) %>% 
  mutate( log_size_t0   = log(tillerN_t0),
          log_size_t0_z = log(tillerN_t0) %>% scale %>% .[,1] )

poar_surv_binned <- poar.surv %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_surv = mean(surv_t1),
            long = unique(long.center),
            bin_n = n())
surv_mean_sizes <- poar_surv_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))


## Growth
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

poar_grow_binned <- poar.grow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_grow = mean(tillerN_t1),
            long = unique(long.center),
            bin_n = n())
grow_mean_sizes <- poar_grow_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

## Flowering
poar.flow <- poar %>% 
  subset( tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          tillerN_t1, flow_t1 ) %>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          unique.block = unique.block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>% 
  rename( sex      = Sex,
          block    = unique.block ) %>% 
  mutate( log_size_t1   = log(tillerN_t1),
          log_size_t1_z = log(tillerN_t1) %>% scale %>% .[,1] )

poar_flow_binned <- poar.flow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_flow = mean(flow_t1),
            long = unique(long.center),
            bin_n = n())
flow_mean_sizes <- poar_flow_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

## Panicles
poar.panic<- poar %>% 
  subset( flowerN_t1 > 0 & tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          tillerN_t1, flowerN_t1 ) %>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          unique.block = unique.block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>% 
  rename( panic_t1 = flowerN_t1,
          sex      = Sex,
          block    = unique.block ) %>% 
  mutate( log_size_t1   = log(tillerN_t1),
          log_size_t1_z = log(tillerN_t1) %>% scale %>% .[,1] )

poar_panic_binned <- poar.panic %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_panic = mean(panic_t1),
            long = unique(long.center),
            bin_n = n())
panic_mean_sizes <- poar_panic_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))


# calculate posterior of sex difference ------------------------------------------------------------------------

# pull out stan coefficients
surv_coef <- rstan::extract(fit_full, pars = quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                                        bsizesex_s, bsizelong_s,blongsex_s,bsizelongsex_s,
                                                        blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s,
                                                        site_tau_s,block_tau_s,source_tau_s))

grow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_g,bsize_g,bsex_g,blong_g,
                                                        bsizesex_g, bsizelong_g,blongsex_g,bsizelongsex_g,
                                                        blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,
                                                        site_tau_g,block_tau_g,source_tau_g,
                                                        phi_g))

flow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_f,bsize_f,bsex_f,blong_f,
                                                        bsizesex_f, bsizelong_f,blongsex_f,bsizelongsex_f,
                                                        blong2_f,bsizelong2_f,blong2sex_f,bsizelong2sex_f,
                                                        site_tau_f,block_tau_f,source_tau_f))

panic_coef <- rstan::extract(fit_full, pars = quote_bare(b0_p,bsize_p,bsex_p,blong_p,
                                                         bsizesex_p, bsizelong_p,blongsex_p,bsizelongsex_p,
                                                         blong2_p,bsizelong2_p,blong2sex_p,bsizelong2sex_p,
                                                         site_tau_p,block_tau_p,source_tau_p))

long_seq <- seq(min(poar_surv_binned$long),max(poar_surv_binned$long),0.1)
n_post_draws <- 500
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
surv_sex_diff_post <- grow_sex_diff_post <- flow_sex_diff_post <- panic_sex_diff_post <- array(NA,dim=c(size_bin_num,length(long_seq),n_post_draws))

for(p in 1:n_post_draws){
  for(i in 1:size_bin_num){
    
    s=1;      
    fem_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                                (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                                (surv_coef$blong_s[post_draws[p]]) * long_seq +
                                (surv_coef$bsizelong_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * long_seq +
                                (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                                (surv_coef$blongsex_s[post_draws[p]]) * long_seq * (s-1) +
                                (surv_coef$bsizelongsex_s[post_draws[p]])  * long_seq * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                (surv_coef$blong2_s[post_draws[p]]) * (long_seq^2) +
                                (surv_coef$bsizelong2_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (long_seq^2) +
                                (surv_coef$blong2sex_s[post_draws[p]]) * (long_seq^2) * (s-1) +
                                (surv_coef$bsizelong2sex_s[post_draws[p]])  * (long_seq^2) * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i]
    )
    fem_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                   (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                   (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                   (grow_coef$blong_g[post_draws[p]]) * long_seq +
                   (grow_coef$bsizelong_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * long_seq +
                   (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                   (grow_coef$blongsex_g[post_draws[p]]) * long_seq * (s-1) +
                   (grow_coef$bsizelongsex_g[post_draws[p]])  * long_seq * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                   (grow_coef$blong2_g[post_draws[p]]) * (long_seq^2) +
                   (grow_coef$bsizelong2_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (long_seq^2) +
                   (grow_coef$blong2sex_g[post_draws[p]]) * (long_seq^2) * (s-1) +
                   (grow_coef$bsizelong2sex_g[post_draws[p]])  * (long_seq^2) * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i]
    )
    fem_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                        (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                        (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                        (flow_coef$blong_f[post_draws[p]]) * long_seq +
                        (flow_coef$bsizelong_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * long_seq +
                        (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                        (flow_coef$blongsex_f[post_draws[p]]) * long_seq * (s-1) +
                        (flow_coef$bsizelongsex_f[post_draws[p]])  * long_seq * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                        (flow_coef$blong2_f[post_draws[p]]) * (long_seq^2) +
                        (flow_coef$bsizelong2_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (long_seq^2) +
                        (flow_coef$blong2sex_f[post_draws[p]]) * (long_seq^2) * (s-1) +
                        (flow_coef$bsizelong2sex_f[post_draws[p]])  * (long_seq^2) * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i]
    )
    fem_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                   (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                   (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                   (panic_coef$blong_p[post_draws[p]]) * long_seq +
                   (panic_coef$bsizelong_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * long_seq +
                   (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                   (panic_coef$blongsex_p[post_draws[p]]) * long_seq * (s-1) +
                   (panic_coef$bsizelongsex_p[post_draws[p]])  * long_seq * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                   (panic_coef$blong2_p[post_draws[p]]) * (long_seq^2) +
                   (panic_coef$bsizelong2_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (long_seq^2) +
                   (panic_coef$blong2sex_p[post_draws[p]]) * (long_seq^2) * (s-1) +
                   (panic_coef$bsizelong2sex_p[post_draws[p]])  * (long_seq^2) * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i]
    )
    
    s=2;       
    male_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                                  (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                  (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                                  (surv_coef$blong_s[post_draws[p]]) * long_seq +
                                  (surv_coef$bsizelong_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * long_seq +
                                  (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                                  (surv_coef$blongsex_s[post_draws[p]]) * long_seq * (s-1) +
                                  (surv_coef$bsizelongsex_s[post_draws[p]])  * long_seq * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                  (surv_coef$blong2_s[post_draws[p]]) * (long_seq^2) +
                                  (surv_coef$bsizelong2_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (long_seq^2) +
                                  (surv_coef$blong2sex_s[post_draws[p]]) * (long_seq^2) * (s-1) +
                                  (surv_coef$bsizelong2sex_s[post_draws[p]])  * (long_seq^2) * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i]
    )
    male_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                   (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                   (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                   (grow_coef$blong_g[post_draws[p]]) * long_seq +
                   (grow_coef$bsizelong_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * long_seq +
                   (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                   (grow_coef$blongsex_g[post_draws[p]]) * long_seq * (s-1) +
                   (grow_coef$bsizelongsex_g[post_draws[p]])  * long_seq * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                   (grow_coef$blong2_g[post_draws[p]]) * (long_seq^2) +
                   (grow_coef$bsizelong2_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (long_seq^2) +
                   (grow_coef$blong2sex_g[post_draws[p]]) * (long_seq^2) * (s-1) +
                   (grow_coef$bsizelong2sex_g[post_draws[p]])  * (long_seq^2) * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i]
    )
    male_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                        (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                        (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                        (flow_coef$blong_f[post_draws[p]]) * long_seq +
                        (flow_coef$bsizelong_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * long_seq +
                        (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                        (flow_coef$blongsex_f[post_draws[p]]) * long_seq * (s-1) +
                        (flow_coef$bsizelongsex_f[post_draws[p]])  * long_seq * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                        (flow_coef$blong2_f[post_draws[p]]) * (long_seq^2) +
                        (flow_coef$bsizelong2_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (long_seq^2) +
                        (flow_coef$blong2sex_f[post_draws[p]]) * (long_seq^2) * (s-1) +
                        (flow_coef$bsizelong2sex_f[post_draws[p]])  * (long_seq^2) * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i]
    )
    male_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                   (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                   (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                   (panic_coef$blong_p[post_draws[p]]) * long_seq +
                   (panic_coef$bsizelong_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * long_seq +
                   (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                   (panic_coef$blongsex_p[post_draws[p]]) * long_seq * (s-1) +
                   (panic_coef$bsizelongsex_p[post_draws[p]])  * long_seq * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                   (panic_coef$blong2_p[post_draws[p]]) * (long_seq^2) +
                   (panic_coef$bsizelong2_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (long_seq^2) +
                   (panic_coef$blong2sex_p[post_draws[p]]) * (long_seq^2) * (s-1) +
                   (panic_coef$bsizelong2sex_p[post_draws[p]])  * (long_seq^2) * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i]
    )
    
    surv_sex_diff_post[i,,p] <-  fem_s-male_s
    grow_sex_diff_post[i,,p] <-  fem_g-male_g
    flow_sex_diff_post[i,,p] <-  fem_f-male_f
    panic_sex_diff_post[i,,p] <-  fem_p-male_p
  }
}

sex_diff_surv_mean <- sex_diff_grow_mean <- sex_diff_flow_mean <- sex_diff_panic_mean <- matrix(NA,size_bin_num,length(long_seq))
sex_diff_surv_95 <- sex_diff_surv_75 <- sex_diff_surv_50 <- sex_diff_surv_25 <- array(NA,dim=c(size_bin_num,length(long_seq),2))
sex_diff_grow_95 <- sex_diff_grow_75 <- sex_diff_grow_50 <- sex_diff_grow_25 <- array(NA,dim=c(size_bin_num,length(long_seq),2))
sex_diff_flow_95 <- sex_diff_flow_75 <- sex_diff_flow_50 <- sex_diff_flow_25 <- array(NA,dim=c(size_bin_num,length(long_seq),2))
sex_diff_panic_95 <- sex_diff_panic_75 <- sex_diff_panic_50 <- sex_diff_panic_25 <- array(NA,dim=c(size_bin_num,length(long_seq),2))

for(s in 1:size_bin_num){
  for(l in 1:length(long_seq)){
    sex_diff_surv_mean[s,l] <- mean(surv_sex_diff_post[s,l,])
    sex_diff_surv_95[s,l,] <- quantile(surv_sex_diff_post[s,l,],probs=c(0.025,0.975))
    sex_diff_surv_75[s,l,] <- quantile(surv_sex_diff_post[s,l,],probs=c(0.125,0.875))
    sex_diff_surv_50[s,l,] <- quantile(surv_sex_diff_post[s,l,],probs=c(0.25,0.75))
    sex_diff_surv_25[s,l,] <- quantile(surv_sex_diff_post[s,l,],probs=c(0.375,0.625))
    
    sex_diff_grow_mean[s,l] <- mean(grow_sex_diff_post[s,l,])
    sex_diff_grow_95[s,l,] <- quantile(grow_sex_diff_post[s,l,],probs=c(0.025,0.975))
    sex_diff_grow_75[s,l,] <- quantile(grow_sex_diff_post[s,l,],probs=c(0.125,0.875))
    sex_diff_grow_50[s,l,] <- quantile(grow_sex_diff_post[s,l,],probs=c(0.25,0.75))
    sex_diff_grow_25[s,l,] <- quantile(grow_sex_diff_post[s,l,],probs=c(0.375,0.625))
    
    sex_diff_flow_mean[s,l] <- mean(flow_sex_diff_post[s,l,])
    sex_diff_flow_95[s,l,] <- quantile(flow_sex_diff_post[s,l,],probs=c(0.025,0.975))
    sex_diff_flow_75[s,l,] <- quantile(flow_sex_diff_post[s,l,],probs=c(0.125,0.875))
    sex_diff_flow_50[s,l,] <- quantile(flow_sex_diff_post[s,l,],probs=c(0.25,0.75))
    sex_diff_flow_25[s,l,] <- quantile(flow_sex_diff_post[s,l,],probs=c(0.375,0.625))
    
    sex_diff_panic_mean[s,l] <- mean(panic_sex_diff_post[s,l,])
    sex_diff_panic_95[s,l,] <- quantile(panic_sex_diff_post[s,l,],probs=c(0.025,0.975))
    sex_diff_panic_75[s,l,] <- quantile(panic_sex_diff_post[s,l,],probs=c(0.125,0.875))
    sex_diff_panic_50[s,l,] <- quantile(panic_sex_diff_post[s,l,],probs=c(0.25,0.75))
    sex_diff_panic_25[s,l,] <- quantile(panic_sex_diff_post[s,l,],probs=c(0.375,0.625))
    
  }
}


# Mega combo figure -------------------------------------------------------
sex_cols <- c("black","white")
sex_lty <- c(1,2)
bin_shapes <- 15:17
diff_col <- "darkgrey"#"dodgerblue"#"black"
diff_alpha <- 0.15

layout.matrix <- rbind(matrix(1:6, nrow = 2, ncol = 3, byrow = F),
                       matrix(7:12, nrow = 2, ncol = 3, byrow = F),
                       matrix(13:18, nrow = 2, ncol = 3, byrow = F),
                       matrix(19:24, nrow = 2, ncol = 3, byrow = F))
pdf("Manuscript/Figures/vital_rates.pdf",height = 10,width = 8,useDingbats = F)
layout(mat = layout.matrix,
       heights = rep(c(2, 1),4), # Heights of the two rows
       widths = c(2, 2, 2))
#layout.show(24)

par(oma=c(3,1,0.5,0.5))
with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,1,0))
    plot(long[size_bin==i] + mean(latlong$Longitude),mean_surv[size_bin==i],type="n",ylim=c(0,1),
         xlab=" ",ylab=" ",xaxt="n");box()
    if(i==1){mtext("Pr(survival)",side=2,line=3)}
    title(LETTERS[i],adj=0)
    for(s in 1:2){
      points(long[sex==s & size_bin==i] + mean(latlong$Longitude),mean_surv[sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(long_seq + mean(latlong$Longitude),invlogit(mean(surv_coef$b0_s) + 
                                mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                mean(surv_coef$bsex_s) * (s-1) +
                                mean(surv_coef$blong_s) * long_seq +
                                mean(surv_coef$bsizelong_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * long_seq +
                                mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                                mean(surv_coef$blongsex_s) * long_seq * (s-1) +
                                mean(surv_coef$bsizelongsex_s)  * long_seq * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                                mean(surv_coef$blong2_s) * (long_seq^2) +
                                mean(surv_coef$bsizelong2_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (long_seq^2) +
                                mean(surv_coef$blong2sex_s) * (long_seq^2) * (s-1) +
                                mean(surv_coef$bsizelong2sex_s)  * (long_seq^2) * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i]
      ),lty=sex_lty[s],lwd=3)
    }
    par(mar=c(2,4,0.5,0))
    plot(long_seq + mean(latlong$Longitude),sex_diff_surv_mean[i,],type="n",lwd=4,
         cex.lab=1.6,xlab=" ",ylab=" ",xaxt="n",ylim=c(min(sex_diff_surv_95[i,,1]),max(sex_diff_surv_95[i,,2])))
    if(i==1){mtext(expression(paste(Delta," (F-M)")),side=2,line=3)}
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_surv_95[i,,1],rev(sex_diff_surv_95[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_surv_75[i,,1],rev(sex_diff_surv_75[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_surv_50[i,,1],rev(sex_diff_surv_50[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_surv_25[i,,1],rev(sex_diff_surv_25[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    lines(long_seq + mean(latlong$Longitude),sex_diff_surv_mean[i,],lwd=2,col=diff_col)
    abline(h=0,lty=2)
  }
})

with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,1,0))
    plot(long[size_bin==i] + mean(latlong$Longitude),mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()    
    if(i==1){mtext("#tillers",side=2,line=3)}
    title(LETTERS[i+3],adj=0)    
    for(s in 1:2){
      points(long[sex==s & size_bin==i] + mean(latlong$Longitude),mean_grow[sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(long_seq + mean(latlong$Longitude),
            exp(mean(grow_coef$b0_g) + 
                  mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                  mean(grow_coef$bsex_g) * (s-1) +
                  mean(grow_coef$blong_g) * long_seq +
                  mean(grow_coef$bsizelong_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * long_seq +
                  mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                  mean(grow_coef$blongsex_g) * long_seq * (s-1) +
                  mean(grow_coef$bsizelongsex_g)  * long_seq * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                  mean(grow_coef$blong2_g) * (long_seq^2) +
                  mean(grow_coef$bsizelong2_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (long_seq^2) +
                  mean(grow_coef$blong2sex_g) * (long_seq^2) * (s-1) +
                  mean(grow_coef$bsizelong2sex_g)  * (long_seq^2) * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i]
            ),lty=sex_lty[s],lwd=3)
    }
    par(mar=c(2,4,0.5,0))
    plot(long_seq + mean(latlong$Longitude),sex_diff_grow_mean[i,],type="n",lwd=4,
         cex.lab=1.6,xlab=" ",ylab=" ",xaxt="n",ylim=c(min(sex_diff_grow_95[i,,1]),max(sex_diff_grow_95[i,,2])))
    if(i==1){mtext(expression(paste(Delta," (F-M)")),side=2,line=3)}
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_grow_95[i,,1],rev(sex_diff_grow_95[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_grow_75[i,,1],rev(sex_diff_grow_75[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_grow_50[i,,1],rev(sex_diff_grow_50[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_grow_25[i,,1],rev(sex_diff_grow_25[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    lines(long_seq + mean(latlong$Longitude),sex_diff_grow_mean[i,],lwd=2,col=diff_col)
    abline(h=0,lty=2)
  }
})


with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,1,0))
    plot(long[size_bin==i] + mean(latlong$Longitude),mean_flow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()
    if(i==1){mtext("Pr(flower)",side=2,line=3)}
    title(LETTERS[i+6],adj=0)    
    for(s in 1:2){
      points(long[sex==s & size_bin==i] + mean(latlong$Longitude),mean_flow[sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(long_seq + mean(latlong$Longitude),
            invlogit(mean(flow_coef$b0_f) + 
                       mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                       mean(flow_coef$bsex_f) * (s-1) +
                       mean(flow_coef$blong_f) * long_seq +
                       mean(flow_coef$bsizelong_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * long_seq +
                       mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                       mean(flow_coef$blongsex_f) * long_seq * (s-1) +
                       mean(flow_coef$bsizelongsex_f)  * long_seq * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                       mean(flow_coef$blong2_f) * (long_seq^2) +
                       mean(flow_coef$bsizelong2_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (long_seq^2) +
                       mean(flow_coef$blong2sex_f) * (long_seq^2) * (s-1) +
                       mean(flow_coef$bsizelong2sex_f)  * (long_seq^2) * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i]
            ),lty=sex_lty[s],lwd=3)
    }
    par(mar=c(2,4,0.5,0))
    plot(long_seq + mean(latlong$Longitude),sex_diff_flow_mean[i,],type="n",lwd=4,
         cex.lab=1.6,xlab=" ",ylab=" ",xaxt="n",ylim=c(min(sex_diff_flow_95[i,,1]),max(sex_diff_flow_95[i,,2])))
    if(i==1){mtext(expression(paste(Delta," (F-M)")),side=2,line=3)}
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_flow_95[i,,1],rev(sex_diff_flow_95[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_flow_75[i,,1],rev(sex_diff_flow_75[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_flow_50[i,,1],rev(sex_diff_flow_50[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_flow_25[i,,1],rev(sex_diff_flow_25[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    lines(long_seq + mean(latlong$Longitude),sex_diff_flow_mean[i,],lwd=2,col=diff_col)
    abline(h=0,lty=2)
  }
})


with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,1,0))
    plot(long[size_bin==i] + mean(latlong$Longitude),mean_panic[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n",ylim=c(0,max(mean_panic[size_bin==i])));box()    
    if(i==1){mtext("#panicles",side=2,line=3)}
    title(LETTERS[i+9],adj=0)    
    for(s in 1:2){
      points(long[sex==s & size_bin==i] + mean(latlong$Longitude),mean_panic[sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=7*(bin_n/max(bin_n)),lwd=2)
      lines(long_seq + mean(latlong$Longitude),
            exp(mean(panic_coef$b0_p) + 
                  mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                  mean(panic_coef$bsex_p) * (s-1) +
                  mean(panic_coef$blong_p) * long_seq +
                  mean(panic_coef$bsizelong_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * long_seq +
                  mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                  mean(panic_coef$blongsex_p) * long_seq * (s-1) +
                  mean(panic_coef$bsizelongsex_p)  * long_seq * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                  mean(panic_coef$blong2_p) * (long_seq^2) +
                  mean(panic_coef$bsizelong2_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (long_seq^2) +
                  mean(panic_coef$blong2sex_p) * (long_seq^2) * (s-1) +
                  mean(panic_coef$bsizelong2sex_p)  * (long_seq^2) * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i]
            ),lty=sex_lty[s],lwd=3)
    }
    par(mar=c(2,4,0.5,0))
    plot(long_seq + mean(latlong$Longitude),sex_diff_panic_mean[i,],type="n",lwd=4,
         cex.lab=1.6,xlab=" ",ylab=" ",ylim=c(min(sex_diff_panic_95[i,,1]),max(sex_diff_panic_95[i,,2])))
    if(i==1){mtext(expression(paste(Delta," (F-M)")),side=2,line=3)}
    mtext("Longitude",side=1,line=3)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_panic_95[i,,1],rev(sex_diff_panic_95[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_panic_75[i,,1],rev(sex_diff_panic_75[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_panic_50[i,,1],rev(sex_diff_panic_50[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
            y=c(sex_diff_panic_25[i,,1],rev(sex_diff_panic_25[i,,2])),
            col=alpha(diff_col,diff_alpha),border=NA)
    lines(long_seq + mean(latlong$Longitude),sex_diff_panic_mean[i,],lwd=2,col=diff_col)
    abline(h=0,lty=2)
  }
})
dev.off()


# OSR trends in common garden ---------------------------------------------

garden_osr <- poar %>% select( year, site, Sex, long.center, flowerN_t1) %>% 
  group_by(year, long.center, site, Sex) %>% 
  summarise(total_panicles = sum(flowerN_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_panicles) %>% 
  rename(fem_pan = `F`,
         mal_pan = `M`) %>% 
  mutate(tot_pan=fem_pan+mal_pan,
         osr = fem_pan/tot_pan,
         year_index = year-2014)

ggplot(garden_osr)+
  geom_point(aes(x=long.center,y=osr,color=as.factor(year)))

garden_sr <- poar %>% select( year, site, Sex, long.center, surv_t1) %>% 
  group_by(year, long.center, site, Sex) %>% 
  summarise(total_plants = sum(surv_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_plants) %>% 
  rename(fem = `F`,
         mal = `M`) %>% 
  mutate(total=fem+mal,
         sr = fem/total)

ggplot(garden_sr)+
  geom_point(aes(x=long.center,y=sr,color=as.factor(year)))

## quick glm fit (might do this in rstan later)
osr_glm <- glm(osr ~ long.center*as.factor(year), weights = tot_pan, family="binomial",data=garden_osr)
sr_glm <- glm(sr ~ long.center*as.factor(year), weights = total, family="binomial",data=garden_sr)


year_cols <- c("#e41a1c","#377eb8","#4daf4a")

pdf("Manuscript/Figures/garden_sex_ratios.pdf",useDingbats = F,height=6,width=12)
par(mfrow=c(1,2),mar=c(5,6,1,1))
plot(garden_osr$long.center + mean(latlong$Longitude),garden_osr$osr,type="n",ylim=c(0,1),cex.lab=1.4,
     xlab="Longitude",ylab="Operational sex ratio\n(proportion female panicles)")
abline(h=0.5,col="gray",lty=2);title("A",adj=0,font=4)
points(garden_osr$long.center[garden_osr$year==2015] + mean(latlong$Longitude),
       garden_osr$osr[garden_osr$year==2015],
       pch=16,cex=log(garden_osr$tot_pan[garden_osr$year==2015]),col=alpha(year_cols[1],0.5))
points(garden_osr$long.center[garden_osr$year==2016] + mean(latlong$Longitude),
       garden_osr$osr[garden_osr$year==2016],
       pch=16,cex=log(garden_osr$tot_pan[garden_osr$year==2016]),col=alpha(year_cols[2],0.5))
points(garden_osr$long.center[garden_osr$year==2017] + mean(latlong$Longitude),
       garden_osr$osr[garden_osr$year==2017],
       pch=16,cex=log(garden_osr$tot_pan[garden_osr$year==2017]),col=alpha(year_cols[3],0.5))
lines(long_seq + mean(latlong$Longitude),invlogit(coef(osr_glm)[1] + coef(osr_glm)[2]*long_seq),col=year_cols[1],lwd=3)
lines(long_seq + mean(latlong$Longitude),invlogit(coef(osr_glm)[1] + coef(osr_glm)[3] + (coef(osr_glm)[2]+coef(osr_glm)[5])*long_seq),col=year_cols[2],lwd=3)
lines(long_seq + mean(latlong$Longitude),invlogit(coef(osr_glm)[1] + coef(osr_glm)[4] + (coef(osr_glm)[2]+coef(osr_glm)[6])*long_seq),col=year_cols[3],lwd=3)
legend("topleft",legend=2015:2017,bty="n",lwd=3,col=year_cols,cex=1.2)

plot(garden_sr$long.center + mean(latlong$Longitude),garden_sr$sr,type="n",ylim=c(0,1),cex.lab=1.4,
     xlab="Longitude",ylab="Sex ratio\n(proportion female plants)")
abline(h=0.5,col="gray",lty=2);title("B",adj=0,font=4)
points(garden_sr$long.center[garden_sr$year==2015] + mean(latlong$Longitude),
       garden_sr$sr[garden_sr$year==2015],
       pch=16,cex=log(garden_sr$total[garden_sr$year==2015]),col=alpha(year_cols[1],0.5))
points(garden_sr$long.center[garden_sr$year==2016] + mean(latlong$Longitude),
       garden_sr$sr[garden_sr$year==2016],
       pch=16,cex=log(garden_sr$total[garden_sr$year==2016]),col=alpha(year_cols[2],0.5))
points(garden_sr$long.center[garden_sr$year==2017] + mean(latlong$Longitude),
       garden_sr$sr[garden_sr$year==2017],
       pch=16,cex=log(garden_sr$total[garden_sr$year==2017]),col=alpha(year_cols[3],0.5))
lines(long_seq + mean(latlong$Longitude),invlogit(coef(sr_glm)[1] + coef(sr_glm)[2]*long_seq),col=year_cols[1],lwd=3)
lines(long_seq + mean(latlong$Longitude),invlogit(coef(sr_glm)[1] + coef(sr_glm)[3] + (coef(sr_glm)[2]+coef(sr_glm)[5])*long_seq),col=year_cols[2],lwd=3)
lines(long_seq + mean(latlong$Longitude),invlogit(coef(sr_glm)[1] + coef(sr_glm)[4] + (coef(sr_glm)[2]+coef(sr_glm)[6])*long_seq),col=year_cols[3],lwd=3)
legend("topleft",legend=2015:2017,bty="n",lwd=3,col=year_cols,cex=1.2)
dev.off()

# Seed viability ----------------------------------------------------------
viab   <- viabVr %>% 
  select( plot, totS, yesMaybe, sr_f ) %>% 
  rename( SR        = sr_f,
          y_viab = yesMaybe,
          tot_seeds_viab = totS) %>% 
  select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit

#viab_pars <- rstan::summary(fit_full, pars=c("v0","a_v"))[[1]][,"mean"]
viab_pars <- rstan::extract(fit_full, pars = quote_bare(v0,a_v,m,lambda_d))

pdf("Manuscript/Figures/seed_viability.pdf",useDingbats = F)
plot(jitter(viab$SR,75),jitter((viab$y_viab / viab$tot_seeds_viab),75),
     type="n",xlab="Operational sex ratio (fraction female panicles)",
     ylab="Seed viability",cex.lab=1.4)
for(p in 1:n_post_draws){
  lines(seq(0,1,0.01),
        viab_pars$v0[post_draws[p]] * (1 - seq(0,1,0.01) ^ viab_pars$a_v[post_draws[p]]),
        col=alpha("darkgrey",0.1))
}
points(jitter(viab$SR,75),jitter((viab$y_viab / viab$tot_seeds_viab),75),
       cex = 5 * (viab$tot_seeds_viab / max(viab$tot_seeds_viab)),lwd=2)
dev.off()

# Site climate variation --------------------------------------------------
library(SPEI)

## read in site coordinates
poar_sites <- read.csv("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/POAR_sites_2014-2017M.csv", stringsAsFactors = F)
poar_sites <- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/POAR_sites_2014-2017M.csv", stringsAsFactors = F)

## calculate SPEI
poar_temp <- poar_sites %>% 
  select(ID1,Year,Latitude,Tave01:Tave12)%>% 
  gather(Tave01:Tave12,key="Month",value="Tave") %>% 
  mutate(Month = as.integer(str_sub(Month,5,6)))%>% 
  arrange(ID1,Year) %>% 
  mutate(transition_year = ifelse(Month>=5,Year,Year-1)) %>% 
  filter(transition_year %in% 2014:2017)
  
poar_prcp <- poar_sites %>% 
  select(ID1,Year,PPT01:PPT12)%>% 
  gather(PPT01:PPT12,key="Month",value="PPT") %>% 
  mutate(Month = as.integer(str_sub(Month,4,5)))%>% 
  arrange(ID1,Year) %>% 
  mutate(transition_year = ifelse(Month>=5,Year,Year-1)) %>% 
  filter(transition_year %in% 2014:2017)

poar_clim <- full_join(poar_temp,poar_prcp,by=c("ID1","Month", "Year","transition_year"))%>% 
  group_by(ID1) %>% 
  mutate(PET = thornthwaite(Tave,unique(Latitude)),
         BAL = PPT-PET) %>% 
  select(ID1,Month,transition_year,BAL)

## set up SPEI matrix
sites <- poar_clim %>% select(ID1) %>% unique()
poar_spei <- matrix(NA,nrow=nrow(sites),ncol=3)
for(i in 1:nrow(poar_spei)){
  poar_spei[i,] <- spei(ts(poar_clim[poar_clim$ID1==sites$ID1[i],'BAL']),12)$fitted[c(12,24,36)]
}
allsites_spei <- data.frame(poar_spei); names(allsites_spei)<-c("2014","2015","2016")
allsites_spei$site <- sites$ID1
allsites_spei <- allsites_spei %>% gather(`2014`:`2016`,key=transition_year,value=spei) %>% 
  mutate(transition_year=as.integer(transition_year))

allsites_prcp <- poar_prcp %>% group_by(ID1,transition_year) %>% summarise(ppt = sum(PPT)) %>% 
  rename(site=ID1) %>% filter(transition_year!=2017)

allsites_temp <- poar_temp %>% group_by(ID1,transition_year) %>% summarise(temp = mean(Tave)) %>% 
  rename(site=ID1) %>% filter(transition_year!=2017)

allsites <- full_join(
  full_join(
    full_join(allsites_prcp,allsites_temp,by=c("site","transition_year")),
    allsites_spei,by=c("site","transition_year")),
                      poar_sites %>% group_by(ID1) %>% summarise(lat = unique(Latitude),
                                           long = unique(Longitude)) %>% rename(site=ID1),
  by=c("site"))

ggplot(allsites)+
  geom_point(aes(x=long,y=spei,color=as.factor(transition_year)))

pdf("Manuscript/Figures/site_precip.pdf",useDingbats = F)
plot(allsites$long,allsites$ppt,type="n",ylab="Annual precipitation (mm)",xlab="Longitude",cex.lab=1.4)
points(allsites$long[allsites$transition_year==2014],allsites$ppt[allsites$transition_year==2014],cex=2,pch=1)
points(allsites$long[allsites$transition_year==2015],allsites$ppt[allsites$transition_year==2015],cex=2,pch=2)
points(allsites$long[allsites$transition_year==2016],allsites$ppt[allsites$transition_year==2016],cex=2,pch=3)
legend("topleft",legend=c(2014:2016),bty="n",pch=1:3,cex=1.2)
## add 30-year normals
dev.off()

# Natural population survey -----------------------------------------------
## The data manipulation steps here come from the script 'POAR_population_preproposal_figure.R', which is a 
## mashup of code from me, Aldo, and maybe some Emily Begnel.

#create the "big" POAR file including data from both years
POAR12_1=read.csv(file="C:/Users/tm9/Dropbox/POAR--Aldo&Tom/POAR GIS Final Data_2012(1.10.2014).csv")
POAR13=read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/POAR GIS Final Data_2013.csv")
POAR12=read.csv(file="C:/Users/tm9/Dropbox/POAR--Aldo&Tom/POAR GIS Final Data.csv")
#POAR13=POAR13[-1,]   #This excludes data from Clear Bay in 2013 

#NOTE Compare the two
POAR12_1=POAR12_1[order(POAR12_1$Population),]
POAR12=POAR12[order(POAR12$Population),c(1:7)]
#DISCREPANCIES with:
#Clear bay: original file misses polygon data, mine is more current
#Wichita Mountains: original file misses point data, mine is more current
#Quartz Mountain: original file has polygons data that I do not have. I therefore use the POAR12 data)
POAR12_1[10,6:7]=POAR12[10,6:7]  #I substitude POAR12 data to my current file
POAR12=POAR12_1
POAR12$year="2012" ; POAR13$year="2013"   

#COMBINE CLEAR BAY AND THUNDERBIRD LAKE
tmp=rbind(POAR12[14,],POAR13[1,])
newDf=data.frame("Population"="ClearBay-Thunderbird", 
                 "Latitude"=apply(tmp[,2:3],2,FUN=mean)[1],"Longitude"=apply(tmp[,2:3],2,FUN=mean)[2],
                 "Total.Points"=apply(tmp[,4:7],2,FUN=sum)[1],"Total.Polygons"=apply(tmp[,4:7],2,FUN=sum)[2],
                 "Total.Male.Infl."=apply(tmp[,4:7],2,FUN=sum)[3],
                 "Total.Female.Infl."=apply(tmp[,4:7],2,FUN=sum)[4],"year"="2013",row.names=16
)

#TWO GRAPHS
#POAR=rbind(POAR12,POAR13)                     #Separate Clear bay and Thunderbird lake
POAR=rbind(POAR12[-14,],rbind(POAR13[-1,],newDf))        #Clear bay and Thunderbird lake Combined
POAR=POAR[order(as.character(POAR$Population)),]
#POAR=rbind(POAR12[,c(1:7,12)],POAR13) #I only include shared columns 
#POAR=rbind(POAR12[,c(1:7,12)],POAR13) #

## Fit Stan model (should be quick and simple)
survey_dat <- list(y = POAR$Total.Female.Infl.,
                   n_trials = POAR$Total.Male.Infl.+POAR$Total.Female.Infl.,
                   longit = POAR$Longitude - mean(POAR$Longitude),
                   N = nrow(POAR))
# simulation parameters
sim_pars <- list(
  warmup = 5000, 
  iter = 25000, 
  thin = 3, 
  chains = 3
)
# fit the "big" model 
fit_surv <- stan(
  file = 'code/stan/poar_survey.stan',
  data = survey_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

summary(fit_surv)$summary
coef_surv <- rstan::extract(fit_surv,par=c("b0","b_long"))
x_long <- seq(min(survey_dat$longit),max(survey_dat$longit),length = 100)

n_post <-500
sample_post <- sample.int(length(coef_surv$b0),n_post)
pdf("Manuscript/Figures/pop_survey.pdf",useDingbats = F)
plot(survey_dat$longit + mean(POAR$Longitude),
     survey_dat$y/survey_dat$n_trials,type="n",
     xlab="Longitude",ylab="Proportion female inflorescences",cex.lab=1.4)
for(i in 1:n_post){
  lines(x_long + mean(POAR$Longitude),
        invlogit(coef_surv$b0[sample_post[i]] + coef_surv$b_long[sample_post[i]] * x_long),
        col=alpha("darkgrey",0.05))
}
points(survey_dat$longit + mean(POAR$Longitude),
       survey_dat$y/survey_dat$n_trials,
       cex=4*(survey_dat$n_trials/max(survey_dat$n_trials))+1,pch=16)
dev.off()

# Lambda-Longitude mean --------------------------------------------------------
## load functions
source("code/twosexMPM.R")

# find the seedling survival rate from our POAU data
POAU <- read.csv(paste0(dir,"/Experiment/Demography/POAR-range-limits/data/POAU.csv"))
sdlg_surv <- POAU %>% filter(year_recruit==year_t,year_t %in% 2014:2016) %>% summarise(sdlg_surv = mean(spring_survival_t1,na.rm=T))

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

F_params <- M_params <- c()
## survival
F_params$surv_mu <- mean_coef$b0_s
F_params$surv_size <- mean_coef$bsize_s
F_params$surv_long <- mean_coef$blong_s 
F_params$surv_size_long <- mean_coef$bsizelong_s 
F_params$surv_long2 <- mean_coef$blong2_s 
F_params$surv_size_long2 <- mean_coef$bsizelong2_s  #0
M_params$surv_mu <- mean_coef$b0_s + mean_coef$bsex_s  
M_params$surv_size <- mean_coef$bsize_s + mean_coef$bsizesex_s 
M_params$surv_long <- mean_coef$blong_s + mean_coef$blongsex_s 
M_params$surv_size_long <- mean_coef$bsizelong_s + mean_coef$bsizelongsex_s 
M_params$surv_long2 <- mean_coef$blong2_s + mean_coef$blong2sex_s #
M_params$surv_size_long2 <- mean_coef$bsizelong2_s + mean_coef$bsizelong2sex_s  #0#
## growth
F_params$grow_mu <- mean_coef$b0_g #+14
F_params$grow_size <- mean_coef$bsize_g 
F_params$grow_long <- mean_coef$blong_g 
F_params$grow_size_long <- mean_coef$bsizelong_g 
F_params$grow_long2 <- mean_coef$blong2_g 
F_params$grow_size_long2 <- mean_coef$bsizelong2_g #0#
F_params$phi_g <- mean_coef$phi_g 
M_params$grow_mu <- mean_coef$b0_g + mean_coef$bsex_g #+14
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
## use POAU seedling survival for females and males
F_params$sdlg_surv <- M_params$sdlg_surv <- sdlg_surv$sdlg_surv
## set max size equal between the sexes
F_params$max_size <- M_params$max_size <- quantile(na.omit(poar$tillerN_t1),probs=0.95) #max(na.omit(poar$tillerN_t0)); 

## long_seq <- seq(min(poar_surv_binned$long),max(poar_surv_binned$long),0.3)
## extend the longitude to bounds and beyond of observed distribution
long_seq <- seq((-104 - mean(latlong$Longitude)),(-94.5 - mean(latlong$Longitude)),0.2)
lambda_long_mean<-lambda_long_mean_2sex<-SR_long_mean<-OSR_long_mean<-c()
n0<-matrix(NA,(F_params$max_size+1)*2,length(long_seq))
max_yrs <- 50

for(l in 1:length(long_seq)){
  #lambda_run <- lambdaSim_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)
  #lambda_long_mean[l] <- lambda_run$lambdatracker[max_yrs]
  #SR_long_mean[l] <- lambda_run$SRtracker[max_yrs]
  #OSR_long_mean[l] <- lambda_run$OSRtracker[max_yrs]
  #n0[,l] <- lambda_run$n0
  
  mat <- megamatrix_delay( F_params,
                           M_params,
                           twosex=F,
                           long=long_seq[l],
                           rfx=rfx_fun())$MEGAmat

  lambda_long_mean[l] <- lambda(mat)
  
  lambda_long_mean_2sex[l] <- lambdaSim_delay( F_params=F_params,
                                               M_params=M_params,
                                               long = long_seq[l],
                                               rfx = rfx_fun(),
                                               max.yrs = max_yrs)$lambdatracker[max_yrs]
  
  #ssd <- stable.stage(mat)
  #SR_long_mean[l] <- sum(ssd[1:(F_params$max_size+1)])
  #n0[,l] <- ssd
  #fem_panic <- ssd[2:(F_params$max_size+1)] %*% 
  #  (pfx(x=1:F_params$max_size,param=F_params,long=long_seq[l],rfx=rfx_fun()) *
  #     nfx(x=1:F_params$max_size,param=F_params,long=long_seq[l],rfx=rfx_fun()))
  #male_panic <- ssd[(M_params$max_size+3):((M_params$max_size+1)*2)] %*% 
  #  (pfx(x=1:M_params$max_size,param=M_params,long=long_seq[l],rfx=rfx_fun()) *
  #     nfx(x=1:M_params$max_size,param=M_params,long=long_seq[l],rfx=rfx_fun()))
  #OSR_long_mean[l] <- fem_panic / (fem_panic + male_panic)
}

par(mfrow=c(3,1))
plot(long_seq
     + mean(latlong$Longitude),
     lambda_long_mean,
     type="l",
     main=F_params$max_size)
lines(long_seq+mean(latlong$Longitude),lambda_long_mean_2sex,type="l",lty=2)
abline(v=c(-103.252677, -95.445907)) # brewster and brazoria county 
plot(long_seq,SR_long_mean,type="b",ylim=c(0,1));abline(h=0.5)
plot(long_seq,OSR_long_mean,type="b",ylim=c(0,1));abline(h=0.5)

## compare 2sex and Fdom models
pdf("Manuscript/Figures/lambda_long_2sex_Fdom.pdf",useDingbats = F)
par(mar=c(5,5,1,1))
plot(long_seq+mean(latlong$Longitude),lambda_long_mean,type="l",lwd=2,
     xlab="Longitude",ylab=expression(paste(lambda)),cex.lab=1.5,ylim=c(1,2.2))
lines(long_seq+mean(latlong$Longitude),lambda_long_mean_2sex,type="l",lwd=2,lty=2)
legend("topright",legend=c("Two-sex","Female-dominant"),lwd=2,lty=c(2,1),bty="n",cex=1.2)
dev.off()


# obtain stable stage distribution for comparison to observed data ----

# stable stage distribution at each longitude
ssd <- matrix(NA, 56, length(long_seq) )

for(l in 1:length(long_seq)){
  
  # stable stage distribution
  ssd[,l] <- lambdaSim_delay( F_params = F_params,
                              M_params = M_params,
                              long     = long_seq[l],
                              rfx      = rfx_fun(),
                              max.yrs  = max_yrs)$n0[,1]
  
}

# let's assume the population is really big
sim_size_distrib <- function( ii ){

  n_round <- c( round( ssd[,ii] * 10000 )[1:28 ],
                round( ssd[,ii] * 10000 )[29:56] )
  
  Map( function(x,y) rep(x,y),
          rep(log(1:28), 2),
          n_round ) %>% 
      unlist %>% 
      data.frame( log_tiller_n = . )

}

size_df <- lapply(1:48, sim_size_distrib ) %>% bind_rows

write.csv(size_df, 'results/ssd.csv', row.names=F)


## LTRE analysis using mean params
## LTRE parameters:
# these are the indices of the intercepts and slopes
LTRE_params <- c(1,2,7,8,14,15,20,21)
# these are the effects of longitude on the intercepts and slopes
betas_long <- c(F_params[c("surv_long","surv_size_long","grow_long","grow_size_long","flow_long","flow_size_long","panic_long","panic_size_long")],
                M_params[c("surv_long","surv_size_long","grow_long","grow_size_long","flow_long","flow_size_long","panic_long","panic_size_long")])
# these are the effects of longitude^2 on the intercepts and slopes
betas_long2 <- c(F_params[c("surv_long2","surv_size_long2","grow_long2","grow_size_long2","flow_long2","flow_size_long2","panic_long2","panic_size_long2")],
                M_params[c("surv_long2","surv_size_long2","grow_long2","grow_size_long2","flow_long2","flow_size_long2","panic_long2","panic_size_long2")])

## Sensitivity of parameters to longitude -- this is just the derivative of a second-order polynomial
dp_dlong_fun <- function(beta_long,beta_long2,long){return(beta_long + 2*beta_long2*long)}

## loop over longitudes and calculated dp_dlong and dlambda_dp
dp_dlong <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(long_seq))
lambda_long <- c()
perturbation <- 0.01
for(l in 1:length(long_seq)){
  #lambda_long[l] <- lambda(megamatrix_delay(F_params,M_params,twosex=F,long=long_seq[l],rfx=rfx_fun())$MEGAmat)
  lambda_long[l] <- lambdaSim_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dlong[,l] <- dp_dlong_fun(beta_long = unlist(betas_long), 
                             beta_long2 = unlist(betas_long2), 
                             long = long_seq[l])
    #cannot vectorize this part unfortunately
    for(p in 1:8){
      F_params_perturb <- F_params; M_params_perturb <- M_params;  
      F_params_perturb[LTRE_params[p]] <- unlist(F_params[LTRE_params[p]]) + perturbation
      #lambda_perturb <- lambda(megamatrix_delay(F_params_perturb,M_params_perturb,twosex=F,long=long_seq[l],rfx=rfx_fun())$MEGAmat)
      lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
      dlambda_dp[p,l] <- (lambda_perturb - lambda_long[l]) / perturbation
    }
    for(p in 9:16){
      F_params_perturb <- F_params; M_params_perturb <- M_params;  
      M_params_perturb[LTRE_params[p-8]] <- unlist(M_params[LTRE_params[p-8]]) + perturbation
      #lambda_perturb <- lambda(megamatrix_delay(F_params_perturb,M_params_perturb,twosex=F,long=long_seq[l],rfx=rfx_fun())$MEGAmat)
      lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
      dlambda_dp[p,l] <- (lambda_perturb - lambda_long[l]) / perturbation
    }
  print(l)
}

# put it all together 
LTRE_out <- dp_dlong * dlambda_dp

par(mfrow=c(2,1))
plot(long_seq+mean(latlong$Longitude),lambda_long_mean,type="b",main=F_params$max_size)

plot(long_seq+mean(latlong$Longitude),colSums(LTRE_out[1:8,]),type="l",lwd=4);abline(h=0)
lines(long_seq[1:47]+mean(latlong$Longitude),
      (lambda_long_mean[2:48]-lambda_long_mean[1:47])/(long_seq[2] - long_seq[1]),col="red",main=F_params$max_size)

# Lambda-Longitude posterior sampling --------------------------------------------------------

## set up output matrices
lambda_long_post <- Fdom_lambda_long <- SR_long_post <- OSR_long_post <- matrix(NA,nrow=n_post_draws,ncol=length(long_seq))
lambda_long_post_rfx <- Fdom_lambda_long_rfx <- SR_long_post_rfx <- OSR_long_post_rfx <- matrix(NA,nrow=n_post_draws,ncol=length(long_seq))

F_params <- M_params <- list()
for(p in 1:n_post_draws){
  #set up param vectors
  ## survival
  F_params$surv_mu <- surv_coef$b0_s[post_draws[p]] 
  F_params$surv_size <- surv_coef$bsize_s[post_draws[p]] 
  F_params$surv_long <- surv_coef$blong_s[post_draws[p]] 
  F_params$surv_size_long <- surv_coef$bsizelong_s[post_draws[p]] 
  F_params$surv_long2 <- surv_coef$blong2_s[post_draws[p]] 
  F_params$surv_size_long2 <- surv_coef$bsizelong2_s[post_draws[p]] 
  M_params$surv_mu <- surv_coef$b0_s[post_draws[p]] + surv_coef$bsex_s[post_draws[p]]
  M_params$surv_size <- surv_coef$bsize_s[post_draws[p]] + surv_coef$bsizesex_s[post_draws[p]]
  M_params$surv_long <- surv_coef$blong_s[post_draws[p]] + surv_coef$blongsex_s[post_draws[p]]
  M_params$surv_size_long <- surv_coef$bsizelong_s[post_draws[p]] + surv_coef$bsizelongsex_s[post_draws[p]]
  M_params$surv_long2 <- surv_coef$blong2_s[post_draws[p]] + surv_coef$blong2sex_s[post_draws[p]]
  M_params$surv_size_long2 <- surv_coef$bsizelong2_s[post_draws[p]] + surv_coef$bsizelong2sex_s[post_draws[p]]
  ## growth
  F_params$grow_mu <- grow_coef$b0_g[post_draws[p]] 
  F_params$grow_size <- grow_coef$bsize_g[post_draws[p]] 
  F_params$grow_long <- grow_coef$blong_g[post_draws[p]] 
  F_params$grow_size_long <- grow_coef$bsizelong_g[post_draws[p]] 
  F_params$grow_long2 <- grow_coef$blong2_g[post_draws[p]] 
  F_params$grow_size_long2 <- grow_coef$bsizelong2_g[post_draws[p]] 
  F_params$phi_g <- grow_coef$phi_g[post_draws[p]] 
  M_params$grow_mu <- grow_coef$b0_g[post_draws[p]] + grow_coef$bsex_g[post_draws[p]]
  M_params$grow_size <- grow_coef$bsize_g[post_draws[p]] + grow_coef$bsizesex_g[post_draws[p]]
  M_params$grow_long <- grow_coef$blong_g[post_draws[p]] + grow_coef$blongsex_g[post_draws[p]]
  M_params$grow_size_long <- grow_coef$bsizelong_g[post_draws[p]] + grow_coef$bsizelongsex_g[post_draws[p]]
  M_params$grow_long2 <- grow_coef$blong2_g[post_draws[p]] + grow_coef$blong2sex_g[post_draws[p]]
  M_params$grow_size_long2 <- grow_coef$bsizelong2_g[post_draws[p]] + grow_coef$bsizelong2sex_g[post_draws[p]]
  M_params$phi_g <- grow_coef$phi_g[post_draws[p]] 
  ## flowering
  F_params$flow_mu <- flow_coef$b0_f[post_draws[p]] 
  F_params$flow_size <- flow_coef$bsize_f[post_draws[p]] 
  F_params$flow_long <- flow_coef$blong_f[post_draws[p]] 
  F_params$flow_size_long <- flow_coef$bsizelong_f[post_draws[p]] 
  F_params$flow_long2 <- flow_coef$blong2_f[post_draws[p]] 
  F_params$flow_size_long2 <- flow_coef$bsizelong2_f[post_draws[p]] 
  M_params$flow_mu <- flow_coef$b0_f[post_draws[p]] + flow_coef$bsex_f[post_draws[p]]
  M_params$flow_size <- flow_coef$bsize_f[post_draws[p]] + flow_coef$bsizesex_f[post_draws[p]]
  M_params$flow_long <- flow_coef$blong_f[post_draws[p]] + flow_coef$blongsex_f[post_draws[p]]
  M_params$flow_size_long <- flow_coef$bsizelong_f[post_draws[p]] + flow_coef$bsizelongsex_f[post_draws[p]]
  M_params$flow_long2 <- flow_coef$blong2_f[post_draws[p]] + flow_coef$blong2sex_f[post_draws[p]]
  M_params$flow_size_long2 <- flow_coef$bsizelong2_f[post_draws[p]] + flow_coef$bsizelong2sex_f[post_draws[p]]
  ## panicles
  F_params$panic_mu <- panic_coef$b0_p[post_draws[p]] 
  F_params$panic_size <- panic_coef$bsize_p[post_draws[p]] 
  F_params$panic_long <- panic_coef$blong_p[post_draws[p]] 
  F_params$panic_size_long <- panic_coef$bsizelong_p[post_draws[p]] 
  F_params$panic_long2 <- panic_coef$blong2_p[post_draws[p]] 
  F_params$panic_size_long2 <- panic_coef$bsizelong2_p[post_draws[p]] 
  M_params$panic_mu <- panic_coef$b0_p[post_draws[p]] + panic_coef$bsex_p[post_draws[p]]
  M_params$panic_size <- panic_coef$bsize_p[post_draws[p]] + panic_coef$bsizesex_p[post_draws[p]]
  M_params$panic_long <- panic_coef$blong_p[post_draws[p]] + panic_coef$blongsex_p[post_draws[p]]
  M_params$panic_size_long <- panic_coef$bsizelong_p[post_draws[p]] + panic_coef$bsizelongsex_p[post_draws[p]]
  M_params$panic_long2 <- panic_coef$blong2_p[post_draws[p]] + panic_coef$blong2sex_p[post_draws[p]]
  M_params$panic_size_long2 <- panic_coef$bsizelong2_p[post_draws[p]] + panic_coef$bsizelong2sex_p[post_draws[p]]
  ## seed viability and misc fertility params
  F_params$v0 <- viab_pars$v0[post_draws[p]] 
  F_params$a_v <- viab_pars$a_v[post_draws[p]] 
  F_params$ov_per_inf <- viab_pars$lambda_d[post_draws[p]] 
  F_params$germ <- viab_pars$m[post_draws[p]] 
  F_params$PSR <- 0.5
  ## use POAU seedling survival for females and males
  F_params$sdlg_surv <- M_params$sdlg_surv <- sdlg_surv$sdlg_surv
  ## set max size equal between the sexes
  F_params$max_size <- M_params$max_size <- quantile(na.omit(poar$tillerN_t0),probs=0.99) #max(na.omit(poar$tillerN_t0)); 

    ## pull out the rfx variances
  rfx <- rfx_fun(site_tau_s = surv_coef$site_tau_s[post_draws[p]],
                   block_tau_s = surv_coef$block_tau_s[post_draws[p]],
                   source_tau_s = surv_coef$source_tau_s[post_draws[p]],
                   site_tau_g = grow_coef$site_tau_g[post_draws[p]],
                   block_tau_g = grow_coef$block_tau_g[post_draws[p]],
                   source_tau_g = grow_coef$source_tau_g[post_draws[p]],
                   site_tau_f = flow_coef$site_tau_f[post_draws[p]],
                   block_tau_f = flow_coef$block_tau_f[post_draws[p]],
                   source_tau_f = flow_coef$source_tau_f[post_draws[p]],
                   site_tau_p = panic_coef$site_tau_p[post_draws[p]],
                   block_tau_p = panic_coef$block_tau_p[post_draws[p]],
                   source_tau_p = panic_coef$source_tau_p[post_draws[p]])
  
  for(l in 1:length(long_seq)){
    Fdom_lambda_long[p,l] <- lambda(A = megamatrix_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),twosex=F)$MEGAmat)
    lambda_run <- lambdaSim_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)
    lambda_long_post[p,l] <- lambda_run$lambdatracker[max_yrs]
    SR_long_post[p,l] <- lambda_run$SRtracker[max_yrs]
    OSR_long_post[p,l] <- lambda_run$OSRtracker[max_yrs]

    #pick up here
    Fdom_lambda_long_rfx[p,l] <- lambda(A = megamatrix_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx,twosex=F)$MEGAmat)
    lambda_run_rfx <- lambdaSim_delay(F_params=F_params,M_params=M_params,long=long_seq[l],rfx=rfx,max.yrs=max_yrs)
    lambda_long_post_rfx[p,l] <- lambda_run_rfx$lambdatracker[max_yrs]
    SR_long_post_rfx[p,l] <- lambda_run_rfx$SRtracker[max_yrs]
    OSR_long_post_rfx[p,l] <- lambda_run_rfx$OSRtracker[max_yrs]
    print(c(p,l))
  }
}

## write output as a list
lambda_long_post_list <- list(Fdom_lambda_long=Fdom_lambda_long,lambda_long_post=lambda_long_post,SR_long_post=SR_long_post,OSR_long_post=OSR_long_post,
                              Fdom_lambda_long_rfx=Fdom_lambda_long_rfx,lambda_long_post_rfx=lambda_long_post_rfx,SR_long_post_rfx=SR_long_post_rfx,OSR_long_post_rfx=OSR_long_post_rfx)
#write_rds(lambda_long_post_list,paste0(dir,"/Experiment/Demography/POAR-range-limits/results/lambda_long_post_list.rds"))
lambda_long_post_list <- read_rds(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/lambda_long_post_list.rds"))

## visualize results
lambda_long_q95 <- lambda_long_q75 <- lambda_long_q50 <- lambda_long_q25 <- matrix(NA,2,length(long_seq))
lambda_long_rfx_q95 <- lambda_long_rfx_q75 <- lambda_long_rfx_q50 <- lambda_long_rfx_q25 <- matrix(NA,2,length(long_seq))
lambda_long_mean <- lambda_long_rfx_mean <- c()
  for(l in 1:length(long_seq)){
    lambda_long_q95[,l] <- quantile(lambda_long_post_list$lambda_long_post[,l],probs=c(0.025,0.975),na.rm = T)
    lambda_long_q75[,l] <- quantile(lambda_long_post_list$lambda_long_post[,l],probs=c(0.125,0.875),na.rm = T)
    lambda_long_q50[,l] <- quantile(lambda_long_post_list$lambda_long_post[,l],probs=c(0.25,0.75),na.rm = T)
    lambda_long_q25[,l] <- quantile(lambda_long_post_list$lambda_long_post[,l],probs=c(0.375,0.625),na.rm = T)
    lambda_long_mean[l] <- mean(lambda_long_post_list$lambda_long_post[,l],na.rm = T)
    
    lambda_long_rfx_q95[,l] <- quantile(lambda_long_post_list$lambda_long_post_rfx[,l],probs=c(0.025,0.975),na.rm = T)
    lambda_long_rfx_q75[,l] <- quantile(lambda_long_post_list$lambda_long_post_rfx[,l],probs=c(0.125,0.875),na.rm = T)
    lambda_long_rfx_q50[,l] <- quantile(lambda_long_post_list$lambda_long_post_rfx[,l],probs=c(0.25,0.75),na.rm = T)
    lambda_long_rfx_q25[,l] <- quantile(lambda_long_post_list$lambda_long_post_rfx[,l],probs=c(0.375,0.625),na.rm = T)
    lambda_long_rfx_mean[l] <- mean(lambda_long_post_list$lambda_long_post_rfx[,l],na.rm = T)
  }

polygon_col <- "darkgrey"
polygon_alpha <- 0.35
ltre_cols <- c("#e41a1c",NA,"#377eb8",NA,"#4daf4a",NA,"#984ea3",NA)
ltre_lty <- c(1,NA,2,NA,3,NA,5,NA)

pdf("Manuscript/Figures/lambda_long_ltre.pdf",useDingbats = F,height=10,width=4)
par(mar=c(5,5,1,1),mfrow=c(3,1))
plot(long_seq + mean(latlong$Longitude),lambda_long,type="n",ylim=c(0,4),lwd=4,
     xlab="Longitude",ylab=expression(paste(lambda)),cex.lab=1.5)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_q95[1,],rev(lambda_long_q95[2,])),
        col=alpha(polygon_col,polygon_alpha),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_q75[1,],rev(lambda_long_q75[2,])),
        col=alpha(polygon_col,polygon_alpha),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_q50[1,],rev(lambda_long_q50[2,])),
        col=alpha(polygon_col,polygon_alpha),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_q25[1,],rev(lambda_long_q25[2,])),
        col=alpha(polygon_col,polygon_alpha),border=NA)
abline(h=1,lty=3)
#lines(long_seq + mean(latlong$Longitude),lambda_long_mean)
abline(v=c(-103.252677,-95.445907),lwd=1) # brewster and brazoria county
title(main="A",adj=0)

plot(long_seq + mean(latlong$Longitude),colSums(LTRE_out[1:8,]),type="l",lwd=4,ylim=c(-0.4,0.5),
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"Longitude")),cex.lab=1.5)
abline(h=0,col="gray")
for(i in seq(1,8,by=2)){
  lines(long_seq + mean(latlong$Longitude),colSums(LTRE_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}
title(main="B",adj=0)

plot(long_seq + mean(latlong$Longitude),colSums(LTRE_out[9:16,]),type="l",lwd=4,ylim=c(-0.4,0.5),
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"Longitude")),cex.lab=1.5)
abline(h=0,col="gray")
for(i in seq(9,16,by=2)){
  lines(long_seq + mean(latlong$Longitude),colSums(LTRE_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
title(main="C",adj=0)
legend("topright",bty="n",legend=c("Survival","Growth","Flowering","Panicles","Total"),lwd=c(2,2,2,2,4),
       lty=c(na.omit(ltre_lty),1),col=c(na.omit(ltre_cols),"black"),cex=1.5)
dev.off()

plot(long_seq + mean(latlong$Longitude),lambda_long_rfx_mean,type="l",ylim=c(0,3),lwd=4)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_rfx_q95[1,],rev(lambda_long_rfx_q95[2,])),
        col=alpha("red",0.2),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_rfx_q75[1,],rev(lambda_long_rfx_q75[2,])),
        col=alpha("red",0.2),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_rfx_q50[1,],rev(lambda_long_rfx_q50[2,])),
        col=alpha("red",0.2),border=NA)
polygon(x=c(long_seq + mean(latlong$Longitude),rev(long_seq + mean(latlong$Longitude))),
        y=c(lambda_long_rfx_q25[1,],rev(lambda_long_rfx_q25[2,])),
        col=alpha("red",0.2),border=NA)
abline(h=1,col="darkgrey")

## draw random effects for site, block, source
#rfx <- data.frame(site = rnorm(4,0,c(surv_coef$site_tau_s,grow_coef$site_tau_g,flow_coef$site_tau_f,panic_coef$site_tau_p)),
#                  block = rnorm(4,0,c(surv_coef$block_tau_s,grow_coef$block_tau_g,flow_coef$block_tau_f,panic_coef$block_tau_p)),
#                  source = rnorm(4,0,c(surv_coef$source_tau_s,grow_coef$source_tau_g,flow_coef$source_tau_f,panic_coef$source_tau_p)))
#rownames(rfx) <- c("surv","grow","flow","panic")


plot(long_seq,Fdom_lambda_long[1,],type="n",ylim=c(0,10))
for(i in 1:44){
  lines(long_seq,Fdom_lambda_long[i,],col=alpha("black",0.2))
} 

plot(long_seq,lambda_long[1,],type="n",ylim=c(0,10))
for(i in 1:44){
  lines(long_seq,lambda_long[i,],col=alpha("black",0.2))
} 

plot((poar$tillerN_t0),poar$tillerN_t1)


F_params$max_size <- quantile(na.omit(poar$tillerN_t0),probs=0.95) #max(na.omit(poar$tillerN_t0)); 
M_params$max_size <- F_params$max_size
rfx <- data.frame(site = rep(0,4),
                  block = rep(0,4),
                  source = rep(0,4))
rownames(rfx) <- c("surv","grow","flow","panic")
lambda(megamatrix(F_params=F_params,M_params=M_params,long=long_seq[l],
                  rfx=rfx,twosex=F)$MEGAmat)


# MS quantities -----------------------------------------------------------

n_survey_pops <- POAR %>% select(Population,year) %>% group_by(year) %>% summarise(n=n())

survey_site_table <- POAR %>% select(Population,Latitude,Longitude,year) %>% remove_rownames() %>% 
  rename(Year_visited = year) %>% 
  mutate(Experimental_source = ifelse(Population=="CooperWMA"|
                                        Population=="Copper Breaks"|
                                        Population=="Quartz_Mountain"|
                                        Population=="Horn_Hill_Cemetery"|
                                        Population=="Lake Arrowhead"|
                                        Population=="South_Llano"|
                                        Population=="Sulfur_Springs","yes","no"))

n_sites <- poar %>% select(site) %>% n_distinct()

n_sources <- poar %>% select(Code) %>% n_distinct()

source_plants <- poar %>% select(Code,ID) %>% unique() %>% group_by(Code) %>% summarise(n = n_distinct(ID))

n_block <- poar %>% select(Block) %>% n_distinct()

n_plants_per_block <- poar %>% filter(year==2015) %>% select(unique.block,ID) %>% group_by(unique.block) %>% summarise(n=n())

poau_sdlg_surv <- round(sdlg_surv$sdlg_surv,2)

survey_range <- range(survey_dat$n_trials)

viab_seeds <- viab %>% mutate(viab_seeds = y_viab==0) %>% filter(SR==1) %>% summarise(mean(viab_seeds))

viab_n <- range(viabVr$totS,na.rm=T)

garden_osr_range <- range(garden_osr$tot_pan[garden_osr$tot_pan>0])

garden_sr_range <- range(garden_sr$total[garden_sr$total>0])

poar_ms_quantities <- list(
  n_survey_pops=n_survey_pops,
  survey_site_table=survey_site_table,
  n_sites=n_sites,
  n_sources=n_sources,
  source_plants=source_plants,
  n_block=n_block,
  n_plants_per_block=n_plants_per_block$n[1],
  poau_sdlg_surv=poau_sdlg_surv,
  survey_range=survey_range,
  viab_seeds=round(viab_seeds$`mean(viab_seeds)`,2),
  viab_n=viab_n,
  garden_osr_range=garden_osr_range,
  garden_sr_range=garden_sr_range
)

write_rds(poar_ms_quantities,"Manuscript/poar_ms_quantities.rds")
