# Fit bayesian model containing all data
#library(scales)
#library(dplyr)
library(rstan)
# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
#library(shinystan)
library(tidyverse)
#library(loo)
library(bayesplot)
#library(gtools)
library(countreg)
#library(actuar)
library(rmutil)
library(actuar)
#options( stringsAsFactors = T)
#source('code/format/plot_binned_prop.R')

# quote a series of bare names
quote_bare <- function( ... ){
    substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
invlogit<-function(x){exp(x)/(1+exp(x))}

# read demographic data
#poar    <- read.csv('data/demography.csv', stringsAsFactors = F)
# Tom's laptop
poar <- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/demography.csv", stringsAsFactors = F)
#viabVr  <- read.csv('data/viability.csv')
viabVr <- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/viability.csv")

# Data formatting -------------------------------------------------

# Survival
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

# growth 
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

# flowering
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


# Panicules 
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


# viability data
viab   <- viabVr %>% 
  select( plot, totS, yesMaybe, sr_f ) %>% 
  rename( SR        = sr_f,
          y_viab = yesMaybe,
          tot_seeds_viab = totS) %>% 
  select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit

# germination data
germ   <- viabVr %>% 
  select( plot, germTot, germFail, sr_f ) %>% 
  rename( SR        = sr_f,
          y_germ    = germTot ) %>% 
  mutate(tot_seeds_germ = y_germ + germFail ) %>% 
  select(y_germ, tot_seeds_germ, SR ) %>% 
  na.omit

# seeds per panicle
seeds   <- viabVr %>% 
  select(SeedN)  %>% 
  na.omit


# # Viab - Germ test --------------------------------------------------------
# viab   <- viabVr %>% 
#   select( plot, totS, yesMaybe, sr_f ) %>% 
#   rename( SR        = sr_f,
#           y_viab = yesMaybe,
#           tot_seeds_viab = totS) %>% 
#   select(y_viab, tot_seeds_viab, SR ) %>% 
#   na.omit
# 
# germ   <- viabVr %>% 
#   select( plot, germTot, germFail, sr_f ) %>% 
#   rename( SR        = sr_f,
#           y_germ    = germTot ) %>% 
#   mutate(tot_seeds_germ = y_germ + germFail ) %>% 
#   select(y_germ, tot_seeds_germ, SR ) %>% 
#   na.omit
# 
# data_viab_germ <- list(
#   n_v       = nrow(viab),
#   y_v       = viab$y_viab,
#   tot_seeds_v = viab$tot_seeds_viab,
#   SR_v        = viab$SR,
#   
#   n_g       = nrow(germ),
#   y_g       = germ$y_germ,
#   tot_seeds_g = germ$tot_seeds_germ,
#   SR_g        = germ$SR  
# )
# 
# # simulation parameters
# sim_pars <- list(
#   warmup = 1000, 
#   iter = 4000, 
#   thin = 2, 
#   chains = 4
# )
# 
# # fit the "big" model 
# fit_viab <- stan(
#   file = 'code/stan/viab_germ_test_betabin.stan',
#   data = data_viab_germ,
#   #pars = quote_bare( v0, a_v, g ),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = 4 )
# 
# mcmc_dens_overlay(fit_viab,par=quote_bare( v0, a_v, g ))
# 
# # PPC
# predV <- rstan::extract(fit_viab, pars = c("predV"))$predV
# phi_v <- rstan::extract(fit_viab, pars = c("phi_v"))$phi_v
# predG <- rstan::extract(fit_viab, pars = c("predG"))$predG
# phi_g <- rstan::extract(fit_viab, pars = c("phi_g"))$phi_g
# 
# n_post_draws <- 500
# post_draws <- sample.int(dim(predV)[1], n_post_draws)
# y_v_sim <- matrix(NA,n_post_draws,length(data_viab_germ$y_v))
# y_g_sim <- matrix(NA,n_post_draws,length(data_viab_germ$y_g))
# 
# for(i in 1:n_post_draws){
#   ## sample viability data (binomial)
#   #y_v_sim[i,] <- rbinom(n=length(data_viab_germ$y_v), size=data_viab_germ$tot_seeds_v, prob = predV[i,])
#   y_v_sim[i,] <- rbetabinom(n=length(data_viab_germ$y_v), size=data_viab_germ$tot_seeds_v, m = predV[i,], s = phi_v[i])
#   ## sample germination data (binomial)
#   #y_g_sim[i,] <- rbinom(n=length(data_viab_germ$y_g), size=data_viab_germ$tot_seeds_g, prob = predG[i,])
#   y_g_sim[i,] <- rbetabinom(n=length(data_viab_germ$y_g), size=data_viab_germ$tot_seeds_g, m = predG[i,], s = phi_g[i])
# }
# ppc_dens_overlay(data_viab_germ$y_v, y_v_sim) ## maybe need beta-binomial?
# ppc_dens_overlay(data_viab_germ$y_g, y_g_sim) ## maybe need beta-binomial?

##############################################################



# data for model
data_all <- list( n_sites    = poar.surv$site %>% n_distinct,
                n_sources  = poar.surv$source %>% n_distinct(),
                
                # survival data
                n_blocks_s = poar.surv$block %>% n_distinct,
                site_s     = poar.surv$site,
                source_s =  poar.surv$source,
                block_s    = poar.surv$block,
                site_block_s = data.frame( site_i  = poar.surv$site,
                                           block_i = poar.surv$block ) %>% 
                                  unique %>% .$site_i,
                male_s     = poar.surv$sex-1,
                long_s     = poar.surv$long.center,
                size_s     = poar.surv$log_size_t0,
                y_s        = poar.surv$surv_t1,
                n_s        = poar.surv$surv_t1 %>% length,
                
                # growth data
                n_blocks_g = poar.grow$block %>% n_distinct,
                site_g     = poar.grow$site,
                source_g =  poar.grow$source,
                block_g    = poar.grow$block,
                site_block_g = data.frame( site_i  = poar.grow$site,
                                           block_i = poar.grow$block ) %>% 
                                  unique %>% .$site_i,
                male_g   = poar.grow$sex-1,
                long_g   = poar.grow$long.center,
                size_g   = poar.grow$log_size_t0,
                y_g      = poar.grow$tillerN_t1,
                n_g      = nrow(poar.grow),
                
                # flowering data
                n_blocks_f = poar.flow$block %>% n_distinct,
                site_f   = poar.flow$site,
                source_f =  poar.flow$source,
                block_f  = poar.flow$block,
                site_block_f = data.frame( site_i  = poar.flow$site,
                                           block_i = poar.flow$block ) %>% 
                                  unique %>% .$site_i,
                male_f   = poar.flow$sex-1,
                long_f   = poar.flow$long.center,
                size_f   = poar.flow$log_size_t1,
                y_f      = poar.flow$flow_t1,
                n_f      = nrow(poar.flow),
                
                # panicle data
                n_blocks_p = poar.panic$block %>% n_distinct,
                site_p   = poar.panic$site,
                source_p =  poar.panic$source,
                block_p  = poar.panic$block,
                site_block_p = data.frame( site_i  = poar.panic$site,
                                           block_i = poar.panic$block ) %>% 
                                  unique %>% .$site_i,
                male_p   = poar.panic$sex-1,
                long_p   = poar.panic$long.center,
                size_p   = poar.panic$log_size_t1,
                y_p      = poar.panic$panic_t1,
                n_p      = nrow(poar.panic),
                
                # viability
                n_v       = nrow(viab),
                y_v       = viab$y_viab,
                tot_seeds_v = viab$tot_seeds_viab,
                SR_v        = viab$SR,
                
                # germination
                n_m       = nrow(germ),
                y_m       = germ$y_germ,
                tot_seeds_m = germ$tot_seeds_germ,
                SR_m        = germ$SR    )


# simulation parameters
sim_pars <- list(
  warmup = 5000, 
  iter = 30000, 
  thin = 3, 
  chains = 3
)

# fit the "big" model 
 fit_full <- stan(
    file = 'code/stan/poar_full.stan',
    data = data_all,
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains )

#saveRDS(fit_full, 'C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_full.rds')
fit_full <- readRDS('C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_full.rds')

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
predS <- rstan::extract(fit_full, pars = c("predS"))$predS
predG <- rstan::extract(fit_full, pars = c("predG"))$predG
phi_G <- rstan::extract(fit_full, pars = c("phi_g"))$phi_g
predF <- rstan::extract(fit_full, pars = c("predF"))$predF
predP <- rstan::extract(fit_full, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_full, pars = c("phi_p"))$phi_p
predV <- rstan::extract(fit_full, pars = c("predV"))$predV
predM <- rstan::extract(fit_full, pars = c("predM"))$predM

n_post_draws <- 500
post_draws <- sample.int(dim(predS)[1], n_post_draws)

y_s_sim <- matrix(NA,n_post_draws,length(data_all$y_s))
y_g_sim <- matrix(NA,n_post_draws,length(data_all$y_g))
y_f_sim <- matrix(NA,n_post_draws,length(data_all$y_f))
y_p_sim <- matrix(NA,n_post_draws,length(data_all$y_p))
y_v_sim <- matrix(NA,n_post_draws,length(data_all$y_v))
y_m_sim <- matrix(NA,n_post_draws,length(data_all$y_m))

for(i in 1:n_post_draws){
  ## sample survival data (bernoulli)
  y_s_sim[i,] <- rbinom(n=length(data_all$y_s), size=1, prob = invlogit(predS[i,]))
  ## sample growth data (zero-truncated NB)
  y_g_sim[i,] <- rztnbinom(n=length(data_all$y_g), mu = exp(predG[i,]), size=phi_G[i])
  ## sample flowering data (bernoulli)
  y_f_sim[i,] <- rbinom(n=length(data_all$y_f), size=1, prob = invlogit(predF[i,]))
  ## sample panicle data (zero-truncated NB)
  y_p_sim[i,] <- rztnbinom(n=length(data_all$y_p), mu = exp(predP[i,]), size=phi_P[i])
  ## sample viability data (binomial)
  y_v_sim[i,] <- rbinom(n=length(data_all$y_v), size=data_all$tot_seeds_v, prob = predV[i,])
  ## sample germination data (binomial)
  y_m_sim[i,] <- rbinom(n=length(data_all$y_m), size=data_all$tot_seeds_m, prob = predM[i,])
}



# Now include "bad" sites -------------------------------------------------

poar_allsites <- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/demography_allsites.csv", stringsAsFactors = F)
poar_allsites <- read.csv("C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/demography_allsites.csv", stringsAsFactors = F)

# Survival
poar_allsites.surv <- poar_allsites %>% 
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

# growth 
poar_allsites.grow <- poar_allsites %>% 
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

# flowering
poar_allsites.flow <- poar_allsites %>% 
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


# Panicules 
poar_allsites.panic<- poar_allsites %>% 
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


# data for model
data_allsites_all <- list( n_sites    = poar_allsites.surv$site %>% n_distinct,
                  n_sources  = poar_allsites.surv$source %>% n_distinct(),
                  
                  # survival data
                  n_blocks_s = poar_allsites.surv$block %>% n_distinct,
                  site_s     = poar_allsites.surv$site,
                  source_s =  poar_allsites.surv$source,
                  block_s    = poar_allsites.surv$block,
                  site_block_s = data.frame( site_i  = poar_allsites.surv$site,
                                             block_i = poar_allsites.surv$block ) %>% 
                    unique %>% .$site_i,
                  male_s     = poar_allsites.surv$sex-1,
                  long_s     = poar_allsites.surv$long.center,
                  size_s     = poar_allsites.surv$log_size_t0,
                  y_s        = poar_allsites.surv$surv_t1,
                  n_s        = poar_allsites.surv$surv_t1 %>% length,
                  
                  # growth data
                  n_blocks_g = poar_allsites.grow$block %>% n_distinct,
                  site_g     = poar_allsites.grow$site,
                  source_g =  poar_allsites.grow$source,
                  block_g    = poar_allsites.grow$block,
                  site_block_g = data.frame( site_i  = poar_allsites.grow$site,
                                             block_i = poar_allsites.grow$block ) %>% 
                    unique %>% .$site_i,
                  male_g   = poar_allsites.grow$sex-1,
                  long_g   = poar_allsites.grow$long.center,
                  size_g   = poar_allsites.grow$log_size_t0,
                  y_g      = poar_allsites.grow$tillerN_t1,
                  n_g      = nrow(poar_allsites.grow),
                  
                  # flowering data
                  n_blocks_f = poar_allsites.flow$block %>% n_distinct,
                  site_f   = poar_allsites.flow$site,
                  source_f =  poar_allsites.flow$source,
                  block_f  = poar_allsites.flow$block,
                  site_block_f = data.frame( site_i  = poar_allsites.flow$site,
                                             block_i = poar_allsites.flow$block ) %>% 
                    unique %>% .$site_i,
                  male_f   = poar_allsites.flow$sex-1,
                  long_f   = poar_allsites.flow$long.center,
                  size_f   = poar_allsites.flow$log_size_t1,
                  y_f      = poar_allsites.flow$flow_t1,
                  n_f      = nrow(poar_allsites.flow),
                  
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
                  n_p      = nrow(poar_allsites.panic),
                  
                  # viability
                  n_v       = nrow(viab),
                  y_v       = viab$y_viab,
                  tot_seeds_v = viab$tot_seeds_viab,
                  SR_v        = viab$SR,
                  
                  # germination
                  n_m       = nrow(germ),
                  y_m       = germ$y_germ,
                  tot_seeds_m = germ$tot_seeds_germ,
                  SR_m        = germ$SR,
                  
                  # seeds per panicle
                  n_d = nrow(seeds),
                  y_d = seeds$SeedN)

# fit the "big" model 
fit_allsites_full <- stan(
  file = 'code/stan/poar_full.stan',
  data = data_allsites_all,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

#saveRDS(fit_allsites_full, 'C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds')
fit_allsites_full <- readRDS('C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds')
#fit_allsites_full <- readRDS('C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsites_full_noLong2intx.rds')
fit_allsites_full <- readRDS('C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsites_full.rds')

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
predS <- rstan::extract(fit_allsites_full, pars = c("predS"))$predS
predG <- rstan::extract(fit_allsites_full, pars = c("predG"))$predG
phi_G <- rstan::extract(fit_allsites_full, pars = c("phi_g"))$phi_g
predF <- rstan::extract(fit_allsites_full, pars = c("predF"))$predF
predP <- rstan::extract(fit_allsites_full, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_allsites_full, pars = c("phi_p"))$phi_p
predV <- rstan::extract(fit_allsites_full, pars = c("predV"))$predV
phi_V <- rstan::extract(fit_allsites_full, pars = c("phi_v"))$phi_v
predM <- rstan::extract(fit_allsites_full, pars = c("predM"))$predM
phi_M <- rstan::extract(fit_allsites_full, pars = c("phi_m"))$phi_m

n_post_draws <- 500
post_draws <- sample.int(dim(predS)[1], n_post_draws)

y_s_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_s))
y_g_sim <- y_g_sim_alt <- matrix(NA,n_post_draws,length(data_allsites_all$y_g))
y_f_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_f))
y_p_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_p))
y_v_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_v))
y_m_sim <- matrix(NA,n_post_draws,length(data_allsites_all$y_m))

for(i in 1:n_post_draws){
  ## sample survival data (bernoulli)
  y_s_sim[i,] <- rbinom(n=length(data_allsites_all$y_s), size=1, prob = invlogit(predS[i,]))
  ## sample flowering data (bernoulli)
  y_f_sim[i,] <- rbinom(n=length(data_allsites_all$y_f), size=1, prob = invlogit(predF[i,]))
  ## sample viability data (beta-binomial)
  y_v_sim[i,] <- rbetabinom(n=length(data_allsites_all$y_v), size=data_allsites_all$tot_seeds_v, m=predV[i,], s=phi_V[i])
  #y_v_sim[i,] <- rbinom(n=length(data_allsites_all$y_v), size=data_allsites_all$tot_seeds_v, prob = predV[i,])
  ## sample germination data (beta-binomial)
  y_m_sim[i,] <- rbetabinom(n=length(data_allsites_all$y_m), size=data_allsites_all$tot_seeds_m, m=predM[i,], s=phi_M[i])
  #y_m_sim[i,] <- rbinom(n=length(data_allsites_all$y_m), size=data_allsites_all$tot_seeds_m, prob = predM[i,])
  
  ## sample growth data (zero-truncated NB)
  for(j in 1:length(data_allsites_all$y_g)){
    y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predG[i,j]), size=phi_G[i]) / (1 - dnbinom(0, mu = exp(predG[i,j]), size=phi_G[i])))
    }
  ## sample panicle data (zero-truncated NB)
  for(j in 1:length(data_allsites_all$y_p)){
    y_p_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predP[i,j]), size=phi_P[i]) / (1 - dnbinom(0, mu = exp(predP[i,j]), size=phi_P[i])))
  }
}

## try a different approach for simulating from ZT-NB
for(i in 1:n_post_draws){
  for(j in 1:length(data_allsites_all$y_g)){
    
  for(k in 1:1000){
    g_draw <- rnbinom(n=1,mu = exp(predG[i,j]), size=phi_G[i])
    if(g_draw>0){y_g_sim_alt[i,j] <- g_draw; break}
  }
  }
}

## survival looks great
ppc_dens_overlay(data_allsites_all$y_s, y_s_sim)

## growth
bayesplot::ppc_dens_overlay(data_allsites_all$y_g, y_g_sim)+xlim(0, 20)
bayesplot::ppc_dens_overlay(data_allsites_all$y_g, y_g_sim_alt)+xlim(0, 20)

## flowering looks great
ppc_dens_overlay(data_allsites_all$y_f, y_f_sim)

## panicles--pretty good
ppc_dens_overlay(data_allsites_all$y_p, y_p_sim)+xlim(0, 20)

## seed viability
ppc_dens_overlay(data_allsites_all$y_v, y_v_sim) 
bayesplot::ppc_intervals(y = data_allsites_all$y_v, yrep = posterior_predict(fit_allsites_full),
                         x = data_allsites_all$SR_v) + ggplot2::xlab("SR")

## seed germination
ppc_dens_overlay(data_allsites_all$y_m, y_m_sim) 

mcmc_dens_overlay(fit_allsites_full,par=quote_bare(blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,phi_g))
mcmc_trace(fit_full,par=quote_bare(blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g))
mcmc_dens_overlay(fit_full,par=quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                          bsizesex_s, bsizelong_s,blongsex_s,bsizelongsex_s,
                                          blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s))


# Growth troubleshooting: Negative Binomial --------------------------------------------

#the PPC looks like the model is under-estimating growth (too many 1-tiller predictions)
#can a different distribution provide a better fit?
library(gamlss)
grow_gamlss_fit <- fitDist(data_allsites_all$y_g,type="counts")
summary(grow_gamlss_fit)
grow_gamlss_fit$fits
## I think this is telling me that a ZA-PIG is best
## ugh, but how to fit that as a regression

## instead, try making the dispersion param of the NB a function of size, sex, and long
# data for growth model
data_grow <- list( n_sites    = poar_allsites.grow$site %>% n_distinct,
                  n_sources  = poar_allsites.grow$source %>% n_distinct(),
                  
                  # growth data
                  n_blocks_g = poar_allsites.grow$block %>% n_distinct,
                  site_g     = poar_allsites.grow$site,
                  source_g =  poar_allsites.grow$source,
                  block_g    = poar_allsites.grow$block,
                  site_block_g = data.frame( site_i  = poar_allsites.grow$site,
                                             block_i = poar_allsites.grow$block ) %>% 
                    unique %>% .$site_i,
                  male_g   = poar_allsites.grow$sex-1,
                  long_g   = poar_allsites.grow$long.center,
                  size_g   = poar_allsites.grow$log_size_t0,
                  y_g      = poar_allsites.grow$tillerN_t1,
                  n_g      = nrow(poar_allsites.grow)
                  )


# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 10000, 
  thin = 3, 
  chains = 3
)

# fit the "big" model 
fit_grow <- stan(
  file = 'code/stan/poar_growth.stan',
  data = data_grow,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

# see how we did
predG <- rstan::extract(fit_grow, pars = c("predG"))$predG
dispG <- rstan::extract(fit_grow, pars = c("dispG"))$dispG

n_post_draws  <- 500
post_draws    <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim       <- y_g_sim_alt <- matrix( NA, 
                                        n_post_draws, 
                                        length(data_grow$y_g) )

# simulate datasets
for(i in 1:n_post_draws){
  ## sample growth data (zero-truncated NB)
  for(j in 1:length(data_allsites_all$y_g)){
    y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predG[i,j]), size=dispG[i,j]) / (1 - dnbinom(0, mu = exp(predG[i,j]), size=dispG[i,j])))
  }
}

# posterior predictive checks
ppc_dens_overlay(data_grow$y_g, y_g_sim) + 
  xlim(0, 50) +
  ggsave( 'nb_regression_truncated.png',
          width = 6.3, height = 6.3 )
## the NB still does not describe the growth data very well,
## even with more flexibility in the dispersion paramete



# Growth troubleshooting: PIG distribution ------------------------------------------------

# fit model
fit_grow_pig <- stan(
                      file = 'code/stan/poar_growth_pig_unequal.stan',
                      data = data_grow,
                      warmup = sim_pars$warmup,
                      iter = sim_pars$iter,
                      thin = sim_pars$thin,
                      chains = sim_pars$chains 
                      )

# OPTIONAL: read already fit model
dir          <- 'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/'
fit_grow_pig <- readRDS( paste0(dir, 'poar_grow_pig_unequal_idiv_cluster.rds') )

# posterior of predictions
predG <- rstan::extract(fit_grow_pig, pars = c("predG"))$predG
theta <- rstan::extract(fit_grow_pig, pars = c("theta"))$theta

# set up simulation of data
n_post_draws  <- 500
post_draws    <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim       <- y_g_sim_alt <- matrix( NA, 
                                        n_post_draws, 
                                        length(data_grow$y_g) )


# fimulate data
for(i in 1:n_post_draws){
  ## sample growth data (zero-truncated PIG)
  for(j in 1:length(data_allsites_all$y_g)){
    
    # probability without truncation
    prob_v <- dpois( 1:1000, 
                     lambda = (predG[i,j] * theta[i,j]) )
    
    # probability for trunctation (denominator)
    prob_t <- (1 - dpois(0, lambda = (predG[i,j] * theta[i,j]) ) )
    
    y_g_sim[i,j] <- sample( x=1:1000, 
                            size = 1, replace = T,
                            prob = prob_v / prob_t )
  }
}

# Posterior predictive check
ppc_dens_overlay(data_grow$y_g, y_g_sim) + 
  xlim(0, 50) +
  ggsave( 'pig_regression_truncated.png',
          width = 6.3, height = 6.3 )

# Posterior predictive checks on final full model with ZT-PIG growth ----
fit_allsites_full <- readRDS("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsite_full_pig_trunc.rds")

predS <- rstan::extract(fit_allsites_full, pars = c("predS"))$predS
predG <- rstan::extract(fit_allsites_full, pars = c("predG"))$predG
sigmaG <- rstan::extract(fit_allsites_full, pars = c("sigma"))$sigma
predF <- rstan::extract(fit_allsites_full, pars = c("predF"))$predF
predP <- rstan::extract(fit_allsites_full, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_allsites_full, pars = c("phi_p"))$phi_p
predV <- rstan::extract(fit_allsites_full, pars = c("predV"))$predV
phi_V <- rstan::extract(fit_allsites_full, pars = c("phi_v"))$phi_v
predM <- rstan::extract(fit_allsites_full, pars = c("predM"))$predM
phi_M <- rstan::extract(fit_allsites_full, pars = c("phi_m"))$phi_m

n_post_draws <- 500
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
  ## sample flowering data (bernoulli)
  y_f_sim[i,] <- rbinom(n=length(data_allsites_all$y_f), size=1, prob = invlogit(predF[i,]))
  ## sample viability data (beta-binomial)
  y_v_sim[i,] <- rbetabinom(n=length(data_allsites_all$y_v), size=data_allsites_all$tot_seeds_v, m=predV[i,], s=phi_V[i])
  ## sample germination data (beta-binomial)
  y_m_sim[i,] <- rbetabinom(n=length(data_allsites_all$y_m), size=data_allsites_all$tot_seeds_m, m=predM[i,], s=phi_M[i])
  ## sample growth data (zero-truncated PIG) -- generates data as legit PIG
  for(i in 1:n_post_draws){
    for(j in 1:length(data_allsites_all$y_g)){
      ## the pig function appears numerically unstable at low probabilities, so here is a hacky solution
      pig<-dpoisinvgauss(0:1000,mean=predG[i,j],shape=(sigmaG[i]*predG[i,j]))
      pig[is.nan(pig) | is.infinite(pig)] <- 0
      pig_trunc_prob <- pig[2:1001] / (1 - ((1 - sum(pig)) + pig[1]))
      y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=pig_trunc_prob)
    } 
  }
  ## sample growth data (zero-truncated PIG) -- generates data as Poissong-IG mixture
  #for(j in 1:length(data_allsites_all$y_g)){
  #  y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dpois(1:1000,lambda=(predG[i,j]*thetaG[i,j]))/(1-dpois(0,lambda=(predG[i,j]*thetaG[i,j]))))
  #}
  ## sample panicle data (zero-truncated NB)
  for(j in 1:length(data_allsites_all$y_p)){
    y_p_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predP[i,j]), size=phi_P[i]) / (1 - dnbinom(0, mu = exp(predP[i,j]), size=phi_P[i])))
  }
}


ppc_dens_overlay(data_allsites_all$y_s, y_s_sim)+
  xlab("Survival status")+ylab("Density")+
  ggtitle(("Survival"))+theme(legend.position = "none")->ppc_surv
ppc_dens_overlay(data_allsites_all$y_g, na.omit(y_g_sim))+xlim(0, 50)+
  xlab("Number of tillers")+ylab("Density")+
  ggtitle(("Growth"))+theme(legend.position = "none")->ppc_grow
ppc_dens_overlay(data_allsites_all$y_f, y_f_sim)+
  xlab("Flowering status")+ylab("Density")+
  ggtitle(("Flowering"))+theme(legend.position = "none")->ppc_flow
ppc_dens_overlay(data_allsites_all$y_p, y_p_sim)+xlim(0, 30)+
  xlab("Number of panicles")+ylab("Density")+
  ggtitle(("Panicles"))+theme(legend.position = "none")->ppc_panic
ppc_dens_overlay(data_allsites_all$y_v, y_v_sim)+
  xlab("Number of viable seeds")+ylab("Density")+
  ggtitle(("Seed viability"))+theme(legend.position = "none")->ppc_viab
ppc_dens_overlay(data_allsites_all$y_m, y_m_sim)+
  xlab("Number of germinants")+ylab("Density")+
  ggtitle(("Seed germination"))+theme(legend.position = "none")->ppc_germ

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf("Manuscript/Figures/PPC.pdf",useDingbats = F,height=10,width=6)
multiplot(ppc_surv,ppc_flow,ppc_viab,
          ppc_grow,ppc_panic,ppc_germ,cols=2)
dev.off()