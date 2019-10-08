# Fit bayesian model containing all data
library(scales)
library(dplyr)
library(rstan)
library(shinystan)
library(tidyverse)
library(loo)
library(bayesplot)
options( stringsAsFactors = T)
source('code/format/plot_binned_prop.R')
# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# quote a series of bare names
quote_bare <- function( ... ){
    substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# read demographic data
#poar    <- read.csv('data/demography.csv', stringsAsFactors = F)
# Tom's Cornell desktop
poar <- read.csv("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/demography.csv", stringsAsFactors = F)
# Tom's laptop
poar <- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/demography.csv", stringsAsFactors = F)
#viabVr  <- read.csv('data/viability.csv')
viabVr <- read.csv("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/viability.csv")
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

# Tom's survival model experiments ----------------------------------------
# data for model
data_l <- list( #response data
                y_s        = poar.surv$surv_t1,
                n_s        = poar.surv$surv_t1 %>% length,
                #covariates
                size_s     = poar.surv$log_size_t0,
                male_s     = poar.surv$sex-1,
                long_s     = poar.surv$long.center,
                n_sites    = poar.surv$site %>% n_distinct,
                site_s     = poar.surv$site,
                n_sources  = poar.surv$source %>% n_distinct(),
                source_s =  poar.surv$source,
                n_blocks = poar.surv$block %>% n_distinct(),
                block_s = poar.surv$block
                
)
# parameters to estimate
params <- quote_bare( b_0, b_size, b_sex, b_long,
                      b_size_sex, b_size_long, b_long_sex, b_size_long_sex,
                      site_tau, block_tau, source_tau)
# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 5000, 
  thin = 3, 
  chains = 4
)
# fit the toy model 
fit_mod <- stan(
  file = 'code/stan/surv_tom.stan',
  data = data_l,
  #pars = params,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

## assess fit and convergence
mcmc_dens_overlay(fit_mod,par=params)
mcmc_trace(fit_mod,par=params)
mcmc_intervals(fit_mod,par=params)

## posterior predictive check
y_s_new <- rstan::extract(fit_mod, pars = c("y_s_new"))
color_scheme_set("brightblue")
ppc_dens_overlay(data_l$y_s, y_s_new[[1]][1:5,])


surv_summary <- summary(fit_mod)
x_surv <- seq(min(poar.surv$log_size_t0),max(poar.surv$log_size_t0),length.out=100)
invlogit <- function(x){exp(x)/(1+exp(x))}

plot(jitter(poar.surv$log_size_t0),jitter(poar.surv$surv_t1))
lines(x_surv,invlogit(surv_summary$summary[,"mean"]["b_0"] + surv_summary$summary[,"mean"]["b_s"] * x_surv),lwd=3)

###########################################################################

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


# Viab - Germ test --------------------------------------------------------
viab   <- viabVr %>% 
  select( plot, totS, yesMaybe, sr_f ) %>% 
  rename( SR        = sr_f,
          y_viab = yesMaybe,
          tot_seeds_viab = totS) %>% 
  select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit

germ   <- viabVr %>% 
  select( plot, germTot, germFail, sr_f ) %>% 
  rename( SR        = sr_f,
          y_germ    = germTot ) %>% 
  mutate(tot_seeds_germ = y_germ + germFail ) %>% 
  select(y_germ, tot_seeds_germ, SR ) %>% 
  na.omit

data_viab_germ <- list(
  n_v       = nrow(viab),
  y_v       = viab$y_viab,
  tot_seeds_v = viab$tot_seeds_viab,
  SR_v        = viab$SR,
  
  n_g       = nrow(germ),
  y_g       = germ$y_germ,
  tot_seeds_g = germ$tot_seeds_germ,
  SR_g        = germ$SR  
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# fit the "big" model 
fit_viab <- stan(
  file = 'code/stan/viab_germ_test.stan',
  data = data_viab_germ,
  pars = quote_bare( v0, a_v, g ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 4 )

mcmc_dens_overlay(fit_viab,par=quote_bare( v0, a_v, g ))

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
  warmup = 2000, 
  iter = 25000, 
  thin = 3, 
  chains = 4
)

# fit the "big" model 
fit_full <- stan(
    file = 'code/stan/poar_full.stan',
    data = data_all,
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = 4 )



saveRDS(fit_mod, 'results/all_vr_site.rds')
