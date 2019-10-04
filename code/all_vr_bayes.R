# Fit bayesian model containing all data
library(scales)
library(dplyr)
library(rstan)
library(shinystan)
library(tidyverse)
library(loo)
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
#viabVr  <- read.csv('data/viability.csv')
viabVr <- read.csv("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/data/viability.csv")

# Data formatting -------------------------------------------------

# Survival
poar.surv <- poar %>% 
                subset(year, code, tillerN_t0>0 ) %>%
                select(site, unique.block, Sex, 
                       long.center, long.scaled, 
                       tillerN_t0, surv_t1) %>% 
                na.omit %>% 
                mutate( site         = site %>% as.factor %>% as.numeric,
                        unique.block = unique.block %>% as.factor %>% as.numeric,
                        Sex          = Sex %>% as.factor %>% as.numeric ) %>% 
                rename( sex   = Sex,
                        block = unique.block ) %>% 
                mutate( log_size_t0   = log(tillerN_t0),
                        log_size_t0_z = log(tillerN_t0) %>% scale %>% .[,1] )

# Tom's survival model experiments ----------------------------------------
# data for model
data_l <- list( y_s        = poar.surv$surv_t1,
                size_s     = poar.surv$log_size_t0,
                n_s        = poar.surv$surv_t1 %>% length
)
# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)
# fit the "big" model 
fit_mod <- stan(
  file = 'code/stan/surv_tom.stan',
  data = data_l,
  pars = quote_bare( b_0, b_s ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 4 )

surv_summary <- summary(fit_mod)
x_surv <- seq(min(poar.surv$log_size_t0),max(poar.surv$log_size_t0),length.out=100)
invlogit <- function(x){exp(x)/(1+exp(x))}

plot(jitter(poar.surv$log_size_t0),jitter(poar.surv$surv_t1))
lines(x_surv,invlogit(surv_summary$summary[,"mean"]["b_0"] + surv_summary$summary[,"mean"]["b_s"] * x_surv),lwd=3)

###########################################################################






# growth 
poar.grow <- poar %>% 
                subset( tillerN_t0 > 0 & tillerN_t1 > 0 ) %>%
                select( site, unique.block, Sex, 
                        long.center, long.scaled, 
                        tillerN_t0, tillerN_t1 ) %>% 
                na.omit %>% 
                mutate( site         = site %>% as.factor %>% as.numeric,
                        unique.block = unique.block %>% as.factor %>% as.numeric,
                        Sex          = Sex %>% as.factor %>% as.numeric ) %>% 
                rename( sex   = Sex,
                        block = unique.block ) %>% 
                mutate( log_size_t0   = log(tillerN_t0),
                        log_size_t0_z = log(tillerN_t0) %>% scale %>% .[,1] )

# flowering
poar.flow <- poar %>% 
                subset( tillerN_t1 > 0 ) %>%
                select( site, unique.block, Sex, 
                        long.center, long.scaled, 
                        tillerN_t1, flow_t1 ) %>% 
                na.omit %>% 
                mutate( site         = site %>% as.factor %>% as.numeric,
                        unique.block = unique.block %>% as.factor %>% as.numeric,
                        Sex          = Sex %>% as.factor %>% as.numeric ) %>% 
                rename( sex      = Sex,
                        block    = unique.block ) %>% 
                mutate( log_size_t1   = log(tillerN_t1),
                        log_size_t1_z = log(tillerN_t1) %>% scale %>% .[,1] )


# Panicules 
poar.panic<- poar %>% 
                subset( flowerN_t1 > 0 & tillerN_t1 > 0 ) %>%
                select( site, unique.block, Sex, 
                        long.center, long.scaled, 
                        tillerN_t1, flowerN_t1 ) %>% 
                na.omit %>% 
                mutate( site         = site %>% as.factor %>% as.numeric,
                        unique.block = unique.block %>% as.factor %>% as.numeric,
                        Sex          = Sex %>% as.factor %>% as.numeric ) %>% 
                rename( panic_t1 = flowerN_t1,
                        sex      = Sex,
                        block    = unique.block ) %>% 
                mutate( log_size_t1   = log(tillerN_t1),
                        log_size_t1_z = log(tillerN_t1) %>% scale %>% .[,1] )


# viability data
viab   <- viabVr %>% 
            select( plot, germTot, germFail, sr_f ) %>% 
            rename( SR        = sr_f,
                    y_germ    = germTot ) %>% 
            mutate( tot_seeds = y_germ + germFail ) %>% 
            select( y_germ, tot_seeds, SR ) %>% 
            na.omit


# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# fit the "big" model 
fit_mod <- stan(
    file = 'code/stan/germ.stan',
    data = data_v,
    pars = quote_bare( v0, a_g ),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = 4 )


# data for model
data_l <- list( n_sites    = poar.surv$site %>% n_distinct,
                n_blocks_s = poar.surv$block %>% n_distinct,
                site_s     = poar.surv$site,
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
                block_f  = poar.flow$block,
                site_block_f = data.frame( site_i  = poar.flow$site,
                                           block_i = poar.flow$block ) %>% 
                                  unique %>% .$site_i,
                male_f   = poar.flow$sex-1,
                long_f   = poar.flow$long.center,
                size_f   = poar.flow$log_size_t1,
                y_f      = poar.flow$flow_t1,
                n_f      = nrow(poar.flow),
                
                # panicule data
                n_blocks_p = poar.panic$block %>% n_distinct,
                site_p   = poar.panic$site,
                block_p  = poar.panic$block,
                site_block_p = data.frame( site_i  = poar.panic$site,
                                           block_i = poar.panic$block ) %>% 
                                  unique %>% .$site_i,
                male_p   = poar.panic$sex-1,
                long_p   = poar.panic$long.center,
                size_p   = poar.panic$log_size_t1,
                y_p      = poar.panic$panic_t1,
                n_p      = nrow(poar.panic),
                
                # germination
                n_v       = nrow(viab),
                y_v       = viab$y_germ,
                tot_seeds = viab$tot_seeds,
                SR        = viab$SR   )




# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# fit the "big" model 
fit_mod <- stan(
    file = 'code/stan/all_vr_s.stan',
    data = data_l,
    pars = quote_bare( site_a_u_s, b_s_s, b_l_s, b_sex_s, b_sex_l_s,
                       site_a_u_g, b_s_g, b_l_g, b_sex_g, b_sex_l_g,
                       site_a_u_f, b_s_f, b_l_f, b_sex_f, b_sex_l_f,
                       site_a_u_p, b_s_p, b_l_p, b_sex_p, b_sex_l_p,
                       v0,         a_v ),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = 4 )

saveRDS(fit_mod, 'results/all_vr_site.rds')
