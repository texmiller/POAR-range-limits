## Fit size-dependent vital rates and build MPM
# library(R2jags)
# library(mcmcplots)
library(scales)
library(dplyr)
library(rstan)
library(shinystan)
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
poar    <- read.csv('data/demography.csv')


# growth ------------------------------------------------
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

# data
data_l <- list( n_sites  = poar.grow$site %>% n_distinct,
                n_blocks = poar.grow$block %>% n_distinct,
                site_s   = poar.grow$site,
                block_s  = poar.grow$block,
                site_block_s = data.frame( site_i  = poar.grow$site,
                                           block_i = poar.grow$block ) %>% 
                                  unique %>% .$site_i,
                male_s   = poar.grow$sex-1,
                long_s   = poar.grow$long.center,
                size_s   = poar.grow$log_size_t0,
                y_s      = poar.grow$tillerN_t1,
                n_s      = nrow(poar.grow) )


# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# general purpose function to fit all models
fit_all_mods <- function(mod_fil){
  
  # choose parameters to store
  if( grepl('_s.rds|_s_nc.rds', mod_fil) ){
    par_list <- c( 'site_a_s', 'b_s', 'b_l', 'b_sex', 'b_sex_l', 'log_lik' )
  } else {
    par_list <- c( 'block_a_s', 'b_s', 'b_l', 'b_sex', 'b_sex_l', 'log_lik' )
  }
    
  obj <- readRDS(mod_fil)
  
  fit_out <- sampling(
      object = obj,
      data = data_l,
      pars = par_list,
      warmup = sim_pars$warmup,
      iter = sim_pars$iter,
      thin = sim_pars$thin,
      chains = sim_pars$chains 
    )
  
   return(fit_out)
   
}


fit_surv <- stan(
  file = 'code/stan/grow/grow_b.stan',
  data = data_l,
  pars = quote_bare( block_a_s, b_s, b_l, b_sex, b_sex_l ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 1 #sim_pars$chains#,
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

all_mods  <- paste0( 'code/stan/negbin/grow_', mod_names, '.stan' )
all_rds   <- gsub('\\.stan','.rds',all_mods)

# fit all models at once
all_fits      <- list()
all_fits[[1]] <- fit_all_mods( all_rds[1] )
all_fits[[2]] <- fit_all_mods( all_rds[2] )
all_fits[[3]] <- fit_all_mods( all_rds[3] )
all_fits[[4]] <- fit_all_mods( all_rds[4] )
all_fits[[5]] <- fit_all_mods( all_rds[5] )
all_fits[[6]] <- fit_all_mods( all_rds[6] )

# Extract log likelihoods
log_liks   <- lapply(all_fits, extract_log_lik)

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                  setNames( paste0('loo_', mod_names) )
loo_df     <- loo::compare(loo_l$loo_s,     loo_l$loo_b,     
                           loo_l$loo_sb_nest,
                           loo_l$loo_s_nc,  loo_l$loo_b_nc,     
                           loo_l$loo_sb_nest_nc ) %>%
                as.data.frame 

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames( paste0('waic_', mod_names) )
waic_df   <- loo::compare(waic_l$waic_s,     waic_l$waic_b,     
                          waic_l$waic_sb_nest,
                          waic_l$waic_s_nc,  waic_l$waic_b_nc,     
                          waic_l$waic_sb_nest_nc ) %>%
                as.data.frame

