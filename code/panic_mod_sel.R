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
poar      <- read.csv('data/demography.csv', stringsAsFactors = F)



# MODEL COMPARISON -----------------------------------------

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
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


# Panicules ------------------------------------------------
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

# data
data_l <- list( n_sites  = poar.panic$site %>% n_distinct,
                n_blocks = poar.panic$block %>% n_distinct,
                site_s   = poar.panic$site,
                block_s  = poar.panic$block,
                site_block_s = data.frame( site_i  = poar.panic$site,
                                           block_i = poar.panic$block ) %>% 
                                  unique %>% .$site_i,
                male_s   = poar.panic$sex-1,
                long_s   = poar.panic$long.center,
                size_s   = poar.panic$log_size_t1,
                y_s      = poar.panic$panic_t1,
                n_s      = nrow(poar.panic) )


# model names and model RDS objects
mod_names <- c( 's','b','sb_nest',
                # non-centered models
                paste0(c('s','b','sb_nest'),'_nc') 
               ) 
all_mods  <- paste0( 'code/stan/negbin/grow_', mod_names, '.stan' )
all_rds   <- gsub('\\.stan','.rds',all_mods)


# fit all models at once
all_fits1     <- lapply(all_rds, fit_all_mods)
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
