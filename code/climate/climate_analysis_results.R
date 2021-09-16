## Author: Tom Miller
## Purpose: fit statistical models to vital rate data using rstan,
## simulate data from fitted models for posterior predictive checks,
## output stan fits for use in graphics and demographic model, generate
## vital rate figures.
## The stan model is not shown here but imported from another file (poar_full.stan).
## Last update: 5/12/2021

# load packages
library(rstan)
# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
library(tidyverse)
library(bayesplot)
library(countreg)
library(rmutil)
library(actuar)
library(SPEI)
library(shinystan)
in_dir  <- "C:/Users/ac22qawo/poar/climate_analysis/"

# read in stan fits
mod_mmatch   <- readRDS( paste0(in_dir, 'poar_climate-7786865-4_4_poar_full_climate_mismatch_nologlik.RDS') )
mod_long     <- readRDS( paste0(in_dir, 'poar_climate-7786866-1_1_poar_full_climate_nologlik.RDS') )
mod_noabio   <- readRDS( paste0(in_dir, 'poar_climate-7786869-3_3_poar_full_noabiotic.RDS') )
mod_clim     <- readRDS( paste0(in_dir, 'poar_climate-7786870-2_2_poar_full_climate_nologlik.RDS') )

# extract parameters
par_mmatch   <- rstan::extract(mod_mmatch) %>% as.data.frame
par_long     <- rstan::extract(mod_long) %>% as.data.frame
par_noabio   <- rstan::extract(mod_noabio) %>% as.data.frame
par_clim     <- rstan::extract(mod_clim) %>% as.data.frame


# effect of precipitation mismatch ----------------------------

# parameter names
par_v       <- c( 'blong_s', 'bmmatch_s',
                  'blong_g', 'bmmatch_g',
                  'blong_f', 'bmmatch_f',
                  'blong_p', 'bmmatch_p' )

# extract kernel density estimation
dens_extr <- function( x, id ){
  x %>% density %>% .[id] %>% .[[1]]
}

# Kernel density estimation
kernel_df <- dplyr::select( par_mmatch,
               par_v) %>% 
  setNames( c( 'Longitude_Survival',  'Mismatch_Survival',
               'Longitude_Growth',    'Mismatch_Growth',
               'Longitude_Flowering', 'Mismatch_Flowering',
               'Longitude_Panicles',  'Mismatch_Panicles' ) ) %>% 
  gather( param, value, 
          Longitude_Survival:Mismatch_Panicules ) %>% 
  group_by( param ) %>% 
  summarise( x = dens_extr(value, 1),   
             y = dens_extr(value, 2) ) %>% 
  ungroup %>% 
  separate( param, c('Parameter','vr'), sep='_' )
  
# plot 
kernel_df %>% 
  subset( Parameter == 'Mismatch' ) %>% 
  ggplot( ) +
  geom_line( aes(x,y),
             lwd = 2) +
  geom_vline( xintercept = 0,
              lty = 2) +
  facet_wrap( ~vr ) +
  theme_minimal() +
  labs( x = 'Effect of precipitation mismatch',
        y = 'Density') +
  scale_color_colorblind() +
  theme( axis.text    = element_text( size = 10),
         axis.title   = element_text( size = 15),
         # legend.text  = element_text( size = 15),
         # legend.title = element_text( size = 15),
         strip.text   = element_text( size = 15) ) +
  ggsave( 'C:/CODE/POAR-range-limits/results/Longitude vs. Mismatch.tiff',
          width = 6.3, height = 6.3, compression = 'lzw')


# effect of precipitation mismatch ----------------------------

par_v     <- c( 'site_tau_s', 'site_tau_g',
                'site_tau_f', 'site_tau_p' )

format_pars <- function(x, mod_lab ){
  
  x %>% 
    setNames( c( 'Survival',  'Growth',
                 'Flowering', 'Panicles' ) ) %>% 
    gather( vr, value, Survival:Panicles ) %>%  
    group_by( vr ) %>% 
    summarise( x = dens_extr(value, 1),   
               y = dens_extr(value, 2) ) %>% 
    ungroup %>% 
    mutate( Model = mod_lab ) 
    
}

long_df   <- dplyr::select(par_long, par_v) %>% format_pars( 'Longitude' )
noabio_df <- dplyr::select(par_noabio, par_v) %>% format_pars( 'No abiotic' ) 
clim_df   <- dplyr::select(par_clim, par_v) %>% format_pars( 'Precipitation' ) 
all_df    <- list( long_df, noabio_df, clim_df ) %>% bind_rows
  
# 
ggplot(all_df) +
  geom_line( aes(x,y, color = Model),
             lwd = 2) +
  facet_wrap( ~vr, scale = 'free' ) +
  scale_color_colorblind() +
  theme_minimal() +
  labs( x = 'Variance of site random effect',
        y = 'Density' ) +
  theme( axis.text    = element_text( size = 10),
         axis.title   = element_text( size = 15),
         legend.text  = element_text( size = 15),
         legend.title = element_text( size = 15),
         strip.text   = element_text( size = 15) ) +
  ggsave( 'C:/CODE/POAR-range-limits/results/Site_variance_by_model.tiff',
          width = 6.3, height = 6.3, compression = 'lzw')


dplyr::select(par_long, par_v) %>% sapply(mean) 
dplyr::select(par_noabio, par_v) %>% sapply(mean)   
dplyr::select(par_clim, par_v) %>% sapply(mean) 


# Read and format data -------------------------------------------------

in_dir  <- "C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/"

#common garden data
# poar_allsites <- read.csv("https://www.dropbox.com/s/bmg6ieiaoz8upri/demography_allsites.csv?dl=1", stringsAsFactors = F)
poar_allsites <- read.csv( paste0(in_dir, "Experiment/Demography/POAR-range-limits/data/demography_allsites.csv"), 
                           stringsAsFactors = F)

# Survival
poar_allsites.surv <- poar_allsites %>% 
  subset(tillerN_t0>0 ) %>%
  select(year, Code, site, unique.block, Sex, 
         long.center, long.scaled, 
         ppt_mismatch_scaled,
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

calc_resid <- function( x, y ){
  
  x2 <- dplyr::select(y,predS.1:predS.2170) %>% 
          sapply(mean) %>% 
          boot::inv.logit()
        
  r <- cbind( x, x2 ) %>% 
        as.data.frame() %>% 
        setNames( c('raw', 'pred') ) %>% 
        mutate( resid = raw - pred ) %>% 
        .$resid

  mean(r^2) %>% sqrt
  
}

calc_resid(poar_allsites.surv$surv_t1, par_long)
calc_resid(poar_allsites.surv$surv_t1, par_clim)
calc_resid(poar_allsites.surv$surv_t1, par_noabio)


# growth 
poar_allsites.grow <- poar_allsites %>% 
  subset( tillerN_t0 > 0 & tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          ppt_mismatch_scaled,
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

calc_resid <- function( x, y ){
  
  x2 <- dplyr::select(y, predG.1:predG.1472) %>% 
          sapply(mean) 
  
  r <- cbind( x, x2 ) %>% 
        as.data.frame() %>% 
        setNames( c('raw', 'pred') ) %>% 
        mutate( resid = raw - pred ) %>% 
        .$resid
  
  mean(r^2) %>% sqrt
  
}

calc_resid(poar_allsites.grow$tillerN_t1, par_long)
calc_resid(poar_allsites.grow$tillerN_t1, par_clim)
calc_resid(poar_allsites.grow$tillerN_t1, par_noabio)


# flowering
poar_allsites.flow <- poar_allsites %>% 
  subset( tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          ppt_mismatch_scaled,
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

calc_resid <- function( x, y ){
  
  x2 <- dplyr::select(y,predF.1:predF.1433) %>% 
    sapply(mean) %>% 
    boot::inv.logit()
  
  r <- cbind( x, x2 ) %>% 
    as.data.frame() %>% 
    setNames( c('raw', 'pred') ) %>% 
    mutate( resid = raw - pred ) %>% 
    .$resid
  
  mean(r^2) %>% sqrt
  
}

calc_resid(poar_allsites.flow$flow_t1, par_long)
calc_resid(poar_allsites.flow$flow_t1, par_clim)
calc_resid(poar_allsites.flow$flow_t1, par_noabio)


# Panicles 
poar_allsites.panic<- poar_allsites %>% 
  subset( flowerN_t1 > 0 & tillerN_t1 > 0 ) %>%
  select( year, Code, site, unique.block, Sex, 
          long.center, long.scaled, 
          ppt_mismatch_scaled,
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


calc_resid <- function( x, y ){
  
  x2 <- dplyr::select(y,predP.1:predP.512) %>% 
    sapply(mean) %>% 
    exp
  
  r <- cbind( x, x2 ) %>% 
    as.data.frame() %>% 
    setNames( c('raw', 'pred') ) %>% 
    mutate( resid = raw - pred ) %>% 
    .$resid
  
  mean(r^2) %>% sqrt
  
}

calc_resid(poar_allsites.panic$panic_t1, par_long)
calc_resid(poar_allsites.panic$panic_t1, par_clim)
calc_resid(poar_allsites.panic$panic_t1, par_noabio)
