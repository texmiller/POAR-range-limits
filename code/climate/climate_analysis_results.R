## Author: Tom Miller
## Purpose: fit statistical models to vital rate data using rstan,
## simulate data from fitted models for posterior predictive checks,
## output stan fits for use in graphics and demographic model, generate
## vital rate figures.
## The stan model is not shown here but imported from another file (poar_full.stan).
## Last update: 5/12/2021

# load packages
library(rstan)
library(loo)
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
library(ggthemes)
in_dir  <- "C:/Users/ac22qawo/poar/climate_analysis/"
in_dir  <- "C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/climate_analyses"

# read in stan fits
# mod_mmatch   <- readRDS( paste0(in_dir, 'poar_climate-7786865-4_4_poar_full_climate_mismatch_nologlik.RDS') )
#mod_mmatch   <- readRDS( paste0(in_dir, 'poar_climate-7885228-5_5_poar_full_climate_mismatch_nologlik.RDS') )
mod_mmatch   <- readRDS(gzcon(url("https://www.dropbox.com/s/27ocb9b3vijzm6t/poar_climate-7885228-5_5_poar_full_climate_mismatch_nologlik.RDS?dl=1")))

#mod_long     <- readRDS( paste0(in_dir, 'poar_climate-7786866-1_1_poar_full_climate_nologlik.RDS') )
mod_long     <- readRDS(gzcon(url("https://www.dropbox.com/s/b1ykbcdvu7u8xtc/poar_climate-7786866-1_1_poar_full_climate_nologlik.RDS?dl=1")))
  
mod_noabio   <- readRDS( paste0(in_dir, 'poar_climate-7786869-3_3_poar_full_noabiotic.RDS') )
mod_noabio   <- readRDS(gzcon(url("https://www.dropbox.com/s/swfubamesi35iev/poar_climate-7786869-3_3_poar_full_noabiotic.RDS?dl=1")))
  
mod_clim     <- readRDS( paste0(in_dir, 'poar_climate-7786870-2_2_poar_full_climate_nologlik.RDS') )

# models with log likelihood for model selection
mod_clim_ll  <- readRDS( paste0(in_dir, 'poar_climate-7814747-2_2_poar_full_climate.RDS') )
mod_long_ll  <- readRDS( paste0(in_dir, 'poar_climate-7668184-1_1full_climate_mod.RDS') )
  
# extract parameters
par_mmatch   <- rstan::extract(mod_mmatch) %>% as.data.frame
par_long     <- rstan::extract(mod_long) %>% as.data.frame
par_noabio   <- rstan::extract(mod_noabio) %>% as.data.frame
par_clim     <- rstan::extract(mod_clim) %>% as.data.frame

## tom's parameter density plots
pdf("Manuscript/Figures/climate_mismatch.pdf",height = 8,width = 8,useDingbats = F)
par(mfrow=c(2,2))
surv_mmdens <- density(par_mmatch$bmmatch_s)
surv_mmCI <- quantile(par_mmatch$bmmatch_s,probs=c(0.05,0.95))
x1 <- surv_mmdens$x[which.min(abs(surv_mmdens$x-surv_mmCI[1]))]
x2 <- surv_mmdens$x[which.min(abs(surv_mmdens$x-surv_mmCI[2]))]
plot(surv_mmdens,lwd=2,main="A) Survival",xlab="Mismatch coefficient",cex.lab=1.4)
abline(v=0,lty=3)
rect(x1,0,x2,0.05*max(surv_mmdens$y),col="gray")

grow_mmdens <- density(par_mmatch$bmmatch_g)
grow_mmCI <- quantile(par_mmatch$bmmatch_g,probs=c(0.05,0.95))
x1 <- grow_mmdens$x[which.min(abs(grow_mmdens$x-grow_mmCI[1]))]
x2 <- grow_mmdens$x[which.min(abs(grow_mmdens$x-grow_mmCI[2]))]
plot(grow_mmdens,lwd=2,main="B) Growth",xlab="Mismatch coefficient",cex.lab=1.4);abline(v=0,lty=3)
rect(x1,0,x2,0.05*max(grow_mmdens$y),col="gray")

flow_mmdens <- density(par_mmatch$bmmatch_f)
flow_mmCI <- quantile(par_mmatch$bmmatch_f,probs=c(0.05,0.95))
x1 <- flow_mmdens$x[which.min(abs(flow_mmdens$x-flow_mmCI[1]))]
x2 <- flow_mmdens$x[which.min(abs(flow_mmdens$x-flow_mmCI[2]))]
plot(flow_mmdens,lwd=2,main="C) Flowering",xlab="Mismatch coefficient",cex.lab=1.4);abline(v=0,lty=3)
rect(x1,0,x2,0.05*max(flow_mmdens$y),col="gray")

pan_mmdens <- density(par_mmatch$bmmatch_p)
pan_mmCI <- quantile(par_mmatch$bmmatch_p,probs=c(0.05,0.95))
x1 <- pan_mmdens$x[which.min(abs(pan_mmdens$x-pan_mmCI[1]))]
x2 <- pan_mmdens$x[which.min(abs(pan_mmdens$x-pan_mmCI[2]))]
plot(pan_mmdens,lwd=2,main="D) Panicles",xlab="Mismatch coefficient",cex.lab=1.4);abline(v=0,lty=3)
rect(x1,0,x2,0.05*max(pan_mmdens$y),col="gray")
dev.off()

# Model selection: longitude vs. climate predictor -----------------------------

mod_names <- data.frame( mod_n = 1:2,
                         model = c('Longitude', 'Climate') ) 

# Do this formatting only once!
form_modsel <- function( x ){
  x %>% 
    as.data.frame %>% 
    tibble::add_column(model = row.names(.), .before = 1) %>% 
    mutate( model = gsub('model','',model) %>% as.numeric ) %>% 
    rename( mod_n = model ) %>%   
    left_join( mod_names )  
}

mod_fit   <- list( mod_long_ll, mod_clim_ll)
log_liks  <- lapply(mod_fit, extract_log_lik)
loo_l     <- lapply(log_liks, loo) %>% 
                setNames( c('loo_long', 'loo_clim') )  
loo_df    <- loo_compare( loo_l$loo_long, 
                          loo_l$loo_clim ) %>% 
                form_modsel %>% 
                dplyr::select( model, elpd_diff:se_looic)

waic_l    <- lapply(log_liks, waic) %>% 
                setNames( c('waic_long', 'waic_clim') )  
waic_df    <- loo_compare( waic_l$waic_long, 
                           waic_l$waic_clim ) %>% 
                form_modsel %>% 
                dplyr::select( model, elpd_diff:se_waic)

out_dir <- 'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/climate_analyses/'

write.csv( loo_df,  paste0(out_dir, 'looic_long_vs_climate.csv'), row.names = F)
write.csv( waic_df, paste0(out_dir, 'waic_long_vs_climate.csv'),  row.names = F)

# effect of precipitation mismatch ----------------------------

# parameter names
par_v       <- c( 'blong_s', 'bmmatch_s',
                  'blong_g', 'bmmatch_g',
                  'blong_f', 'bmmatch_f',
                  'blong_p', 'bmmatch_p' )

# extract kernel density estimation values (both xs and ys)
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
          Longitude_Survival:Mismatch_Panicles ) %>% 
  group_by( param ) %>% 
  summarise( x = dens_extr(value, 1),   
             y = dens_extr(value, 2) ) %>% 
  ungroup %>% 
  separate( param, c('Parameter','vr'), sep='_' )
  
# plot longitude vs. mismatch
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
  ggsave( 'C:/CODE/POAR-range-limits/results/Precipitation Mismatch.tiff',
          width = 6.3, height = 6.3, compression = 'lzw')


# effect of precipitation mismatch ----------------------------

# names of parameter of site-level random effect variance
par_v     <- c( 'site_tau_s', 'site_tau_g',
                'site_tau_f', 'site_tau_p' )

# Format a parameters data frame across models
format_pars <- function(x, mod_lab ){
  
  x %>% 
    setNames( c( 'Survival',  'Growth',
                 'Flowering', 'Panicles' ) ) %>% 
    gather( vr, value, Survival:Panicles ) %>%  
    group_by( vr ) %>% 
    # Calculate kernel density estimation
    summarise( x = dens_extr(value, 1),   
               y = dens_extr(value, 2) ) %>% 
    ungroup %>% 
    mutate( Model = mod_lab ) 
    
}

# format parameters, and put them all in one data frame
long_df   <- dplyr::select(par_long, par_v) %>% format_pars( 'Longitude' )
noabio_df <- dplyr::select(par_noabio, par_v) %>% format_pars( 'No abiotic' ) 
clim_df   <- dplyr::select(par_clim, par_v) %>% format_pars( 'Precipitation' ) 
all_df    <- list( long_df, noabio_df, clim_df ) %>% bind_rows
  
# Plot estimates of variance by model
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

# Calculate RMSE for survival data
calc_rmse_s <- function( x, y ){
  
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

# Store RMSE values for survival data
surv_rmse <- data.frame( rmse = c( 
                calc_rmse_s(poar_allsites.surv$surv_t1, par_long),
                calc_rmse_s(poar_allsites.surv$surv_t1, par_clim),
                calc_rmse_s(poar_allsites.surv$surv_t1, par_noabio) ) 
                        ) %>% 
                mutate( model    = c('Longitude',
                                     'Precipitation',
                                     'No Abiotic'),
                        response = 'Survival' )


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

# calculate rmse for growth data
calc_rmse_g <- function( x, y ){
  
  x2 <- dplyr::select(y, predG.1:predG.1472) %>% 
          sapply(mean) 
  
  r <- cbind( x, x2 ) %>% 
        as.data.frame() %>% 
        setNames( c('raw', 'pred') ) %>% 
        mutate( resid = raw - pred ) %>% 
        .$resid
  
  mean(r^2) %>% sqrt
  
}

# Store RMSE values for survival data
grow_rmse <- data.frame( rmse = c( 
  calc_rmse_g(poar_allsites.grow$tillerN_t1, par_long),
  calc_rmse_g(poar_allsites.grow$tillerN_t1, par_clim),
  calc_rmse_g(poar_allsites.grow$tillerN_t1, par_noabio) )
                        ) %>% 
  mutate( model    = c('Longitude',
                       'Precipitation',
                       'No Abiotic'),
          response = 'Growth' )


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

# calculate rmse for flowering data
calc_resid_f <- function( x, y ){
  
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

# Store RMSE values for survival data
flow_rmse <- data.frame( rmse = c( 
  calc_resid_f(poar_allsites.flow$flow_t1, par_long),
  calc_resid_f(poar_allsites.flow$flow_t1, par_clim),
  calc_resid_f(poar_allsites.flow$flow_t1, par_noabio) )
                        ) %>% 
  mutate( model    = c('Longitude',
                       'Precipitation',
                       'No Abiotic'),
          response = 'Flowering' )


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

# calculate rmse for panicle data
calc_resid_p <- function( x, y ){
  
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

# Store RMSE values for survival data
panic_rmse <- data.frame( rmse = c( 
  calc_resid_p(poar_allsites.panic$panic_t1, par_long),
  calc_resid_p(poar_allsites.panic$panic_t1, par_clim),
  calc_resid_p(poar_allsites.panic$panic_t1, par_noabio) )
                        ) %>% 
  mutate( model    = c('Longitude',
                       'Precipitation',
                       'No Abiotic'),
          response = 'Panicles' )

# store it all in one place!
list( surv_rmse, 
      grow_rmse, 
      flow_rmse, 
      panic_rmse ) %>% 
  bind_rows %>% 
  write.csv( 'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/POAR-range-limits/results/climate_analyses/rmse_climate_mods.csv',
             row.names = F )
