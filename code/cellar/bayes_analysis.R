## Fit size-dependent vital rates and build MPM
# library(R2jags)
# library(mcmcplots)
library(scales)
library(dplyr)
library(rstan)
library(shinystan)
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
viabVr  <- read.csv('data/viability.csv')



#1. Survival -------------------------------------------
# poar.surv <- poar %>% 
#               select( site, unique.block, Sex, 
#                       long.center, long.scaled,
#                       tillerN_t0, surv_t1 ) %>% 
#               # convert factors to numeric
#               mutate_at( c('site', 'unique.block', 'Sex'), 
#                          as.numeric) %>% 
#               na.omit %>% 
#               setNames( c("site","block", "sex",
#                           "long.center", "long.scaled",
#                           "tillerN_t0", "surv_t1") ) %>% 
#               # individuals need be alive at start of transition
#               subset( tillerN_t0>0 ) %>%
#               mutate( log_size_t0 = log(tillerN_t0) )
# 

#1. Survival
poar.surv<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t0,poar$surv_t1)))
names(poar.surv)<-c("site","block","sex","long.center","long.scaled","tillerN_t0","surv_t1")
poar.surv<-subset(poar.surv,tillerN_t0>0)
N.sites<-max(poar.surv$site)
site.surv<-poar.surv$site
N.blocks<-max(poar.surv$block)
block.surv<-poar.surv$block
male.surv<-poar.surv$sex-1
long.surv<-poar.surv$long.center
y.surv<-poar.surv$surv_t1
size.surv<-log(poar.surv$tillerN_t0)
poar.surv <- mutate(poar.surv, log_size_t0 = log(tillerN_t0) )
N.obs.surv<-nrow(poar.surv)


#2. Growth
poar.grow<-na.omit(data.frame(cbind(poar$site,
                                    poar$unique.block,
                                    poar$Sex,
                                    poar$long.center,
                                    poar$long.scaled,
                                    poar$tillerN_t0,
                                    poar$tillerN_t1)))
names(poar.grow)<-c("site","block","sex","long.center","long.scaled","tillerN_t0","tillerN_t1")
max(poar.surv$site)==max(poar.grow$site)
max(poar.surv$block)==max(poar.grow$block)## sites and blocks line up
poar.grow<-subset(poar.grow,tillerN_t0>0 & tillerN_t1>0)
site.grow<-poar.grow$site
block.grow<-poar.grow$block
male.grow<-poar.grow$sex-1
long.grow<-poar.grow$long.center
size.grow<-log(poar.grow$tillerN_t0)
y.grow<-poar.grow$tillerN_t1
N.obs.grow<-nrow(poar.grow)


# 3. Flowering
poar.flow<-na.omit(data.frame(cbind(poar$site,
                                    poar$unique.block,
                                    poar$Sex,
                                    poar$long.center,
                                    poar$long.scaled,
                                    poar$tillerN_t1,
                                    poar$flow_t1))) 
names(poar.flow)<-c("site","block","sex","long.center","long.scaled","tillerN_t1","flow_t1")
max(poar.surv$site)==max(poar.flow$site)
max(poar.surv$block)==max(poar.flow$block)## sites and blocks line up
poar.flow<-subset(poar.flow,tillerN_t1>0)
site.flow<-poar.flow$site
block.flow<-poar.flow$block
male.flow<-poar.flow$sex-1
long.flow<-poar.flow$long.center
size.flow<-log(poar.flow$tillerN_t1)
y.flow<-poar.flow$flow_t1
N.obs.flow<-nrow(poar.flow)
poar.flow <- mutate( poar.flow, log_size_t1 = log(tillerN_t1) )
# plot(size.flow,y.flow)


#4. Panicles
poar.panic<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t1,poar$flowerN_t1)))
names(poar.panic)<-c("site","block","sex","long.center","long.scaled","tillerN_t1","panic_t1")
poar.panic<-subset(poar.panic,panic_t1>0 & tillerN_t1>0)
max(poar.surv$site)==max(poar.panic$site)
max(poar.surv$block)==max(poar.panic$block)## sites and blocks line up
site.panic<-poar.panic$site
block.panic<-poar.panic$block
male.panic<-poar.panic$sex-1
long.panic<-poar.panic$long.center
size.panic<-log(poar.panic$tillerN_t1)
y.panic<-poar.panic$panic_t1
N.obs.panic<-nrow(poar.panic)


# plot all vital rates together!
par(mfrow=c(2,2), mar=c(3,3,0.1,0.1), mgp=c(1.6,0.6,0))
# plot(size.surv,y.surv)
# binned survival function (more intuitive than raw data)
plot_binned_prop(poar.surv, 10, log_size_t0, surv_t1)
plot(size.grow,y.grow)
# binned survival function (more intuitive than raw data)
plot_binned_prop(poar.flow, 10, log_size_t1, flow_t1)
plot(size.panic,y.panic)


## Misc data for prediction
size.pred<-seq(0,max(log(poar$tillerN_t1),na.rm=T),length.out = 30)
W.long<-min(poar$long.center)
E.long<-max(poar$long.center)

#5. Germination
poar.germ    <- na.omit(data.frame(cbind(viabVr$plot,viabVr$germTot,viabVr$germFail,viabVr$sr_f)))
poar.germ    <- setNames( poar.germ, c("plot","germTot","germFail","sr_f") )
y_germ       <- poar.germ$germTot
tot_seeds    <- poar.germ$germTot+poar.germ$germFail
SR           <- poar.germ$sr_f
plot_germ    <- poar.germ$plot
N.germ.plots <- max(poar.germ$plot)
N.obs.germ   <- nrow(poar.germ)



data_l <- list(N.sites = N.sites,
               N.blocks= N.blocks,
               
               site.surv  = site.surv,
               block.surv = block.surv,
               male.surv  = male.surv,
               long.surv  = long.surv,
               size_surv  = size.surv,
               y_surv     = y.surv,
               n_surv     = N.obs.surv,

               site.grow=site.grow,
               block.grow=block.grow,
               male.grow=male.grow,
               long.grow=long.grow,
               size.grow=size.grow,
               y.grow=y.grow,
               N.obs.grow=N.obs.grow,
               
               site.flow=site.flow,
               block.flow=block.flow,
               male.flow=male.flow,
               long.flow=long.flow,
               size.flow=size.flow,
               y.flow=y.flow,
               N.obs.flow=N.obs.flow,
               
               site.panic=site.panic,
               block.panic=block.panic,
               male.panic=male.panic,
               long.panic=long.panic,
               y.panic=y.panic,
               size.panic=size.panic,
               N.obs.panic=N.obs.panic,
               
               size.pred=size.pred,
               N.sizes=length(size.pred),
               W.long=W.long,
               E.long=E.long,
               
               y.germ=y.germ,
               tot.seeds=tot.seeds,
               SR=SR,
               plot.germ=plot.germ,
               N.germ.plots=N.germ.plots,
               N.obs.germ=N.obs.germ)


#1. Survival
poar      <- read.csv('data/demography.csv', stringsAsFactors = F)
poar.surv <- poar %>% 
                subset( tillerN_t0>0 ) %>%
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


N.sites     <- max(poar.surv$site)
site.surv   <- poar.surv$site
N.blocks    <- max(poar.surv$block)
block.surv  <- poar.surv$block
male.surv   <- poar.surv$sex-1
long.surv   <- poar.surv$long.center
y.surv      <- poar.surv$surv_t1
size.surv   <- log(poar.surv$tillerN_t0)
N.obs.surv  <- nrow(poar.surv)


data_l <- list( n_sites  = N.sites,
                n_blocks = block.surv %>% n_distinct,
                site_s   = site.surv,
                block_s  = block.surv,
                site_block_s = data.frame( site_i  = site.surv,
                                           block_i = block.surv ) %>% 
                                  unique %>% .$site_i,
                male_s   = male.surv,
                long_s   = long.surv,
                size_s   = size.surv,
                y_s      = y.surv,
                n_s      = N.obs.surv )

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# fit survival model
fit_surv <- stan(
  file = 'code/stan/surv_s.stan',
  data = data_l,
  pars = quote_bare( site_a_s, b_s, b_l, b_sex, b_sex_l ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = 1 #sim_pars$chains#,
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

shinystan::launch_shinystan(fit_surv)

# DOESN'T CONVERGE
# fit survival model
fit_surv_nc <- stan(
  file = 'code/stan/surv_sb_nest_nc.stan',
  data = data_l,
  pars = quote_bare( site_a_s, block_a_s, b_s, b_l, b_sex, b_sex_l ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


shinystan::launch_shinystan(fit_surv)

library(lme4)
mod <- glmer(surv_t1 ~ log_size_t0 + sex*long.center + (1 | block),
             data=poar.surv, family='binomial')


ranef(mod)$block[,1] %>% hist
plot(mod)


# create predicted values
coef  <- summary(fit_surv)$summary[,'mean']
x_seq <- seq( min(poar.surv$log_size_t0),
              max(poar.surv$log_size_t0),
              length.out = 100)
y_prd <- boot::inv.logit( coef['a'] + coef['b_s'] * x_seq )

# plot pred VS 
par(mfrow=c(2,2), mar=c(3,3,0.1,0.1), mgp=c(1.6,0.6,0))
plot_binned_prop(poar.surv, 10, log_size_t0, surv_t1)
lines(x_seq, y_prd)

# germination model ---------------------------------------

#5. Germination
poar.germ    <- na.omit(data.frame(cbind(viabVr$plot,viabVr$germTot,viabVr$germFail,viabVr$sr_f)))
poar.germ    <- setNames( poar.germ, c("plot","germTot","germFail","sr_f") )
y_germ       <- poar.germ$germTot
tot_seeds    <- poar.germ$germTot+poar.germ$germFail
SR           <- poar.germ$sr_f
plot_germ    <- poar.germ$plot
N.germ.plots <- max(poar.germ$plot)
N.obs.germ   <- nrow(poar.germ)


# cia
data_l <- list( n_sites  = N.sites,
                n_blocks = block.surv %>% n_distinct,
                site_s   = site.surv,
                block_s  = block.surv,
                site_block_s = data.frame( site_i  = site.surv,
                                           block_i = block.surv ) %>% 
                                  unique %>% .$site_i,
                male_s   = male.surv,
                long_s   = long.surv,
                size_s   = size.surv,
                y_s      = y.surv,
                n_s      = N.obs.surv,
                
                n_g      = poar.germ$germTot %>% length,
                y_germ   = poar.germ$germTot,
                tot_seeds= poar.germ$germTot+poar.germ$germFail,
                SR       = poar.germ$sr_f )


plot( data_l$SR, data_l$y_germ)

glm( cbind(data_l$y_germ,data_l$tot_seeds-data_l$y_germ) ~ 1, 
     family='binomial') %>% coef %>% boot::inv.logit()
data_l$y_germ


# fit survival model
fit_germ <- stan(
  file = 'code/stan/germ.stan',
  data = data_l,
  pars = quote_bare( v0, a_g ),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  # init = list( list( v0 = 0.5, a_g = 2),
  #              list( v0 = 0.5, a_g = 2) ),
  chains = 3
  # control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

y_germ ~ binomial(tot_seeds, prob_g)

dbinom(data_l$y_germ[115], 
       data_l$tot_seeds[115], 0.1, log=T)

data_l$y_germ[115]
data_l$tot_seeds[115]

fit_germ <- fit_surv
coef     <- summary(fit_germ)$summary['mean']


fun <- function(x){
  v0*(1 - x^alpha)
}
v0    <- 0.57 ## this should get a beta prior
alpha <- 9.38 ## this should get gamma or lognormal
x     <-seq(0,1,0.01)
plot(data_l$SR, data_l$y_germ/data_l$tot_seeds)
lines(x,fun(x),type="l",xlim=c(0,1),ylim=c(0,1))

## seed_viab ~ binomial(fun(x),N_seeds)

