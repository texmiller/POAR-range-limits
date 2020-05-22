## Purpose: compare predicted and observed size distributions, simulate range limits model with higher growth,
## thus larger size distributions, and hence more opportunity for sex ratio skew
## Are we underestimating the role of males simply because the parameterized MPM does not allow plants to get as big as they do]
## in real populations?
library(dplyr)
library(mgcv)
library(scales)
in_dir <- "C:/Users/tm9/Dropbox/POAR--Aldo&Tom/"

# Size distribution comparison --------------------------------------------
## first generate a figure comparing MPM SSD, common garden size distribution, and natural population size distribution

#1.Natural populations
files_12  <- list.files( paste0(in_dir, 'Data/2012/csv') )
files_13  <- list.files( paste0(in_dir, 'Data/2013/csv') )

# read polygon files
read_poly <- function(x, year){
  read.csv( paste0(in_dir,"Data/",year,"/csv/",x),
            stringsAsFactors = F ) %>% 
    select(Sex, GNSS_Area) %>% 
    mutate( Population = strsplit(x, "_Pol.csv")[[1]][1] )
}

poly_12 <- lapply(files_12, read_poly, 2012) %>% bind_rows
poly_13 <- lapply(files_13, read_poly, 2013) %>% bind_rows

# put polygon data together, and format it!
poly_df <- bind_rows(poly_12, poly_13) %>% 
  rename( area = GNSS_Area ) %>% 
  mutate( id       = as.character( 1:nrow(.) ),
          log_area = log(area),
          type     = 'Observational' ) %>% 
  select( -Population )

## here is a histogram of area pooling all individuals and populations, but what we really want is tiller number (conversion is below)
hist(poly_df$area)

#2.Common gardens
# format experimental data
exp_df <- read.csv( paste0(in_dir,'Range limits/Experiment/Demography/POAR-range-limits/data/demography.csv'),
                    stringsAsFactors = F )
exp_df %>% 
  ggplot()+
  geom_histogram(aes(x=log(tillerN_t1))) +
  facet_grid(.~year)

## common garden data include both tiller counts and area measurements in at least one year, right? check
exp_df %>% 
  select(year,MaxWidth_t1) %>% 
  group_by(year) %>% 
  summarise(sum(!is.na(MaxWidth_t1))) ## we actually did it in both 2016 and 2017

tiller_area_df <- exp_df %>% 
  # area of ellipse is pi*a*b where a is length/2 and b is width/2
  mutate(area_t1 = ifelse(MaxWidth_t1==0 | MaxLength_t1==0,NA,(pi*(MaxWidth_t1/2)*(MaxLength_t1/2))/10000),
         log_area_t1 = log(area_t1)) %>% 
  select(tillerN_t1,area_t1,log_area_t1) %>% 
  arrange(log_area_t1) %>% 
  na.omit

hist(tiller_area_df$area_t1)
plot(tiller_area_df$log_area_t1,log(tiller_area_df$tillerN_t1))

tiller_area_model <- gam(log(tillerN_t1) ~ s(log_area_t1), data=tiller_area_df)
summary(tiller_area_model)
lines(tiller_area_df$log_area_t1,predict(tiller_area_model),lwd=4)

## use fitted gam to predict tiller numbers for natural populations
poly_df$logtillers_pred <- predict.gam(tiller_area_model,data.frame(log_area_t1 = poly_df$log_area))


#3.MPM SSD
## this reproduces the mean parameter analysis in POAR_range_limits_analysis, and also does some more stuff
source("code/twosexMPM.R")
# read in MCMC output and data files
fit_full <- readRDS(paste0(in_dir,"Range limits/Experiment/Demography/POAR-range-limits/results/fit_allsite_full_pig_trunc.rds"))
POAU <- read.csv(paste0(in_dir,"Range limits//Experiment/Demography/POAR-range-limits/data/POAU.csv"))
sdlg_surv <- POAU %>% filter(year_recruit==year_t,year_t %in% 2014:2016) %>% summarise(sdlg_surv = mean(spring_survival_t1,na.rm=T))
latlong <- read.csv(paste0(in_dir,"Range limits/Experiment/Demography/data/SiteLatLong.csv"))
poar <- read.csv(paste0(in_dir,"Range limits//Experiment/Demography/POAR-range-limits/data/demography_allsites.csv"), stringsAsFactors = F)

# pull out stan coefficients
mean_coef <- lapply(rstan::extract(fit_full, pars = quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                                               bsizesex_s, bsizelong_s,blongsex_s,bsizelongsex_s,
                                                               blong2_s,bsizelong2_s,blong2sex_s,bsizelong2sex_s,
                                                               site_tau_s,block_tau_s,source_tau_s,
                                                               
                                                               b0_g,bsize_g,bsex_g,blong_g,
                                                               bsizesex_g, bsizelong_g,blongsex_g,bsizelongsex_g,
                                                               blong2_g,bsizelong2_g,blong2sex_g,bsizelong2sex_g,
                                                               site_tau_g,block_tau_g,source_tau_g,
                                                               sigma,
                                                               
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
F_params$sigma_g <- mean_coef$sigma 
M_params$grow_mu <- mean_coef$b0_g + mean_coef$bsex_g #+14
M_params$grow_size <- mean_coef$bsize_g + mean_coef$bsizesex_g 
M_params$grow_long <- mean_coef$blong_g + mean_coef$blongsex_g 
M_params$grow_size_long <- mean_coef$bsizelong_g + mean_coef$bsizelongsex_g 
M_params$grow_long2 <- mean_coef$blong2_g + mean_coef$blong2sex_g #0#
M_params$grow_size_long2 <- mean_coef$bsizelong2_g + mean_coef$bsizelong2sex_g #0#
M_params$sigma_g <- mean_coef$sigma 
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
F_params$max_size <- M_params$max_size <- round(quantile(na.omit(poar$tillerN_t1),probs=0.99)) #max(na.omit(poar$tillerN_t0)); 

## do these SSD simulations for middle of the range (longitude 0)
growth_factor <- c(1,3) ## <- going with 3x increase in grow_mu
new_mat_size <- round(F_params$max_size*1.5) ## <- increasing max size by 50%
## FRIDAY: use a smaller matrix (less than double standard max) but with 2x growth
ssd<-matrix(NA,new_mat_size,length(growth_factor))
max_yrs <- 20

for(g in 1:length(growth_factor)){
  F_params_new <- F_params;   M_params_new <- M_params
  # adjust growth size slope (so their increase with respect to current size is greater)
  F_params_new$grow_mu <- F_params$grow_mu*growth_factor[g]
  M_params_new$grow_mu <- M_params$grow_mu*growth_factor[g]
  #F_params_new$grow_size <- F_params$grow_size*growth_factor[g]
  #M_params_new$grow_size <- M_params$grow_size*growth_factor[g]
  F_params_new$max_size <- M_params_new$max_size <- new_mat_size
  
  print(g/length(growth_factor))
  lambda_run <- lambdaSim_delay( F_params=F_params_new,
                                 M_params=M_params_new,
                                 long = 0,
                                 rfx = rfx_fun(),
                                 max.yrs = max_yrs)
  ## drop "recruit" stage and combine females and males
  MandF <- lambda_run$n0[2:(F_params_new$max_size+1)] + lambda_run$n0[(F_params_new$max_size+3):((F_params_new$max_size+1)*2)] 
  ssd[,g] <- MandF/sum(MandF)
}

## visualize natural and common garden size distributions and MPM SSD
pdf("Manuscript/Figures/size_dist.pdf",useDingbats = F,height=9,width=5)
par(mfrow=c(3,1),mar=c(5,4,2,1))
# nat pops
hist(poly_df$logtillers_pred,xlim=c(0,10),main=" ",freq = F,col=alpha("blue",0.5),xlab="log(Tillers)",breaks=20)
title("A - Natural populations",adj=0)
# common garden
hist(log(exp_df$tillerN_t1),xlim=c(0,10),main=" ",col=alpha("red",0.5),freq = F,xlab="log(Tillers)",breaks=20)
title("B - Common garden",adj=0)
# MPM
hist(log(sample(1:new_mat_size,size=10000,replace=T,prob=ssd[,1])),xlim=c(0,10),main=" ",freq = F,col=alpha("black",0.5),xlab="log(Tillers)")
hist(log(sample(1:new_mat_size,size=10000,replace=T,prob=ssd[,2])),xlim=c(0,10),main=" ",freq = F,col=alpha("white",0.5),add=T)
legend("topright",legend=c("Field model","Elevated growth"),fill=c(alpha("black",0.5),alpha("white",0.5)),bty="n",cex=1.5)
title("C - MPM stable size distribution",adj=0)
dev.off()

## OK, the simulation experiment will have grow.mu increased by a factor of 2.5 -- now simulate across longitudes
long_seq_extend <- seq((-104 - mean(latlong$Longitude)),(-94.5 - mean(latlong$Longitude)),length.out = 20)
lambda_long_gsim<-lambda_long_gsim_2sex<-SR_long_gsim<-OSR_long_gsim<-matrix(NA,nrow=length(growth_factor),ncol=length(long_seq_extend))

for(g in 1:length(growth_factor)){
  F_params_new <- F_params;   M_params_new <- M_params
  # adjust growth size slope (so their increase with respect to current size is greater)
  F_params_new$grow_mu <- F_params$grow_mu*growth_factor[g]
  M_params_new$grow_mu <- M_params$grow_mu*growth_factor[g]
  #F_params_new$grow_size <- F_params$grow_size*growth_factor[g]
  #M_params_new$grow_size <- M_params$grow_size*growth_factor[g]
  F_params_new$max_size <- M_params_new$max_size <- new_mat_size
  
for(l in 1:length(long_seq_extend)){
  print(c(g,l))
  #linear model for comparison
  mat <- megamatrix_delay( F_params_new,
                           M_params_new,
                           twosex=F,
                           long=long_seq_extend[l],
                           rfx=rfx_fun())$MEGAmat
  lambda_long_gsim[g,l] <- lambda(mat)
  
  #2-sex model
  lambda_run <- lambdaSim_delay( F_params=F_params_new,
                                 M_params=M_params_new,
                                 long = long_seq_extend[l],
                                 rfx = rfx_fun(),
                                 max.yrs = max_yrs)
  lambda_long_gsim_2sex[g,l] <- lambda_run$lambdatracker[max_yrs]
  SR_long_gsim[g,l] <- lambda_run$SRtracker[max_yrs]
  OSR_long_gsim[g,l] <- lambda_run$OSRtracker[max_yrs]
}
}

par(mfrow=c(3,1),mar=c(5,4,1,1))

plot(long_seq_extend+mean(latlong$Longitude),OSR_long_gsim[1,],type="n",ylim=c(0,1),
     xlab="Longitude",ylab="Operational sex ratio (fraction female panicles)")
lines(long_seq_extend+mean(latlong$Longitude),OSR_long_gsim[1,],lwd=4)
lines(long_seq_extend+mean(latlong$Longitude),OSR_long_gsim[2,],lwd=4,lty=2)
abline(h=0.5,lty=3,col="lightgray")

plot(long_seq_extend+mean(latlong$Longitude),viab(params=F_params,twosex=T,OSR=OSR_long_gsim[1,]),
     ylim=c(0,1),lwd=4,type="l",xlab="Longitude",ylab="Predicted seed viability")
lines(long_seq_extend+mean(latlong$Longitude),viab(params=F_params,twosex=T,OSR=OSR_long_gsim[2,]),lwd=4,lty=2)

## Now do LTRE on fake parameter set
# these are the indices of the intercepts and slopes
F_params_new <- F_params;   M_params_new <- M_params
F_params_new$grow_mu <- F_params$grow_mu*growth_factor[2]
M_params_new$grow_mu <- M_params$grow_mu*growth_factor[2]
F_params_new$max_size <- M_params_new$max_size <- new_mat_size

LTRE_params <- c(1,2,7,8,14,15,20,21)
# these are the effects of longitude on the intercepts and slopes
betas_long <- c(F_params_new[c("surv_long","surv_size_long","grow_long","grow_size_long","flow_long","flow_size_long","panic_long","panic_size_long")],
                M_params_new[c("surv_long","surv_size_long","grow_long","grow_size_long","flow_long","flow_size_long","panic_long","panic_size_long")])
# these are the effects of longitude^2 on the intercepts and slopes
betas_long2 <- c(F_params_new[c("surv_long2","surv_size_long2","grow_long2","grow_size_long2","flow_long2","flow_size_long2","panic_long2","panic_size_long2")],
                 M_params_new[c("surv_long2","surv_size_long2","grow_long2","grow_size_long2","flow_long2","flow_size_long2","panic_long2","panic_size_long2")])

## Sensitivity of parameters to longitude -- this is just the derivative of a second-order polynomial
dp_dlong_fun <- function(beta_long,beta_long2,long){return(beta_long + 2*beta_long2*long)}

## loop over longitudes and calculated dp_dlong and dlambda_dp
dp_dlong <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(long_seq_extend))
lambda_long <- lambda_long_Fdom <- c()
perturbation <- 0.01
for(l in 1:length(long_seq_extend)){
  lambda_long_Fdom[l] <- lambda(megamatrix_delay(F_params_new,M_params_new,twosex=F,long=long_seq_extend[l],rfx=rfx_fun())$MEGAmat)
  lambda_long[l] <- lambdaSim_delay(F_params=F_params_new,M_params=M_params_new,
                                    long=long_seq_extend[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dlong[,l] <- dp_dlong_fun(beta_long = unlist(betas_long), 
                               beta_long2 = unlist(betas_long2), 
                               long = long_seq_extend[l])
  #cannot vectorize this part unfortunately
  for(p in 1:8){
    F_params_perturb <- F_params_new; M_params_perturb <- M_params_new;  
    F_params_perturb[LTRE_params[p]] <- unlist(F_params_new[LTRE_params[p]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,long=long_seq_extend[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_long[l]) / perturbation
  }
  for(p in 9:16){
    F_params_perturb <- F_params_new; M_params_perturb <- M_params_new;  
    M_params_perturb[LTRE_params[p-8]] <- unlist(M_params_new[LTRE_params[p-8]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,long=long_seq_extend[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_long[l]) / perturbation
  }
  print(l)
}

## write LTRE output based on posterior mean parameters
write_rds(list(dp_dlong=dp_dlong,dlambda_dp=dlambda_dp),paste0(dir,"/Experiment/Demography/POAR-range-limits/results/POAR_LTRE_growth_experiment.rds"))
POAR_LTRE_growth_experiment <- read_rds(paste0(dir,"/Experiment/Demography/POAR-range-limits/results/POAR_LTRE_growth_experiment.rds"))

# put it all together 
LTRE_out <- POAR_LTRE_growth_experiment$dp_dlong * POAR_LTRE_growth_experiment$dlambda_dp

par(mfrow=c(2,1))

plot(long_seq_extend + mean(latlong$Longitude),colSums(LTRE_out[1:8,]),type="l",lwd=4,#ylim=c(-0.4,0.5),
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"Longitude")),cex.lab=1.5)
abline(h=0,col="gray")
for(i in seq(1,8,by=2)){
  lines(long_seq_extend + mean(latlong$Longitude),colSums(LTRE_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}

plot(long_seq_extend + mean(latlong$Longitude),colSums(LTRE_out[9:16,]),type="l",lwd=4,#ylim=c(-0.4,0.5),
     xlab="Longitude",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"Longitude")),cex.lab=1.5)
abline(h=0,col="gray")
for(i in seq(9,16,by=2)){
  lines(long_seq_extend + mean(latlong$Longitude),colSums(LTRE_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
legend("topright",bty="n",legend=c("Survival","Growth","Flowering","Panicles","Total"),lwd=c(2,2,2,2,4),
       lty=c(na.omit(ltre_lty),1),col=c(na.omit(ltre_cols),"black"),cex=1.5)
