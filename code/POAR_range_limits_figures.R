## Load the Stan output for full vital rate model
source('code/all_vr_bayes.R')


# Posterior predictive checks ---------------------------------------------
ppc_dens_overlay(data_all$y_s, y_s_sim)
ppc_dens_overlay(data_all$y_g, y_g_sim)+xlim(0, 100)
ppc_dens_overlay(data_all$y_f, y_f_sim)
ppc_dens_overlay(data_all$y_p, y_p_sim)+xlim(0, 50)
ppc_dens_overlay(data_all$y_v, y_v_sim) ## maybe need beta-binomial?
ppc_dens_overlay(data_all$y_m, y_m_sim) ## maybe need beta-binomial?


# Core vital rates --------------------------------------------------------
# estimate sex- and longitude-specific vital rates for small, medium, and large plants

# first bin size groups
size_bin_num <- 3
sex_cols <- c("red","blue")
bin_shapes <- 15:17

## Survival
poar_surv_binned <- poar.surv %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_surv = mean(surv_t1),
            long = unique(long.center),
            bin_n = n())
surv_mean_sizes <- poar_surv_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

# pull out stan coefficients
surv_coef <- rstan::extract(fit_full, pars = quote_bare(b0_s,bsize_s,bsex_s,blong_s,
                                                        bsizesex_s, bsizelong_s,blongsex_s,
                                                        bsizelongsex_s))
long_seq <- seq(min(poar_surv_binned$long),max(poar_surv_binned$long),0.1)

win.graph()
par(mfrow=c(4,3))
with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    plot(long[size_bin==i],mean_surv[size_bin==i],type="n",ylim=c(0,1))
    for(s in 1:2){
      points(long[sex==s & size_bin==i],mean_surv[sex==s & size_bin==i],
             col=sex_cols[s],pch=16,cex=2)
      lines(long_seq,
            invlogit(mean(surv_coef$b0_s) + 
                       mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] +
                       mean(surv_coef$bsex_s) * (s-1) +
                       mean(surv_coef$blong_s) * long_seq +
                       mean(surv_coef$bsizelong_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * long_seq +
                       mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                       mean(surv_coef$blongsex_s) * long_seq * (s-1) +
                       mean(surv_coef$bsizelongsex_s)  * long_seq * (s-1) * surv_mean_sizes$size[surv_mean_sizes$sex==s & surv_mean_sizes$size_bin==i]
                       ),
            col=sex_cols[s],lwd=3)
    }
  }
})

## Growth
poar_grow_binned <- poar.grow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_grow = mean(tillerN_t1),
            long = unique(long.center),
            bin_n = n())
grow_mean_sizes <- poar_grow_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

# pull out stan coefficients
grow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_g,bsize_g,bsex_g,blong_g,
                                                        bsizesex_g, bsizelong_g,blongsex_g,
                                                        bsizelongsex_g))

with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    plot(long[size_bin==i],mean_grow[size_bin==i],type="n")
    for(s in 1:2){
      points(long[sex==s & size_bin==i],mean_grow[sex==s & size_bin==i],
             col=sex_cols[s],pch=16,cex=2)
      lines(long_seq,
            exp(mean(grow_coef$b0_g) + 
                       mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] +
                       mean(grow_coef$bsex_g) * (s-1) +
                       mean(grow_coef$blong_g) * long_seq +
                       mean(grow_coef$bsizelong_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * long_seq +
                       mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                       mean(grow_coef$blongsex_g) * long_seq * (s-1) +
                       mean(grow_coef$bsizelongsex_g)  * long_seq * (s-1) * grow_mean_sizes$size[grow_mean_sizes$sex==s & grow_mean_sizes$size_bin==i]
            ),
            col=sex_cols[s],lwd=3)
    }
  }
})

## Flowering
poar_flow_binned <- poar.flow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_flow = mean(flow_t1),
            long = unique(long.center),
            bin_n = n())
flow_mean_sizes <- poar_flow_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

# pull out stan coefficients
flow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_f,bsize_f,bsex_f,blong_f,
                                                        bsizesex_f, bsizelong_f,blongsex_f,
                                                        bsizelongsex_f))

with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    plot(long[size_bin==i],mean_flow[size_bin==i],type="n",ylim=c(0,1))
    for(s in 1:2){
      points(long[sex==s & size_bin==i],mean_flow[sex==s & size_bin==i],
             col=sex_cols[s],pch=16,cex=2)
      lines(long_seq,
            invlogit(mean(flow_coef$b0_f) + 
                       mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] +
                       mean(flow_coef$bsex_f) * (s-1) +
                       mean(flow_coef$blong_f) * long_seq +
                       mean(flow_coef$bsizelong_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * long_seq +
                       mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                       mean(flow_coef$blongsex_f) * long_seq * (s-1) +
                       mean(flow_coef$bsizelongsex_f)  * long_seq * (s-1) * flow_mean_sizes$size[flow_mean_sizes$sex==s & flow_mean_sizes$size_bin==i]
            ),
            col=sex_cols[s],lwd=3)
    }
  }
})


## Panicles
poar_panic_binned <- poar.panic %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_panic = mean(panic_t1),
            long = unique(long.center),
            bin_n = n())
panic_mean_sizes <- poar_panic_binned %>% group_by(sex,size_bin) %>% summarise(size = mean(mean_size))

# pull out stan coefficients
panic_coef <- rstan::extract(fit_full, pars = quote_bare(b0_p,bsize_p,bsex_p,blong_p,
                                                        bsizesex_p, bsizelong_p,blongsex_p,
                                                        bsizelongsex_p))

with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    plot(long[size_bin==i],mean_panic[size_bin==i],type="n")
    for(s in 1:2){
      points(long[sex==s & size_bin==i],mean_panic[sex==s & size_bin==i],
             col=sex_cols[s],pch=16,cex=2)
      lines(long_seq,
            exp(mean(panic_coef$b0_p) + 
                  mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] +
                  mean(panic_coef$bsex_p) * (s-1) +
                  mean(panic_coef$blong_p) * long_seq +
                  mean(panic_coef$bsizelong_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * long_seq +
                  mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                  mean(panic_coef$blongsex_p) * long_seq * (s-1) +
                  mean(panic_coef$bsizelongsex_p)  * long_seq * (s-1) * panic_mean_sizes$size[panic_mean_sizes$sex==s & panic_mean_sizes$size_bin==i]
            ),
            col=sex_cols[s],lwd=3)
    }
  }
})

# Seed viability ----------------------------------------------------------
viab_pars <- rstan::summary(fit_full, pars=c("v0","a_v"))[[1]][,"mean"]

plot(data_all$SR_v,(data_all$y_v / data_all$tot_seeds_v),
     cex = 3 * (data_all$tot_seeds_v / max(data_all$tot_seeds_v)))
lines(seq(0,1,0.01),
      viab_pars["v0"] * (1 - seq(0,1,0.01) ^ viab_pars["a_v"]))


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


plot(allsites$long,allsites$ppt,type="n",ylab="Annual precipitation (mm)",xlab="Longitude",cex.lab=1.4)
points(allsites$long[allsites$transition_year==2014],allsites$ppt[allsites$transition_year==2014],cex=2,pch=1)
points(allsites$long[allsites$transition_year==2015],allsites$ppt[allsites$transition_year==2015],cex=2,pch=2)
points(allsites$long[allsites$transition_year==2016],allsites$ppt[allsites$transition_year==2016],cex=2,pch=3)
legend("topleft",legend=c(2014:2016),bty="n",pch=1:3,cex=1.2)
## add 30-year normals