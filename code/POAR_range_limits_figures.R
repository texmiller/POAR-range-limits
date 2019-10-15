library(tidyverse)
source('all_vr_bayes.R')


# Posterior predictive checks ---------------------------------------------
ppc_dens_overlay(data_all$y_s, y_s_sim)
ppc_dens_overlay(data_all$y_g, y_g_sim)+xlim(0, 100)
ppc_dens_overlay(data_all$y_f, y_f_sim)
ppc_dens_overlay(data_all$y_p, y_p_sim)+xlim(0, 100)
ppc_dens_overlay(data_all$y_v, y_v_sim) ## maybe need beta-binomial?
ppc_dens_overlay(data_all$y_m, y_m_sim) ## maybe need beta-binomial?


# Core vital rates --------------------------------------------------------
# estimate sex- and longitude-specific vital rates for small, medium, and large plants


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

## calculate SPEI
ozona_temp <- poar_sites %>% filter(ID1 == "ozona") %>% 
  select(ID1,Year,Latitude,Tave01:Tave12)%>% 
  gather(Tave01:Tave12,key="Month",value="Tave") %>% 
  mutate(Month = as.integer(str_sub(Month,5,6)))%>% 
  arrange(Year) %>% 
  mutate(transition_year = ifelse(Month>=5,Year,Year-1)) %>% 
  filter(transition_year %in% 2014:2017)
  
ozona_prcp <- poar_sites %>% filter(ID1 == "ozona") %>% 
  select(ID1,Year,PPT01:PPT12)%>% 
  gather(PPT01:PPT12,key="Month",value="PPT") %>% 
  mutate(Month = as.integer(str_sub(Month,4,5)))%>% 
  arrange(Year) %>% 
  mutate(transition_year = ifelse(Month>=5,Year,Year-1)) %>% 
  filter(transition_year %in% 2014:2017)

ozona <- full_join(ozona_temp,ozona_prcp,by=c("ID1","Month", "Year","transition_year")) %>% 
  mutate(PET = thornthwaite(Tave,unique(Latitude)),
         BAL = PPT-PET)

spei12 <- spei(ts(ozona[,'BAL']),12)

