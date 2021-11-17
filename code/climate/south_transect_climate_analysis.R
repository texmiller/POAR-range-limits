## Purpose quantify how well we can connect longitude to climate using the south transect sites
## This is for an NSF proposal (proof of concept that long. can be converted to climate effects)
## Author: Tom Miller
## Last update: 26 Oct 2021

library(tidyverse)
library(lme4)

# Read demography data -------------------------------------------------
poar_Stransect <- read.csv("https://www.dropbox.com/s/xk4225mn8btqhbm/demography_allsites.csv?dl=1", stringsAsFactors = F) %>% 
  ## S transect sites in 2016 and 2017 (2015 had less than one transition year)
  subset(site%in%c("shsu","Katy","lostP","llano","ozona")) %>% 
  subset(year%in%2016:2017)

## write out lat-longs for climate data
poar_Stransect %>% select(site,Latitude,Longitude) %>% unique() %>% 
  mutate(ID2=".",
         el=".") %>% 
  rename(ID1=site,lat=Latitude,long=Longitude) %>% 
  select(ID1,ID2,lat,long,el) -> poar_Stransect_latlong
#write.csv(poar_Stransect_latlong,"code/climate/Stransect_lat_long.csv",row.names=F)

#read in monthly climate variables from climateNA
clim <- read.csv("code/climate/Stransect_lat_long_2015-2017MP.csv") 

clim %>% 
  select(Year,ID1,Latitude,Longitude,PPT01:PPT12) %>% 
  subset(Year==2015) %>% 
  pivot_longer(PPT01:PPT12,names_to="month",values_to="PPT") -> clim15
clim %>% 
  select(Year,ID1,Latitude,Longitude,PPT01:PPT12) %>% 
  subset(Year==2016) %>% 
  pivot_longer(PPT01:PPT12,names_to="month",values_to="PPT") -> clim16
clim %>% 
  select(Year,ID1,Latitude,Longitude,PPT01:PPT12) %>% 
  subset(Year==2017) %>% 
  pivot_longer(PPT01:PPT12,names_to="month",values_to="PPT") -> clim17
## bind years and calculate growing season precip by year
bind_rows(clim15,clim16,clim17) %>% 
  mutate(site=ID1,
         Month=as.integer(substr(month,4,5)),
         clim_year=ifelse(Month>=5,Year+1,Year),
         growing_season=Month%in%c(2,3,4)) %>% 
         #growing_season=Month%in%c(10,11,12,1,2,3,4)) %>% 
  filter(clim_year%in%2016:2017,growing_season==T)%>% 
  group_by(site,clim_year) %>% 
  summarise(precip=sum(PPT))-> poar_Stransect_clim

## merge climate and demography data
## I double checked aldo's data formatting and the "year" variable is year t1 of the t0-t1 transition
demog_clim <- left_join(poar_Stransect,poar_Stransect_clim,by=c("site","year"="clim_year"))
plot(demog_clim$Longitude,demog_clim$precip)

## Fit mixed models
grow_dat<-demog_clim %>% 
  select(log_ratio,tillerN_t0,Sex,site,unique.block,Code,precip) %>% 
  drop_na()

grow_site <- lmer(log_ratio~log(tillerN_t0)+Sex+(1|site)+(1|unique.block)+(1|Code),data=grow_dat)
grow_dat$resid<-residuals(grow_site)
plot(grow_dat$precip,grow_dat$resid)


grow_clim <- lmer(log_ratio~log(tillerN_t0)+Sex+precip+(1|site)+(1|unique.block)+(1|Code),data=grow_dat)

plot(demog_clim$Longitude,demog_clim$log_ratio)
plot(demog_clim$precip,demog_clim$log_ratio)
cor.test(demog_clim$precip,demog_clim$log_ratio)

## flowering
flow <- glmer(flowerN_t1>0~log(tillerN_t1)+(1|site)+(1|unique.block)+(1|Code),family="binomial",data=demog_clim)
flow_precip <- glmer(flow_t1~log(tillerN_t1)+precip+(1|site),family="binomial",data=demog_clim)
flow_precip2 <- glmer(flow_t1~log(tillerN_t0)+I(precip)^2+(1|site),family="binomial",data=demog_clim)

demog_clim %>% 
  group_by(site,year,Latitude,Longitude,precip) %>% 
  summarise(flow=mean(flow_t1,na.rm=T))->flow_dat
plot(flow_dat$Longitude,flow_dat$flow)
plot(flow_dat$precip,flow_dat$flow)
