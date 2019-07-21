## Author: Tom
## Purpose: Analysis of POAR range limits data
## last update: 8/18/2017
## Notes: Data formatting and manipulation in Demography/Data/f14_s17_demo_data.R
## We know there were some merge issues with f15 and s15 -- as of now, this issue is not corrected
library(lme4)
library(bbmle)
library(dplyr)
library(tidyr)
library(xlsx)

setwd("D:/Dropbox/poar--Aldo&Tom/Range limits/Experiment/Demography")
setwd("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography")

# useful functions ----------------------------------------------------
length2 <- function(x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}
invlogit<-function(x){exp(x)/(1+exp(x))}
# Symbols for the sexes 
sex_symbol <- function(x){
  out <- x %>% 
    mutate(symb = factor(Sex, labels=c("21","24")) ) %>%
    mutate(symb = as.character(symb) ) %>%
    mutate(symb = as.integer(symb) ) %>%
    mutate(col = factor(Sex, labels=c("blue","red")) ) %>%
    mutate(col = as.character(col) )
  return(out)
}


# read data -----------------------------------------------------------------
poar    <- read.csv("data/f14_s17_data.csv")
latlong <- read.csv("data/SiteLatLong.csv")

str(latlong)

hist(latlong$Longitude)


# format ---------------------------------------------------------------

## rescale longitude
latlong$long.scaled<-(latlong$Longitude-mean(latlong$Longitude))/sd(latlong$Longitude)
## or centering may be easier for converting back to real long
latlong$long.center<-latlong$Longitude-mean(latlong$Longitude)

# add lat/long to data frame
poar    <- merge(poar,latlong,by="site")
poar    <- poar %>% 
  # add a unique block ID to fit as random effect
  mutate(unique.block = interaction(site,Block) ) %>%
  mutate(flow_t1 = as.numeric(flowerN_t1>0) ) %>% # add flowering (bernoulli)
  mutate(pc_fn = flowerN_t1 / tillerN_t1) %>%
  mutate(pc_fn = replace(pc_fn, pc_fn == Inf | is.nan(pc_fn), NA) ) %>%
  mutate(surv_t1 = as.numeric(tillerN_t1>0) ) %>% # add survival (bernoulli)
  mutate(delta.tiller = tillerN_t1 - tillerN_t0) %>% # change in size
  mutate(log_ratio = log(tillerN_t1 / tillerN_t0) ) %>% # log ratio
  mutate(log_ratio = replace(log_ratio, log_ratio == -Inf | log_ratio == Inf | is.nan(log_ratio), NA) ) %>%
  ## if clone area is zero, make it arbitrarily small
  mutate(clone_area_t0 = replace(clone_area_t0, clone_area_t0==0, 1) ) %>%
  mutate(clone_area_t1 = replace(clone_area_t1, clone_area_t1==0, 1) ) %>%
  # flag ploidy of provenance
  mutate(ploidy = "low") %>%
  mutate(ploidy = replace(ploidy, Code == "HHC" | Code == "LLELA", "high") ) %>%
  sex_symbol()


# Drop the 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
poar  <- subset(poar, site != "llela" & site != "lubbock" & site != "wf" )
unique(poar$site) ##this works
#str(poar.drop)

hist(log(poar$clone_area_t1[poar$year==2017]))
hist(log(poar$tillerN_t1[poar$year==2017]))

## Visualize trends in raw data with respect to longitude

## overall change in size
plot(poar$tillerN_t0,poar$tillerN_t1,ylim=c(0,100))
abline(0,1)

## final size in area
plot(poar$Longitude[poar$year==2017 & poar$Sex=="M"],jitter(log(poar$clone_area_t1[poar$year==2017 & poar$Sex=="M"]),factor=50),
     xlab="Longitude",ylab="log final size (clone area)")
points(poar$Longitude[poar$year==2017 & poar$Sex=="F"]+0.1,jitter(log(poar$clone_area_t1[poar$year==2017 & poar$Sex=="F"]),factor=50),col="red")

## final size in tillers
plot(poar$Longitude[poar$year==2017 & poar$Sex=="M"],jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="M"]),factor=50),
     xlab="Longitude",ylab="log final size (# Tillers)")
points(poar$Longitude[poar$year==2017 & poar$Sex=="F"]+0.1,jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="F"]),factor=50),col="red")

plot(poar$Latitude[poar$year==2017 & poar$Sex=="M"],jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="M"]),factor=50),
     xlab="Longitude",ylab="log final size (# Tillers)")
points(poar$Latitude[poar$year==2017 & poar$Sex=="F"]+0.1,jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="F"]),factor=50),col="red")


grow.models<-list()
grow.models[[1]]<-lmer(log(tillerN_t1)~(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[2]]<-lmer(log(tillerN_t1)~Longitude+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[3]]<-lmer(log(tillerN_t1)~Sex+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[4]]<-lmer(log(tillerN_t1)~Sex+Longitude+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[5]]<-lmer(log(tillerN_t1)~Sex*Longitude+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
AICtab(grow.models, weights=T)

long.seq<-seq(min(poar$Longitude),max(poar$Longitude),0.1)
lines(long.seq,fixef(grow.models[[2]])[1]+fixef(grow.models[[2]])[2]*long.seq,col="black",lwd=3)

## same but rescaled
grow.models<-list()
grow.models[[1]]<-lmer(log(tillerN_t1)~(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[2]]<-lmer(log(tillerN_t1)~long.scaled+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[3]]<-lmer(log(tillerN_t1)~Sex+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[4]]<-lmer(log(tillerN_t1)~Sex+long.scaled+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
grow.models[[5]]<-lmer(log(tillerN_t1)~Sex*long.scaled+(1|site) + (1|unique.block),data=subset(poar,year==2017 & tillerN_t1>0))
AICtab(grow.models, weights=T)

plot(poar$long.scaled[poar$year==2017 & poar$Sex=="M"],jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="M"]),factor=50),
     xlab="Longitude scaled",ylab="log final size (# Tillers)")
points(poar$long.scaled[poar$year==2017 & poar$Sex=="F"]+0.1,jitter(log(poar$tillerN_t1[poar$year==2017 & poar$Sex=="F"]),factor=50),col="red")
long.seq<-seq(min(poar$long.scaled),max(poar$long.scaled),0.1)
lines(long.seq,fixef(grow.models[[2]])[1]+fixef(grow.models[[2]])[2]*long.seq,col="black",lwd=3)

## Now annual growth, pooling years, a log ratio
## annual growth rate (all years)
ann.grow.models<-list()
ann.grow.models[[1]]<-lmer(log_ratio~(1|site) + (1|unique.block),data=poar)
ann.grow.models[[2]]<-lmer(log_ratio~long.scaled+(1|site) + (1|unique.block),data=poar)
ann.grow.models[[3]]<-lmer(log_ratio~Sex+(1|site) + (1|unique.block),data=poar)
ann.grow.models[[4]]<-lmer(log_ratio~Sex+long.scaled+(1|site) + (1|unique.block),data=poar)
ann.grow.models[[5]]<-lmer(log_ratio~Sex*long.scaled+(1|site) + (1|unique.block),data=poar)
AICtab(ann.grow.models, weights=T)

plot(poar$long.scaled[poar$Sex=="M"],jitter(poar$log_ratio[poar$Sex=="M"],factor=50),
     xlab="Longitude",ylab="Annual growth (log ratio)")
points(poar$long.scaled[poar$Sex=="F"]+0.05,jitter(poar$log_ratio[poar$Sex=="F"],factor=50),col="red")
long.seq<-seq(min(poar$long.scaled),max(poar$long.scaled),0.1)
lines(long.seq,fixef(ann.grow.models[[2]])[1]+fixef(ann.grow.models[[2]])[2]*long.seq,col="black",lwd=3)
abline(h=0,col="gray")

plot(poar$long.scaled,jitter(poar$log_ratio,factor=50),
     xlab="Longitude",ylab="Annual growth (log ratio)")
long.seq<-seq(min(poar$long.scaled),max(poar$long.scaled),0.1)
lines(long.seq,fixef(ann.grow.models[[2]])[1]+fixef(ann.grow.models[[2]])[2]*long.seq,col="black",lwd=3)
abline(h=0,col="gray")

abline(fixef(ann.grow.models[[5]])[1],
       fixef(ann.grow.models[[5]])[3],col="red")
abline((fixef(ann.grow.models[[5]])[1]+fixef(ann.grow.models[[5]])[2]),
       (fixef(ann.grow.models[[5]])[3]+fixef(ann.grow.models[[5]])[4]))

## are growth, survival, and reproduction size-dependent?
plot(log(poar$tillerN_t0),jitter(poar$surv_t1))
plot(log(poar$tillerN_t0),jitter(poar$log_ratio,factor=500))
plot(log(poar$tillerN_t1),jitter(poar$flow_t1))
plot(log(poar$tillerN_t1),jitter(poar$flowerN_t1))
## These will be great for MPM

# survival
# format data
surv <- poar %>% 
  select(site,Sex,year,surv_t1) %>%
  subset(!is.na(surv_t1)) %>%
  group_by(site,Sex,year) %>%
  summarise( site_surv = sum(surv_t1,na.rm=T),
             total_ind = 42) %>%  ### I replaced n() with 42 because I want to see total survival from the initial population size (42 planted per site)
  mutate(surv_pr = site_surv / total_ind) %>%
  as.data.frame() %>%
  inner_join(latlong) %>%
  sex_symbol()


## annual survival probability, pooling years
surv.ann.models<-list()
surv.ann.models[[1]]<-glmer(surv_t1~(1|site)+(1|unique.block),family="binomial",data=poar)
surv.ann.models[[2]]<-glmer(surv_t1~long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
surv.ann.models[[3]]<-glmer(surv_t1~Sex+(1|site)+(1|unique.block),family="binomial",data=poar)
surv.ann.models[[4]]<-glmer(surv_t1~Sex+long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
surv.ann.models[[5]]<-glmer(surv_t1~Sex*long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
AICtab(surv.ann.models, weights=T)

surv.cum.models<-list()
surv.cum.models[[1]]<-glmer(surv_pr~(1|site),weights=total_ind,family="binomial",data=subset(surv,year==2017))
surv.cum.models[[2]]<-glmer(surv_pr~long.scaled+(1|site),weights=total_ind,family="binomial",data=subset(surv,year==2017))
surv.cum.models[[3]]<-glmer(surv_pr~Sex+(1|site),weights=total_ind,family="binomial",data=subset(surv,year==2017))
surv.cum.models[[4]]<-glmer(surv_pr~Sex+long.scaled+(1|site),weights=total_ind,family="binomial",data=subset(surv,year==2017))
surv.cum.models[[5]]<-glmer(surv_pr~Sex*long.scaled+(1|site),weights=total_ind,family="binomial",data=subset(surv,year==2017))
AICtab(surv.cum.models, weights=T)

## final survival
plot(surv$Longitude[surv$year==2017],surv$surv_t1[surv$year==2017],pch=surv$symb[surv$year==2017],col=surv$col[surv$year==2017])
abline(h=invlogit(fixef(surv.cum.models[[3]])[1]),col="blue")
abline(h=invlogit(fixef(surv.cum.models[[3]])[1]+fixef(surv.cum.models[[3]])[2]),col="red")

## annual survival
plot(jitter(poar$long.scaled[poar$Sex=="M"],factor=5),jitter(poar$surv_t1[poar$Sex=="M"],factor=0.05),
     xlab="Longitude scaled",ylab="# Flowers")
points(jitter(poar$long.scaled[poar$Sex=="F"],factor=5),jitter(poar$surv_t1[poar$Sex=="F"],factor=0.05)-0.025,col="red")
abline(h=invlogit(fixef(surv.ann.models[[1]])))

### Probability of flowering

prflow.models<-list()
prflow.models[[1]]<-glmer(flowerN_t1>0~(1|site)+(1|unique.block),family="binomial",data=poar)
prflow.models[[2]]<-glmer(flowerN_t1>0~long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
prflow.models[[3]]<-glmer(flowerN_t1>0~Sex+(1|site)+(1|unique.block),family="binomial",data=poar)
prflow.models[[4]]<-glmer(flowerN_t1>0~Sex+long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
prflow.models[[5]]<-glmer(flowerN_t1>0~Sex*long.scaled+(1|site)+(1|unique.block),family="binomial",data=poar)
AICtab(prflow.models, weights=T)

plot(jitter(poar$long.scaled[poar$Sex=="M"],factor=5),jitter(as.numeric(poar$flowerN_t1[poar$Sex=="M"]>0),factor=0.05),
     xlab="Longitude scaled",ylab="Probability of flowering")
points(jitter(poar$long.scaled[poar$Sex=="F"],factor=5),jitter(as.numeric(poar$flowerN_t1[poar$Sex=="F"]>0),factor=0.05)-0.025,col="red")
lines(long.seq,invlogit(fixef(prflow.models[[4]])[1]+fixef(prflow.models[[4]])[3]*long.seq),col="red",lwd=3)
lines(long.seq,invlogit(fixef(prflow.models[[4]])[1]+fixef(prflow.models[[4]])[2]+fixef(prflow.models[[4]])[3]*long.seq),col="black",lwd=3)

##model 5 is a close second...see what it looks like
lines(long.seq,invlogit(fixef(prflow.models[[5]])[1]+fixef(prflow.models[[5]])[3]*long.seq),col="red",lwd=1,lty=2)
lines(long.seq,invlogit(fixef(prflow.models[[5]])[1]+fixef(prflow.models[[5]])[2]+(fixef(prflow.models[[5]])[3]+fixef(prflow.models[[5]])[4])*long.seq),col="black",lwd=1,lty=2)


##### Number of flowers
plot(poar$long.scaled[poar$Sex=="M" & poar$flowerN_t1>0],(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0]),
     xlab="Longitude scaled",ylab="# Flowers")
points(poar$long.scaled[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0]),col="red")

plot(poar$Latitude[poar$Sex=="M" & poar$flowerN_t1>0],(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0]),
     xlab="Longitude scaled",ylab="# Flowers",ylim=c(0,50))
points(poar$Latitude[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0]),col="red")



plot(poar$long.scaled[poar$Sex=="M" & poar$flowerN_t1>0],log(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0]),
     xlab="Longitude scaled",ylab="# Flowers")
points(poar$long.scaled[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,log(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0]),col="red")

##add individual random effect to approximate NegBinom
poar$unique.ID<-interaction(poar$site,poar$ID,poar$Code,poar$Sex)

Nflow.models<-list()
Nflow.models[[1]]<-glmer(flowerN_t1~(1|unique.ID)+(1|site)+(1|unique.block),family="poisson",data=subset(poar,flowerN_t1>0))
Nflow.models[[2]]<-glmer(flowerN_t1~long.scaled+(1|unique.ID)+(1|site)+(1|unique.block),family="poisson",data=subset(poar,flowerN_t1>0))
Nflow.models[[3]]<-glmer(flowerN_t1~Sex+(1|site)+(1|unique.ID)+(1|unique.block),family="poisson",data=subset(poar,flowerN_t1>0))
Nflow.models[[4]]<-glmer(flowerN_t1~Sex+long.scaled+(1|unique.ID)+(1|site)+(1|unique.block),family="poisson",data=subset(poar,flowerN_t1>0))
Nflow.models[[5]]<-glmer(flowerN_t1~Sex*long.scaled+(1|unique.ID)+(1|site)+(1|unique.block),family="poisson",data=subset(poar,flowerN_t1>0))
AICtab(Nflow.models, weights=T)

plot(poar$long.scaled[poar$Sex=="M" & poar$flowerN_t1>0],jitter(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0],factor=5),
     xlab="Longitude scaled",ylab="# Flowers",ylim=c(0,20))
points(poar$long.scaled[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,jitter(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0],factor=5),col="red")
lines(long.seq,exp(fixef(Nflow.models[[2]])[1] + fixef(Nflow.models[[2]])[2]*long.seq),lwd=3)


### same but log flowers
logNflow.models<-list()
logNflow.models[[1]]<-lmer(log(flowerN_t1)~(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))
logNflow.models[[2]]<-lmer(log(flowerN_t1)~long.scaled+(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))
logNflow.models[[3]]<-lmer(log(flowerN_t1)~Sex+(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))
logNflow.models[[4]]<-lmer(log(flowerN_t1)~Sex+long.scaled+(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))
logNflow.models[[5]]<-lmer(log(flowerN_t1)~Sex*long.scaled+(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))
AICtab(logNflow.models, weights=T)

plot(poar$long.scaled[poar$Sex=="M" & poar$flowerN_t1>0],jitter(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0],factor=5),
     xlab="Longitude scaled",ylab="# Flowers",ylim=c(0,20))
points(poar$long.scaled[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,jitter(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0],factor=5),col="red")
lines(long.seq,exp(fixef(Nflow.models[[5]])[1] + fixef(Nflow.models[[5]])[3]*long.seq),lwd=3,col="red")
lines(long.seq,exp(fixef(Nflow.models[[5]])[1] + fixef(Nflow.models[[5]])[2] + 
                     (fixef(Nflow.models[[5]])[3]+fixef(Nflow.models[[5]])[4])*long.seq),lwd=3)

## can these results generate a cline in flower sr? apparently yes
f.flowers<-c()
m.flowers<-c()
sr.flowers<-c()
for(i in 1:length(long.seq)){
  f.flowers[i]<-invlogit(fixef(prflow.models[[4]])[1]+fixef(prflow.models[[4]])[3]*long.seq[i])*exp(fixef(Nflow.models[[2]])[1] + fixef(Nflow.models[[2]])[2]*long.seq[i])
  m.flowers[i]<-invlogit(fixef(prflow.models[[4]])[1]+fixef(prflow.models[[4]])[2]+fixef(prflow.models[[4]])[3]*long.seq[i])*exp(fixef(Nflow.models[[2]])[1] + fixef(Nflow.models[[2]])[2]*long.seq[i])
  sr.flowers[i]<-f.flowers[i]/(f.flowers[i]+m.flowers[i])
}
plot(long.seq,sr.flowers)


## what about flower length?
poar$mean.flowerlength_t1<-rowMeans(cbind(poar$flowerLength1_t1,poar$flowerLength2_t1,poar$flowerLength3_t1),na.rm=T)
plot(poar$long.scaled[poar$Sex=="M" & poar$flowerN_t1>0],(poar$mean.flowerlength_t1[poar$Sex=="M" & poar$flowerN_t1>0]),
     xlab="Longitude scaled",ylab="Flower length",ylim=c(0,20))
points(poar$long.scaled[poar$Sex=="F" & poar$flowerN_t1>0]+0.05,(poar$mean.flowerlength_t1[poar$Sex=="F" & poar$flowerN_t1>0]),col="red")
range(poar$mean.flowerlength_t1)
poar$mean.flowerlength_t1[is.nan(poar$mean.flowerlength_t1)]<-NA
hist(log(poar$mean.flowerlength_t1))
range(poar$mean.flowerlength_t1,na.rm=T)
##why are there zeroes?
poar[which(poar$mean.flowerlength_t1==0),]
#...because there were zeroes entered. Need to check raw data sheets from 2015
flowlength.models<-list()
flowlength.models[[1]]<-lmer(log(mean.flowerlength_t1)~(1|site)+(1|unique.block),data=subset(poar,flowerN_t1>0))



## Calculate raw panicle sex ratio by site*year
flowers<-aggregate(flowerN_t1~site*year,data=poar,FUN=sum)
female.flowers<-aggregate(flowerN_t1~site*year,data=subset(poar, Sex=="F"),FUN=sum)
panicle.sr<-data.frame(merge(flowers,female.flowers,by=c("site","year")))
names(panicle.sr)[3]<-"total.panicles"
names(panicle.sr)[4]<-"female.panicles"
panicle.sr$sr<-panicle.sr$female.panicles / panicle.sr$total.panicles
panicle.sr    <- merge(panicle.sr,latlong,by="site")

plot(panicle.sr$long.scaled,panicle.sr$sr,cex=log(panicle.sr$total.panicles)+1)
plot(panicle.sr$long.scaled,panicle.sr$sr,pch=16)

sr.models <- list()
sr.models[[1]]<-glmer(sr ~ (1|site),weights = total.panicles, family="binomial",data=panicle.sr)
sr.models[[2]]<-glmer(sr ~ long.center + (1|site),weights = total.panicles,family="binomial",data=panicle.sr)
AICtab(sr.models,weights=T)
anova(sr.models[[1]],sr.models[[2]])

long.seq<-seq(min(poar$long.center),max(poar$long.center),0.1)
win.graph()
plot(panicle.sr$long.center+mean(poar$Longitude),panicle.sr$sr,cex=log(panicle.sr$total.panicles)+1,
     pch=21,bg="tomato",cex.lab=1.4,
     ylab="Operational sex ratio (proportion female)",xlab="Longitude")
lines(long.seq+mean(poar$Longitude),invlogit(fixef(sr.models[[2]])[1]+fixef(sr.models[[2]])[2]*long.seq),lwd=4)

## compare to Aldo's
# sex ratios
sex_r <- poar %>% 
  group_by(site,Sex,year) %>%
  summarise( site_flow = sum(flowerN_t1,na.rm=T) ) %>%
  as.data.frame() %>%
  inner_join(latlong) %>%
  select(site, Sex, year, site_flow,Longitude,Latitude) %>%
  spread(Sex, site_flow) %>%
  mutate( sr = F / (F + M) )

### save poar data frame
write.csv(poar,"C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/analysis/poar.csv")
