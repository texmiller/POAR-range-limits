### Author: Aldo
### Purpose: 2017 range limits experiment analysis
### Last update: 4/18/2017
setwd("C:/Users/ac79/Downloads/Dropbox/poar--Aldo&Tom/Range limits/Experiment/Demography")
library(lme4)
library(bbmle)
library(dplyr)

# useful functions ----------------------------------------------------
length2 <- function(x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}
invlogit<-function(x){exp(x)/(1+exp(x))}

# read data -----------------------------------------------------------------
poar    <- read.csv("data/f14_s17_data.csv")
latlong <- read.csv("data/SiteLatLong.csv")

# format ---------------------------------------------------------------

# add lat/long to data frame
poar    <- merge(poar,latlong,by="site")
poar    <- poar %>% 
            # add a unique block ID to fit as random effect
            mutate(unique.block = interaction(site,Block) ) %>%
            mutate(flow_t1 = flowerN_t1>0) %>% # add flowering (bernoulli)
            mutate(delta.tiller = tillerN_t1 - tillerN_t0)    # change in size
            mutate(log_ratio = log(tillerN_t1 / tillerN_t0) ) # log ratio


# Drop the 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
poar.drop  <- subset(poar, site != "llela" & site != "lubbock" & site != "wf" )
unique(poar.drop$site) ##this works
str(poar.drop)


# plot data -----------------------------------------------------------

# survival
poar

# growth
boxplot(log_ratio ~ sex, data = poar)


# Analyses 

# growth -------------------------------------------------------------- 

gr      <- list()
gr[[1]] <- lmer(delta.tiller ~ (1|site)+(1|unique.block), data = poar.drop)
gr[[2]] <- lmer(delta.tiller ~ Sex + (1|site)+(1|unique.block), data = poar.drop)
gr[[3]] <- lmer(delta.tiller ~ Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[4]] <- lmer(delta.tiller ~ Sex + Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[5]] <- lmer(delta.tiller ~ Sex * Longitude + (1|site)+(1|unique.block), data = poar.drop)
AICtab(gr, weights=T)


## survival analysis -- test sex*longitude with random block effect   ***convergence problems with fixed longitude and random site!
surv0<-glmer(tillerN_t1>0 ~ (1|site)+(1|unique.block),family="binomial",data=poar.drop)
surv1<-glmer(tillerN_t1>0 ~ Sex + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
surv2<-glmer(tillerN_t1>0 ~ Longitude + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
surv3<-glmer(tillerN_t1>0 ~ Sex + Longitude + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
surv4<-glmer(tillerN_t1>0 ~ Sex * Longitude + (1|site)+ (1|unique.block),family="binomial",data=poar.drop)
AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)
## FUCK!

## switch to glms
surv0<-glm(tillerN_t1>0 ~ 1,family="binomial",data=poar.drop)
surv1<-glm(tillerN_t1>0 ~ Sex ,family="binomial",data=poar.drop)
surv2<-glm(tillerN_t1>0 ~ Longitude ,family="binomial",data=poar.drop)
surv3<-glm(tillerN_t1>0 ~ Sex + Longitude ,family="binomial",data=poar.drop)
surv4<-glm(tillerN_t1>0 ~ Sex * Longitude ,family="binomial",data=poar.drop)
AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)


## flowering analysis
flow0<-glmer(flowerN_t1>0 ~ (1|unique.block),family="binomial",data=poar.drop)
flow1<-glmer(flowerN_t1>0 ~ Sex + (1|unique.block),family="binomial",data=poar.drop)
flow2<-glmer(flowerN_t1>0 ~ Longitude + (1|unique.block),family="binomial",data=poar.drop)
flow3<-glmer(flowerN_t1>0 ~ Sex + Longitude + (1|unique.block),family="binomial",data=poar.drop)
flow4<-glmer(flowerN_t1>0 ~ Sex * Longitude + (1|unique.block),family="binomial",data=poar.drop)
AICtab(flow0,flow1,flow2,flow3,flow4,weights=T,sort=F)

flow0<-glm(flowerN_t1>0 ~ 1,family="binomial",data=poar)
flow1<-glm(flowerN_t1>0 ~ Sex ,family="binomial",data=poar)
flow2<-glm(flowerN_t1>0 ~ Longitude ,family="binomial",data=poar)
flow3<-glm(flowerN_t1>0 ~ Sex + Longitude ,family="binomial",data=poar)
flow4<-glm(flowerN_t1>0 ~ Sex * Longitude ,family="binomial",data=poar)
flow5<-glm(flowerN_t1>0 ~ Longitude + I(Longitude^2),family="binomial",data=poar)
flow6<-glm(flowerN_t1>0 ~ Sex + Longitude + I(Longitude^2),family="binomial",data=poar)
flow7<-glm(flowerN_t1>0 ~ Sex * Longitude + Sex*I(Longitude^2),family="binomial",data=poar)


AICtab(flow0,flow1,flow2,flow3,flow4,flow5,flow6,flow7,weights=T,sort=F)
anova(flow0,flow4,test="Chisq")


aggregate.data<-poar.drop
P.flower.aggregate<-merge(data.frame(cbind(aggregate(flow_t1 ~ site*Sex,sum,na.rm=T,data=aggregate.data)[,1:2],
aggregate(flow_t1 ~ site*Sex,sum,na.rm=T,data=aggregate.data)[,3]/
  aggregate(flow_t1 ~ site*Sex,length2,na.rm=T,data=aggregate.data)[3])),latlong,by="site")

win.graph()
plot(P.flower.aggregate$Longitude,P.flower.aggregate$flow_t1,type="n",ylim=c(0,1),
     xlab="Longitude",ylab="Probability of flowering")
points(P.flower.aggregate$Longitude[P.flower.aggregate$Sex=="F"],
       P.flower.aggregate$flow_t1[P.flower.aggregate$Sex=="F"],pch=21,bg="tomato",cex=1.4)
points(P.flower.aggregate$Longitude[P.flower.aggregate$Sex=="M"],
       P.flower.aggregate$flow_t1[P.flower.aggregate$Sex=="M"],pch=1,,cex=1.4)

long.dum<-seq(min(P.flower.aggregate$Longitude),max(P.flower.aggregate$Longitude),0.0001)
lines(long.dum,invlogit(coef(flow4)[1]+coef(flow4)[3]*long.dum),lwd=3)
lines(long.dum,invlogit(coef(flow4)[1]+coef(flow4)[2]+
                          (coef(flow4)[3]+coef(flow4)[4])*long.dum),lwd=3,lty=2)
legend("topleft",legend=c("Female","Male"),pch=c(21,1),lty=c(1,2),pt.bg="tomato",lwd=2,
       bty="n",cex=1.1,pt.cex=1.4)



lines(long.dum,invlogit(fixef(flow4)[1]+fixef(flow4)[3]*long.dum),lwd=3,col="tomato")
lines(long.dum,invlogit(fixef(flow4)[1]+fixef(flow4)[2]+
                          (fixef(flow4)[3]+fixef(flow4)[4])*long.dum),lwd=3,lty=2,col="tomato")


## how many plants in total?
aggregate(tillerN_t0 ~ Sex,length2,na.rm=T,data=poar)
