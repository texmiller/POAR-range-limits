### Author: Tom
### Purpose: preliminary analysis of range limits experiment for NSF proposal
### Last update: 7/14/2015
library(lme4);library(bbmle)
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}
invlogit<-function(x){exp(x)/(1+exp(x))}

## read in data frame that Aldo constructed (see "fall14_spr15.R and associated metadata)
POAR<-read.csv("D:/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/f14-s15DemoData.csv")
POAR<-read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/f14-s15DemoData.csv")

## add lat/long to data frame
latlong<-read.csv("D:/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/SiteLatLong.csv")
latlong<-read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/SiteLatLong.csv")
POAR<-merge(POAR,latlong,by="site")
## add a unique block ID to fit as random effect
POAR$unique.block<-interaction(POAR$site,POAR$Block)
## add flowering (bernoulli)
POAR$flow_t1<-POAR$flowerN_t1>0
## change in size
POAR$delta.tiller<-POAR$tillerN_t1-POAR$tillerN_t0

boxplot(delta.tiller~Sex*Longitude,data=POAR,col=c("gray","tomato"))

## here is a version of the data frame that drops 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
POAR.drop<-POAR[!(POAR$site=="llela" | POAR$site=="lubbock" | POAR$site=="wf"),]
unique(POAR.drop$site) ##this works
str(POAR.drop)

## survival analysis -- test sex*longitude with random block effect   ***convergence problems with fixed longitude and random site!
surv0<-glmer(tillerN_t1>0 ~ (1|site)+(1|unique.block),family="binomial",data=POAR.drop)
surv1<-glmer(tillerN_t1>0 ~ Sex + (1|site)+(1|unique.block),family="binomial",data=POAR.drop)
surv2<-glmer(tillerN_t1>0 ~ Longitude + (1|site)+(1|unique.block),family="binomial",data=POAR.drop)
surv3<-glmer(tillerN_t1>0 ~ Sex + Longitude + (1|site)+(1|unique.block),family="binomial",data=POAR.drop)
surv4<-glmer(tillerN_t1>0 ~ Sex * Longitude + (1|site)+ (1|unique.block),family="binomial",data=POAR.drop)
AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)
## FUCK!

## switch to glms
surv0<-glm(tillerN_t1>0 ~ 1,family="binomial",data=POAR.drop)
surv1<-glm(tillerN_t1>0 ~ Sex ,family="binomial",data=POAR.drop)
surv2<-glm(tillerN_t1>0 ~ Longitude ,family="binomial",data=POAR.drop)
surv3<-glm(tillerN_t1>0 ~ Sex + Longitude ,family="binomial",data=POAR.drop)
surv4<-glm(tillerN_t1>0 ~ Sex * Longitude ,family="binomial",data=POAR.drop)
AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)


## flowering analysis
flow0<-glmer(flowerN_t1>0 ~ (1|unique.block),family="binomial",data=POAR.drop)
flow1<-glmer(flowerN_t1>0 ~ Sex + (1|unique.block),family="binomial",data=POAR.drop)
flow2<-glmer(flowerN_t1>0 ~ Longitude + (1|unique.block),family="binomial",data=POAR.drop)
flow3<-glmer(flowerN_t1>0 ~ Sex + Longitude + (1|unique.block),family="binomial",data=POAR.drop)
flow4<-glmer(flowerN_t1>0 ~ Sex * Longitude + (1|unique.block),family="binomial",data=POAR.drop)
AICtab(flow0,flow1,flow2,flow3,flow4,weights=T,sort=F)

flow0<-glm(flowerN_t1>0 ~ 1,family="binomial",data=POAR)
flow1<-glm(flowerN_t1>0 ~ Sex ,family="binomial",data=POAR)
flow2<-glm(flowerN_t1>0 ~ Longitude ,family="binomial",data=POAR)
flow3<-glm(flowerN_t1>0 ~ Sex + Longitude ,family="binomial",data=POAR)
flow4<-glm(flowerN_t1>0 ~ Sex * Longitude ,family="binomial",data=POAR)
flow5<-glm(flowerN_t1>0 ~ Longitude + I(Longitude^2),family="binomial",data=POAR)
flow6<-glm(flowerN_t1>0 ~ Sex + Longitude + I(Longitude^2),family="binomial",data=POAR)
flow7<-glm(flowerN_t1>0 ~ Sex * Longitude + Sex*I(Longitude^2),family="binomial",data=POAR)


AICtab(flow0,flow1,flow2,flow3,flow4,flow5,flow6,flow7,weights=T,sort=F)
anova(flow0,flow4,test="Chisq")


aggregate.data<-POAR.drop
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
aggregate(tillerN_t0 ~ Sex,length2,na.rm=T,data=POAR)
