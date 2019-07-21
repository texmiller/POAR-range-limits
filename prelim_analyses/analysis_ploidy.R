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
# Symbols for the sexes 
ploidy_symbol <- function(x){
  out <- x %>% 
    mutate(symb = factor(ploidy, labels=c("21","24")) ) %>%
    mutate(symb = as.character(symb) ) %>%
    mutate(symb = as.integer(symb) ) %>%
    mutate(col = factor(ploidy, labels=c("blue","red")) ) %>%
    mutate(col = as.character(col) )
  return(out)
}

# read data -----------------------------------------------------------------
poar    <- read.csv("data/f14_s17_data.csv")
latlong <- read.csv("data/SiteLatLong.csv")

# format ---------------------------------------------------------------

# add lat/long to data frame
poar    <- merge(poar,latlong,by="site")
poar    <- poar %>% 
            # add a unique block ID to fit as random effect
            mutate(unique.block = interaction(site,Block) ) %>%
            mutate(flow_t1 = as.numeric(flowerN_t1>0) ) %>% # add flowering (bernoulli)
            mutate(surv_t1 = as.numeric(tillerN_t1>0) ) %>% # add survival (bernoulli)
            mutate(delta.tiller = tillerN_t1 - tillerN_t0) %>% # change in size
            mutate(log_ratio = log(tillerN_t1 / tillerN_t0) ) %>% # log ratio
            mutate(log_ratio = replace(log_ratio, log_ratio == -Inf | log_ratio == Inf | is.nan(log_ratio), NA) ) %>%
            # flag ploidy of provenance
            mutate(ploidy = "low") %>%
            mutate(ploidy = replace(ploidy, Code == "HHC" | Code == "LLELA", "high") ) %>%
            ploidy_symbol()
  
# Drop the 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
poar.drop  <- subset(poar, site != "llela" & site != "lubbock" & site != "wf" )
unique(poar.drop$site) ##this works
#str(poar.drop)


# plot data -----------------------------------------------------------

# survival
plot(jitter(surv_t1) ~ tillerN_t0, data=poar)

# flowering
plot(jitter(flow_t1) ~ tillerN_t0, data=poar)

# fertility (n. of flowers)
plot(jitter(flowerN_t1) ~ tillerN_t0, data=poar)

# growth
plot(tillerN_t1 ~ tillerN_t0, data=poar)
plot(clone_area_t1 ~ clone_area_t0, data=poar)

# 
plot(tillerN_t1 ~ tillerN_t0, data = poar)

# growth measurement
boxplot(log_ratio ~ Sex, data = poar)
boxplot(delta.tiller ~ Sex, data = poar)

summary(lm(log_ratio ~ ploidy*Sex, data = poar))

# Analyses --------------------------------------------------------------

# growth  
gr      <- list()
gr[[1]] <- lmer(delta.tiller ~ (1|site)+(1|unique.block), data = poar.drop)
gr[[2]] <- lmer(delta.tiller ~ ploidy + (1|site)+(1|unique.block), data = poar.drop)
gr[[3]] <- lmer(delta.tiller ~ Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[4]] <- lmer(delta.tiller ~ ploidy + Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[5]] <- lmer(delta.tiller ~ ploidy * Longitude + (1|site)+(1|unique.block), data = poar.drop)
AICtab(gr, weights=T)

gr      <- list()
gr[[1]] <- lmer(log_ratio ~ (1|site)+(1|unique.block), data = poar.drop)
gr[[2]] <- lmer(log_ratio ~ ploidy + (1|site)+(1|unique.block), data = poar.drop)
gr[[3]] <- lmer(log_ratio ~ Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[4]] <- lmer(log_ratio ~ ploidy + Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[5]] <- lmer(log_ratio ~ ploidy * Longitude + (1|site)+(1|unique.block), data = poar.drop)
AICtab(gr, weights=T)


gr      <- list()
gr[[1]] <- glm(tillerN_t1 ~ tillerN_t0 , family="poisson",data = poar.drop)
gr[[2]] <- glm(tillerN_t1 ~ tillerN_t0 + ploidy , family="poisson",data = poar.drop)
gr[[3]] <- glm(tillerN_t1 ~ tillerN_t0 + Longitude , family="poisson",data = poar.drop)
gr[[4]] <- glm(tillerN_t1 ~ tillerN_t0 + ploidy + Longitude , family="poisson",data = poar.drop)
gr[[5]] <- glm(tillerN_t1 ~ tillerN_t0 + ploidy * Longitude , family="poisson",data = poar.drop)
AICtab(gr, weights=T)


plot()

# survival
sr      <- list()
sr[[1]] <- glm(tillerN_t1>0 ~ tillerN_t0, family="binomial",data = poar.drop)
sr[[2]] <- glm(tillerN_t1>0 ~ tillerN_t0 , family="binomial",data = poar.drop)
sr[[3]] <- glm(tillerN_t1>0 ~ tillerN_t0 + ploidy , family="binomial",data = poar.drop)
sr[[4]] <- glm(tillerN_t1>0 ~ tillerN_t0 * ploidy , family="binomial",data = poar.drop)
sr[[5]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude , family="binomial",data = poar.drop)
sr[[6]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude + ploidy, family="binomial",data = poar.drop)
sr[[7]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude * ploidy, family="binomial",data = poar.drop)
AICtab(sr, weights=T)


# flowering
fl      <- list()
fl[[1]] <- glm(tillerN_t1>0 ~ tillerN_t0, family="binomial",data = poar.drop)
fl[[2]] <- glm(tillerN_t1>0 ~ tillerN_t0 , family="binomial",data = poar.drop)
fl[[3]] <- glm(tillerN_t1>0 ~ tillerN_t0 + ploidy , family="binomial",data = poar.drop)
fl[[4]] <- glm(tillerN_t1>0 ~ tillerN_t0 * ploidy , family="binomial",data = poar.drop)
fl[[5]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude , family="binomial",data = poar.drop)
fl[[6]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude + ploidy, family="binomial",data = poar.drop)
fl[[7]] <- glm(tillerN_t1>0 ~ tillerN_t0 + Longitude * ploidy, family="binomial",data = poar.drop)
AICtab(fl, weights=T)


gr      <- list()
gr[[1]] <- glm(tillerN_t1 ~ 1 , family="poisson",data = poar.drop)
gr[[2]] <- glm(tillerN_t1 ~ tillerN_t0 , family="poisson",data = poar.drop)
gr[[3]] <- glm(tillerN_t1 ~ tillerN_t0 + ploidy , family="poisson",data = poar.drop)
gr[[4]] <- glm(tillerN_t1 ~ tillerN_t0 * ploidy , family="poisson",data = poar.drop)
AICtab(gr, weights=T)


plot(tillerN_t1 ~ tillerN_t0, data = poar.drop,
     pch = symb, bg = col)

xSeq <- seq(0,90)
xSeq1 <- seq(0,30)
beta <- coef(gr[[4]])
y_h  <- beta[1] + beta[2]*xSeq1
y_l  <- beta[1] + beta[2]*xSeq + beta[3] + beta[4]*xSeq
lines(xSeq1, exp(y_h), col = "blue", lwd=2) 
lines(xSeq, exp(y_l), col = "red", lwd=2)

# survival
surv0<-glm(tillerN_t1>0 ~ tillerN_t0 + 1,family="binomial",data=poar.drop)
surv1<-glm(tillerN_t1>0 ~ tillerN_t0 + ploidy ,family="binomial",data=poar.drop)
surv2<-glm(tillerN_t1>0 ~ tillerN_t0 + Longitude ,family="binomial",data=poar.drop)
surv3<-glm(tillerN_t1>0 ~ tillerN_t0 + ploidy + Longitude ,family="binomial",data=poar.drop)
surv4<-glm(tillerN_t1>0 ~ tillerN_t0 + ploidy * Longitude ,family="binomial",data=poar.drop)
AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)

# flowering prob
flow0<-glm(flowerN_t1>0 ~ 1,family="binomial",data=poar)
flow1<-glm(flowerN_t1>0 ~ ploidy ,family="binomial",data=poar)
flow2<-glm(flowerN_t1>0 ~ Longitude ,family="binomial",data=poar)
flow3<-glm(flowerN_t1>0 ~ ploidy + Longitude ,family="binomial",data=poar)
flow4<-glm(flowerN_t1>0 ~ ploidy * Longitude ,family="binomial",data=poar)
AICtab(flow0,flow1,flow2,flow3,flow4,weights=T,sort=F)

# n of flowers
n_flow0<-glm(flowerN_t1 ~ 1,family="poisson",data=poar)
n_flow1<-glm(flowerN_t1 ~ ploidy ,family="poisson",data=poar)
n_flow2<-glm(flowerN_t1 ~ Longitude ,family="poisson",data=poar)
n_flow3<-glm(flowerN_t1 ~ ploidy + Longitude ,family="poisson",data=poar)
n_flow4<-glm(flowerN_t1 ~ ploidy * Longitude ,family="poisson",data=poar)
AICtab(n_flow0,n_flow1,n_flow2,n_flow3,n_flow4,weights=T,sort=F)

