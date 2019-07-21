### Author: Aldo
### Purpose: 2017 range limits experiment analysis
### Last update: 4/18/2017
setwd("C:/Users/ac79/Downloads/Dropbox/poar--Aldo&Tom/Range limits/Experiment/Demography")
library(lme4)
library(bbmle)
library(dplyr)
library(tidyr)

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

# format ---------------------------------------------------------------

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
            # flag ploidy of provenance
            mutate(ploidy = "low") %>%
            mutate(ploidy = replace(ploidy, Code == "HHC" | Code == "LLELA", "high") ) %>%
            sex_symbol()
        

# Drop the 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
poar.drop  <- subset(poar, site != "llela" & site != "lubbock" & site != "wf" )
unique(poar.drop$site) ##this works
#str(poar.drop)


# plot data -----------------------------------------------------------

# survival
plot(jitter(surv_t1) ~ tillerN_t0, data=poar, pch = 16)

# flowering
plot(jitter(flow_t1) ~ tillerN_t0, data=poar, pch = 16)

# fertility (n. of flowers)
plot(jitter(flowerN_t1) ~ tillerN_t0, data=poar, pch = 16)

# growth
plot(tillerN_t1 ~ tillerN_t0, data=poar, pch = 16)
plot(clone_area_t1 ~ clone_area_t0, data=poar, pch = 16)

# growth measurement
boxplot(log_ratio ~ Sex, data = poar)
boxplot(delta.tiller ~ Sex, data = poar)

summary(lm(log_ratio ~ ploidy*Sex, data = poar))


# Analyses --------------------------------------------------------------

# growth  
gr      <- list()
gr[[1]] <- lmer(log_ratio ~ (1|site)+(1|unique.block), data = poar.drop)
gr[[2]] <- lmer(log_ratio ~ Sex + (1|site)+(1|unique.block), data = poar.drop)
gr[[3]] <- lmer(log_ratio ~ Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[4]] <- lmer(log_ratio ~ Sex + Longitude + (1|site)+(1|unique.block), data = poar.drop)
gr[[5]] <- lmer(log_ratio ~ Sex * Longitude + (1|site)+(1|unique.block), data = poar.drop)
AICtab(gr, weights=T)

## survival analysis -- test sex*longitude with random block effect   ***convergence problems with fixed longitude and random site!
#surv0<-glmer(tillerN_t1>0 ~ (1|site)+(1|unique.block),family="binomial",data=poar.drop)
#surv1<-glmer(tillerN_t1>0 ~ Sex + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
#surv2<-glmer(tillerN_t1>0 ~ Longitude + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
#surv3<-glmer(tillerN_t1>0 ~ Sex + Longitude + (1|site)+(1|unique.block),family="binomial",data=poar.drop)
# surv4<-glmer(tillerN_t1>0 ~ Sex * Longitude + (1|site),family="binomial",data=poar.drop)
# AICtab(surv0,surv1,surv2,surv3,surv4,weights=T,sort=F)
# ## FUCK!

# survival (glms)
sr   <- list()
sr[[1]]<-glm(tillerN_t1>0 ~ tillerN_t0 + 1,family="binomial",data=poar.drop)
sr[[2]]<-glm(tillerN_t1>0 ~ tillerN_t0 + Sex ,family="binomial",data=poar.drop)
sr[[3]]<-glm(tillerN_t1>0 ~ tillerN_t0 + Longitude ,family="binomial",data=poar.drop)
sr[[4]]<-glm(tillerN_t1>0 ~ tillerN_t0 + Sex + Longitude ,family="binomial",data=poar.drop)
sr[[5]]<-glm(tillerN_t1>0 ~ tillerN_t0 + Sex * Longitude ,family="binomial",data=poar.drop)
AICtab(sr,weights=T)

# flowering probability
# fl <- list()
# fl[[1]]<-glmer(flowerN_t1>0 ~ (1|unique.block),family="binomial",data=poar.drop)
# fl[[2]]<-glmer(flowerN_t1>0 ~ Sex ,family="binomial",data=poar.drop)
# fl[[3]]<-glmer(flowerN_t1>0 ~ Longitude ,family="binomial",data=poar.drop)
# fl[[4]]<-glmer(flowerN_t1>0 ~ Sex + Longitude ,family="binomial",data=poar.drop)
# fl[[5]]<-glmer(flowerN_t1>0 ~ Sex * Longitude ,family="binomial",data=poar.drop)
# AICtab(fl,weights=T)
fl <- list()
fl[[1]]<-glm(flowerN_t1>0 ~ tillerN_t1,family="binomial",data=poar.drop)
fl[[2]]<-glm(flowerN_t1>0 ~ tillerN_t1 + Sex ,family="binomial",data=poar.drop)
fl[[3]]<-glm(flowerN_t1>0 ~ tillerN_t1 + Longitude ,family="binomial",data=poar.drop)
fl[[4]]<-glm(flowerN_t1>0 ~ tillerN_t1 + Sex + Longitude ,family="binomial",data=poar.drop)
fl[[5]]<-glm(flowerN_t1>0 ~ tillerN_t1 + Sex * Longitude ,family="binomial",data=poar.drop)
AICtab(fl,weights=T)


# sex ratios
sex_r <- poar.drop %>% 
          group_by(site,Sex,year) %>%
          summarise( site_flow = sum(flowerN_t1,na.rm=T) ) %>%
          as.data.frame() %>%
          inner_join(latlong) %>%
          select(site, Sex, year, site_flow,Longitude,Latitude) %>%
          spread(Sex, site_flow) %>%
          mutate( sr = F / (F + M) )

# sex ratios
fe <- list()
fe[[1]]<-glm(sr ~ 1,weights = M+F, family="binomial",data=sex_r)
fe[[2]]<-glm(sr ~ Longitude,weights = M+F,family="binomial",data=sex_r)
fe[[3]]<-glm(sr ~ Latitude ,weights = M+F,family="binomial",data=sex_r)
fe[[4]]<-glm(sr ~ Longitude + Latitude,weights = M+F,family="binomial",data=sex_r)
fe[[5]]<-glm(sr ~ Longitude * Latitude,weights = M+F,family="binomial",data=sex_r)
AICtab(fe,weights=T)


# plots ------------------------------------------------------------------
tiff("results/vr_models.tiff",unit="in", width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(3,3,0.1,0.4),mgp=c(1.6,0.5,0),cex.lab=1.4,cex.axis=1,
    oma=c(0,0,0.2,0))

# write "xSeq" one only
xSeq <- seq(min(latlong$Longitude), max(latlong$Longitude), by = 1)


# survival
# format data
surv <- poar.drop %>% 
  select(site,Sex,year,surv_t1) %>%
  subset(!is.na(surv_t1)) %>%
  group_by(site,Sex,year) %>%
  summarise( site_surv = sum(surv_t1,na.rm=T),
             total_ind = n()) %>%
  mutate(surv_pr = site_surv / total_ind) %>%
  as.data.frame() %>%
  inner_join(latlong) %>%
  sex_symbol()

m_s  <- mean(poar.drop$tillerN_t0, na.rm=T)
beta <- coef(sr[[3]])
y  <- beta[1] + beta[2]*m_s + beta[3]*xSeq
plot(surv_pr ~ Longitude, pch = symb, data = surv, bg = col,
     ylab= "Probability of survival", ylim = c(0,1))
lines(xSeq, invlogit(y), lwd = 2)
#lines(xSeq, invlogit(y_f), col = "blue")

legend("bottomleft", c("male", "female"), pch = c(17,16), col = c("red","blue"), 
       bty="n")

# growth (log ratio)
plot(log_ratio ~ Longitude, pch = symb, bg = col, data = poar.drop,
     ylab = "Tiller log ratio (growth)")
beta  <- fixef(gr[[3]])
y     <- beta[1] + beta[2]*xSeq
lines(xSeq, y, lwd=2)

# flowering 
# format data
flow <- poar.drop %>% 
          group_by(site,Sex,year) %>%
          summarise( site_flow = sum(flow_t1,na.rm=T),
                     total_ind = n()) %>%
          mutate(flow_pr = site_flow / total_ind) %>%
          as.data.frame() %>%
          inner_join(latlong) %>%
          sex_symbol()

m_s  <- mean(poar.drop$tillerN_t1,na.rm=T)
beta <- coef(fl[[5]])
y_m  <- beta[1] + beta[2]*m_s + beta[3] + beta[4]*xSeq + beta[5]*xSeq
y_f  <- beta[1] + beta[2]*m_s + beta[4]*xSeq
plot(flow_pr ~ Longitude, pch = symb, data = flow, bg = col,
     ylab= "Probability of flowering", ylim = c(0,1))
lines(xSeq, invlogit(y_m), col = "red", lty = 2, lwd=2)
lines(xSeq, invlogit(y_f), col = "blue", lwd=2)

legend("topleft", c("male", "female"), pch = c(17,16), , lwd=2,col = c("red","blue"), 
       lty = c(2,1), bty="n")

# Sex ratios
plot(sr ~ Longitude, pch = 16, data=sex_r,
     ylab = "Panicle sex ratio")
beta <- coef(fe[[2]])
y  <- invlogit(beta[1] + beta[2]*xSeq)
lines(xSeq, y, lwd = 2)

dev.off()

