### Bayesian vital rate models for POAR range limits
library(R2jags)
library(mcmcplots)
setwd("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/analysis")
poar<-read.csv("poar.csv")

## Model
## Notes: The Gamma-Poisson mixture appears very unstable for panicle production.
## I am therefore going with Poisson model, with normally distributed individual ran effect

sink("POAR_vital_rates_unstructured.txt")
cat("
    model{
    
    ## Priors
    #1. Survival
    alpha.surv~dnorm(0,0.001)   ##hyperprior for survival
    sigma.surv.site~dunif(0,1000) ##site variance
    tau.surv.site<-1/(sigma.surv.site*sigma.surv.site)
    sigma.surv.block~dunif(0,1000) ##site variance
    tau.surv.block<-1/(sigma.surv.block*sigma.surv.block)
    beta.surv.sex~dnorm(0,0.001)
    beta.surv.long~dnorm(0,0.001)
    beta.surv.sex.long~dnorm(0,0.001)

    #2. Growth
    alpha.grow~dnorm(0,0.001)   ##hyperprior for growth
    sigma.grow.site~dunif(0,1000) ##site variance
    tau.grow.site<-1/(sigma.grow.site*sigma.grow.site)
    sigma.grow.block~dunif(0,1000) ##site variance
    tau.grow.block<-1/(sigma.grow.block*sigma.grow.block)
    sigma.grow~dunif(0,1000) ##residual growth variance
    tau.grow<-1/(sigma.grow*sigma.grow)    
    beta.grow.sex~dnorm(0,0.001)
    beta.grow.long~dunif(-1,1)
    beta.grow.sex.long~dnorm(0,0.001)

    #3. Flowering
    alpha.flow~dnorm(0,0.001)   ##hyperprior for flowival
    sigma.flow.site~dunif(0,1000) ##site variance
    tau.flow.site<-1/(sigma.flow.site*sigma.flow.site)
    sigma.flow.block~dunif(0,1000) ##site variance
    tau.flow.block<-1/(sigma.flow.block*sigma.flow.block)
    beta.flow.sex~dnorm(0,0.001)
    beta.flow.long~dnorm(0,0.001)
    beta.flow.sex.long~dnorm(0,0.001)

    #4. Panicles
    alpha.panic~dnorm(0,0.001)   ##hyperprior for panicival
    sigma.panic.site~dunif(0,1000) ##site variance
    tau.panic.site<-1/(sigma.panic.site*sigma.panic.site)
    sigma.panic.block~dunif(0,1000) ##site variance
    tau.panic.block<-1/(sigma.panic.block*sigma.panic.block)
    beta.panic.sex~dnorm(0,0.001)
    beta.panic.long~dnorm(0,0.001)
    beta.panic.sex.long~dnorm(0,0.001)
    sigma.panic.ind~dunif(0,1000) ##individual random effect for overdispersion
    tau.panic.ind<-1/(sigma.panic.ind*sigma.panic.ind)
    #alpha.gamma.panic~dunif(0,100)

    # Assign random effects for sites and blocks
    for(i in 1:N.sites){
    alpha.surv.site[i]~dnorm(alpha.surv,tau.surv.site)
    alpha.grow.site[i]~dnorm(alpha.grow,tau.grow.site)
    alpha.flow.site[i]~dnorm(alpha.flow,tau.flow.site)
    alpha.panic.site[i]~dnorm(alpha.panic,tau.panic.site)
    for(j in 1:N.blocks){
    alpha.surv.site.block[i,j]~dnorm(alpha.surv.site[i],tau.surv.block)
    alpha.grow.site.block[i,j]~dnorm(alpha.grow.site[i],tau.grow.block)
    alpha.flow.site.block[i,j]~dnorm(alpha.flow.site[i],tau.flow.block)
    alpha.panic.site.block[i,j]~dnorm(alpha.panic.site[i],tau.panic.block)
    }
    }
    
    #1. Survival Likelihood 
    for(i in 1:N.obs.surv){
    y.surv[i]~dbern(p.surv[i])
    logit(p.surv[i])<-alpha.surv.site.block[site.surv[i],block.surv[i]] + 
                 beta.surv.sex*male.surv[i]+
                 beta.surv.long*long.surv[i]+
                 beta.surv.sex.long*male.surv[i]*long.surv[i]
    }
    #2. Growth Likelihood 
    for(i in 1:N.obs.grow){
    y.grow[i]~dnorm(mu[i],tau.grow)
    mu[i]<-alpha.grow.site.block[site.grow[i],block.grow[i]] + 
            beta.grow.sex*male.grow[i]+
            beta.grow.long*long.grow[i]+
            beta.grow.sex.long*male.grow[i]*long.grow[i]
    }
    #3. Flowering Likelihood 
    for(i in 1:N.obs.flow){
    y.flow[i]~dbern(p.flow[i])
    logit(p.flow[i])<-alpha.flow.site.block[site.flow[i],block.flow[i]] + 
                 beta.flow.sex*male.flow[i]+
                 beta.flow.long*long.flow[i]+
                 beta.flow.sex.long*male.flow[i]*long.flow[i]
    }
    #4. Panicles Likelihood 
    for(i in 1:N.obs.panic){
    y.panic[i]~dpois(mu.panic[i])T(1,)
    #mustar[i]<-rho[i]*mu.panic[i]
    log(mu.panic[i])<-alpha.panic.site.block[site.panic[i],block.panic[i]] + 
    beta.panic.sex*male.panic[i]+
    beta.panic.long*long.panic[i]+
    beta.panic.sex.long*male.panic[i]*long.panic[i]
    #rho[i]~dnorm(0,tau.panic.ind)
    #rho[i]~dgamma(alpha.gamma.panic,alpha.gamma.panic)
    }

    ## Prediction
    for(i in 1:N.longitudes){
    surv.F.pred[i]<-exp(alpha.surv+beta.surv.long*longitudes[i])/(1+exp(alpha.surv+beta.surv.long*longitudes[i]))
    surv.M.pred[i]<-exp(alpha.surv+beta.surv.sex+(beta.surv.long+beta.surv.sex.long)*longitudes[i])/(1+exp(alpha.surv+beta.surv.sex+(beta.surv.long+beta.surv.sex.long)*longitudes[i]))
    grow.F.pred[i]<-alpha.grow+beta.grow.long*longitudes[i]
    grow.M.pred[i]<-alpha.grow+beta.grow.sex+(beta.grow.long+beta.grow.sex.long)*longitudes[i]
    flow.F.pred[i]<-exp(alpha.flow+beta.flow.long*longitudes[i])/(1+exp(alpha.flow+beta.flow.long*longitudes[i]))
    flow.M.pred[i]<-exp(alpha.flow+beta.flow.sex+(beta.flow.long+beta.flow.sex.long)*longitudes[i])/(1+exp(alpha.flow+beta.flow.sex+(beta.flow.long+beta.flow.sex.long)*longitudes[i]))
    panic.F.pred[i]<-exp(alpha.panic+beta.panic.long*longitudes[i])
    panic.M.pred[i]<-exp(alpha.panic+beta.panic.sex+(beta.panic.long+beta.panic.sex.long)*longitudes[i])
    }

    ## Posterior predictive check
    # Panicle
    for(i in 1:N.obs.panic){
    panic.resi[i] <- (y.panic[i] - mu.panic[i])/sqrt(mu.panic[i])
    y.panic.new[i] ~ dpois(mu.panic[i])T(1,)
    panic.resi.new[i] <- (y.panic.new[i] - mu.panic[i])/sqrt(mu.panic[i])
    panic.D[i] <- pow(panic.resi[i],2)
    panic.D.new[i] <- pow(panic.resi.new[i],2)
    }
    panic.fit <- sum(panic.D[])
    panic.fit.new <- sum(panic.D.new[])

    }##end model
    ",fill=T)
sink()


#### Bundle data
#1. Survival
poar.surv<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$surv_t1)))
names(poar.surv)<-c("site","block","sex","long.center","long.scaled","surv_t1")
N.sites<-max(poar.surv$site)
site.surv<-poar.surv$site
N.blocks<-max(poar.surv$block)
block.surv<-poar.surv$block
male.surv<-poar.surv$sex-1
long.surv<-poar.surv$long.center
y.surv<-poar.surv$surv_t1
N.obs.surv<-nrow(poar.surv)
#2. Growth
poar.grow<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$log_ratio)))
names(poar.grow)<-c("site","block","sex","long.center","long.scaled","logratio")
max(poar.surv$site)==max(poar.grow$site)
max(poar.surv$block)==max(poar.grow$block)## sites and blocks line up
site.grow<-poar.grow$site
block.grow<-poar.grow$block
male.grow<-poar.grow$sex-1
long.grow<-poar.grow$long.center
y.grow<-poar.grow$logratio
N.obs.grow<-nrow(poar.grow)
#3. Flowering
poar.flow<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$flow_t1)))
names(poar.flow)<-c("site","block","sex","long.center","long.scaled","flow_t1")
max(poar.surv$site)==max(poar.flow$site)
max(poar.surv$block)==max(poar.flow$block)## sites and blocks line up
site.flow<-poar.flow$site
block.flow<-poar.flow$block
male.flow<-poar.flow$sex-1
long.flow<-poar.flow$long.center
y.flow<-poar.flow$flow_t1
N.obs.flow<-nrow(poar.flow)
#4. Panicles
poar.panic<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$flowerN_t1)))
names(poar.panic)<-c("site","block","sex","long.center","long.scaled","panic_t1")
poar.panic<-subset(poar.panic,panic_t1>0)
max(poar.surv$site)==max(poar.panic$site)
max(poar.surv$block)==max(poar.panic$block)## sites and blocks line up
site.panic<-poar.panic$site
block.panic<-poar.panic$block
male.panic<-poar.panic$sex-1
long.panic<-poar.panic$long.center
y.panic<-poar.panic$panic_t1
N.obs.panic<-nrow(poar.panic)

long.seq<-seq(min(poar$long.center),max(poar$long.center),length.out = 30)

jag.data<-list(N.sites=N.sites,
               N.blocks=N.blocks,
               
               site.surv=site.surv,
               block.surv=block.surv,
               male.surv=male.surv,
               long.surv=long.surv,
               y.surv=y.surv,
               N.obs.surv=N.obs.surv,
               
               site.grow=site.grow,
               block.grow=block.grow,
               male.grow=male.grow,
               long.grow=long.grow,
               y.grow=y.grow,
               N.obs.grow=N.obs.grow,
               
               site.flow=site.flow,
               block.flow=block.flow,
               male.flow=male.flow,
               long.flow=long.flow,
               y.flow=y.flow,
               N.obs.flow=N.obs.flow,
               
               site.panic=site.panic,
               block.panic=block.panic,
               male.panic=male.panic,
               long.panic=long.panic,
               y.panic=y.panic,
               N.obs.panic=N.obs.panic,
               
               longitudes=long.seq,
               N.longitudes=length(long.seq))

inits<-function(){list(alpha.surv=rnorm(1,0,2),
                       sigma.surv.site=rlnorm(1),
                       sigma.surv.block=rlnorm(1),
                       beta.surv.sex=rnorm(1,0,2),
                       beta.surv.long=rnorm(1,0,2),
                       beta.surv.sex.long=rnorm(1,0,2),
                       
                       alpha.grow=rnorm(1,0,2),
                       sigma.grow.site=rlnorm(1),
                       sigma.grow.block=rlnorm(1),
                       beta.grow.sex=rnorm(1,0,2),
                       beta.grow.long=runif(1,-1,1),
                       beta.grow.sex.long=rnorm(1,0,2),
                       sigma.grow=rlnorm(1),
                       
                       alpha.flow=rnorm(1,0,2),
                       sigma.flow.site=rlnorm(1),
                       sigma.flow.block=rlnorm(1),
                       beta.flow.sex=rnorm(1,0,2),
                       beta.flow.long=rnorm(1,0,2),
                       beta.flow.sex.long=rnorm(1,0,2),
                       
                       alpha.panic=rnorm(1,0,2),
                       sigma.panic.site=rlnorm(1),
                       sigma.panic.block=rlnorm(1),
                       beta.panic.sex=rnorm(1,0,2),
                       beta.panic.long=rnorm(1,0,2),
                       beta.panic.sex.long=rnorm(1,0,2),
                       sigma.panic.ind=rlnorm(1))}

## Params to estimate
parameters<-c("alpha.surv","beta.surv.sex","beta.surv.long","beta.surv.sex.long",
              "alpha.grow","beta.grow.sex","beta.grow.long","beta.grow.sex.long",
              "alpha.flow","beta.flow.sex","beta.flow.long","beta.flow.sex.long",
              "alpha.panic","beta.panic.sex","beta.panic.long","beta.panic.sex.long",
              "surv.F.pred","surv.M.pred","grow.F.pred","grow.M.pred",
              "flow.F.pred","flow.M.pred","panic.F.pred","panic.M.pred",
              "panic.fit","panic.fit.new","sigma.panic.ind","alpha.panic.site.block")

## MCMC settings
ni<-20000
nb<-2000
nt<-10
nc<-3

## run JAGS
poar.out<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,model.file="POAR_vital_rates_unstructured.txt",
                      n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())
mcmcplot(poar.out,c("alpha.panic","beta.panic.sex","beta.panic.long","beta.panic.sex.long",
                    "panic.F.pred","panic.M.pred"))#,"sigma.panic.ind"))

## posterior predictive check for panicles
bpvalue <- mean(poar.out$BUGSoutput$sims.list$panic.fit.new >
                  poar.out$BUGSoutput$sims.list$panic.fit)
plot(poar.out$BUGSoutput$sims.list$panic.fit,
     poar.out$BUGSoutput$sims.list$panic.fit.new)
abline(0,1)
## graphics
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(poar$long.center[poar$Sex=="M"]+mean(poar$Longitude),jitter(poar$log_ratio[poar$Sex=="M"],factor=50),
     xlab="Longitude",ylab="Annual growth (log ratio)")
points(poar$long.center[poar$Sex=="F"]+mean(poar$Longitude)+0.15,jitter(poar$log_ratio[poar$Sex=="F"],factor=50),col="red")
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$grow.F.pred,col="red",lwd=3)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$grow.M.pred,lwd=3)

#plot(jitter(poar$long.scaled[poar$Sex=="M"],factor=5)+mean(poar$Longitude),jitter(poar$surv_t1[poar$Sex=="M"],factor=0.05),
#     xlab="Longitude scaled",ylab="# Flowers")
#points(jitter(poar$long.scaled[poar$Sex=="F"],factor=5)+mean(poar$Longitude),jitter(poar$surv_t1[poar$Sex=="F"],factor=0.05)-0.025,col="red")
surv.plot<-aggregate(y.surv~sex*long.center,data=poar.surv,FUN=mean)
plot(surv.plot$long.center[surv.plot$sex==1]+mean(poar$Longitude),
     surv.plot$y.surv[surv.plot$sex==1],pch=16,col="red",cex=2,
     xlab="Longitude",ylab="Pr(Survival)",ylim=c(0,1))
points(surv.plot$long.center[surv.plot$sex==2]+mean(poar$Longitude),
       surv.plot$y.surv[surv.plot$sex==2],pch=16,cex=2)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$surv.F.pred,col="red",lwd=3)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$surv.M.pred,lwd=3)

#plot(jitter(poar$long.scaled[poar$Sex=="M"],factor=5)+mean(poar$Longitude),jitter(poar$flow_t1[poar$Sex=="M"],factor=0.05),
#     xlab="Longitude scaled",ylab="# Flowers")
#points(jitter(poar$long.scaled[poar$Sex=="F"],factor=5)+mean(poar$Longitude),jitter(poar$flow_t1[poar$Sex=="F"],factor=0.05)-0.025,col="red")
flow.plot<-aggregate(y.flow~sex*long.center,data=poar.flow,FUN=mean)
plot(flow.plot$long.center[flow.plot$sex==1]+mean(poar$Longitude),
     flow.plot$y.flow[flow.plot$sex==1],pch=16,col="red",cex=2,
     xlab="Longitude",ylab="Pr(Flowering)",ylim=c(0,1))
points(flow.plot$long.center[flow.plot$sex==2]+mean(poar$Longitude),
       flow.plot$y.flow[flow.plot$sex==2],pch=16,cex=2)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$flow.F.pred,col="red",lwd=3)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$flow.M.pred,lwd=3)

plot(poar$long.center[poar$Sex=="M" & poar$flowerN_t1>0]+mean(poar$Longitude),
     jitter(poar$flowerN_t1[poar$Sex=="M" & poar$flowerN_t1>0],factor=5),
     xlab="Longitude",ylab="# Flowers",ylim=c(0,20))
points(poar$long.center[poar$Sex=="F" & poar$flowerN_t1>0]+mean(poar$Longitude)+0.15,
       jitter(poar$flowerN_t1[poar$Sex=="F" & poar$flowerN_t1>0],factor=5),col="red")
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$panic.F.pred,col="red",lwd=3)
lines(long.seq+mean(poar$Longitude),poar.out$BUGSoutput$mean$panic.M.pred,lwd=3)
