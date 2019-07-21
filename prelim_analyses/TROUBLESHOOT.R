## Fit size-dependent vital rates and build MPM
library(R2jags)
library(mcmcplots)
setwd("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/analysis")
poar<-read.csv("poar.csv")


sink("POAR_vital_rates_sizestructured_trouble.txt")
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
    beta.surv.size~dnorm(0,0.001)

    #2. Growth
    alpha.grow~dnorm(0,0.001)   ##hyperprior for growth
    sigma.grow.site~dunif(0,1000) ##site variance
    tau.grow.site<-1/(sigma.grow.site*sigma.grow.site)
    sigma.grow.block~dunif(0,1000) ##site variance
    tau.grow.block<-1/(sigma.grow.block*sigma.grow.block)
    beta.grow.sex~dnorm(0,0.001)
    beta.grow.long~dnorm(0,0.001)
    beta.grow.sex.long~dnorm(0,0.001)
    beta.grow.size~dnorm(0,0.001)
    alpha.gamma.grow~dunif(0,100)

    #3. Flowering
    alpha.flow~dnorm(0,0.001)   ##hyperprior for flowival
    sigma.flow.site~dunif(0,1000) ##site variance
    tau.flow.site<-1/(sigma.flow.site*sigma.flow.site)
    sigma.flow.block~dunif(0,1000) ##site variance
    tau.flow.block<-1/(sigma.flow.block*sigma.flow.block)
    beta.flow.sex~dnorm(0,0.001)
    beta.flow.long~dnorm(0,0.001)
    beta.flow.sex.long~dnorm(0,0.001)
    beta.flow.size~dnorm(0,0.001)
    
    #4. Panicles
    alpha.panic~dnorm(0,0.001)   ##hyperprior for panicival
    sigma.panic.site~dunif(0,1000) ##site variance
    tau.panic.site<-1/(sigma.panic.site*sigma.panic.site)
    sigma.panic.block~dunif(0,1000) ##site variance
    tau.panic.block<-1/(sigma.panic.block*sigma.panic.block)
    beta.panic.sex~dnorm(0,0.001)
    beta.panic.long~dnorm(0,0.001)
    beta.panic.sex.long~dnorm(0,0.001)
    beta.panic.size~dnorm(0,0.001)
    sigma.panic.ind~dunif(0,1000) ##individual random effect for overdispersion
    tau.panic.ind<-1/(sigma.panic.ind*sigma.panic.ind)

    # Assign random effects for sites and blocks - for panicles just block effects
    for(i in 1:N.sites){
    alpha.surv.site[i]~dnorm(alpha.surv,tau.surv.site)
    alpha.flow.site[i]~dnorm(alpha.flow,tau.flow.site)
    for(j in 1:N.blocks){
    alpha.surv.site.block[i,j]~dnorm(alpha.surv.site[i],tau.surv.block)
    alpha.flow.site.block[i,j]~dnorm(alpha.flow.site[i],tau.flow.block)
    }
    }
    ##for panicles just block effects
    for(j in 1:N.blocks){
      alpha.grow.block[j]~dnorm(alpha.grow,tau.grow.block)
      alpha.panic.block[j]~dnorm(alpha.panic,tau.panic.block)
    }

    #1. Survival Likelihood 
    for(i in 1:N.obs.surv){
    y.surv[i]~dbern(p.surv[i])
    logit(p.surv[i])<-alpha.surv.site.block[site.surv[i],block.surv[i]] + 
    beta.surv.sex*male.surv[i]+
    beta.surv.long*long.surv[i]+
    beta.surv.sex.long*male.surv[i]*long.surv[i]+
    beta.surv.size*size.surv[i]
    }
    #2. Growth Likelihood 
    for(i in 1:N.obs.grow){
    y.grow[i]~dpois(mustar.grow[i])T(1,)
    mustar.grow[i]<-mu.grow[i]*rho.grow[i]
    rho.grow[i]~dgamma(alpha.gamma.grow,alpha.gamma.grow)
    log(mu.grow[i])<-alpha.grow.block[block.grow[i]] + 
    beta.grow.sex*male.grow[i]+
    beta.grow.long*long.grow[i]+
    beta.grow.sex.long*male.grow[i]*long.grow[i]+
    beta.grow.size*size.grow[i]
    }
    #3. Flowering Likelihood 
    for(i in 1:N.obs.flow){
    y.flow[i]~dbern(p.flow[i])
    logit(p.flow[i])<-alpha.flow.site.block[site.flow[i],block.flow[i]] + 
    beta.flow.sex*male.flow[i]+
    beta.flow.long*long.flow[i]+
    beta.flow.sex.long*male.flow[i]*long.flow[i]+
    beta.flow.size*size.flow[i]
    }
    #4. Panicles Likelihood 
    for(i in 1:N.obs.panic){
    y.panic[i]~dpois(mu.panic[i])T(1,)
    ##try dropping site effects here, just block effects
    log(mu.panic[i])<-alpha.panic.block[block.panic[i]] + 
    beta.panic.sex*male.panic[i]+
    beta.panic.long*long.panic[i]+
    beta.panic.sex.long*male.panic[i]*long.panic[i]+
    beta.panic.size*size.panic[i]+rho[i]
    rho[i]~dnorm(0,tau.panic.ind)
    }

    ## Prediction
    for(i in 1:N.sizes){
    logit(surv.F.W.pred[i])<-alpha.surv+beta.surv.size*size.pred[i]+beta.surv.long*W.long
    logit(surv.M.W.pred[i])<-alpha.surv+beta.surv.sex+beta.surv.size*size.pred[i]+(beta.surv.long+beta.surv.sex.long)*W.long
    logit(surv.F.E.pred[i])<-alpha.surv+beta.surv.size*size.pred[i]+beta.surv.long*E.long
    logit(surv.M.E.pred[i])<-alpha.surv+beta.surv.sex+beta.surv.size*size.pred[i]+(beta.surv.long+beta.surv.sex.long)*E.long

    grow.F.W.pred[i]<-exp(alpha.grow+beta.grow.size*size.pred[i]+beta.grow.long*W.long)
    grow.M.W.pred[i]<-exp(alpha.grow+beta.grow.sex+beta.grow.size*size.pred[i]+(beta.grow.long+beta.grow.sex.long)*W.long)
    grow.F.E.pred[i]<-exp(alpha.grow+beta.grow.size*size.pred[i]+beta.grow.long*E.long)
    grow.M.E.pred[i]<-exp(alpha.grow+beta.grow.sex+beta.grow.size*size.pred[i]+(beta.grow.long+beta.grow.sex.long)*E.long)

    logit(flow.F.W.pred[i])<-alpha.flow+beta.flow.size*size.pred[i]+beta.flow.long*W.long
    logit(flow.M.W.pred[i])<-alpha.flow+beta.flow.sex+beta.flow.size*size.pred[i]+(beta.flow.long+beta.flow.sex.long)*W.long
    logit(flow.F.E.pred[i])<-alpha.flow+beta.flow.size*size.pred[i]+beta.flow.long*E.long
    logit(flow.M.E.pred[i])<-alpha.flow+beta.flow.sex+beta.flow.size*size.pred[i]+(beta.flow.long+beta.flow.sex.long)*E.long
    
    log(panic.F.W.pred[i])<-alpha.panic+beta.panic.size*size.pred[i]+beta.panic.long*W.long
    log(panic.M.W.pred[i])<-alpha.panic+beta.panic.sex+beta.panic.size*size.pred[i]+(beta.panic.long+beta.panic.sex.long)*W.long
    log(panic.F.E.pred[i])<-alpha.panic+beta.panic.size*size.pred[i]+beta.panic.long*E.long
    log(panic.M.E.pred[i])<-alpha.panic+beta.panic.sex+beta.panic.size*size.pred[i]+(beta.panic.long+beta.panic.sex.long)*E.long
    }
    }##end model
    ",fill=T)
sink()

## Posterior predictive check
# Growth
for(i in 1:N.obs.grow){
  grow.resi[i] <- (y.grow[i] - mustar.grow[i])/sqrt(mustar.grow[i])
  y.grow.new[i] ~ dpois(mustar.grow[i])
  grow.resi.new[i] <- (y.grow.new[i] - mustar.grow[i])/sqrt(mustar.grow[i])
  grow.D[i] <- pow(grow.resi[i],2)
  grow.D.new[i] <- pow(grow.resi.new[i],2)
}
grow.fit <- sum(grow.D[])
grow.fit.new <- sum(grow.D.new[])

# Panicle
for(i in 1:N.obs.panic){
  panic.resi[i] <- (y.panic[i] - mu.panic[i])/sqrt(mu.panic[i])
  y.panic.new[i] ~ dpois(mu.panic[i])
  panic.resi.new[i] <- (y.panic.new[i] - mu.panic[i])/sqrt(mu.panic[i])
  panic.D[i] <- pow(panic.resi[i],2)
  panic.D.new[i] <- pow(panic.resi.new[i],2)
}
panic.fit <- sum(panic.D[])
panic.fit.new <- sum(panic.D.new[])

#### Bundle data
#1. Survival
poar.surv<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t0,poar$surv_t1)))
names(poar.surv)<-c("site","block","sex","long.center","long.scaled","tillerN_t0","surv_t1")
poar.surv<-subset(poar.surv,tillerN_t0>0)
N.sites<-max(poar.surv$site)
site.surv<-poar.surv$site
N.blocks<-max(poar.surv$block)
block.surv<-poar.surv$block
male.surv<-poar.surv$sex-1
long.surv<-poar.surv$long.center
y.surv<-poar.surv$surv_t1
size.surv<-log(poar.surv$tillerN_t0)
N.obs.surv<-nrow(poar.surv)
par(mfrow=c(2,2));plot(size.surv,y.surv)
#2. Growth
poar.grow<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t0,poar$tillerN_t1)))
names(poar.grow)<-c("site","block","sex","long.center","long.scaled","tillerN_t0","tillerN_t1")
max(poar.surv$site)==max(poar.grow$site)
max(poar.surv$block)==max(poar.grow$block)## sites and blocks line up
poar.grow<-subset(poar.grow,tillerN_t0>0 & tillerN_t1>0)
site.grow<-poar.grow$site
block.grow<-poar.grow$block
male.grow<-poar.grow$sex-1
long.grow<-poar.grow$long.center
size.grow<-log(poar.grow$tillerN_t0)
y.grow<-poar.grow$tillerN_t1
N.obs.grow<-nrow(poar.grow)
plot(size.grow,y.grow)
#3. Flowering
poar.flow<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t1,poar$flow_t1)))
names(poar.flow)<-c("site","block","sex","long.center","long.scaled","tillerN_t1","flow_t1")
max(poar.surv$site)==max(poar.flow$site)
max(poar.surv$block)==max(poar.flow$block)## sites and blocks line up
poar.flow<-subset(poar.flow,tillerN_t1>0)
site.flow<-poar.flow$site
block.flow<-poar.flow$block
male.flow<-poar.flow$sex-1
long.flow<-poar.flow$long.center
size.flow<-log(poar.flow$tillerN_t1)
y.flow<-poar.flow$flow_t1
N.obs.flow<-nrow(poar.flow)
plot(size.flow,y.flow)
#4. Panicles
poar.panic<-na.omit(data.frame(cbind(poar$site,poar$unique.block,poar$Sex,poar$long.center,poar$long.scaled,poar$tillerN_t1,poar$flowerN_t1)))
names(poar.panic)<-c("site","block","sex","long.center","long.scaled","tillerN_t1","panic_t1")
poar.panic<-subset(poar.panic,panic_t1>0 & tillerN_t1>0)
max(poar.surv$site)==max(poar.panic$site)
max(poar.surv$block)==max(poar.panic$block)## sites and blocks line up
site.panic<-poar.panic$site
block.panic<-poar.panic$block
male.panic<-poar.panic$sex-1
long.panic<-poar.panic$long.center
size.panic<-log(poar.panic$tillerN_t1)
y.panic<-poar.panic$panic_t1
N.obs.panic<-nrow(poar.panic)
plot(size.panic,y.panic)
## Misc data for prediction
size.pred<-seq(0,max(log(poar$tillerN_t1),na.rm=T),length.out = 30)
W.long<-min(poar$long.center)
E.long<-max(poar$long.center)


jag.data<-list(N.sites=N.sites,
               N.blocks=N.blocks,
               
               site.surv=site.surv,
               block.surv=block.surv,
               male.surv=male.surv,
               long.surv=long.surv,
               size.surv=size.surv,
               y.surv=y.surv,
               N.obs.surv=N.obs.surv,

               site.grow=site.grow,
               block.grow=block.grow,
               male.grow=male.grow,
               long.grow=long.grow,
               size.grow=size.grow,
               y.grow=y.grow,
               N.obs.grow=N.obs.grow,
               
               site.flow=site.flow,
               block.flow=block.flow,
               male.flow=male.flow,
               long.flow=long.flow,
               size.flow=size.flow,
               y.flow=y.flow,
               N.obs.flow=N.obs.flow,
               
               site.panic=site.panic,
               block.panic=block.panic,
               male.panic=male.panic,
               long.panic=long.panic,
               y.panic=y.panic,
               size.panic=size.panic,
               N.obs.panic=N.obs.panic,
               
               size.pred=size.pred,
               N.sizes=length(size.pred),
               W.long=W.long,
               E.long=E.long)

inits<-function(){list(alpha.surv=rnorm(1,0,1),
                       sigma.surv.site=rlnorm(1),
                       sigma.surv.block=rlnorm(1),
                       beta.surv.sex=rnorm(1,0,1),
                       beta.surv.long=rnorm(1,0,1),
                       beta.surv.sex.long=rnorm(1,0,1),
                       beta.surv.size=rnorm(1,0,1),
                       
                       alpha.grow=rnorm(1,0,1),
                       sigma.grow.site=rlnorm(1),
                       sigma.grow.block=rlnorm(1),
                       beta.grow.sex=rnorm(1,0,1),
                       beta.grow.long=rnorm(1,0,1),
                       beta.grow.sex.long=rnorm(1,0,1),
                       beta.grow.size=rnorm(1,0,1),
                       alpha.gamma.grow=rlnorm(1),
                       
                       alpha.flow=rnorm(1,0,1),
                       sigma.flow.site=rlnorm(1),
                       sigma.flow.block=rlnorm(1),
                       beta.flow.sex=rnorm(1,0,1),
                       beta.flow.long=rnorm(1,0,1),
                       beta.flow.sex.long=rnorm(1,0,1),
                       beta.flow.size=rnorm(1,0,1),
                       
                       alpha.panic=rnorm(1,0,1),
                       sigma.panic.site=rlnorm(1),
                       sigma.panic.block=rlnorm(1),
                       beta.panic.sex=rnorm(1,0,1),
                       beta.panic.long=rnorm(1,0,1),
                       beta.panic.sex.long=rnorm(1,0,1),
                       beta.panic.size=rnorm(1,0,1),
                       sigma.panic.ind=rlnorm(1))}

## Params to estimate
parameters<-c("alpha.surv","beta.surv.sex","beta.surv.long","beta.surv.sex.long","beta.surv.size",
              "alpha.grow","beta.grow.sex","beta.grow.long","beta.grow.sex.long","beta.grow.size","alpha.gamma.grow",
              "alpha.flow","beta.flow.sex","beta.flow.long","beta.flow.sex.long","beta.flow.size",
              "alpha.panic","beta.panic.sex","beta.panic.long","beta.panic.sex.long","beta.panic.size",
              "surv.F.W.pred","surv.M.W.pred","surv.F.E.pred","surv.M.E.pred",
              "grow.F.W.pred","grow.M.W.pred","grow.F.E.pred","grow.M.E.pred",
              "flow.F.W.pred","flow.M.W.pred","flow.F.E.pred","flow.M.E.pred",
              "panic.F.W.pred","panic.M.W.pred","panic.F.E.pred","panic.M.E.pred")
#,"panic.fit","panic.fit.new","grow.fit","grow.fit.new"

## MCMC settings
ni<-10000
nb<-1000
nt<-10
nc<-3

## run JAGS
poar.trouble.out<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,model.file="POAR_vital_rates_sizestructured_trouble.txt",
               n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())
mcmcplot(poar.trouble.out,c("alpha.panic","beta.panic.sex","beta.panic.long","beta.panic.sex.long","beta.panic.size",
                            "alpha.grow","beta.grow.sex","beta.grow.long","beta.grow.sex.long","beta.grow.size","alpha.gamma.grow"))

par(mfrow=c(2,2))
plot(size.surv,y.surv)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$surv.F.W.pred,col="red")
lines(size.pred,poar.trouble.out$BUGSoutput$mean$surv.M.W.pred)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$surv.F.E.pred,col="red",lty=2)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$surv.M.E.pred,lty=2)

plot(size.grow,y.grow)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$grow.F.W.pred,col="red")
lines(size.pred,poar.trouble.out$BUGSoutput$mean$grow.M.W.pred)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$grow.F.E.pred,col="red",lty=2)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$grow.M.E.pred,lty=2)

plot(size.flow,y.flow)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$flow.F.W.pred,col="red")
lines(size.pred,poar.trouble.out$BUGSoutput$mean$flow.M.W.pred)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$flow.F.E.pred,col="red",lty=2)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$flow.M.E.pred,lty=2)

plot(size.panic,y.panic)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$panic.F.W.pred,col="red")
lines(size.pred,poar.trouble.out$BUGSoutput$mean$panic.M.W.pred)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$panic.F.E.pred,col="red",lty=2)
lines(size.pred,poar.trouble.out$BUGSoutput$mean$panic.M.E.pred,lty=2)

