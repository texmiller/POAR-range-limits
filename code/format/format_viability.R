# Format viability data
library(scales)
rl_dir <- "C:/cloud/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/"

# germination data from response surface
germ   <- read.csv( paste0( rs_dir, "Data/Spring 2014/viability/germination_data_LLELA.csv") )
viab   <- read.csv( paste0( rs_dir, "Data/Spring 2014/viability/viabilityData.csv") )
vr     <- read.csv( paste0( rs_dir, "Analysis1/fallSpringLLELA.csv") )

#1. format Germination 
germ$germTot=apply(germ[,grep("Census",names(germ))],1,sum)
#sum up regrowth as a check - this need be less than germTot
germ$germRegrowthTot=apply(germ[,grep("Re.growth",names(germ))],1,sum) 
sum(germ$germRegrowthTot > germ$germTot,na.rm=T) #this should be 0
germ$germFail=25-germ$germTot
names(germ)[1]="plot"

#2. Format tetrazolium data
## First, remove rows that contain no data 
rI=which(is.na(viab$Yes) & is.na(viab$No) & is.na(viab$Maybe) & is.na(viab$Brown))
viab=viab[-rI,]
viab$Brown[is.na(viab$Brown)]=0
viab$totS=viab$Yes + viab$No + viab$Maybe + viab$Brown  #Total seeds (for binomial model)
viab$yesMaybe=viab$Yes + viab$Maybe               #Conservative way to identify "viable seeds"
viab$fail=viab$totS-viab$Yes                      #"Fail" column for the glm() analyses.
viab$failMaybe=viab$totS-viab$yesMaybe
names(viab)[1]="plot" 

#3. Format plot data
plotD=unique(vr[,c("plot","F","M","F_flow","M_flow")])
#Calculate total number of flowers in each plot 
focal_f_flow=aggregate(flower ~ plot, sum, data=subset(vr,sex=="f")) ; names(focal_f_flow)[2]="F_flow_f"
focal_m_flow=aggregate(flower ~ plot, sum, data=subset(vr,sex=="m")) ; names(focal_m_flow)[2]="M_flow_f"
ch=merge(plotD,focal_f_flow,all=T) #Merge data 
ch=merge(ch,focal_m_flow,all=T)
ch$F_flow_f[is.na(ch$F_flow_f)]=0 #Substitute NAs with zeros
ch$M_flow_f[is.na(ch$M_flow_f)]=0
#IF #1. F or M are <=5, then #2. Substitute #_flow with #_flow_f  
ch[which(ch$F<=5),]$F_flow = ch[which(ch$F<=5),]$F_flow_f 
ch[which(ch$M<=5),]$M_flow = ch[which(ch$M<=5),]$M_flow_f 
#Calculate sex ratios
ch$totFlow=ch$M_flow + ch$F_flow #Tot number of flowers
ch$sr_f=ch$F_flow/ch$totFlow 
ch$sr_f2=ch$sr_f ^ 2

#Merge the three datasets
tmp=merge(viab,germ[,c("plot","F_ID","P_ID","germTot","germFail")],all=T)
tmp$focalI=matrix(unlist(strsplit(as.character(tmp$F_ID),"F-")),nrow(tmp),2,byrow=T)[,2]
tmp$focalI=paste0("f",as.numeric(tmp$focalI))
viabVr=merge(ch,tmp) #Final file

#Calculate viability/germination ratios
viabVr$tetra_ratio        <- viabVr$Yes / viabVr$totS 
viabVr$tetra_maybe_ratio  <- viabVr$yesMaybe / viabVr$totS
viabVr$germ_ratio         <- viabVr$germTot / (viabVr$germTot + viabVr$germFail)

FINAL.germ<-glmer(cbind(germTot,germFail) ~ sr_f + (1|plot),family="binomial", data=viabVr)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,
     xlab="Sex ratio",ylab="Germination rate")
yMean=invlogit(fixef(FINAL.germ)[1]+fixef(FINAL.germ)[2]*xSeq)
lines(xSeq,yMean,lwd=2,col="black")


# write result out -------------------------------------------

write.csv(viabVr, 
          'C:/CODE/POAR-range-limits/data/viability.csv', 
          row.names=F)
