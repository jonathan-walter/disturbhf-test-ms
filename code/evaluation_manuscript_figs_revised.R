## Interpret evaluations of all simulation tests

rm(list=ls())

library(here)
library(viridis)

## <------------------------------- fix file paths

## load simulation output
sim.flat.negwedge<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_negwedge_20210323_t5.rds")
sim.flat.poswedge<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_poswedge_20210323_t5.rds")
sim.flat.negs<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_negs_20210323_t5.rds")
sim.flat.poss<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_poss_20210323_t5.rds")
sim.sine.negwedge<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_20210323_t5.rds")
sim.sine.poswedge<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_poswedge_20210323_t5.rds")
sim.sine.negs<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negs_20210323_t5.rds")
sim.sine.poss<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_poss_20210323_t5.rds")

sim.flat.negwedge2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_negwedge_adapt_20210323_t5.rds")
sim.flat.poswedge2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_poswedge_adapt_20210323_t5.rds")
sim.flat.negs2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_negs_adapt_20210323_t5.rds")
sim.flat.poss2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_flat_poss_adapt_20210323_t5.rds")
sim.sine.negwedge2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_adapt_20210323_t5.rds")
sim.sine.poswedge2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_poswedge_adapt_20210323_t5.rds")
sim.sine.negs2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negs_adapt_20210323_t5.rds")
sim.sine.poss2<-readRDS("~/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_poss_adapt_20210323_t5.rds")

##-----------------------------------------------------------------------------
## helper functions

## Extract evaluation components 

getEval<-function(sims, which.eval){
  keys<-NULL
  evals<-NULL
  for(ii in 1:length(sims)){
    key.ii<-sims[[ii]]$key
    if(which.eval=="normal"){
      eval.ii<-sims[[ii]]$eval
      evals<-rbind(evals, eval.ii)
    }
    if(which.eval=="filtered"){
      eval.ii<-sims[[ii]]$filteval
      evals<-rbind(evals, eval.ii)
    }
    key.ii<-cbind(rep(ii, nrow(eval.ii)),
                  matrix(key.ii, nrow=nrow(eval.ii), ncol=length(key.ii),byrow=T))
    key.ii<-as.data.frame(key.ii)
    colnames(key.ii)<-c("simrep","dt","ref.end","dday","rday","recovtime","severity")
    keys<-rbind(keys, key.ii)
  }
  return(cbind(keys,evals))
}

## Function for cleaning evals

cleanEval<-function(evals){
  sim<-unique(evals$simrep)
  recovtime<-rep(NA, length(sim))
  severity<-rep(NA, length(sim))
  detect<-rep(NA, length(sim))
  
  for(ii in sim){
    tmp<-evals[evals$simrep == ii,]
    recovtime[ii]<-tmp$recovtime[1]
    severity[ii]<-tmp$severity[1]
    detect[ii]<-max(tmp$detect)
  }
  
  return(data.frame(simrep=sim, recovtime=recovtime, severity=severity, detect=detect))
  
}


## Function to create and plot success as a function of disturbance 
## severity and duration

tprSurface<-function(evals, along.s, along.r, plotit=TRUE, title=NULL){
  out<-matrix(NA, nrow=length(along.s)-1, ncol=length(along.r)-1)
  for(ii in 1:(length(along.s)-1)){
    for(jj in 1:(length(along.r)-1)){
      out[ii,jj]<-mean(evals$detect[
        evals$severity >= along.s[ii] & evals$severity < along.s[ii+1] &
          evals$recovtime >= along.r[jj] & evals$recovtime < along.r[jj+1]], 
        na.rm=T
      )
    }
  }
  if(plotit){
    
    layout(matrix(c(1,2),ncol=2),widths=c(0.8,0.2))
    par(mar=c(4.1,4.1,2.1,0.1))
    image(x=along.s, y=along.r, z=out, xlab="Severity", ylab="Duration",zlim=c(0,1),col=viridis(25))
    if(!is.null(title)){
      mtext(title)
    }
    par(mar=c(4.1,4.1,2.1,1))
    image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",ylab="True positive rate")
    
    layout(matrix(1))
  }
  return(out)
}

## Function to create and plot number of false alarms as function of disturbance 
## severity and duration

falmSurface<-function(evals, along.s, along.r, plotit=TRUE, title=NULL, zmax=NULL, denom=1){
  out<-matrix(NA, nrow=length(along.s)-1, ncol=length(along.r)-1)
  
  nfails<-rep(NA, max(evals$simrep))
  severity<-rep(NA, max(evals$simrep))
  recovtime<-rep(NA, max(evals$simrep))
  for(nn in 1:max(evals$simrep)){
    nfails[nn]<-sum(evals$simrep==nn & evals$detect==0)/denom
    severity[nn]<-evals$severity[evals$simrep==nn][1]
    recovtime[nn]<-evals$recovtime[evals$simrep==nn][1]
  }
  
  for(ii in 1:(length(along.s)-1)){
    for(jj in 1:(length(along.r)-1)){
      
      out[ii,jj]<-mean(nfails[
        severity >= along.s[ii] & severity < along.s[ii+1] &
          recovtime >= along.r[jj] & recovtime < along.r[jj+1]], 
        na.rm=T
      )
    }
  }
  if(plotit){
    
    if(is.null(zmax)){
      zmax<-ceiling(max(out))
    }
    
    
    layout(matrix(c(1,2),ncol=2),widths=c(0.8,0.2))
    par(mar=c(4.1,4.1,2.1,0.1))
    image(x=along.s, y=along.r, z=out, xlab="Severity", ylab="Duration",zlim=c(0,zmax),col=viridis(25))
    if(!is.null(title)){
      mtext(title)
    }
    par(mar=c(4.1,4.1,2.1,1))
    image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",yaxt="n", ylab="False positive rate")
    axis(2,at=seq(0,1,length.out=4),labels=round(seq(0,zmax,length.out=4),1))
    
    layout(matrix(1))
  }
  return(out)
}


errorSurface<-function(evals, along.s, along.r, plotit=TRUE, title=NULL){
  out<-matrix(NA, nrow=length(along.s)-1, ncol=length(along.r)-1)
  for(ii in 1:(length(along.s)-1)){
    for(jj in 1:(length(along.r)-1)){
      out[ii,jj]<-median(evals$delta.recov[
        evals$severity >= along.s[ii] & evals$severity < along.s[ii+1] &
          evals$recovtime >= along.r[jj] & evals$recovtime < along.r[jj+1]], 
        na.rm=T
      )
    }
  }
  
  if(plotit){
    zmax<-ceiling(max(out))
    layout(matrix(c(1,2),ncol=2),widths=c(0.8,0.2))
    par(mar=c(4.1,4.1,2.1,0.1))
    image(x=along.s, y=along.r, z=out, xlab="Severity", ylab="Duration",col=viridis(25))
    if(!is.null(title)){
      mtext(title)
    }
    par(mar=c(4.1,4.1,2.1,1))
    image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",yaxt="n",ylab="Error in recovery date (true-estimated)")
    axis(2,at=seq(0,1,length.out=4),labels=round(seq(floor(min(out)),zmax,length.out=4),1))
    layout(matrix(1))
  }
  return(out)
}


eval.flat.negwedge<-getEval(sim.flat.negwedge, which.eval = "normal")
eval.flat.poswedge<-getEval(sim.flat.poswedge, which.eval = "normal")
eval.flat.negs<-getEval(sim.flat.negs, which.eval="normal")
eval.flat.poss<-getEval(sim.flat.poss, which.eval="normal")
eval.sine.negwedge<-getEval(sim.sine.negwedge, which.eval = "normal")
eval.sine.poswedge<-getEval(sim.sine.poswedge, which.eval = "normal")
eval.sine.negs<-getEval(sim.sine.negs, which.eval="normal")
eval.sine.poss<-getEval(sim.sine.poss, which.eval="normal")

eval.flat.negwedge.adapt<-getEval(sim.flat.negwedge2, which.eval = "normal")
eval.flat.poswedge.adapt<-getEval(sim.flat.poswedge2, which.eval = "normal")
eval.flat.negs.adapt<-getEval(sim.flat.negs2, which.eval="normal")
eval.flat.poss.adapt<-getEval(sim.flat.poss2, which.eval="normal")
eval.sine.negwedge.adapt<-getEval(sim.sine.negwedge2, which.eval = "normal")
eval.sine.poswedge.adapt<-getEval(sim.sine.poswedge2, which.eval = "normal")
eval.sine.negs.adapt<-getEval(sim.sine.negs2, which.eval="normal")
eval.sine.poss.adapt<-getEval(sim.sine.poss2, which.eval="normal")


cln.flat.negwedge<-cleanEval(eval.flat.negwedge)
cln.flat.poswedge<-cleanEval(eval.flat.poswedge)
cln.flat.negs<-cleanEval(eval.flat.negs)
cln.flat.poss<-cleanEval(eval.flat.poss)
cln.sine.negwedge<-cleanEval(eval.sine.negwedge)
cln.sine.poswedge<-cleanEval(eval.sine.poswedge)
cln.sine.negs<-cleanEval(eval.sine.negs)
cln.sine.poss<-cleanEval(eval.sine.poss)

cln.flat.negwedge.adapt<-cleanEval(eval.flat.negwedge.adapt)
cln.flat.poswedge.adapt<-cleanEval(eval.flat.poswedge.adapt)
cln.flat.negs.adapt<-cleanEval(eval.flat.negs.adapt)
cln.flat.poss.adapt<-cleanEval(eval.flat.poss.adapt)
cln.sine.negwedge.adapt<-cleanEval(eval.sine.negwedge.adapt)
cln.sine.poswedge.adapt<-cleanEval(eval.sine.poswedge.adapt)
cln.sine.negs.adapt<-cleanEval(eval.sine.negs.adapt)
cln.sine.poss.adapt<-cleanEval(eval.sine.poss.adapt)



tpr.flat.negwedge<-tprSurface(cln.flat.negwedge,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.poswedge<-tprSurface(cln.flat.poswedge,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.negs<-tprSurface(cln.flat.negs,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.poss<-tprSurface(cln.flat.poss,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.negwedge<-tprSurface(cln.sine.negwedge,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.poswedge<-tprSurface(cln.sine.poswedge,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.negs<-tprSurface(cln.sine.negs,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.poss<-tprSurface(cln.sine.poss,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))

tpr.flat.negwedge.adapt<-tprSurface(cln.flat.negwedge.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.poswedge.adapt<-tprSurface(cln.flat.poswedge.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.negs.adapt<-tprSurface(cln.flat.negs.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.flat.poss.adapt<-tprSurface(cln.flat.poss.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.negwedge.adapt<-tprSurface(cln.sine.negwedge.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.poswedge.adapt<-tprSurface(cln.sine.poswedge.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.negs.adapt<-tprSurface(cln.sine.negs.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
tpr.sine.poss.adapt<-tprSurface(cln.sine.poss.adapt,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))


denom <- sum(!is.na(sim.sine.negwedge[[1]]$empdiff))

fpr.flat.negwedge<-falmSurface(eval.flat.negwedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.poswedge<-falmSurface(eval.flat.poswedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.negs<-falmSurface(eval.flat.negs, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.poss<-falmSurface(eval.flat.poss, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.negwedge<-falmSurface(eval.sine.negwedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.poswedge<-falmSurface(eval.sine.poswedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.negs<-falmSurface(eval.sine.negs, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.poss<-falmSurface(eval.sine.poss, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

fpr.flat.negwedge.adapt<-falmSurface(eval.flat.negwedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.poswedge.adapt<-falmSurface(eval.flat.poswedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.negs.adapt<-falmSurface(eval.flat.negs.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.flat.poss.adapt<-falmSurface(eval.flat.poss.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.negwedge.adapt<-falmSurface(eval.sine.negwedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.poswedge.adapt<-falmSurface(eval.sine.poswedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.negs.adapt<-falmSurface(eval.sine.negs.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)
fpr.sine.poss.adapt<-falmSurface(eval.sine.poss.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)


hits.flat.negwedge<-eval.flat.negwedge[eval.flat.negwedge$detect==1,]
hits.flat.negs<-eval.flat.negs[eval.flat.negs$detect==1,]
hits.flat.negwedge.adapt<-eval.flat.negwedge.adapt[eval.flat.negwedge.adapt$detect==1,]
hits.flat.negs.adapt<-eval.flat.negs.adapt[eval.flat.negs.adapt$detect==1,]
hits.sine.negwedge<-eval.sine.negwedge[eval.sine.negwedge$detect==1,]
hits.sine.negs<-eval.sine.negs[eval.sine.negs$detect==1,]
hits.sine.negwedge.adapt<-eval.sine.negwedge.adapt[eval.sine.negwedge.adapt$detect==1,]
hits.sine.negs.adapt<-eval.sine.negs.adapt[eval.sine.negs.adapt$detect==1,]

hits.flat.poswedge<-eval.flat.poswedge[eval.flat.poswedge$detect==1,]
hits.flat.poss<-eval.flat.poss[eval.flat.poss$detect==1,]
hits.flat.poswedge.adapt<-eval.flat.poswedge.adapt[eval.flat.poswedge.adapt$detect==1,]
hits.flat.poss.adapt<-eval.flat.poss.adapt[eval.flat.poss.adapt$detect==1,]
hits.sine.poswedge<-eval.sine.poswedge[eval.sine.poswedge$detect==1,]
hits.sine.poss<-eval.sine.poss[eval.sine.poss$detect==1,]
hits.sine.poswedge.adapt<-eval.sine.poswedge.adapt[eval.sine.poswedge.adapt$detect==1,]
hits.sine.poss.adapt<-eval.sine.poss.adapt[eval.sine.poss.adapt$detect==1,]


error.flat.negwedge<-errorSurface(hits.flat.negwedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.poswedge<-errorSurface(hits.flat.poswedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.negs<-errorSurface(hits.flat.negs, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.poss<-errorSurface(hits.flat.poss, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.negwedge<-errorSurface(hits.sine.negwedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.poswedge<-errorSurface(hits.sine.poswedge, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.negs<-errorSurface(hits.sine.negs, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.poss<-errorSurface(hits.sine.poss, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))

error.flat.negwedge.adapt<-errorSurface(hits.flat.negwedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.poswedge.adapt<-errorSurface(hits.flat.poswedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.negs.adapt<-errorSurface(hits.flat.negs.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.flat.poss.adapt<-errorSurface(hits.flat.poss.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.negwedge.adapt<-errorSurface(hits.sine.negwedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.poswedge.adapt<-errorSurface(hits.sine.poswedge.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.negs.adapt<-errorSurface(hits.sine.negs.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
error.sine.poss.adapt<-errorSurface(hits.sine.poss.adapt, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))

pal.tpr<-colorRampPalette(c("#e0ecf4","#9ebcda","#8856a7"))
pal.fpr<-colorRampPalette(c("#e0f3db","#a8ddb5","#43a2ca"))

laymat4<-matrix(c(1:8,10,10,9,9),4,3)
along.s<-seq(0,5,by=0.5)
along.r<-seq(3,60,by=5)




#Figure 3: TPR for negative disturbances
summary(c(tpr.flat.negwedge))
summary(c(tpr.flat.negs))
summary(c(tpr.sine.negwedge))
summary(c(tpr.sine.negs))
summary(c(tpr.flat.negwedge.adapt))
summary(c(tpr.flat.negs.adapt))
summary(c(tpr.sine.negwedge.adapt))
summary(c(tpr.sine.negs.adapt))

summary(c(tpr.sine.negs.adapt-tpr.sine.negs))
summary(c(tpr.sine.negwedge.adapt-tpr.sine.negwedge))

# jointmax<-ceiling(max(c(fpr.flat.negwedge,fpr.flat.negs,fpr.flat.negwedge.adapt,fpr.flat.negs.adapt)))

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/fig3_tpr_negative.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
image(along.s, along.r, tpr.flat.negwedge, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("Wedge disturbance\nFlat background",2,cex=0.9,line=0.25,at=7/8,outer=T)
image(along.s, along.r, tpr.sine.negwedge, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("Wedge disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=5/8,outer=T)
mtext("S disturbance\nFlat background",2,cex=0.9,line=0.25,at=3/8,outer=T)
image(along.s, along.r, tpr.flat.negs, col=rev(viridis(25)), zlim=c(0,1))
mtext("S disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=1/8,outer=T)
image(along.s, along.r, tpr.sine.negs, col=rev(viridis(25)), zlim=c(0,1))

image(along.s, along.r, tpr.flat.negwedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
image(along.s, along.r, tpr.flat.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))
image(along.s, along.r, tpr.sine.negwedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
image(along.s, along.r, tpr.sine.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=rev(viridis(25)),xaxt="n",ylab="True detection rate")

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=3.5)
mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()


#Figure 4: error for recovery date surfaces-------------------------------------------------------

joint.min<-min(c(error.flat.negwedge, error.flat.negs, error.sine.negs, error.sine.negwedge,
                 error.flat.negwedge.adapt, error.flat.negs.adapt, error.sine.negs.adapt, error.sine.negwedge.adapt))
joint.max<-max(c(error.flat.negwedge, error.flat.negs, error.sine.negs, error.sine.negwedge,
                 error.flat.negwedge.adapt, error.flat.negs.adapt, error.sine.negs.adapt, error.sine.negwedge.adapt))

n=25
# errorRamp <- colorRampPalette(colors=c("red","lightgrey","blue"))
# errorPal <- errorRamp(n)
errorPal <- turbo(n)
plotScale <- function(x, joint.max){
  xx <- x + joint.max
  xx <- xx/max(xx)
  xx <- round(xx*n)
  return(xx)
}

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/fig4_recoverror_negative.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
image(along.s, along.r, plotScale(error.flat.negwedge,joint.max), col=errorPal, 
      ylab="", zlim=c(1,n))
mtext("Wedge disturbance\nFlat background",2,cex=0.9,line=0.25,at=7/8,outer=T)
image(along.s, along.r, plotScale(error.sine.negwedge,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
mtext("Wedge disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=5/8,outer=T)
mtext("S disturbance\nFlat background",2,cex=0.9,line=0.25,at=3/8,outer=T)
image(along.s, along.r, plotScale(error.flat.negs,joint.max), col=errorPal, 
      zlim=c(1,n))
mtext("S disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=1/8,outer=T)
image(along.s, along.r, plotScale(error.sine.negs,joint.max), col=errorPal, zlim=c(1,n))

image(along.s, along.r, plotScale(error.flat.negwedge.adapt,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
image(along.s, along.r, plotScale(error.flat.negs.adapt,joint.max), col=errorPal, 
      zlim=c(1,n))
image(along.s, along.r, plotScale(error.sine.negwedge.adapt,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
image(along.s, along.r, plotScale(error.sine.negs.adapt,joint.max), col=errorPal,
      zlim=c(1,n))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:51)),col=errorPal,xaxt="n",yaxt="n",ylab="Recovery date error (true-estimated)")
axis(2,at=seq(0,1,length.out=5),labels=round(seq(floor(-1*joint.max),ceiling(joint.max),length.out=5),1))

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=3.5)
mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()


#Figure 3b: error for recovery date histograms ----------------------------------------------------


png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/hist_recoverror_negative.png",
    units="in", res=300, width=6.5,height=8)

par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0),mfcol=c(4,2))
hist(hits.flat.negwedge$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("Wedge disturbance",2,cex=0.9,line=0.25,at=7/8,outer=T)
hist(hits.flat.negwedge.adapt$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("Wedge disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=5/8, outer=T)
hist(hits.flat.negs$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("S disturbance",2,cex=0.9,line=0.25,at=3/8,outer=T)
hist(hits.flat.negs.adapt$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("S disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=1/8, outer=T)

hist(hits.sine.negwedge$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("Wedge disturbance",2,cex=0.9,line=0.25,at=7/8,outer=T)
hist(hits.sine.negwedge.adapt$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("Wedge disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=5/8, outer=T)
hist(hits.sine.negs$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("S disturbance",2,cex=0.9,line=0.25,at=3/8,outer=T)
hist(hits.sine.negs.adapt$delta.recov, xlab="", ylab="", col="lightgrey", main="")
mtext("S disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=1/8, outer=T)

mtext("Error in recovery date",1,outer=T,line=0.5,at=0.9/2)
mtext("Frequency",2,outer=T,line=3.5)
mtext("Flat background",outer=T,at=0.9/4,line=-0.25,cex=0.8)
mtext("Seasonal background", outer=T, at=0.9/4*3, line=-0.25,cex=0.8)

dev.off()



#Figure SX: false positiver rate surfaces-------------------------------------------------------
joint.min<-min(c(fpr.flat.negwedge, fpr.flat.negs, fpr.sine.negs, fpr.sine.negwedge,
                 fpr.flat.negwedge.adapt, fpr.flat.negs.adapt, fpr.sine.negs.adapt, fpr.sine.negwedge.adapt))
joint.max<-max(c(fpr.flat.negwedge, fpr.flat.negs, fpr.sine.negs, fpr.sine.negwedge,
                 fpr.flat.negwedge.adapt, fpr.flat.negs.adapt, fpr.sine.negs.adapt, fpr.sine.negwedge.adapt))

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figSX_fpr_negative.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
image(along.s, along.r, fpr.flat.negwedge, col=viridis(25), zlim=c(joint.min,joint.max), ylab="")
mtext("Wedge disturbance\nFlat background",2,cex=0.9,line=0.25,at=7/8,outer=T)
image(along.s, along.r, fpr.sine.negwedge, col=viridis(25), zlim=c(joint.min,joint.max), ylab="")
mtext("Wedge disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=5/8,outer=T)
mtext("S disturbance\nFlat background",2,cex=0.9,line=0.25,at=3/8,outer=T)
image(along.s, along.r, fpr.flat.negs, col=viridis(25), zlim=c(joint.min,joint.max))
mtext("S disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=1/8,outer=T)
image(along.s, along.r, fpr.sine.negs, col=viridis(25), zlim=c(joint.min,joint.max))

image(along.s, along.r, fpr.flat.negwedge.adapt, col=viridis(25), zlim=c(joint.min,joint.max), ylab="")
image(along.s, along.r, fpr.flat.negs.adapt, col=viridis(25), zlim=c(joint.min,joint.max))
image(along.s, along.r, fpr.sine.negwedge.adapt, col=viridis(25), zlim=c(joint.min,joint.max), ylab="")
image(along.s, along.r, fpr.sine.negs.adapt, col=viridis(25), zlim=c(joint.min,joint.max))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",yaxt="n",ylab="False positive rate")
axis(2,at=seq(0,1,length.out=4),labels=round(seq(joint.min,joint.max,length.out=4),3))

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=3.5)
mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()


### for positive disturbances --------------------------------------------------------


#Figure SX: TPR for positive disturbances
summary(c(tpr.flat.poswedge))
summary(c(tpr.flat.poss))
summary(c(tpr.sine.poswedge))
summary(c(tpr.sine.poss))
summary(c(tpr.flat.poswedge.adapt))
summary(c(tpr.flat.poss.adapt))
summary(c(tpr.sine.poswedge.adapt))
summary(c(tpr.sine.poss.adapt))

summary(c(tpr.sine.poss.adapt-tpr.sine.poss))
summary(c(tpr.sine.poswedge.adapt-tpr.sine.poswedge))

# jointmax<-ceiling(max(c(fpr.flat.poswedge,fpr.flat.poss,fpr.flat.poswedge.adapt,fpr.flat.poss.adapt)))

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figSX_tpr_positive.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
image(along.s, along.r, tpr.flat.poswedge, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("Wedge disturbance\nFlat background",2,cex=0.9,line=0.25,at=7/8,outer=T)
image(along.s, along.r, tpr.sine.poswedge, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("Wedge disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=5/8,outer=T)
mtext("S disturbance\nFlat background",2,cex=0.9,line=0.25,at=3/8,outer=T)
image(along.s, along.r, tpr.flat.poss, col=rev(viridis(25)), zlim=c(0,1))
mtext("S disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=1/8,outer=T)
image(along.s, along.r, tpr.sine.poss, col=rev(viridis(25)), zlim=c(0,1))

image(along.s, along.r, tpr.flat.poswedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
image(along.s, along.r, tpr.flat.poss.adapt, col=rev(viridis(25)), zlim=c(0,1))
image(along.s, along.r, tpr.sine.poswedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
image(along.s, along.r, tpr.sine.poss.adapt, col=rev(viridis(25)), zlim=c(0,1))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=rev(viridis(25)),xaxt="n",ylab="True detection rate")

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=3.5)
mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()


#Figure SX: error for recovery date surfaces-------------------------------------------------------

joint.min<-min(c(error.flat.poswedge, error.flat.poss, error.sine.poss, error.sine.poswedge,
                 error.flat.poswedge.adapt, error.flat.poss.adapt, error.sine.poss.adapt, error.sine.poswedge.adapt))
joint.max<-max(c(error.flat.poswedge, error.flat.poss, error.sine.poss, error.sine.poswedge,
                 error.flat.poswedge.adapt, error.flat.poss.adapt, error.sine.poss.adapt, error.sine.poswedge.adapt))

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figSX_recoverror_positive.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
image(along.s, along.r, plotScale(error.flat.poswedge,joint.max), col=errorPal, 
      ylab="", zlim=c(1,n))
mtext("Wedge disturbance\nFlat background",2,cex=0.9,line=0.25,at=7/8,outer=T)
image(along.s, along.r, plotScale(error.sine.poswedge,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
mtext("Wedge disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=5/8,outer=T)
mtext("S disturbance\nFlat background",2,cex=0.9,line=0.25,at=3/8,outer=T)
image(along.s, along.r, plotScale(error.flat.poss,joint.max), col=errorPal, 
      zlim=c(1,n))
mtext("S disturbance\nSeasonal background",2,cex=0.9,line=0.25,at=1/8,outer=T)
image(along.s, along.r, plotScale(error.sine.poss,joint.max), col=errorPal, zlim=c(1,n))

image(along.s, along.r, plotScale(error.flat.poswedge.adapt,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
image(along.s, along.r, plotScale(error.flat.poss.adapt,joint.max), col=errorPal, 
      zlim=c(1,n))
image(along.s, along.r, plotScale(error.sine.poswedge.adapt,joint.max), col=errorPal, 
      zlim=c(1,n), ylab="")
image(along.s, along.r, plotScale(error.sine.poss.adapt,joint.max), col=errorPal,
      zlim=c(1,n))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:51)),col=errorPal,xaxt="n",yaxt="n",ylab="Recovery date error (true-estimated)")
axis(2,at=seq(0,1,length.out=5),labels=round(seq(floor(-1*joint.max),ceiling(joint.max),length.out=5),1))

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=3.5)
mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()






#Figure 3: TPR/FPR for seasonal negwedge, seasonal negs, seasonal negwedge adaptive, flat negs adaptive
# 
# sum(fpr.flat.negs)/sum(fpr.flat.negs.adapt)
# sum(fpr.flat.negwedge)/sum(fpr.flat.negwedge.adapt)
# sum(fpr.sine.negs)/sum(fpr.sine.negs.adapt)
# sum(fpr.sine.negwedge)/sum(fpr.sine.negwedge.adapt)
# 
# jointmax<-ceiling(max(c(fpr.sine.negwedge,fpr.sine.negs,fpr.sine.negwedge.adapt,fpr.sine.negs.adapt)))

# 
# png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figS_tprfpr_seasonalbackground.png",
#     units="in", res=300, width=6.5,height=8)
# 
# layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
# # par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
# # image(along.s, along.r, tpr.sine.negwedge, col=pal.tpr(25), zlim=c(0,1), ylab="")
# # mtext("Wedge disturbance",2,cex=0.9,line=0.25,at=7/8,outer=T)
# # image(along.s, along.r, tpr.sine.negwedge.adapt, col=pal.tpr(25), zlim=c(0,1), ylab="")
# # mtext("Wedge disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=5/8, outer=T)
# # image(along.s, along.r, tpr.sine.negs, col=pal.tpr(25), zlim=c(0,1))
# # mtext("S disturbance",2,cex=0.9,line=0.25,at=3/8,outer=T)
# # image(along.s, along.r, tpr.sine.negs.adapt, col=pal.tpr(25), zlim=c(0,1))
# # mtext("S disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=1/8, outer=T)
# 
# # image(along.s, along.r, fpr.flat.negwedge, col=pal.fpr(25), zlim=c(0,jointmax))
# # image(along.s, along.r, fpr.flat.negwedge.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# # image(along.s, along.r, fpr.flat.negs, col=pal.fpr(25), zlim=c(0,jointmax))
# # image(along.s, along.r, fpr.flat.negs.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# 
# image(along.s, along.r, fpr.sine.negwedge, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.negwedge.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.negs, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.negs.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# 
# par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
# image(z=t(matrix(1:25)),col=pal.tpr(25),xaxt="n",ylab="True detection rate")
# image(z=t(matrix(1:25)),col=pal.fpr(25),xaxt="n",yaxt="n", ylab="False detection rate")
# axis(2,at=seq(0,1,length.out=5),labels=round(seq(0,jointmax,length.out=5),1))
# 
# mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
# mtext("Duration (days)",2,outer=T,line=3.5)
# mtext("Seasonal background variation",outer=T,at=0.9/2,line=-0.25)
# 
# dev.off()
# 
# 
# #Figure S1: TPR/FPR for flat poswedge, flat poss, flat poswedge adaptive, flat poss adaptive
# 
# jointmax<-ceiling(max(c(fpr.flat.poswedge,fpr.flat.poss,fpr.flat.poswedge.adapt,fpr.flat.poss.adapt)))
# 
# png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figS1_tprfpr_flatbackground_positive.png",
#     units="in", res=300, width=6.5,height=8)
# 
# layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
# par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
# image(along.s, along.r, tpr.flat.poswedge, col=pal.tpr(25), zlim=c(0,1), ylab="")
# mtext("Wedge disturbance",2,cex=0.9,line=0.25,at=7/8,outer=T)
# image(along.s, along.r, tpr.flat.poswedge.adapt, col=pal.tpr(25), zlim=c(0,1), ylab="")
# mtext("Wedge disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=5/8, outer=T)
# image(along.s, along.r, tpr.flat.poss, col=pal.tpr(25), zlim=c(0,1))
# mtext("S disturbance",2,cex=0.9,line=0.25,at=3/8,outer=T)
# image(along.s, along.r, tpr.flat.poss.adapt, col=pal.tpr(25), zlim=c(0,1))
# mtext("S disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=1/8, outer=T)
# 
# image(along.s, along.r, fpr.flat.poswedge, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.flat.poswedge.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.flat.poss, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.flat.poss.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# 
# par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
# image(z=t(matrix(1:25)),col=pal.tpr(25),xaxt="n",ylab="True detetion rate")
# image(z=t(matrix(1:25)),col=pal.fpr(25),xaxt="n",yaxt="n", ylab="False detection rate")
# axis(2,at=seq(0,1,length.out=5),labels=round(seq(0,jointmax,length.out=5),1))
# 
# mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
# mtext("Duration (days)",2,outer=T,line=3.5)
# mtext("Flat background variation",outer=T,at=0.9/2,line=-0.25)
# 
# dev.off()
# 
# 
# #Figure S2: TPR/FPR for seasonal poswedge, seasonal poss, seasonal poswedge adaptive, seasonal poss adaptive
# 
# jointmax<-ceiling(max(c(fpr.sine.poswedge,fpr.sine.poss,fpr.sine.poswedge.adapt,fpr.sine.poss.adapt)))
# 
# png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figS2_tprfpr_seasonalbackground_pos.png",
#     units="in", res=300, width=6.5,height=8)
# 
# layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
# par(mar=c(2.1,2.1,1.1,1.1),oma=c(2.1,5.1,1.3,0))
# image(along.s, along.r, tpr.sine.poswedge, col=pal.tpr(25), zlim=c(0,1), ylab="")
# mtext("Wedge disturbance",2,cex=0.9,line=0.25,at=7/8,outer=T)
# image(along.s, along.r, tpr.sine.poswedge.adapt, col=pal.tpr(25), zlim=c(0,1), ylab="")
# mtext("Wedge disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=5/8, outer=T)
# image(along.s, along.r, tpr.sine.poss, col=pal.tpr(25), zlim=c(0,1))
# mtext("S disturbance",2,cex=0.9,line=0.25,at=3/8,outer=T)
# image(along.s, along.r, tpr.sine.poss.adapt, col=pal.tpr(25), zlim=c(0,1))
# mtext("S disturbance\nadaptive reference",2, cex=0.9, line=0.25, at=1/8, outer=T)
# 
# image(along.s, along.r, fpr.sine.poswedge, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.poswedge.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.poss, col=pal.fpr(25), zlim=c(0,jointmax))
# image(along.s, along.r, fpr.sine.poss.adapt, col=pal.fpr(25), zlim=c(0,jointmax))
# 
# par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
# image(z=t(matrix(1:25)),col=pal.tpr(25),xaxt="n",ylab="True detection rate")
# image(z=t(matrix(1:25)),col=pal.fpr(25),xaxt="n",yaxt="n", ylab="False detection rate")
# axis(2,at=seq(0,1,length.out=5),labels=round(seq(0,jointmax,length.out=5),1))
# 
# mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
# mtext("Duration (days)",2,outer=T,line=3.5)
# mtext("Seasonal background variation",outer=T,at=0.9/2,line=-0.25)
# 
# dev.off()


## Figure for ERI poster

laymat=matrix(1:4, nrow=1, byrow=TRUE)

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/posterFig_tpr_error.png",
    units="in", res=300, width=8,height=4)

layout(laymat, widths=rep(c(0.4,0.1),2))
par(mar=c(3.5,3.5,1.1,0.5), mgp=c(2,0.8,0), cex.axis=1.1, cex.lab=1.2)
image(along.s, along.r, tpr.sine.negwedge.adapt, col=rev(viridis(25)), zlim=c(0,1),
      xlab="Severity", ylab="Duration (days)")
par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2,0.8,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=rev(viridis(25)),xaxt="n",ylab="True detection rate")

par(mar=c(3.5,3.5,1.1,0.5), mgp=c(2,0.8,0), cex.axis=1.1, cex.lab=1.2)
image(along.s, along.r, error.sine.negwedge.adapt, col=viridis(25), ylab="Duration (days)",
      xlab="Severity")
par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2,0.8,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",yaxt="n",ylab="Recovery date error (true-estimated)")
axis(2,at=seq(0,1,length.out=4),labels=round(seq(floor(min(error.sine.negwedge.adapt)),
                                                 ceiling(max(error.sine.negwedge.adapt)),length.out=4),1))

dev.off()