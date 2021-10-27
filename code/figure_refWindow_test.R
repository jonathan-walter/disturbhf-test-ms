## Interpret evaluations of all simulation tests

rm(list=ls())

library(here)
library(viridis)

## ----------------------------------------------------------------------------
## Load data

res_w7.5 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest75d.rds")
res_w15 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest15d.rds")
res_w30 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest30d.rds")
res_w60 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest60dx.rds")
res_w120 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest120d.rds")
res_w240 <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest240d.rds")


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

## ------------------------------------------------------------------------------------------------
## Extract stats from results

denom <- sum(!is.na(res_w60[[1]]$empdiff$ddiff))

eval.w7.5<-getEval(res_w7.5, which.eval = "normal")
cln.w7.5<-cleanEval(eval.w7.5)
tpr.w7.5<-tprSurface(cln.w7.5,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w7.5<-falmSurface(eval.w7.5, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

eval.w15<-getEval(res_w15, which.eval = "normal")
cln.w15<-cleanEval(eval.w15)
tpr.w15<-tprSurface(cln.w15,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w15<-falmSurface(eval.w15, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

eval.w30<-getEval(res_w30, which.eval = "normal")
cln.w30<-cleanEval(eval.w30)
tpr.w30<-tprSurface(cln.w30,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w30<-falmSurface(eval.w30, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

eval.w60<-getEval(res_w60, which.eval = "normal")
cln.w60<-cleanEval(eval.w60)
tpr.w60<-tprSurface(cln.w60,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w60<-falmSurface(eval.w60, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

eval.w120<-getEval(res_w120, which.eval = "normal")
cln.w120<-cleanEval(eval.w120)
tpr.w120<-tprSurface(cln.w120,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w120<-falmSurface(eval.w120, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)

eval.w240<-getEval(res_w240, which.eval = "normal")
cln.w240<-cleanEval(eval.w240)
tpr.w240<-tprSurface(cln.w240,along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5))
fpr.w240<-falmSurface(eval.w240, along.s=seq(0,5,by=0.5),along.r=seq(3,60,by=5), denom=denom)




## -------------------------------------------------------------------------------------
## Nice plots

laymat4<-matrix(c(1:6,8,7,7),3,3)
along.s<-seq(0,5,by=0.5)
along.r<-seq(3,60,by=5)


png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figX_test_refWidth_tpr.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,2.1,1.1),oma=c(2.1,3.1,1.3,0))
image(along.s, along.r, tpr.w7.5, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("refwidth = 7.5 days", cex=0.8, line=0.25)
image(along.s, along.r, tpr.w15, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("refwidth = 15 days", cex=0.8, line=0.25)
image(along.s, along.r, tpr.w30, col=rev(viridis(25)), zlim=c(0,1))
mtext("refwidth = 30 days", cex=0.8, line=0.25)
image(along.s, along.r, tpr.w60, col=rev(viridis(25)), zlim=c(0,1))
mtext("refwidth = 60 days", cex=0.8, line=0.25)
image(along.s, along.r, tpr.w120, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("refwidth = 120 days", cex=0.8, line=0.25)
image(along.s, along.r, tpr.w240, col=rev(viridis(25)), zlim=c(0,1), ylab="")
mtext("refwidth = 240 days", cex=0.8, line=0.25)
# image(along.s, along.r, tpr.flat.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))
# image(along.s, along.r, tpr.sine.negwedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
# image(along.s, along.r, tpr.sine.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=rev(viridis(25)),xaxt="n",ylab="True detection rate")

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=1.5)
# mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
# mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()



jointmax=max(c(fpr.w7.5, fpr.w15, fpr.w30, fpr.w60, fpr.w120, fpr.w240))

png("~/Box Sync/EstuaryStormResilience/AlgorithmManuscript/figX_test_refWidth_fpr.png",
    units="in", res=300, width=6.5,height=8)

layout(laymat4, widths=c(0.87/2,0.87/2,0.13))
par(mar=c(2.1,2.1,2.1,1.1),oma=c(2.1,3.1,1.3,0))
image(along.s, along.r, fpr.w7.5, col=viridis(25), zlim=c(0,jointmax), ylab="")
mtext("refwidth = 7.5 days", cex=0.8, line=0.25)
image(along.s, along.r, fpr.w15, col=viridis(25), zlim=c(0,jointmax), ylab="")
mtext("refwidth = 15 days", cex=0.8, line=0.25)
image(along.s, along.r, fpr.w30, col=viridis(25), zlim=c(0,jointmax))
mtext("refwidth = 30 days", cex=0.8, line=0.25)
image(along.s, along.r, fpr.w60, col=viridis(25), zlim=c(0,jointmax))
mtext("refwidth = 60 days", cex=0.8, line=0.25)
image(along.s, along.r, fpr.w120, col=viridis(25), zlim=c(0,jointmax), ylab="")
mtext("refwidth = 120 days", cex=0.8, line=0.25)
image(along.s, along.r, fpr.w240, col=viridis(25), zlim=c(0,jointmax), ylab="")
mtext("refwidth = 240 days", cex=0.8, line=0.25)
# image(along.s, along.r, tpr.flat.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))
# image(along.s, along.r, tpr.sine.negwedge.adapt, col=rev(viridis(25)), zlim=c(0,1), ylab="")
# image(along.s, along.r, tpr.sine.negs.adapt, col=rev(viridis(25)), zlim=c(0,1))

par(mar=c(2.1,3.6,1.1,1.1), mgp=c(2.25,1,0), cex.axis=1.1, cex.lab=1.2)
image(z=t(matrix(1:25)),col=viridis(25),xaxt="n",ylab="False positive rate",yaxt="n")
axis(2, at=seq(0,1,length.out=4), labels=round(seq(0, jointmax,length.out=4), digits=3))

mtext("Severity",1,outer=T,line=0.5,at=0.9/2)
mtext("Duration (days)",2,outer=T,line=1.5)
# mtext("Standard reference",outer=T,at=0.9/4,line=-0.5,cex=0.8)
# mtext("Adaptive reference", outer=T, at=0.9/4*3, line=-0.5,cex=0.8)

dev.off()

