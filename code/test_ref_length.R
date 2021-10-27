## test effect of changing length of adaptive reference window

rm(list=ls())

library(here)
library(disturbhf)

source(here("./code/simdisturb.R"))
source(here("./code/evalstats.R"))

sims<-simDisturb(bkgrnd="sinusoidal",disturbtype="neg.wedge",nreps=5000,dist.start=366+120/2)


## 120 days -----------------------------------------------------------------------------------------------

refwidth=120 #days

# Full tests: sinusoidal background, negative wedge
sims120<-sims
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sims120)){
  dt<-sims120[[ii]]$key["dt"]
  sims120[[ii]]$key["dday"]<-sims120[[ii]]$key["dday"]-365
  sims120[[ii]]$key["rday"]<-sims120[[ii]]$key["rday"]-365
  
  refy<-sims120[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims120[[ii]]$simts[(366/dt):length(sims120[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims120[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims120[[ii]]$empalarm<-disturbalarm(sims120[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims120[[ii]]$filtalarm<-alarmfilter(sims120[[ii]]$empalarm,dmin.dist=3)
  sims120[[ii]]$eval<-evalStats(sims120[[ii]]$key, sims120[[ii]]$empalarm, tol=5)
  sims120[[ii]]$filteval<-evalStats(sims120[[ii]]$key, sims120[[ii]]$filtalarm, tol=5)
}

saveRDS(sims120, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest120d.rds")

## 60 days -----------------------------------------------------------------------------------------------

refwidth=60 #days

# Full tests: sinusoidal background, negative wedge

sims60 <- sims

for(ii in 1:length(sims60)){
  dt<-sims60[[ii]]$key["dt"]
  sims60[[ii]]$key["dday"]<-sims60[[ii]]$key["dday"]-365
  sims60[[ii]]$key["rday"]<-sims60[[ii]]$key["rday"]-365
  
  refy<-sims60[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims60[[ii]]$simts[(366/dt):length(sims60[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims60[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims60[[ii]]$empalarm<-disturbalarm(sims60[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims60[[ii]]$filtalarm<-alarmfilter(sims60[[ii]]$empalarm,dmin.dist=3)
  sims60[[ii]]$eval<-evalStats(sims60[[ii]]$key, sims60[[ii]]$empalarm, tol=5)
  sims60[[ii]]$filteval<-evalStats(sims60[[ii]]$key, sims60[[ii]]$filtalarm, tol=5)
}

saveRDS(sims60, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest60d.rds")


## 30 days -----------------------------------------------------------------------------------------------

refwidth=30 #days

# Full tests: sinusoidal background, negative wedge
sims30<-sims
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sims30)){
  dt<-sims30[[ii]]$key["dt"]
  sims30[[ii]]$key["dday"]<-sims30[[ii]]$key["dday"]-365
  sims30[[ii]]$key["rday"]<-sims30[[ii]]$key["rday"]-365
  
  refy<-sims30[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims30[[ii]]$simts[(366/dt):length(sims30[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims30[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims30[[ii]]$empalarm<-disturbalarm(sims30[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims30[[ii]]$filtalarm<-alarmfilter(sims30[[ii]]$empalarm,dmin.dist=3)
  sims30[[ii]]$eval<-evalStats(sims30[[ii]]$key, sims30[[ii]]$empalarm, tol=5)
  sims30[[ii]]$filteval<-evalStats(sims30[[ii]]$key, sims30[[ii]]$filtalarm, tol=5)
}

saveRDS(sims30, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest30d.rds")



## 15 days -----------------------------------------------------------------------------------------------

refwidth=15 #days

# Full tests: sinusoidal background, negative wedge
sims15<-sims
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sims15)){
  dt<-sims15[[ii]]$key["dt"]
  sims15[[ii]]$key["dday"]<-sims15[[ii]]$key["dday"]-365
  sims15[[ii]]$key["rday"]<-sims15[[ii]]$key["rday"]-365
  
  refy<-sims15[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims15[[ii]]$simts[(366/dt):length(sims15[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims15[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims15[[ii]]$empalarm<-disturbalarm(sims15[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims15[[ii]]$filtalarm<-alarmfilter(sims15[[ii]]$empalarm,dmin.dist=3)
  sims15[[ii]]$eval<-evalStats(sims15[[ii]]$key, sims15[[ii]]$empalarm, tol=5)
  sims15[[ii]]$filteval<-evalStats(sims15[[ii]]$key, sims15[[ii]]$filtalarm, tol=5)
}

saveRDS(sims15, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest15d.rds")



## 7.5 days -----------------------------------------------------------------------------------------------

refwidth=7.5 #days

# Full tests: sinusoidal background, negative wedge
sims7.5<-sims
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sims7.5)){
  dt<-sims7.5[[ii]]$key["dt"]
  sims7.5[[ii]]$key["dday"]<-sims7.5[[ii]]$key["dday"]-365
  sims7.5[[ii]]$key["rday"]<-sims7.5[[ii]]$key["rday"]-365
  
  refy<-sims7.5[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims7.5[[ii]]$simts[(366/dt):length(sims7.5[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims7.5[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims7.5[[ii]]$empalarm<-disturbalarm(sims7.5[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims7.5[[ii]]$filtalarm<-alarmfilter(sims7.5[[ii]]$empalarm,dmin.dist=3)
  sims7.5[[ii]]$eval<-evalStats(sims7.5[[ii]]$key, sims7.5[[ii]]$empalarm, tol=5)
  sims7.5[[ii]]$filteval<-evalStats(sims7.5[[ii]]$key, sims7.5[[ii]]$filtalarm, tol=5)
}

saveRDS(sims7.5, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest75d.rds")



## 240 days -----------------------------------------------------------------------------------------------

refwidth=240 #days


sims <- readRDS("/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest75d.rds")

# Full tests: sinusoidal background, negative wedge
sims240<-sims
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sims240)){
  dt<-sims240[[ii]]$key["dt"]
  # sims240[[ii]]$key["dday"]<-sims240[[ii]]$key["dday"]-365
  # sims240[[ii]]$key["rday"]<-sims240[[ii]]$key["rday"]-365
  # 
  refy<-sims240[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sims240[[ii]]$simts[(366/dt):length(sims240[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sims240[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6)
  sims240[[ii]]$empalarm<-disturbalarm(sims240[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sims240[[ii]]$filtalarm<-alarmfilter(sims240[[ii]]$empalarm,dmin.dist=3)
  sims240[[ii]]$eval<-evalStats(sims240[[ii]]$key, sims240[[ii]]$empalarm, tol=5)
  sims240[[ii]]$filteval<-evalStats(sims240[[ii]]$key, sims240[[ii]]$filtalarm, tol=5)
}

saveRDS(sims240, file="/Users/jonathanwalter/Box Sync/EstuaryStormResilience/AlgorithmTestOutput/sim_sin_negwedge_reftest240d.rds")
