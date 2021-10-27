## Simulation tests of disturbance and recovery detection algorithm

rm(list=ls())

#library(devtools) #only run for new/updated install of disturbhf
#install_github("jonathan-walter/disturbhf")

library(here)
library(disturbhf)

source(here("./code/simdisturb.R"))
source(here("./code/evalstats.R"))


# ## Demo: flat background, negative wedge 
# sim.demo<-simDisturb(bkgrnd="flat",disturbtype="neg.wedge",nreps=5)
# plot(sim.demo[[1]]$simts, type="l")
# 
# for(ii in 1:length(sim.demo)){
#   simts<-sim.demo[[ii]]$simts
#   dt<-sim.demo[[ii]]$key["dt"]
#   tt<-(1:length(simts))*dt
#   refy<-simts[1:(sim.demo[[ii]]$key["ref.end"]/dt)]
#   sim.demo[[ii]]$empdiff<-mwdistdiffz(yy=simts,tt=tt,refy=refy,wwidth=5/dt)
#   sim.demo[[ii]]$empalarm<-disturbalarm(sim.demo[[ii]]$empdiff, dthresh=2, rthresh=0.5)
#   sim.demo[[ii]]$eval<-evalStats(sim.demo[[ii]]$key, sim.demo[[ii]]$empalarm, tol=5)
#   print(sim.demo[[ii]]$eval)
# }

# # Full tests: flat background, negative wedge
sttime<-Sys.time()

sim.flat.negwedge<-simDisturb(bkgrnd="flat",disturbtype="neg.wedge",nreps=5000)
#plot(sim.flat.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.negwedge)){
  dt<-sim.flat.negwedge[[ii]]$key["dt"]
  sim.flat.negwedge[[ii]]$key["dday"]<-sim.flat.negwedge[[ii]]$key["dday"]-sim.flat.negwedge[[ii]]$key["ref.end"]
  sim.flat.negwedge[[ii]]$key["rday"]<-sim.flat.negwedge[[ii]]$key["rday"]-sim.flat.negwedge[[ii]]$key["ref.end"]
  refy<-sim.flat.negwedge[[ii]]$simts[1:(sim.flat.negwedge[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.negwedge[[ii]]$simts[(sim.flat.negwedge[[ii]]$key["ref.end"]/dt+1):length(sim.flat.negwedge[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  #tt<-(1:length(yy))*dt

  sim.flat.negwedge[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.negwedge[[ii]]$empalarm<-disturbalarm(sim.flat.negwedge[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.negwedge[[ii]]$filtalarm<-alarmfilter(sim.flat.negwedge[[ii]]$empalarm,dmin.dist=3)
  sim.flat.negwedge[[ii]]$eval<-evalStats(sim.flat.negwedge[[ii]]$key, sim.flat.negwedge[[ii]]$empalarm, tol=5)
  sim.flat.negwedge[[ii]]$filteval<-evalStats(sim.flat.negwedge[[ii]]$key, sim.flat.negwedge[[ii]]$filtalarm, tol=5)
  #print(sim.flat.negwedge[[ii]]$eval)
}

print(Sys.time()-sttime)
saveRDS(sim.flat.negwedge, file="sim_flat_negwedge_integral_20211011_t5.rds")

## Full tests: flat background, positive wedge
sim.flat.poswedge<-simDisturb(bkgrnd="flat",disturbtype="pos.wedge",nreps=5000)
#plot(sim.flat.poswedge[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.poswedge)){
  dt<-sim.flat.poswedge[[ii]]$key["dt"]
  sim.flat.poswedge[[ii]]$key["dday"]<-sim.flat.poswedge[[ii]]$key["dday"]-sim.flat.poswedge[[ii]]$key["ref.end"]
  sim.flat.poswedge[[ii]]$key["rday"]<-sim.flat.poswedge[[ii]]$key["rday"]-sim.flat.poswedge[[ii]]$key["ref.end"]
  refy<-sim.flat.poswedge[[ii]]$simts[1:(sim.flat.poswedge[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.poswedge[[ii]]$simts[(sim.flat.poswedge[[ii]]$key["ref.end"]/dt+1):length(sim.flat.poswedge[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.poswedge[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.poswedge[[ii]]$empalarm<-disturbalarm(sim.flat.poswedge[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.poswedge[[ii]]$filtalarm<-alarmfilter(sim.flat.poswedge[[ii]]$empalarm,dmin.dist=3)
  sim.flat.poswedge[[ii]]$eval<-evalStats(sim.flat.poswedge[[ii]]$key, sim.flat.poswedge[[ii]]$empalarm, tol=2.5)
  sim.flat.poswedge[[ii]]$filteval<-evalStats(sim.flat.poswedge[[ii]]$key, sim.flat.poswedge[[ii]]$filtalarm, tol=2.5)
  #print(sim.flat.poswedge[[ii]]$eval)
}
saveRDS(sim.flat.poswedge, file="sim_flat_poswedge_integral_20211011_t5.rds")

# Full tests: flat background, negative S
sim.flat.negs<-simDisturb(bkgrnd="flat",disturbtype="neg.s",nreps=5000)
#plot(sim.flat.negs[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.negs)){
  dt<-sim.flat.negs[[ii]]$key["dt"]
  sim.flat.negs[[ii]]$key["dday"]<-sim.flat.negs[[ii]]$key["dday"]-sim.flat.negs[[ii]]$key["ref.end"]
  sim.flat.negs[[ii]]$key["rday"]<-sim.flat.negs[[ii]]$key["rday"]-sim.flat.negs[[ii]]$key["ref.end"]
  refy<-sim.flat.negs[[ii]]$simts[1:(sim.flat.negs[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.negs[[ii]]$simts[(sim.flat.negs[[ii]]$key["ref.end"]/dt+1):length(sim.flat.negs[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)

  sim.flat.negs[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.negs[[ii]]$empalarm<-disturbalarm(sim.flat.negs[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.negs[[ii]]$filtalarm<-alarmfilter(sim.flat.negs[[ii]]$empalarm,dmin.dist=3)
  sim.flat.negs[[ii]]$eval<-evalStats(sim.flat.negs[[ii]]$key, sim.flat.negs[[ii]]$empalarm, tol=5)
  sim.flat.negs[[ii]]$filteval<-evalStats(sim.flat.negs[[ii]]$key, sim.flat.negs[[ii]]$filtalarm, tol=5)
  #print(sim.flat.negs[[ii]]$eval)
}
saveRDS(sim.flat.negs, file="sim_flat_negs_integral_20211011_t5.rds")

## Full tests: flat background, positive s
sim.flat.poss<-simDisturb(bkgrnd="flat",disturbtype="pos.s",nreps=5000)
#plot(sim.flat.poss[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.poss)){
  dt<-sim.flat.poss[[ii]]$key["dt"]
  sim.flat.poss[[ii]]$key["dday"]<-sim.flat.poss[[ii]]$key["dday"]-sim.flat.poss[[ii]]$key["ref.end"]
  sim.flat.poss[[ii]]$key["rday"]<-sim.flat.poss[[ii]]$key["rday"]-sim.flat.poss[[ii]]$key["ref.end"]
  refy<-sim.flat.poss[[ii]]$simts[1:(sim.flat.poss[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.poss[[ii]]$simts[(sim.flat.poss[[ii]]$key["ref.end"]/dt+1):length(sim.flat.poss[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.poss[[ii]]$empdiff<-mwdistdiffz(testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.poss[[ii]]$empalarm<-disturbalarm(sim.flat.poss[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.poss[[ii]]$filtalarm<-alarmfilter(sim.flat.poss[[ii]]$empalarm,dmin.dist=3)
  sim.flat.poss[[ii]]$eval<-evalStats(sim.flat.poss[[ii]]$key, sim.flat.poss[[ii]]$empalarm, tol=5)
  sim.flat.poss[[ii]]$filteval<-evalStats(sim.flat.poss[[ii]]$key, sim.flat.poss[[ii]]$filtalarm, tol=5)
  #print(sim.flat.poss[[ii]]$eval)
}
saveRDS(sim.flat.poss, file="sim_flat_poss_integral_20211011_t5.rds")

## Full tests: sinusoidal background, negative wedge
sim.sin.negwedge<-simDisturb(bkgrnd="sinusoidal",disturbtype="neg.wedge",nreps=5000)
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.negwedge)){
  dt<-sim.sin.negwedge[[ii]]$key["dt"]
  sim.sin.negwedge[[ii]]$key["dday"]<-sim.sin.negwedge[[ii]]$key["dday"]-sim.sin.negwedge[[ii]]$key["ref.end"]
  sim.sin.negwedge[[ii]]$key["rday"]<-sim.sin.negwedge[[ii]]$key["rday"]-sim.sin.negwedge[[ii]]$key["ref.end"]
  refy<-sim.sin.negwedge[[ii]]$simts[1:(sim.sin.negwedge[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.negwedge[[ii]]$simts[(sim.sin.negwedge[[ii]]$key["ref.end"]/dt+1):length(sim.sin.negwedge[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)

  sim.sin.negwedge[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.negwedge[[ii]]$empalarm<-disturbalarm(sim.sin.negwedge[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.negwedge[[ii]]$filtalarm<-alarmfilter(sim.sin.negwedge[[ii]]$empalarm,dmin.dist=3)
  sim.sin.negwedge[[ii]]$eval<-evalStats(sim.sin.negwedge[[ii]]$key, sim.sin.negwedge[[ii]]$empalarm, tol=5)
  sim.sin.negwedge[[ii]]$filteval<-evalStats(sim.sin.negwedge[[ii]]$key, sim.sin.negwedge[[ii]]$filtalarm, tol=5)
  #print(sim.sin.negwedge[[ii]]$eval)
}
saveRDS(sim.sin.negwedge, file="sim_sin_negwedge_integral_20211011_t5.rds")

## Full tests: sinusoidal background, positive wedge
sim.sin.poswedge<-simDisturb(bkgrnd="sinusoidal",disturbtype="pos.wedge",nreps=5000)
#plot(sim.sin.poswedge[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.poswedge)){
  dt<-sim.sin.poswedge[[ii]]$key["dt"]
  sim.sin.poswedge[[ii]]$key["dday"]<-sim.sin.poswedge[[ii]]$key["dday"]-sim.sin.poswedge[[ii]]$key["ref.end"]
  sim.sin.poswedge[[ii]]$key["rday"]<-sim.sin.poswedge[[ii]]$key["rday"]-sim.sin.poswedge[[ii]]$key["ref.end"]
  refy<-sim.sin.poswedge[[ii]]$simts[1:(sim.sin.poswedge[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.poswedge[[ii]]$simts[(sim.sin.poswedge[[ii]]$key["ref.end"]/dt+1):length(sim.sin.poswedge[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.sin.poswedge[[ii]]$empdiff<-mwdistdiffz(testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.poswedge[[ii]]$empalarm<-disturbalarm(sim.sin.poswedge[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.poswedge[[ii]]$filtalarm<-alarmfilter(sim.sin.poswedge[[ii]]$empalarm,dmin.dist=3)
  sim.sin.poswedge[[ii]]$eval<-evalStats(sim.sin.poswedge[[ii]]$key, sim.sin.poswedge[[ii]]$empalarm, tol=5)
  sim.sin.poswedge[[ii]]$filteval<-evalStats(sim.sin.poswedge[[ii]]$key, sim.sin.poswedge[[ii]]$filtalarm, tol=5)
  #print(sim.sin.poswedge[[ii]]$eval)
}
saveRDS(sim.sin.poswedge, file="sim_sin_poswedge_integral_20211011_t5.rds")

## Full tests: sinusoidal background, negative s
sim.sin.negs<-simDisturb(bkgrnd="sinusoidal",disturbtype="neg.s",nreps=5000)
#plot(sim.sin.negs[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.negs)){
  dt<-sim.sin.negs[[ii]]$key["dt"]
  sim.sin.negs[[ii]]$key["dday"]<-sim.sin.negs[[ii]]$key["dday"]-sim.sin.negs[[ii]]$key["ref.end"]
  sim.sin.negs[[ii]]$key["rday"]<-sim.sin.negs[[ii]]$key["rday"]-sim.sin.negs[[ii]]$key["ref.end"]
  refy<-sim.sin.negs[[ii]]$simts[1:(sim.sin.negs[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.negs[[ii]]$simts[(sim.sin.negs[[ii]]$key["ref.end"]/dt+1):length(sim.sin.negs[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
 
  sim.sin.negs[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.negs[[ii]]$empalarm<-disturbalarm(sim.sin.negs[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.negs[[ii]]$filtalarm<-alarmfilter(sim.sin.negs[[ii]]$empalarm,dmin.dist=3)
  sim.sin.negs[[ii]]$eval<-evalStats(sim.sin.negs[[ii]]$key, sim.sin.negs[[ii]]$empalarm, tol=5)
  sim.sin.negs[[ii]]$filteval<-evalStats(sim.sin.negs[[ii]]$key, sim.sin.negs[[ii]]$filtalarm, tol=5)
  #print(sim.sin.negs[[ii]]$eval)
}
saveRDS(sim.sin.negs, file="sim_sin_negs_integral_20211011_t5.rds")

## Full tests: sinusoidal background, positive s
sim.sin.poss<-simDisturb(bkgrnd="sinusoidal",disturbtype="pos.s",nreps=5000)
#plot(sim.sin.poss[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.poss)){
  dt<-sim.sin.poss[[ii]]$key["dt"]
  sim.sin.poss[[ii]]$key["dday"]<-sim.sin.poss[[ii]]$key["dday"]-sim.sin.poss[[ii]]$key["ref.end"]
  sim.sin.poss[[ii]]$key["rday"]<-sim.sin.poss[[ii]]$key["rday"]-sim.sin.poss[[ii]]$key["ref.end"]
  refy<-sim.sin.poss[[ii]]$simts[1:(sim.sin.poss[[ii]]$key["ref.end"]/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.poss[[ii]]$simts[(sim.sin.poss[[ii]]$key["ref.end"]/dt+1):length(sim.sin.poss[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.sin.poss[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.poss[[ii]]$empalarm<-disturbalarm(sim.sin.poss[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.poss[[ii]]$filtalarm<-alarmfilter(sim.sin.poss[[ii]]$empalarm,dmin.dist=3)
  sim.sin.poss[[ii]]$eval<-evalStats(sim.sin.poss[[ii]]$key, sim.sin.poss[[ii]]$empalarm, tol=5)
  sim.sin.poss[[ii]]$filteval<-evalStats(sim.sin.poss[[ii]]$key, sim.sin.poss[[ii]]$filtalarm, tol=5)
  #print(sim.sin.poss[[ii]]$eval)
}
saveRDS(sim.sin.poss, file="sim_sin_poss_integral_20211011_t5.rds")

## ADD an adaptive reference window --------------------------------------------------------------------------------

refwidth=60 #days

# Full tests: sinusoidal background, negative wedge
sim.sin.negwedge2<-simDisturb(bkgrnd="sinusoidal",disturbtype="neg.wedge",nreps=5000,dist.start=366+refwidth/2)
#plot(sim.sin.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.negwedge2)){
  dt<-sim.sin.negwedge2[[ii]]$key["dt"]
  sim.sin.negwedge2[[ii]]$key["dday"]<-sim.sin.negwedge2[[ii]]$key["dday"]-365
  sim.sin.negwedge2[[ii]]$key["rday"]<-sim.sin.negwedge2[[ii]]$key["rday"]-365

  refy<-sim.sin.negwedge2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.negwedge2[[ii]]$simts[(366/dt):length(sim.sin.negwedge2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)

  sim.sin.negwedge2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.negwedge2[[ii]]$empalarm<-disturbalarm(sim.sin.negwedge2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.negwedge2[[ii]]$filtalarm<-alarmfilter(sim.sin.negwedge2[[ii]]$empalarm,dmin.dist=3)
  sim.sin.negwedge2[[ii]]$eval<-evalStats(sim.sin.negwedge2[[ii]]$key, sim.sin.negwedge2[[ii]]$empalarm, tol=5)
  sim.sin.negwedge2[[ii]]$filteval<-evalStats(sim.sin.negwedge2[[ii]]$key, sim.sin.negwedge2[[ii]]$filtalarm, tol=5)
}

saveRDS(sim.sin.negwedge2, file="sim_sin_negwedge_adapt_integral_20211011_t5.rds")


## Full tests: sinusoidal background, positive wedge
sim.sin.poswedge2<-simDisturb(bkgrnd="sinusoidal",disturbtype="pos.wedge",nreps=5000,dist.start=366+30)
#plot(sim.sin.poswedge2[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.poswedge2)){
  dt<-sim.sin.poswedge2[[ii]]$key["dt"]
  sim.sin.poswedge2[[ii]]$key["dday"]<-sim.sin.poswedge2[[ii]]$key["dday"]-365
  sim.sin.poswedge2[[ii]]$key["rday"]<-sim.sin.poswedge2[[ii]]$key["rday"]-365
  
  refy<-sim.sin.poswedge2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.poswedge2[[ii]]$simts[(366/dt):length(sim.sin.poswedge2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.sin.poswedge2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.poswedge2[[ii]]$empalarm<-disturbalarm(sim.sin.poswedge2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.poswedge2[[ii]]$filtalarm<-alarmfilter(sim.sin.poswedge2[[ii]]$empalarm,dmin.dist=3)
  sim.sin.poswedge2[[ii]]$eval<-evalStats(sim.sin.poswedge2[[ii]]$key, sim.sin.poswedge2[[ii]]$empalarm, tol=5)
  sim.sin.poswedge2[[ii]]$filteval<-evalStats(sim.sin.poswedge2[[ii]]$key, sim.sin.poswedge2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.sin.poswedge2, file="sim_sin_poswedge_adapt_integral_20211011_t5.rds")

## Full tests: sinusoidal background, negative s
sim.sin.negs2<-simDisturb(bkgrnd="sinusoidal",disturbtype="neg.s",nreps=5000,dist.start=366+30)
#plot(sim.sin.negs2[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.negs2)){
  dt<-sim.sin.negs2[[ii]]$key["dt"]
  sim.sin.negs2[[ii]]$key["dday"]<-sim.sin.negs2[[ii]]$key["dday"]-365
  sim.sin.negs2[[ii]]$key["rday"]<-sim.sin.negs2[[ii]]$key["rday"]-365
  
  refy<-sim.sin.negs2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.negs2[[ii]]$simts[(366/dt):length(sim.sin.negs2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.sin.negs2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.negs2[[ii]]$empalarm<-disturbalarm(sim.sin.negs2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.negs2[[ii]]$filtalarm<-alarmfilter(sim.sin.negs2[[ii]]$empalarm,dmin.dist=3)
  sim.sin.negs2[[ii]]$eval<-evalStats(sim.sin.negs2[[ii]]$key, sim.sin.negs2[[ii]]$empalarm, tol=5)
  sim.sin.negs2[[ii]]$filteval<-evalStats(sim.sin.negs2[[ii]]$key, sim.sin.negs2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.sin.negs2, file="sim_sin_negs_adapt_integral_20211011_t5.rds")

## Full tests: sinusoidal background, positive s
sim.sin.poss2<-simDisturb(bkgrnd="sinusoidal",disturbtype="pos.s",nreps=5000,dist.start=366+30)
#plot(sim.sin.poss2[[1]]$simts, type="l")

for(ii in 1:length(sim.sin.poss2)){
  dt<-sim.sin.poss2[[ii]]$key["dt"]
  sim.sin.poss2[[ii]]$key["dday"]<-sim.sin.poss2[[ii]]$key["dday"]-365
  sim.sin.poss2[[ii]]$key["rday"]<-sim.sin.poss2[[ii]]$key["rday"]-365
  
  refy<-sim.sin.poss2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.sin.poss2[[ii]]$simts[(366/dt):length(sim.sin.poss2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.sin.poss2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.sin.poss2[[ii]]$empalarm<-disturbalarm(sim.sin.poss2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.sin.poss2[[ii]]$filtalarm<-alarmfilter(sim.sin.poss2[[ii]]$empalarm,dmin.dist=3)
  sim.sin.poss2[[ii]]$eval<-evalStats(sim.sin.poss2[[ii]]$key, sim.sin.poss2[[ii]]$empalarm, tol=5)
  sim.sin.poss2[[ii]]$filteval<-evalStats(sim.sin.poss2[[ii]]$key, sim.sin.poss2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.sin.poss2, file="sim_sin_poss_adapt_integral_20211011_t5.rds")



## Full tests: flat background, negative wedge
sim.flat.negwedge2<-simDisturb(bkgrnd="flat",disturbtype="neg.wedge",nreps=5000,dist.start=366+refwidth/2) ##<---- REDO!!!!!!!
#plot(sim.flat.negwedge[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.negwedge2)){
  dt<-sim.flat.negwedge2[[ii]]$key["dt"]
  sim.flat.negwedge2[[ii]]$key["dday"]<-sim.flat.negwedge2[[ii]]$key["dday"]-365
  sim.flat.negwedge2[[ii]]$key["rday"]<-sim.flat.negwedge2[[ii]]$key["rday"]-365
  
  refy<-sim.flat.negwedge2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.negwedge2[[ii]]$simts[(366/dt):length(sim.flat.negwedge2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.negwedge2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.negwedge2[[ii]]$empalarm<-disturbalarm(sim.flat.negwedge2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.negwedge2[[ii]]$filtalarm<-alarmfilter(sim.flat.negwedge2[[ii]]$empalarm,dmin.dist=3)
  sim.flat.negwedge2[[ii]]$eval<-evalStats(sim.flat.negwedge2[[ii]]$key, sim.flat.negwedge2[[ii]]$empalarm, tol=5)
  sim.flat.negwedge2[[ii]]$filteval<-evalStats(sim.flat.negwedge2[[ii]]$key, sim.flat.negwedge2[[ii]]$filtalarm, tol=5)
}

saveRDS(sim.flat.negwedge2, file="sim_flat_negwedge_adapt_integral_20211011_t5.rds")


## Full tests: flat background, positive wedge
sim.flat.poswedge2<-simDisturb(bkgrnd="flat",disturbtype="pos.wedge",nreps=5000,dist.start=366+30)
#plot(sim.flat.poswedge2[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.poswedge2)){
  dt<-sim.flat.poswedge2[[ii]]$key["dt"]
  sim.flat.poswedge2[[ii]]$key["dday"]<-sim.flat.poswedge2[[ii]]$key["dday"]-365
  sim.flat.poswedge2[[ii]]$key["rday"]<-sim.flat.poswedge2[[ii]]$key["rday"]-365
  
  refy<-sim.flat.poswedge2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.poswedge2[[ii]]$simts[(366/dt):length(sim.flat.poswedge2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.poswedge2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.poswedge2[[ii]]$empalarm<-disturbalarm(sim.flat.poswedge2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.poswedge2[[ii]]$filtalarm<-alarmfilter(sim.flat.poswedge2[[ii]]$empalarm,dmin.dist=3)
  sim.flat.poswedge2[[ii]]$eval<-evalStats(sim.flat.poswedge2[[ii]]$key, sim.flat.poswedge2[[ii]]$empalarm, tol=5)
  sim.flat.poswedge2[[ii]]$filteval<-evalStats(sim.flat.poswedge2[[ii]]$key, sim.flat.poswedge2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.flat.poswedge2, file="sim_flat_poswedge_adapt_integral_20211011_t5.rds")

## Full tests: flat background, negative s
sim.flat.negs2<-simDisturb(bkgrnd="flat",disturbtype="neg.s",nreps=5000,dist.start=366+30)
#plot(sim.flat.negs2[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.negs2)){
  dt<-sim.flat.negs2[[ii]]$key["dt"]
  sim.flat.negs2[[ii]]$key["dday"]<-sim.flat.negs2[[ii]]$key["dday"]-365
  sim.flat.negs2[[ii]]$key["rday"]<-sim.flat.negs2[[ii]]$key["rday"]-365
  
  refy<-sim.flat.negs2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.negs2[[ii]]$simts[(366/dt):length(sim.flat.negs2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.negs2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.negs2[[ii]]$empalarm<-disturbalarm(sim.flat.negs2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.negs2[[ii]]$filtalarm<-alarmfilter(sim.flat.negs2[[ii]]$empalarm,dmin.dist=3)
  sim.flat.negs2[[ii]]$eval<-evalStats(sim.flat.negs2[[ii]]$key, sim.flat.negs2[[ii]]$empalarm, tol=5)
  sim.flat.negs2[[ii]]$filteval<-evalStats(sim.flat.negs2[[ii]]$key, sim.flat.negs2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.flat.negs2, file="sim_flat_negs_adapt_integral_20211011_t5.rds")

## Full tests: flat background, positive s
sim.flat.poss2<-simDisturb(bkgrnd="flat",disturbtype="pos.s",nreps=5000,dist.start=366+30)
#plot(sim.flat.poss2[[1]]$simts, type="l")

for(ii in 1:length(sim.flat.poss2)){
  dt<-sim.flat.poss2[[ii]]$key["dt"]
  sim.flat.poss2[[ii]]$key["dday"]<-sim.flat.poss2[[ii]]$key["dday"]-365
  sim.flat.poss2[[ii]]$key["rday"]<-sim.flat.poss2[[ii]]$key["rday"]-365
  
  refy<-sim.flat.poss2[[ii]]$simts[1:(365/dt)]
  refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
  testy<-sim.flat.poss2[[ii]]$simts[(366/dt):length(sim.flat.poss2[[ii]]$simts)]
  testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
  
  sim.flat.poss2[[ii]]$empdiff<-mwdistdiffz(testy=testy,refy=refy,refwidth=refwidth/dt,wwidth=5/dt, stride=6, ddiff_method="integral")
  sim.flat.poss2[[ii]]$empalarm<-disturbalarm(sim.flat.poss2[[ii]]$empdiff, dthresh=2, rthresh=0.5)
  sim.flat.poss2[[ii]]$filtalarm<-alarmfilter(sim.flat.poss2[[ii]]$empalarm,dmin.dist=3)
  sim.flat.poss2[[ii]]$eval<-evalStats(sim.flat.poss2[[ii]]$key, sim.flat.poss2[[ii]]$empalarm, tol=5)
  sim.flat.poss2[[ii]]$filteval<-evalStats(sim.flat.poss2[[ii]]$key, sim.flat.poss2[[ii]]$filtalarm, tol=5)
}
saveRDS(sim.flat.poss2, file="sim_flat_poss_adapt_integral_20211011_t5.rds")

