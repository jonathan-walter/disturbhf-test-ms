## Simulation tests of disturbance and recovery detection algorithm

rm(list=ls())

#library(devtools) #only run for new/updated install of disturbhf
#install_github("jonathan-walter/disturbhf")

library(here)
library(disturbhf)
library(heatwaveR)

source(here("./code/simdisturb.R"))
source(here("./code/evalstats.R"))
source(here("./code/doCUSUM.R"))


# ## Demo: flat background, negative wedge 
# sim.demo<-simDisturb(bkgrnd="flat",disturbtype="neg.wedge",nreps=5)
# plot(sim.demo[[1]]$simts, type="l")
# 
# for(ii in 1:length(sim.demo)){
#   dt<-sim.demo[[ii]]$key["dt"]
#   sim.demo[[ii]]$key["dday"]<-sim.demo[[ii]]$key["dday"]-sim.demo[[ii]]$key["ref.end"]
#   sim.demo[[ii]]$key["rday"]<-sim.demo[[ii]]$key["rday"]-sim.demo[[ii]]$key["ref.end"]
#   refy<-sim.demo[[ii]]$simts[1:(sim.demo[[ii]]$key["ref.end"]/dt)]
#   refy<-data.frame(tt=(1:length(refy))*dt, yy=refy)
#   testy<-sim.demo[[ii]]$simts[(sim.demo[[ii]]$key["ref.end"]/dt+1):length(sim.demo[[ii]]$simts)]
#   testy<-data.frame(tt=(1:length(testy))*dt, yy=testy)
#   
#   sim.demo[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
#   sim.demo[[ii]]$eval<-evalStats(sim.demo[[ii]]$key, sim.demo[[ii]]$empalarm, tol=5)
#   
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

  sim.flat.negwedge[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.flat.negwedge[[ii]]$eval<-evalStats(sim.flat.negwedge[[ii]]$key, sim.flat.negwedge[[ii]]$empalarm, tol=5)
}

print(Sys.time()-sttime)
saveRDS(sim.flat.negwedge, file="sim_flat_negwedge_cusum_20211016_t5.rds")

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
  
  sim.flat.poswedge[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.flat.poswedge[[ii]]$eval<-evalStats(sim.flat.poswedge[[ii]]$key, sim.flat.poswedge[[ii]]$empalarm, tol=5)
}
saveRDS(sim.flat.poswedge, file="sim_flat_poswedge_cusum_20211016_t5.rds")

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

  sim.flat.negs[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.flat.negs[[ii]]$eval<-evalStats(sim.flat.negs[[ii]]$key, sim.flat.negs[[ii]]$empalarm, tol=5)
}
saveRDS(sim.flat.negs, file="sim_flat_negs_cusum_20211016_t5.rds")

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
  
  sim.flat.poss[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.flat.poss[[ii]]$eval<-evalStats(sim.flat.poss[[ii]]$key, sim.flat.poss[[ii]]$empalarm, tol=5)
}
saveRDS(sim.flat.poss, file="sim_flat_poss_cusum_20211016_t5.rds")

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

  sim.sin.negwedge[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.sin.negwedge[[ii]]$eval<-evalStats(sim.sin.negwedge[[ii]]$key, sim.sin.negwedge[[ii]]$empalarm, tol=5)
}
saveRDS(sim.sin.negwedge, file="sim_sin_negwedge_cusum_20211016_t5.rds")

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
  
  sim.sin.poswedge[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.sin.poswedge[[ii]]$eval<-evalStats(sim.sin.poswedge[[ii]]$key, sim.sin.poswedge[[ii]]$empalarm, tol=5)
}
saveRDS(sim.sin.poswedge, file="sim_sin_poswedge_cusum_20211016_t5.rds")

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
 
  sim.sin.negs[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.sin.negs[[ii]]$eval<-evalStats(sim.sin.negs[[ii]]$key, sim.sin.negs[[ii]]$empalarm, tol=5)
}
saveRDS(sim.sin.negs, file="sim_sin_negs_cusum_20211016_t5.rds")

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
  
  sim.sin.poss[[ii]]$empalarm <- doCUSUM(testy$tt, testy$yy)
  sim.sin.poss[[ii]]$eval<-evalStats(sim.sin.poss[[ii]]$key, sim.sin.poss[[ii]]$empalarm, tol=5)
}
saveRDS(sim.sin.poss, file="sim_sin_poss_cusum_20211016_t5.rds")

