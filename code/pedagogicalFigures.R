## Pedagogical figures

rm(list=ls())

#library(devtools) #only run for new/updated install of disturbhf
#install_github("jonathan-walter/disturbhf")

library(here)
library(viridis)
library(disturbhf)

#Create 4 simulated time series

source(here("./code/simdisturb.R"))

set.seed(11)
ex.flat.negwedge<-simDisturb(bkgrnd="flat", disturbtype="neg.wedge", recov.lim=c(44,45), sev.lim=c(3,3), nreps=1)[[1]]$simts
tt=1:length(ex.flat.negwedge)*(1/24)
ex.flat.negs<-simDisturb(bkgrnd="flat", disturbtype="neg.s", recov.lim=c(44,45), sev.lim=c(3,3), nreps=1)[[1]]$simts
ex.sine.negwedge<-simDisturb(bkgrnd="sinusoidal", disturbtype="neg.wedge", recov.lim=c(44,45), sev.lim=c(3,3), nreps=1)[[1]]$simts
ex.sine.negs<-simDisturb(bkgrnd="sinusoidal", disturbtype="neg.s", recov.lim=c(44,45), sev.lim=c(3,3), nreps=1)[[1]]$simts


#Illustrate how method works
laymat<-matrix(c(1,1,2,3,4,4), nrow=3, ncol=2, byrow=TRUE)

eval<-seq(-4,2.5,0.05)
ref.ecdf<-ecdf(ex.flat.negwedge[tt<366])
tnorm<-c(400,405)
tdist<-c(620,625)
norm.ecdf<-ecdf(ex.flat.negwedge[tt>=tnorm[1] & tt < tnorm[2]])
dist.ecdf<-ecdf(ex.flat.negwedge[tt>=tdist[1] & tt < tdist[2]])

max(norm.ecdf(eval)-ref.ecdf(eval))
max(dist.ecdf(eval)-ref.ecdf(eval))
xline<-which.max(dist.ecdf(eval)-ref.ecdf(eval))

result<-mwdistdiffz(testy = data.frame(tt=tt[tt>=366], yy=ex.flat.negwedge[tt>=366]),
                    refy = data.frame(tt=tt[tt<366], yy=ex.flat.negwedge[tt<366]),
                    wwidth=5/(1/24),
                    stride=6
                    )
alarm<-disturbalarm(result)

png("fig1_methodIllustration.png", units="in", res=300, height=4.5, width=6)
par(mar=c(3.1,3.1,0.5,0.5), mgp=c(1.8,0.5,0), tcl=-0.2)
layout(laymat)

plot(tt, ex.flat.negwedge, type="l", ylim=c(-4.2,4.2), xlab="Time", ylab="x")
abline(v=366, col="grey")
rect(xleft=tnorm[1], xright=tnorm[2], ybottom=-5, ytop=5, density=15, col="blue")
rect(xleft=tdist[1], xright=tdist[2], ybottom=-5, ytop=5, density=15, col="red")
text(par("usr")[1]+0.025*diff(par("usr")[1:2]),
     par("usr")[4]-0.075*diff(par("usr")[3:4]),
     "a)")

#mtext("Flat background, wedge disturbance", cex=0.8)
plot(eval, ref.ecdf(eval), type="l", xlab="x", ylab="ECDF")
lines(eval, norm.ecdf(eval), col="blue")
text(-2,0.9,expression(paste(d[w]," = 0.052")))
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.075*diff(par("usr")[3:4]),
     "b)")

plot(eval, ref.ecdf(eval), type="l", xlab="x", ylab="ECDF")
lines(eval, dist.ecdf(eval), col="red")
text(-2,0.9,expression(paste(d[w]," = 0.760")))
segments(x0=eval[xline], x1=eval[xline], y0=ref.ecdf(eval)[xline], y1=dist.ecdf(eval)[xline],
         col="grey",lty=2)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.075*diff(par("usr")[3:4]),
     "c)")

plot(result$wleft+2.5,result$zz, type="l", xlab="Time", ylab="Z-score")
abline(v=612.9583, col="grey", lty=2)
abline(v=660, col="grey", lty=2)
text(par("usr")[1]+0.025*diff(par("usr")[1:2]),
     par("usr")[4]-0.075*diff(par("usr")[3:4]),
     "d)")

dev.off()

png("fig2_simulationExamples.png", units="in", res=300, height=6.5, width=6.5)
par(mfrow=c(4,1), mar=c(3.1,2.9,1.5,1.1), oma=c(1.1,0,0,0))

plot(tt, ex.flat.negwedge, type="l", ylim=c(-4.2,4.2))
abline(v=366, col="grey")
mtext("Flat background, wedge disturbance", cex=0.8)
plot(tt, ex.flat.negs, type="l", ylim=c(-4.2,4.2))
abline(v=366, col="grey")
mtext("Flat background, s-shaped disturbance", cex=0.8)
plot(tt, ex.sine.negwedge, type="l", ylim=c(-4.2,4.2))
abline(v=366, col="grey")
mtext("Sinusoidal background, wedge disturbance", cex=0.8)
plot(tt, ex.sine.negs, type="l", ylim=c(-4.2,4.2))
abline(v=366, col="grey")
mtext("Sinusoidal background, s-shaped disturbance", cex=0.8)

mtext("Time",1,outer=T,line=-0.3)

dev.off()


