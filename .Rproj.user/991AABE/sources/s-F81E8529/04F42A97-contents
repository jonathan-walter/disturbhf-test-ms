## empirical examples of disturbhf methods

rm(list=ls())

#library(devtools) #only run for new/updated install of disturbhf
#install_github("jonathan-walter/disturbhf")

library(here)
library(disturbhf)
library(lubridate)



##-----------------------------------------------------------------------------
## NERRS data -- Apalachicola -- Real Date Oct. 10, 2018

apa.dat<-read.csv(here("./data/apa_michael.csv"), stringsAsFactors = F)

apa.dat$DateTimeFormatted<-as.POSIXct(apa.dat$DateTimeFormatted)

apa.dat<-apa.dat[!is.na(apa.dat$DateTimeFormatted),]

# par(mfrow=c(2,1))
# plot(apa.dat$DateTimeFormatted, apa.dat$Sal, type="l")
# plot(apa.dat$DateTimeFormatted, apa.dat$DO_pct, type="l")

dt=0.25/24 #15-minute data

sum(duplicated(apa.dat$DateTimeFormatted))
apa.dat[apa.dat$DateTimeFormatted %in% apa.dat$DateTimeFormatted[duplicated(apa.dat$DateTimeFormatted)],]

test.start<-as.POSIXct("2018-01-01")

# testy<-data.frame(tt=apa.dat$DateTimeFormatted[apa.dat$DateTimeFormatted>=test.start], 
#                   yy=apa.dat$Sal[apa.dat$DateTimeFormatted>=test.start])
# refy<-data.frame(tt=apa.dat$DateTimeFormatted[apa.dat$DateTimeFormatted<test.start], 
#                  yy=apa.dat$Sal[apa.dat$DateTimeFormatted<test.start])
# 
# refy<-refy[!is.na(refy$tt),]
# 
# dd.apa.Sal<-mwdistdiffz(testy, refy, wwidth=5/dt, refwidth=30/dt, dx=1, stride=12, dmin=0.5)
# alarm.apa.Sal<-disturbalarm(dd.apa.Sal, dthresh=1.5, rthresh=0.5)

testy<-data.frame(tt=apa.dat$DateTimeFormatted[apa.dat$DateTimeFormatted>=test.start], 
                  yy=apa.dat$DO_pct[apa.dat$DateTimeFormatted>=test.start])
refy<-data.frame(tt=apa.dat$DateTimeFormatted[apa.dat$DateTimeFormatted<test.start], 
                 yy=apa.dat$DO_pct[apa.dat$DateTimeFormatted<test.start])

dd.apa.DO_pct<-mwdistdiffz(testy, refy, wwidth=5/dt, refwidth=60/dt, dx=1, stride=12, dmin=0.5)
alarm.apa.DO_pct<-disturbalarm(dd.apa.DO_pct, dthresh=3, rthresh=0.5)

# par(mfrow=c(4,1), mar=c(4.1,4.1,1.1,1.1))
# 
# plot(apa.dat$DateTimeFormatted, apa.dat$Sal, type="l")
# abline(v=alarm.apa.Sal$dist.date, col="red")
# abline(v=alarm.apa.Sal$recov.date, col="blue")
# 
# plot(dd.apa.Sal$wleft, dd.apa.Sal$zz, type="l", xlim=range(apa.dat$DateTimeFormatted))
# 
# plot(apa.dat$DateTimeFormatted, apa.dat$DO_pct, type="l")
# abline(v=alarm.apa.DO_pct$dist.date, col="red")
# abline(v=alarm.apa.DO_pct$recov.date, col="blue")
# 
# plot(dd.apa.DO_pct$wleft, dd.apa.DO_pct$zz, type="l", xlim=range(apa.dat$DateTimeFormatted))

## nice plot

png("fig4_nerrs_casestudy.png",
     units="in", res=300, width=8.5, height=5)

par(mar=c(3,3,1.1,1.1), mfrow=c(2,1), mgp=c(1.7,0.5,0), tcl=-0.3)

plot(apa.dat$DateTimeFormatted, apa.dat$DO_pct, xlim=c(as.POSIXct("2017-01-01"),as.POSIXct("2018-01-01")),
     ylab="DO", xlab="2017", type="l", col="black")
mtext("Reference period", line=0.1)

plot(apa.dat$DateTimeFormatted, apa.dat$DO_pct, xlim=c(as.POSIXct("2018-01-01"),as.POSIXct("2019-01-01")),
     ylab="DO", xlab="2018", type="l", col="black")
mtext("Evaluation period", line=0.1)
rect(alarm.apa.DO_pct$dist.date, -10, alarm.apa.DO_pct$recov.date, 120, col="red", density=9)
abline(v=as.POSIXct("2018-10-10"), lwd=2.5, col="blue")

dev.off()


##-----------------------------------------------------------------------------
## Ameriflux data -- Shark River Slough

srs.dat<-read.csv(here("./data/AMF_US-Skr_BASE_HH_1-1.csv")
                  , skip=2, na.strings="-9999", stringsAsFactors=F)

srs.dat<-srs.dat[,!apply(srs.dat,2,function(x){all(is.na(x))})]

apply(srs.dat,2,function(x){sum(is.na(x))})

srs.dat<-srs.dat[ ,colnames(srs.dat) %in% c("TIMESTAMP_START","NEE_PI","LE","TA")]

srs.dat$Date<-as.POSIXct(paste(substr(srs.dat$TIMESTAMP_START,1,4)
                                 ,substr(srs.dat$TIMESTAMP_START,5,6)
                                 ,substr(srs.dat$TIMESTAMP_START,7,8)
                                 ,substr(srs.dat$TIMESTAMP_START,9,10)
                                 ,substr(srs.dat$TIMESTAMP_START,11,12)
                                 ,sep="-")
                           ,format="%Y-%m-%d-%H-%M")
srs.dat<-srs.dat[,colnames(srs.dat)!="TIMESTAMP_START"]

srs.dat<-srs.dat[-c(1:288),] #drop leading NAs of focal variables


#srs.dat<-srs.dat[year(srs.dat$Date)<=2014,]
srs.dat<-srs.dat[!is.na(srs.dat$Date),]
srs.dat<-srs.dat[srs.dat$Date >= as.POSIXct("2006-01-01"),]
srs.dat<-srs.dat[srs.dat$Date < as.POSIXct("2010-07-01"),]



# par(mfrow=c(3,1))
# plot(srs.dat$Date, srs.dat$LE, type="l")
# plot(srs.dat$Date, srs.dat$NEE_PI, type="l")
# plot(srs.dat$Date, srs.dat$TA, type="l")

dt=0.5/24 #half-hourly data
# testy<-data.frame(tt=srs.dat$Date, yy=srs.dat$LE)
# refy<-data.frame(tt=srs.dat$Date, yy=srs.dat$LE)
# 
# dd.srs.LE<-mwdistdiffz(testy, refy, wwidth=5/dt, refwidth=30/dt, dx=1, stride=6, dmin=0.5)
# alarm.srs.LE<-disturbalarm(dd.srs.LE, dthresh=2, rthresh=0.5)

test.start<-as.POSIXct("2009-07-01")


testy<-data.frame(tt=srs.dat$Date[srs.dat$Date >= test.start], yy=srs.dat$NEE_PI[srs.dat$Date >= test.start])
refy<-data.frame(tt=srs.dat$Date[srs.dat$Date < test.start], yy=srs.dat$NEE_PI[srs.dat$Date < test.start])

dd.srs.NEE<-mwdistdiffz(testy, refy, wwidth=5/dt, refwidth=60/dt, dx=1, stride=6, dmin=0.5)
alarm.srs.NEE<-disturbalarm(dd.srs.NEE, dthresh=3, rthresh=0.5)

testy<-data.frame(tt=srs.dat$Date[srs.dat$Date >= test.start], yy=srs.dat$TA[srs.dat$Date >= test.start])
refy<-data.frame(tt=srs.dat$Date[srs.dat$Date < test.start], yy=srs.dat$TA[srs.dat$Date < test.start])

dd.srs.TA<-mwdistdiffz(testy, refy, wwidth=5/dt, refwidth=60/dt, dx=1, stride=6, dmin=0.5)
alarm.srs.TA<-disturbalarm(dd.srs.TA, dthresh=3, rthresh=0.5)

par(mfrow=c(2,1))
# plot(srs.dat$Date, srs.dat$LE, type="l")
# abline(v=alarm.srs.LE$dist.date, col="red")
# abline(v=alarm.srs.LE$recov.date, col="blue")

plot(srs.dat$Date, srs.dat$NEE_PI, type="l")
abline(v=alarm.srs.NEE$dist.date, col="red")
abline(v=alarm.srs.NEE$recov.date, col="blue")

plot(srs.dat$Date, srs.dat$TA, type="l")
abline(v=alarm.srs.TA$dist.date, col="red")
abline(v=alarm.srs.TA$recov.date, col="blue")


srs.ref<-srs.dat[srs.dat$Date < test.start,]
srs.ref$DateAdj<-rep(NA, nrow(srs.ref))

flag<-NULL

for(ii in 1:nrow(srs.ref)){
  if(month(srs.ref$Date[ii])==2 & day(srs.ref$Date[ii])==29){next}
  if(yday(srs.ref$Date[ii])>=yday(as.POSIXct("2009-07-01"))){
    tmp<-as.character(srs.ref$Date[ii])
    substr(tmp,1,4)<-"2009"
    srs.ref$DateAdj[ii]<-as.POSIXct(tmp)
  }
  if(yday(srs.ref$Date[ii])<yday(as.POSIXct("2009-07-01"))){
    tmp<-as.character(srs.ref$Date[ii])
    substr(tmp,1,4)<-"2010"
    srs.ref$DateAdj[ii]<-as.POSIXct(tmp)
  }
  
}

srs.ref<-srs.ref[!is.na(srs.ref$DateAdj),]

# srs.ref$doy<-yday(srs.ref$Date)
# srs.ref<-srs.ref[srs.ref$doy != 366,]
# 
srs.ref.NEE.agg<-aggregate(NEE_PI ~ DateAdj, data=srs.ref, FUN=quantile, probs=c(0.25,0.75), type=4)
srs.ref.TA.agg<-aggregate(TA ~ DateAdj, data=srs.ref, FUN=quantile, probs=c(0.25,0.75), type=4)




# plot(NA,NA, xlim=c(as.POSIXct("2009-07-01"),as.POSIXct("2010-07-01")), 
#      ylim=c(min(srs.ref.NEE.agg$NEE_PI),max(srs.ref.NEE.agg$NEE_PI)))
# polygon(x=c(srs.ref.NEE.agg$DateAdj, rev(srs.ref.NEE.agg$DateAdj), srs.ref.NEE.agg$DateAdj[1]),
#         y=c(srs.ref.NEE.agg$NEE_PI[,2], rev(srs.ref.NEE.agg$NEE_PI[,1]), srs.ref.NEE.agg$NEE_PI[1,1]),
#         col="grey", border="grey"
# )



tiff("fig5_srs_casestudy.tif",
     units="in", res=300, width=8.5, height=5)

par(mar=c(3,3,1.1,1.1), mfrow=c(2,1), mgp=c(1.7,0.5,0), tcl=-0.3)

plot(srs.dat$Date, srs.dat$TA, xlim=c(as.POSIXct("2009-07-01"),as.POSIXct("2010-07-01")),
     ylab="Air temperature", xlab="2009-2010", type="l", col="black", xaxs="i")
polygon(x=c(srs.ref.TA.agg$DateAdj, rev(srs.ref.TA.agg$DateAdj), srs.ref.TA.agg$DateAdj[1]),
        y=c(srs.ref.TA.agg$TA[,2], rev(srs.ref.TA.agg$TA[,1]), srs.ref.TA.agg$TA[1,1]),
        col="grey", border="grey")
lines(srs.dat$Date, srs.dat$TA)

rect(alarm.srs.TA$dist.date, -10, alarm.srs.TA$recov.date, 120, col="red", density=9)


plot(srs.dat$Date, srs.dat$NEE_PI, xlim=c(as.POSIXct("2009-07-01"),as.POSIXct("2010-07-01")),
     ylab="NEE", xlab="2009-2010", type="l", col="black", xaxs="i")
polygon(x=c(srs.ref.NEE.agg$DateAdj, rev(srs.ref.NEE.agg$DateAdj), srs.ref.NEE.agg$DateAdj[1]),
        y=c(srs.ref.NEE.agg$NEE_PI[,2], rev(srs.ref.NEE.agg$NEE_PI[,1]), srs.ref.NEE.agg$NEE_PI[1,1]),
        col="grey", border="grey")
lines(srs.dat$Date, srs.dat$NEE_PI)
rect(alarm.srs.NEE$dist.date, -40, alarm.srs.NEE$recov.date, 120, col="red", density=9)


dev.off()



##-----------------------------------------------------------------------------
## Hudson river chlorophyll and rotifers

hdr.dat<-read.csv(here("./data/Hudson_Kingston_1987thru2015_reduced.csv")
                  , stringsAsFactors = F)

fix.date<-function(baddate){
  
  reform<-function(x){
    paste(x[1],
          x[2],
          paste(x[3]),sep="-")
  }
  
  tmp<-strsplit(baddate,"/")
  
  for(ii in 1:length(tmp)){
    
    if(as.numeric(tmp[[ii]][2]) < 15){
      tmp[[ii]][2]<-"1"
    }
    else{
      tmp[[ii]][2]<-"15"
    }
    if(tmp[[ii]][3] > 50){ #fix years to 4 digits
      tmp[[ii]][3] <- paste0("19",tmp[[ii]][3])
    }
    else{
      tmp[[ii]][3] <- paste0("20",tmp[[ii]][3])
    }
  }
  
  tmp<-unlist(lapply(tmp,reform))
  tmp<-as.POSIXct(tmp, format="%m-%d-%Y")

  return(tmp)
}


hdr.dat$DATE<-fix.date(hdr.dat$DATE)
hdr.dat<-hdr.dat[!duplicated(hdr.dat$DATE),]

par(mfrow=c(2,1), mar=c(4.1,4.1,1.1,1.1))
plot(hdr.dat$DATE, hdr.dat$CHL, type="l")
plot(hdr.dat$DATE, hdr.dat$ROT, type="l")


# testy<-data.frame(tt=hdr.dat$DATE, yy=hdr.dat$CHL)
# refy<-data.frame(tt=hdr.dat$DATE, yy=hdr.dat$CHL)
# 
# dd.hdr.Chl<-mwdistdiffz(testy, refy, wwidth=36, dx=1, stride=1, dmin=0.33)
# alarm.hdr.Chl<-disturbalarm(dd.hdr.Chl, dthresh=1, rthresh=0.5)

# par(mfrow=c(2,1), mar=c(4.1,4.1,1.1,1.1))
# plot(hdr.dat$DATE, hdr.dat$CHL, type="l")
# plot(dd.hdr.Chl$wleft, dd.hdr.Chl$zz)

test.start<-as.POSIXct("1992-01-01")


testy<-data.frame(tt=hdr.dat$DATE[hdr.dat$DATE>=test.start], yy=hdr.dat$ROT[hdr.dat$DATE>=test.start])
refy<-data.frame(tt=hdr.dat$DATE[hdr.dat$DATE<test.start], yy=hdr.dat$ROT[hdr.dat$DATE<test.start])

dd.hdr.Rot<-mwdistdiffz(testy, refy, wwidth=36, dx=1, stride=1, dmin=0.5)
alarm.hdr.Rot<-disturbalarm(dd.hdr.Rot, dthresh=2, rthresh=0.5)

par(mfrow=c(2,1), mar=c(4.1,4.1,1.1,1.1))
plot(hdr.dat$DATE, hdr.dat$ROT, type="l")
abline(v=alarm.hdr.Rot$dist.date, col="red")
abline(v=alarm.hdr.Rot$recov.date, col="blue")
plot(dd.hdr.Rot$wleft, dd.hdr.Rot$zz)


tiff("fig6_hudrot_casestudy.tif",
     units="in", res=300, width=8.5, height=2.5)

par(mar=c(4.1,4.1,1.1,1.1))

plot(hdr.dat$DATE, hdr.dat$ROT,
     ylab="Rotifer abundance", xlab="Date", type="l", col="black")

rect(alarm.hdr.Rot$dist.date, -1000, alarm.hdr.Rot$recov.date, 3000, col="red", density=9)

abline(v=as.POSIXct("1992-10-01"), lwd=2.5, col="blue")

dev.off()

