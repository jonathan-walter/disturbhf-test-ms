#' Simulates time series including with disturbances
#' 
#' @param tmax number of days in simulation, defaults to 730 days (2 years)
#' @param dt time interval between observations, defaults to 1/24 (hourly data)
#' @param bkgrnd background variability type. one of "flat" or "sinusoidal"
#' @param disturbtype disturbance type. one of "pos.wedge","neg.wedge,"pos.s","neg.s"
#' @param dist.start at what timestep can random disturbances start. Defaults to day 366
#' @param recov.lim lower and upper limits on recovery time. Defaults to c(3,60)
#' @param sev.lim lower and upper limits on disturbance severity; 
#' relative to variance of background time series
#' @param k.var constants controlling variance of the seasonal, diurnal, and noise components
#' @param nreps number of simulated time series to produce
#' 
#' @return \code{simDisturb} returns a list with \code{nreps} elements, each of which is
#' itself a list having the elements:
#' \item{key} a vector containing \code{dt}, the true disturbance start date, recovery date, recovery time, and severity
#' \item{simts} a vector containing the simulated time series
#' 
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#' 
#' @export

## Function that generates simulated time series to apply algorithms to
simDisturb<-function(tmax=730,dt=1/24,bkgrnd,disturbtype,dist.start=366,recov.lim=c(3,60),sev.lim=c(0,5),k.var=c(0.5,1,0.5),nreps=100){

  #some basic error handling
  if(!bkgrnd %in% c("flat","sinusoidal")){
    stop("bkgrnd must be one of 'flat' or 'sinusoidal'.")
  }
  if(!disturbtype %in% c("pos.wedge","neg.wedge","pos.s","neg.s")){
    stop("disturbtype must be one of 'pos.wedge','neg.wedge','pos.s', or 'neg.s'.")
  }
  
  tt=1:(tmax/dt) #vector of timesteps
  ref.end<-dist.start-1
  out<-list()
  
  for(rep in 1:nreps){
    
    dday<-sample(dist.start:(tmax-recov.lim[2]-1), 1) #random day of disturbance in 2nd year
    ddur<-sample(recov.lim[1]:recov.lim[2], 1) #random duration of disturbance
    sev<-runif(1, sev.lim[1], sev.lim[2])
    
    if(bkgrnd=="flat"){
      #make background time series
      bts<-k.var[2]*cos(seq(0, 2*pi*tmax, length.out=length(tt)))
      bts<-bts + k.var[3]*rnorm(length(tt))
      
      #make disturbance time series
      dts<-rep(0,length(tt))
      if(disturbtype=="pos.wedge"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-seq(sev*var(bts),0,length.out=ddur/dt)
      }
      if(disturbtype=="neg.wedge"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-seq(-1*sev*var(bts),0,length.out=ddur/dt)
      }
      if(disturbtype=="pos.s"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-sev*var(bts)*sin(seq(0,2*pi,length.out=ddur/dt))  
      }
      if(disturbtype=="neg.s"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<--1*sev*var(bts)*sin(seq(0,2*pi,length.out=ddur/dt))  
      }
      tts<-bts+dts
    }
    if(bkgrnd=="sinusoidal"){
      #make background time series
      bts<-k.var[1]*sin(seq(0,4*pi,length.out=length(tt)))
      bts<-bts + k.var[2]*cos(seq(0, 2*pi*tmax, length.out=length(tt)))
      bts<-bts + k.var[3]*rnorm(length(tt))
      
      #make disturbance time series
      dts<-rep(0,length(tt))
      if(disturbtype=="pos.wedge"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-seq(sev*var(bts),0,length.out=ddur/dt)
      }
      if(disturbtype=="neg.wedge"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-seq(-1*sev*var(bts),0,length.out=ddur/dt)
      }
      if(disturbtype=="pos.s"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<-sev*var(bts)*sin(seq(0,2*pi,length.out=ddur/dt))  
      }
      if(disturbtype=="neg.s"){
        dts[tt >= dday/dt & tt < (dday+ddur)/dt]<--1*sev*var(bts)*sin(seq(0,2*pi,length.out=ddur/dt))  
      }
      tts<-bts+dts
    }
    key<-c(dt, ref.end, dday, dday+ddur, ddur, sev)
    names(key)<-c("dt","ref.end","dday","rday","recovtime","severity")
    out[[rep]]<-list(key=key, simts=tts)
    
  }#end for loop
  return(out)
}#close function