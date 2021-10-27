doCUSUM <- function(tt, simts, alpha=0.05){
  
  require(strucchange)
  
  cusum <- efp(simts~1, type="OLS-CUSUM")
  p.cusum<-sctest(cusum)
  
  if(p.cusum$p.value >= alpha){
    out <- data.frame(dist.date=NA,
                      recov.date=NA,
                      tdiff=NA,
                      peakz=NA,
                      peak.date=NA)
  }
  else{
    if(abs(max(cusum$process)) > abs(min(cusum$process))){
      bdy<-boundary(cusum)[1]
      d.index <- min(which(cusum$process>bdy))
      postdist <- cusum$process[d.index:length(cusum$process)]
      r.index <- min(which(postdist <= bdy/2)) + d.index
      
      out <- data.frame(dist.date=tt[d.index],
                        recov.date=tt[r.index],
                        tdiff=tt[r.index]-tt[d.index],
                        peakz=max(cusum$process),
                        peak.date=tt[which.max(cusum$process)])
    }
    if(abs(max(cusum$process)) < abs(min(cusum$process))){
      bdy <- boundary(cusum)[1]*-1
      d.index <- min(which(cusum$process<bdy))
      postdist <- cusum$process[d.index:length(cusum$process)]
      r.index <- min(which(postdist >= bdy/2)) + d.index
      
      out <- data.frame(dist.date=tt[d.index],
                        recov.date=tt[r.index],
                        tdiff=tt[r.index]-tt[d.index],
                        peakz=max(cusum$process),
                        peak.date=tt[which.max(cusum$process)])
      
      
    }

  }
  return(out)
}