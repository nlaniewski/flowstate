peaks.valleys <- function(x,trim=TRUE,plot=FALSE){
  d <- stats::density(x)
  ##
  diff.1 <- diff(d$y)
  diff.2 <- diff(sign(diff.1))
  ##
  peaks <- which(diff.2 == -2) + 1
  valleys <- which(diff.2 == 2) + 1
  ##
  if(trim){
    peaks <- peaks[d$y[peaks] > mean(d$y)]
    valleys <- valleys[d$y[valleys] > mean(d$y)]
  }
  ##
  if(plot){
    plot(d)
    graphics::abline(v = d$x[peaks],col="blue",lty="dashed")
    graphics::abline(v = d$x[valleys],col="red",lty="dashed")
  }else{
    return(list(d = d,peaks = peaks,valleys = valleys))
  }
}

peaks.cut <- function(peaks.valleys.result,threshold = 0.001,plot=FALSE){
  if(length(peaks.valleys.result$peaks)==2){
    cuts <- c("left"=1,"right"=2)
  }
  ##
  cuts <- sapply(names(cuts),function(side){
    i <- sort(peaks.valleys.result$peaks)[cuts[[side]]]
    while(peaks.valleys.result$d$y[i] > threshold){
      i <- get(if(side=='left'){"-"}else if(side=='right'){"+"})(i,1)
    }
    return(i)
  })
  ##
  if(plot){
    plot(peaks.valleys.result$d)
    graphics::abline(v = peaks.valleys.result$d$x[cuts],col="red",lty="dashed")
  }
  ##
  return(peaks.valleys.result$d$x[cuts])
}

peak.height.cut<-function(j,height.value=.5,plot=F){
  ##density
  d<-stats::density(j)
  ##peak values
  i <- which.max(d$y)
  peak.max.y <- d$y[i]
  peak.max.x <- d$x[i]
  ##height adjust
  i.height <- which.min(abs(d$y - peak.max.y*height.value))
  i.adjust<-abs(i-i.height)
  bins<-c(i-i.adjust,i+i.adjust)
  if(plot){
    ##plot derivatives
    plot(d)
    graphics::abline(v=peak.max.x)
    graphics::abline(v=d$x[bins],col='red',lty='dotted')
    graphics::abline(h=d$y[i.height],col = 'blue',lty = 'dashed')
  }
  ##return bin values
  d$x[bins]
}

peak.values<-function(j,return.n.peaks=2,plot=F){
  d<-stats::density(j)
  peaks.i<-which(diff(sign(diff(d$y)))==-2)+1
  peaks.i<-peaks.i[order(d$y[peaks.i],decreasing = T)][seq(return.n.peaks)]
  peaks.values<-d$x[peaks.i]
  if(plot){
    plot(d)
    graphics::abline(v=peaks.values)
  }
  return(peaks.values)
}

peak.height.cuts<-function(j,peak.trim=TRUE,which.peak=NULL,height.value=.5,plot=FALSE){
  ##density
  d<-stats::density(j)
  ##peaks
  peaks.i<-which(diff(sign(diff(d$y)))==-2)+1
  if(peak.trim){
    peaks.i <- peaks.i[which(d$y[peaks.i] > mean(d$y))]
  }
  ##peak-based cuts; store values in a list
  cuts<-lapply(peaks.i,function(peak){
    peak.height.adjust <- d$y[peak]*height.value
    i<-peak
    while(!d$y[i]<peak.height.adjust){
      i<-i-1
    }
    cut.left<-d$x[i]
    ##
    i<-peak
    while(d$y[i]>peak.height.adjust){
      i<-i+1
    }
    cut.right<-d$x[i]
    ##
    return(c(cut.left,cut.right))
  })
  ##
  if(!is.null(which.peak)){
    cuts<-cuts[which.peak]
  }
  ##plot results
  if(plot){
    plot(d)
    graphics::abline(v=unlist(cuts),col = 'red', lty = 'dotted')
  }
  ##
  return(cuts)
}
