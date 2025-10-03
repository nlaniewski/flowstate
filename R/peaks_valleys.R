peaks.valley <- function(x,n.peaks=2,plot=FALSE,...){
  d <- stats::density(x)
  ##
  diff.1 <- diff(d$y)
  diff.2 <- diff(sign(diff.1))
  ##
  peaks <- which(diff.2 == -2) + 1
  valleys <- which(diff.2 == 2) + 1
  ##
  peaks <- peaks[order(d$y[peaks],decreasing = T)[seq(n.peaks)]]
  if(n.peaks == 2){
    valley <- valleys[data.table::`%between%`(valleys,peaks)]
  }
  ##
  if(plot){
    plot(d,...)
    graphics::abline(v = d$x[peaks],col="blue",lty="dashed")
    if(n.peaks == 2){
      graphics::abline(v = d$x[valley],col="red",lty="dashed")
    }
  }
  ##
  res <- list(d = d,peaks = peaks)
  if(n.peaks == 2){
    res$valley <- valley
  }
  return(res)
}
peak.height.cut <- function(x,which.peak=1,height.cut = 0.25,bisect.peak=NULL,plot=FALSE,...){
  res <- peaks.valley(x)
  if(which.peak==1){
    res$peaks <- res$peaks[1]
  }else if(which.peak == 2){
    res$peaks <- res$peaks[2]
  }
  threshold.value <- res$d$y[res$peaks] * height.cut
  cuts <- sapply(c('left','right'),function(side){
    i <- res$peaks
    while(res$d$y[i] >= threshold.value){
      i <- get(if(side=='left'){"-"}else if(side=='right'){"+"})(i,1)
    }
    return(i)
  })
  if(!is.null(bisect.peak)){
    side <- match.arg(bisect.peak,choices = c('left','right'))
    if(side == 'left'){
      cuts[['right']] <- res$peaks[1]
    }else if(side == 'right'){
      cuts[['left']] <- res$peaks[1]
    }
  }
  if(plot){
    plot(res$d,...)
    graphics::abline(v = res$d$x[res$peaks],col="blue",lty="dashed")
    graphics::abline(v = res$d$x[cuts],col="red",lty="dashed")
  }
  return(res$d$x[cuts])
}

