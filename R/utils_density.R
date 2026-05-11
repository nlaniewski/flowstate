peak.bounds <- function(x, adjust = 1, height.threshold = 0.01, plot = F){
  d <- stats::density(x, adjust = adjust)
  pv <- diff(sign(diff(d$y)))
  peaks.i <- which(pv == -2) + 1
  peak.major.i <- peaks.i[which.max(d$y[peaks.i])]
  ##
  threshold.value <- d$y[peak.major.i] * height.threshold
  ##
  bounds.i <- sapply(c('left','right'), function(side){
    i <- peak.major.i
    while(d$y[i] > threshold.value){
      i <- get(if(side=='left'){"-"}else if(side=='right'){"+"})(i,1)
      if(pv[i] != 2){
        i
      }else{
        break
      }
    }
    return(i)
  })
  ##
  if(plot){
    plot(d)
    graphics::abline(v = d$x[peak.major.i], col = "red")
    graphics::abline(v = d$x[bounds.i], col = "red", lty = "dashed")
  }
  ##
  return(d$x[bounds.i])
}

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

kde.contour <- function (
    x,
    y,
    bandwidth.adjust = 4,
    grid.points = 50,
    threshold = 0.75,
    plot = F,
    ...
)
{
  ## pre-calculate bandwidth
  h <- sapply(list(x, y), function(j) MASS::bandwidth.nrd(j))
  ## adjust bandwidth
  if (!is.null(bandwidth.adjust)) {
    h <- h * bandwidth.adjust
  }
  ## kernel density estimate
  dens <- MASS::kde2d(x = x, y = y, h = h, n = grid.points)
  ## contour lines
  cl <- grDevices::contourLines(dens)
  ## set an index; start with outermost contour
  i <- 1
  ## calculate the fraction of total points within the indexed contour;
  ## repeat (move inward) until defined 'threshold' is met and return the index
  repeat{
    in.contour <- sp::point.in.polygon(
      point.x = x,
      point.y = y,
      pol.x = cl[[i]]$x,
      pol.y = cl[[i]]$y
    )
    fraction.in.contour <- sum(in.contour) / length(in.contour)
    if(fraction.in.contour >= threshold){
      i <- i + 1
    }else{
      i <- i - 1
      break
    }
  }
  ## indexed contour line
  cl.i <- cl[[i]]
  ## visualize
  if (plot) {
    p <- ggplot2::ggplot(data = NULL) +
      ggplot2::geom_hex(ggplot2::aes(x = x, y = y), bins = 200) +
      viridis::scale_fill_viridis(
        option = 'plasma',
        limits = c(0, 5),
        oob = scales::squish) +
      ggplot2::geom_path(ggplot2::aes(x = cl.i$x, y = cl.i$y), linewidth = 1, color = "red") +
      ggplot2::guides(fill = 'none') +
      ggplot2::labs(
        x = sub("\\[top.i\\]", "",deparse(substitute(x))),
        y = sub("\\[top.i\\]", "",deparse(substitute(y)))
      ) +
      ggplot2::xlim(0, 4E6) +
      ggplot2::ylim(0, 4E6)
    suppressWarnings(print(p))
  }
  ## return the indexed contour line
  return(cl.i[c("x", "y")])
}
