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
spill.matrix.from.controls<-function(
    batch.specific.control.paths,
    sample.id = '$FIL',
    pattern.sub = NULL,#" |-|Stained|Control",
    scatter.select.ignore = NULL,#list(VioletG=c("FSC_A"))
    fluor.pattern = "Blue|Violet|Red|Green",
    universal.unstained = TRUE,
    internal.negative = NULL,
    plot = FALSE,
    return.spill = FALSE
){
  ##flowstate objects; single-stained controls; concatenated; no transformation
  fs<-flowstate::read.flowstate(
    batch.specific.control.paths,
    colnames.type = "N",
    cofactor = NULL,
    sample.id = sample.id,
    concatenate = TRUE
  )
  ##pattern substitution for 'sample.id'
  if(!is.null(pattern.sub)){
    data.table::set(
      x = fs$data,
      j = "sample.id",
      value = as.factor(
        fs$keywords[
          ,
          rep(gsub(pattern.sub, "", j), as.numeric(`$TOT`)),
          env = list(j = sample.id)
        ]
      )
    )
    if(!fs$data[,all(grep("Unstained",levels(sample.id),value = T,invert = T) %in% names(.SD))]){
      stop("sample.id/control names do not fully match to data column names;
           modify 'pattern.sub' to acheive equivalence.")
    }
  }
  ##density distributions to find singlet bead population;
  ##assumes singlet beads are represented by the highest/most dense peak;
  ##cut FSC and SSC peaks (at % height = x) to isolate singlet beads using 'peak.height.cuts' function
  ##works for all but 'ARC LIVEDEAD (Violet G)' beads
  fs$data[
    # i = sample.id != "VioletG550/40",
    j = scatter.select :=
      data.table::`%between%`(FSC_A,peak.height.cuts(FSC_A,height.value = 0.25,which.peak = 1)[[1]]) &
      data.table::`%between%`(SSC_A,peak.height.cuts(SSC_A,height.value = 0.25,which.peak = 1)[[1]]),
    by = sample.id
  ]
  ##scatter select override for the case where a 'dominant/singular' bead peak cannot be defined
  if(!is.null(scatter.select.ignore)){
    for(detector in names(scatter.select.ignore)){
      ##reset scatter.select to FALSE
      fs$data[grepl(detector,sample.id),scatter.select := FALSE]
      ##re-index, ignoring the defined scatter parameter
      for(scatter in c('FSC_A','SSC_A')[!c('FSC_A','SSC_A') %in% scatter.select.ignore[[detector]]]){
        data.table::set(
          x = fs$data,
          i = fs$data[,.I[grepl(detector,sample.id)]],
          j = 'scatter.select',
          value = fs$data[
            i = grepl(detector,sample.id),
            j = data.table::`%between%`(j,peak.height.cuts(j,height.value = 0.25,which.peak = 1)[[1]]),
            env = list(j = scatter)
          ]
        )
      }
    }
  }
  ##plot results for 'peak.height.cuts'
  if(plot){
    dt.melt <- data.table::melt(
      fs$data[,.(sample.id,scatter.select,FSC_A,SSC_A)],
      measure.vars = c('FSC_A','SSC_A')
    )
    ##
    p.list <- sapply(c('FSC_A','SSC_A'),function(scatter){
      ggplot2::ggplot(
        data = dt.melt[variable == scatter],
        mapping = ggplot2::aes(value)
      ) +
        ggplot2::geom_density() +
        ggplot2::facet_wrap(~sample.id,scales = 'free') +
        ggplot2::theme(
          axis.text = ggplot2::element_blank()
        ) +
        ggplot2::geom_vline(
          data = dt.melt[
            i = variable == scatter & (scatter.select),
            j = range(value),
            by=sample.id
          ],
          mapping = ggplot2::aes(xintercept = V1),
          linetype = 'dashed',
          color = "red"
        ) +
        ggplot2::labs(x = scatter)
    },simplify = F)
    print(p.list)
  }
  ##fluorescence detectors (PMTs) used to measure fluorescence values
  fluors<-grep(fluor.pattern,names(fs$data),value = T)
  ##match the order of the controls/files
  fluors<-fluors[stats::na.omit(match(fs$data[,levels(sample.id)],fluors))]
  ##mean to find the peak detector for each fluor/control
  ##test for equivalency between control name and peak detector (highest mean value)
  res<-fs$data[
    i = sample.id %in% fluors & (scatter.select),
    j = .(detector.peak = names(which.max(sapply(.SD,mean)))),
    .SDcols = fluors,
    by = sample.id
  ][,all(sample.id==detector.peak)]
  if(!res){stop("Control names do not match peak detector names.")}
  ##density distribution of peak detector using selected scatter events (singlets);
  ##assumes two major peaks (negative/positive);
  ##bisect positive peak using return value from 'peak.height.cuts' function; second peak
  fs$data[,detector.select := FALSE]
  for(fluor in fluors){
    data.table::set(
      x = fs$data,
      i = fs$data[,.I[sample.id %in% fluor & (scatter.select)]],
      j = 'detector.select',
      value =fs$data[
        i = sample.id %in% fluor & (scatter.select),
        j = data.table::`%between%`(j,peak.height.cuts(j,which.peak = 2)[[1]]),
        env = list(j = fluor)
      ]
    )
  }
  ##plot results
  if(plot){
    dt.melt <- data.table::rbindlist(
      lapply(split(fs$data[sample.id %in% fluors],by='sample.id',drop = T),function(dt){
        dt[
          i = (scatter.select),
          j = .(value = j,detector.select,sample.id),
          env = list(j = as.character(unique(dt[['sample.id']])))
        ]
      })
    )
    ##
    p <- ggplot2::ggplot(
      data = dt.melt,
      mapping = ggplot2::aes(value)
    ) +
      ggplot2::geom_density() +
      ggplot2::facet_wrap(~sample.id,scales = 'free') +
      ggplot2::theme(
        axis.text = ggplot2::element_blank()
      ) +
      ggplot2::geom_vline(
        data = dt.melt[(detector.select),range(value),by=sample.id],
        mapping = ggplot2::aes(xintercept = V1),
        linetype = 'dashed',
        color = "red"
      )
    print(p)
  }
  ##detector medians using selected detector positive events
  ##converted to matrix -> spill
  spill <- as.matrix(fs$data[
    i = (detector.select),
    j = lapply(.SD,stats::median),
    .SDcols = fluors,
    by=sample.id
  ][,!'sample.id'])
  ##background/unstained medians
  ##converted to matrix
  if(universal.unstained){
    ##background/unstained fluorescence
    unstained.medians <- as.matrix(fs$data[
      i = (scatter.select) & grepl("Unstained",sample.id),
      j = lapply(.SD,stats::median),
      .SDcols = fluors,
      by = sample.id
    ][,!'sample.id'])
    ##subtract background/unstained median fluorescence
    ##negative values zeroed out
    spill <- pmax(sweep(spill,2,unstained.medians),0)
  }
  ##normalize; row-wise; x/max(x)
  spill <- t(apply(spill,1,function(x){x/max(x)}))
  ##matrix to data.table; assign spill to fs$spill list element
  fs$spill <- data.table::as.data.table(spill)
  ##return
  if(return.spill){
    data.table::setnames(fs$spill,new=fs$parameters[,N[match(names(fs$spill),N.alias)]])
    return(fs$spill)
  }else{
    fs[]
  }
}
