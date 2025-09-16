raw.fluorescence.check<-function(fcs.file.path){
  pars<-readFCStext(fcs.file.path) |>
    parameters.to.data.table()
  if("TYPE" %in% names(pars)){
    if(!"Raw_Fluorescence" %in% pars[['TYPE']]){
      stop(
        paste(
          fcs.file.path,
          "'Raw_Fluorescence' value not found in 'TYPE' keyword; is this a reference control?",
          sep = "\n"
        )
      )
    }
  }else{
    stop(
      paste(
        fcs.file.path,
        "'TYPE' keyword not found; is this a full spectrum .fcs file (version 3.1)?",
        sep = "\n"
      )
    )
  }
}

spectral.signature.medians <- function(fcs.file.path,subtract.internal.negative=FALSE){
  ##test;
  ##stop if not FCS 3.1 (has a 'TYPE' keyword);
  ##stop if not a reference control ('TYPE' == 'Raw_Fluorescence')
  raw.fluorescence.check(fcs.file.path)
  ##flowstate object; no transformation
  fs<-read.flowstate(
    fcs.file.path,
    colnames.type = "N"
  )
  ##density distributions to find singlet bead population;
  ##assumes singlet beads are represented by the highest/most dense peak;
  ##cut FSC and SSC peaks (at half-height) to isolate singlet beads using 'peak.height.cut' function
  fs$data[
    ,
    j = scatter.select :=
      data.table::`%between%`(FSC_A,peak.height.cut(FSC_A)) &
      data.table::`%between%`(SSC_A,peak.height.cut(SSC_A)),
    by = sample.id
  ]
  ##fluorescence detectors used to measure raw fluorescence values
  detectors.fluors<-fs$parameters[TYPE=="Raw_Fluorescence",N.alias]
  ##mean to find the peak detector
  detector.peak<-fs$data[
    i = (scatter.select),
    j = names(which.max(sapply(.SD,mean))),
    .SDcols = detectors.fluors
  ]
  ##density distribution of peak detector using selected scatter events (singlets);
  ##assumes two major peaks (negative/positive);
  ##bisect positive peak using return value from 'peak.values' function; 'top' n events
  fs$data[
    ,
    detector.select := scatter.select == TRUE & j>ifelse(
      grepl("Unstained",levels(sample.id)),
      min(peak.values(j)),
      max(peak.values(j))
    ),
    env = list(j = detector.peak)
  ]
  ##detector medians for 'top' n events -- positive
  detector.medians.positive<-fs$data[
    detector.select==TRUE,
    lapply(.SD,stats::median),
    .SDcols = detectors.fluors,
    by ='sample.id'
  ] |>
    data.table::melt(
      data=_,
      id.vars = 'sample.id',
      value.name = 'emission',
      variable.name = 'detector'
    )
  if(subtract.internal.negative){
    ##negative peak using return value from 'peak.values' function; 'bottom' n events
    fs$data[
      ,
      j = detector.select := scatter.select == TRUE & j<=min(peak.values(j)),
      env = list(j = detector.peak)
    ]
    ##detector medians for 'bottom' n events -- negative
    detector.medians.negative<-fs$data[
      detector.select==TRUE,
      lapply(.SD,stats::median),
      .SDcols = detectors.fluors,
      by ='sample.id'
    ] |>
      data.table::melt(
        data=_,
        id.vars = 'sample.id',
        value.name = 'emission',
        variable.name = 'detector'
      )
    detector.medians.positive[,emission := emission-detector.medians.negative[['emission']]]
  }
  ##
  detector.medians.positive[,emission.normalized := emission/max(emission)]
  ##
  detector.medians.positive[
    ,
    c('marker','fluor','tissue.type') := data.table::tstrsplit(sample.id," ",type.convert = factor)
  ]
  ##
  detector.medians.positive[]
}

spectral.signature.plot<-function(raw.fluor.control.path,subtract.internal.negative=FALSE){
  ##
  detector.medians <- spectral.signature.medians(raw.fluor.control.path,subtract.internal.negative)
  ##plot
  ggplot2::ggplot(
    data=detector.medians,
    mapping = ggplot2::aes(detector,emission.normalized)
  ) +
    ggplot2::geom_line(linewidth=0.5,group=1) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      xintercept = detector.medians[,which.max(emission.normalized)],
      linetype = 'dashed'
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        color = "black",
        angle = 90,
        vjust = 0.5,
        hjust=1)
    ) +
    ggplot2::labs(
      title = paste0(
        "Spectral Signature: Normalized Emission",
        paste0(rep(" ",10),collapse = ""),
        detector.medians[,levels(sample.id)]
      ),
      subtitle = paste(
        paste("Fluorophore:",detector.medians[,levels(fluor)]),
        paste("Marker:", detector.medians[,levels(marker)]),
        paste("Tissue Type:",detector.medians[,levels(tissue.type)]),
        # paste("Project:",unique(dts.medians$parameters.non$`$PROJ`)),
        sep = "\n"
      ),
      x = "Detector",
      y = "Normalized Emission",
      caption = paste("Peak Detector:", detector.medians[,detector[which.max(emission.normalized)]])
    )
  ##
}
# ggplot2::ggplot(
#   data = spectral.signature.medians(
#     "./inst/extdata/TNFa PE (Beads).fcs",
#     subtract.internal.negative = TRUE
#   ),
#   mapping = ggplot2::aes(detector,emission.normalized)) +
#   geom_line(group=1)
# ##
# pos<-spectral.signature.medians(
#   "./inst/extdata/TNFa PE (Beads).fcs",
#   subtract.internal.negative = FALSE
# )
# neg<-spectral.signature.medians(
#   "./inst/extdata/Unstained (Beads).fcs",
#   subtract.internal.negative = FALSE
# )
# pos[,emission := emission-neg[['emission']]]
# pos[,emission.normalized := emission/max(emission)]
#
# ggplot2::ggplot(
#   data = pos,
#   mapping = ggplot2::aes(detector,emission.normalized)) +
#   geom_line(group=1)
