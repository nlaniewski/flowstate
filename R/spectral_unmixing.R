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
reference.group.keywords <- function(
    flowstate.object,
    internal.negative = NULL
)
{
  ##column used to match 'sample.id' in [['data']] to keyword/value in [['keywords']]
  col.identifier <- j.match.keyword.to.data.sample.id(flowstate.object)
  ##add keywords for eventual use in data.table by/group operations
  flowstate.object$keywords[
    ,
    group.type := stringr::str_extract(j,"Cells|Beads") |> tolower(),
    env = list(j = col.identifier)
  ]
  flowstate.object$keywords[
    i = grep("Unstained",i),
    j = c('reference.type','subtract.type') := list("universal negative","subtrahend"),
    env = list(i = col.identifier)
  ]
  flowstate.object$keywords[is.na(reference.type),reference.type := "spectra"]
  flowstate.object$keywords[reference.type == 'spectra',subtract.type := "external"]
  if(!is.null(internal.negative)){
    flowstate.object$keywords[
      i = i %in% internal.negative,
      j = subtract.type := "internal",
      env = list(i = col.identifier)
    ]
  }
  flowstate.object$keywords[,autofluorescence := FALSE]
  flowstate.object$keywords[
    i = grep("Unstained.*Cells",i),
    j = autofluorescence := TRUE,
    env = list(i = col.identifier)
  ]
  flowstate.object$keywords[
    ,
    j = scatter.population := ifelse(group.type == 'beads','beads','lymphocytes')
  ]
  ##add the new keywords to [['data']]; conversion to factor using 'type.covert'
  cols.factor <- c('group.type','reference.type','subtract.type',
                   'autofluorescence','scatter.population')
  flowstate::add.keywords.to.data(flowstate.object,cols.factor,type.convert = T)
}
##helper function for events selection based on density distribution/peak height
##selects for expected bead, cellular/PBMC scatter populations using scatter dimensions as input
##returns a vector (logical) that is assigned using ':='
##can also return (base) plots of density distributions/selection bounds
select.scatter.population <- function(
    ...,
    population = c("bead","lymphocyte","monocyte"),
    plot = FALSE,
    sample.id = NULL
)
{
  dots.names <- sapply(substitute(list(...)), deparse)[-1]
  if(plot){
    graphics::par(mfcol = c(length(dots.names), 1), oma = c(0, 0, 1, 0))
    on.exit(graphics::par(mfcol = c(1, 1), oma = c(0, 0, 0, 0)))
  }
  for(i in seq_along(dots.names)){
    cuts <- peak.height.cut(
      x = if(!exists('vec')){
        list(...)[[i]]
      }else{
        list(...)[[i]][vec]
      },
      height.cut = 0.5,
      which.peak = ifelse(population == "monocyte" & dots.names[i] == "SSC_A",2,1),
      plot = plot,
      main = paste(dots.names[i],"Select")
    )
    if(!exists('vec')){
      vec <- data.table::`%between%`(list(...)[[i]],cuts)
    }else{
      vec[which(vec)] <- data.table::`%between%`(list(...)[[i]][vec],cuts)
    }
  }
  ##
  if(!is.null(sample.id) & plot){
    graphics::mtext(sample.id, side = 3, line = -1, outer = TRUE)
  }
  return(vec)
}
##helper function for peak detector events selection based on variance/mean
##selects for an optimal distribution of expressing events in the peak detector channel
##returns vectors (character and logical) that are assigned using ':='
select.detector <- function(
    .SD,
    method = c('beads','cells'),#as.character(.BY$group.type)
    subtract.type,#.BY$subtract.type
    n=100,
    plot = FALSE,
    sample.id#.BY$sample.id
)
{
  if(grepl("Unstained",sample.id)){
    detector.peak <-names(which.max(sapply(.SD,mean)))
    cuts <- peak.height.cut(
      x = .SD[[detector.peak]],
      which.peak = 1,
      height.cut = 0.10,
      plot = F
    )
    vec <- data.table::`%between%`(.SD[[detector.peak]],cuts)
  }else{
    detector.peak.var <- names(which.max(sapply(.SD,stats::var)))
    ord <- order(.SD[[detector.peak.var]],decreasing = T)[1:n]
    detector.peak.mean.ord <-names(which.max(sapply(.SD[ord],mean)))
    if(detector.peak.var == detector.peak.mean.ord){
      detector.peak <- detector.peak.var
    }else{
      detector.peak <- detector.peak.mean.ord
      ord <- order(.SD[[detector.peak]],decreasing = T)[1:n]
    }
    if(subtract.type == 'internal'){
      # ord <- c(ord,order(.SD[[detector.peak]])[1:n])
      cuts <- peak.height.cut(
        x = .SD[[detector.peak]],
        which.peak = 1,
        height.cut = 0.80,
        plot = F
      )
      vec.neg <- data.table::`%between%`(.SD[[detector.peak]],cuts)
    }
    method <- match.arg(method)
    if(method == 'beads'){
      cuts <- peak.height.cut(
        x = .SD[[detector.peak]],
        which.peak = 2,
        height.cut = 0.99,
        plot = F
      )
      vec <- data.table::`%between%`(.SD[[detector.peak]],cuts)
    }else if(method == 'cells'){
      vec <- vector(length = length(.SD[[detector.peak]]))
      vec[ord] <- TRUE
      if(subtract.type == 'internal'){
        vec[which(vec.neg)] <- TRUE
      }
    }
  }
  ##
  if(plot){
    p <- ggplot2::ggplot() +
      ggplot2::geom_hex(
        data = NULL,
        mapping = ggplot2::aes(
          x = seq_along(.SD[[detector.peak]])[which(!vec)],
          y = .SD[[detector.peak]][which(!vec)]),
        bins = 100,
        alpha = 0.6,
        show.legend = F
      ) +
      ggplot2::scale_fill_gradient(low = 'gray',high = 'black') +
      ggplot2::geom_hex(
        data = NULL,
        mapping = ggplot2::aes(
          x = seq_along(.SD[[detector.peak]])[which(vec)],
          y = .SD[[detector.peak]][which(vec)]),
        fill = 'red',
        bins = 100,
        show.legend = F
      ) +
      ggplot2::labs(
        title = sprintf("%s\nDetector (Peak): %s",sample.id,detector.peak),
        caption = sprintf("%d 'spectral' events selected in %s for calculating per-detector (n = %d) medians.",
                          length(which(vec)), detector.peak, length(.SD)),
        x = 'Index',
        y = detector.peak
      )
    print(p)
  }
  ##
  list(detector.peak,vec)
}
#' @title A flowstate containing spectrally-associated bead/cellular events.
#' @description
#' This function is entirely dependent on the following:
#' \itemize{
#'   \item Spectroflo software raw reference controls  -- tested using SpectroFlo 3.3.0.
#'   \item Spectroflo software naming convention:
#'   \itemize{
#'     \item Directory -- './Raw/Reference Group/...'
#'     \item File names -- 'MARKER FLUOR (Beads|Cells).fcs' (literal space characters separating)
#'     \item Quality reference controls -- well resolved, dominant scatter populations
#'   }
#' }
#'
#' @param raw.reference.group.directory Character vector; a file path to a SpectroFlo Raw Reference Group directory.
#' @param internal.negative Character vector; spectral events associated with an internal negative population will be generated for any named reference controls.
#' @param n Numeric; default 500. Orders the top `n` expressing events and aids in auto-detection of peak detector; used by an internal helper function.
#' @param plot.select Logical; default `FALSE`. Plots diagnostic/QC output for evaluating function performance.
#' @param plot.select.dir Character vector; a file path for storing plot output.
#'
#' @returns A flowstate containing spectrally-associated bead/cellular events. `[['data']]` will contain appended columns/metadata.
#' @export
#'
#' @examples
#' raw.ref.dir <- system.file("extdata", package = "flowstate") |>
#' list.dirs() |>
#' grep('Reference Group', x = _, value = TRUE)
#'
#' ref <- reference.group.spectral.events(
#' raw.reference.group.directory = raw.ref.dir,
#' internal.negative = NULL,
#' n = 500,
#' plot.select = TRUE,
#' plot.select.dir = tempdir()
#' )
#'
#' #plot output
#' list.files(tempdir(),pattern = "select.*.pdf")
#'
#' #spectral events
#' ref$data
#'
reference.group.spectral.events <- function(
    raw.reference.group.directory,
    internal.negative = NULL,
    n = 500,
    plot.select = FALSE,
    plot.select.dir = tempdir()
)
{
  ##paths to raw .fcs files; reference group controls
  ref.paths <- list.files(raw.reference.group.directory,full.names = T)
  ##test for 'Raw' fluorescence
  sapply(ref.paths,raw.fluorescence.check) |> invisible()
  ##read multiple files; concatenate
  ref <- flowstate::read.flowstate(
    ref.paths,
    colnames.type = 'N',
    concatenate = TRUE
  ) |> suppressMessages()
  ##update with reference control-specific keywords
  reference.group.keywords(
    flowstate.object = ref,
    internal.negative = internal.negative
  )
  ##alias column
  j.match <- j.match.parameters.to.data(ref)
  ##detector columns
  cols.detector <- ref$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##select against saturating events -- scatter and fluors; mark saturating events as FALSE
  ref$data[,select.nonsaturating := TRUE]
  for(j in c(grep("[FS]SC",names(ref$data),value = T),cols.detector)){
    data.table::set(
      x = ref$data,
      i = which(ref$data[[j]]>=4194304),
      j = 'select.nonsaturating',
      value = FALSE
    )
  }
  ##select for expected bead, cellular/PBMC scatter populations
  ##capture plot output
  if(plot.select){
    file.out.name <- sprintf("%s_select_population.pdf", unique(ref$keywords$`$PROJ`))
    grDevices::pdf(file.path(plot.select.dir,file.out.name),width = 6,height = 8)
  }
  ref$data[
    (select.nonsaturating),
    j = select.population := select.scatter.population(
      ... = FSC_A,SSC_A,SSCB_A,
      population = .BY$scatter.population,
      plot = ifelse(plot.select,T,F),
      sample.id = .BY$sample.id
    ),
    by = .(sample.id,scatter.population)#c(ref$data[,names(.SD),.SDcols = !is.numeric])
  ]
  if(plot.select) grDevices::dev.off()
  # list.files(tempdir(),full.names = T, pattern = file.out.name) |> normalizePath() |> shell.exec()
  ##transform
  flowstate::flowstate.transform(
    ref,
    .j = cols.detector,
    transform.type = "asinh",
    cofactor = 5000
  ) |> suppressMessages()
  ##peak detector events selection based on variance/mean
  ##selects for an optimal distribution of expressing events in the peak detector channel
  ##returns vectors (character and logical) that are assigned using ':='
  ##capture plot output
  if(plot.select){
    file.out.name <- sprintf("%s_select_detector.pdf",unique(ref$keywords$`$PROJ`))
    grDevices::pdf(file.path(plot.select.dir,file.out.name),width = 7,height = 7)
  }
  ref$data[
    (select.population),
    j = c('detector.peak', 'select.detector') :=
      select.detector(
        .SD[(select.population)],
        method = as.character(.BY$group.type),
        subtract.type = .BY$subtract.type,
        n = n,
        plot = ifelse(plot.select,T,F),
        sample.id = .BY$sample.id
      ),
    .SDcols = cols.detector,
    by = .(sample.id,group.type,subtract.type)
  ]
  if(plot.select) grDevices::dev.off()
  # list.files(tempdir(),full.names = T, pattern = file.out.name) |> normalizePath() |> shell.exec()
  ##subset the flowstate.object to include only 'select.detector' events
  ref <- subset(ref,select.detector == TRUE)
  ##update keywords with peak detectors; factor
  ref$keywords[
    ,
    detector.peak := ref$data[
      ,
      j = .(detector.peak=unique(detector.peak)),
      by = .(sample.id)][['detector.peak']]
  ]
  ref$keywords[,detector.peak := factor(detector.peak,levels = cols.detector)]
  flowstate::add.keywords.to.data(ref,'detector.peak')
  ##order [['data']] by 'detector.peak'
  data.table::setorder(ref$data,detector.peak)
  ##transform inverse; [['data']] to linear/raw
  flowstate.transform.inverse(ref)
  ##NULL-out some now-redundanct columns
  ref$data[,grep('select',names(.SD),value = T) := NULL]
  ##for 'internal' negatives
  for(.i in internal.negative){
    detector <- ref$data[
      i = sample.id == .i,
      j = as.character(unique(detector.peak))
    ]
    x <- ref$data[sample.id == .i][[detector]] |> sort()
    i.neg <- x |> diff() |> which.max()
    val <- x[i.neg]
    data.table::set(
      x = ref$data,
      i = ref$data[,.I[sample.id == .i & j <= val],env = list(j = detector)],
      j = 'subtract.type',
      value = 'subtrahend'
    )
  }
  ##return the flowstate.object
  ref[]
}
plot_spectral.ridges <- function(flowstate.object.reference,plot.dir=tempdir()){
  j.match <- j.match.parameters.to.data(flowstate.object.reference)
  cols.detector <- flowstate.object.reference$parameters[TYPE == 'Raw_Fluorescence'][[j.match]]
  cols.by <- flowstate.object.reference$data[,names(.SD),.SDcols = is.factor]
  # lims <- flowstate.object.reference$data[
  #   ,
  #   sapply(.SD,function(j){range(asinh(j/5000))}) |> range(),
  #   .SDcols = cols.detector
  # ]
  ##
  file.out.name <- sprintf("%s_spectral_ridges.pdf",
                           unique(flowstate.object.reference$keywords$`$PROJ`))
  grDevices::pdf(file.path(plot.dir,file.out.name),width = 6,height = 10)
  on.exit(grDevices::dev.off())
  message("Printing 'Spectral Ridges' to device (.pdf)...")
  ##
  p <- flowstate.object.reference$data[
    ,
    j = {
      message(.BY$sample.id)
      x <- NULL
      p <- ggplot2::ggplot(
        data = data.table::melt(.SD,measure.vars = cols.detector),
        mapping = ggplot2::aes(x = asinh(value/5000), y = variable, fill = ggplot2::after_stat(x))
      ) +
        ggridges::geom_density_ridges_gradient(show.legend = F) +
        # ggplot2::xlim(lims) +
        ggplot2::scale_y_discrete(sec.axis = ggplot2::dup_axis()) +
        viridis::scale_fill_viridis(option = 'plasma') +
        ggplot2::labs(
          x = "Expression (Density)",
          y = "Detector"
        ) +
        ggplot2::labs(
          title = sprintf("%s\nDetector (Peak): %s",.BY$sample.id,.BY$detector.peak),
          subtitle = "Spectral Ridges",
          caption = sprintf("%s events (n = %d); transformed using asinh(x/5000).",
                            ifelse(.BY$subtract.type == 'subtrahend','Negative','Positive'),
                            .N)
        ) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank()
        )
      print(p) |> suppressMessages()
      data.table::data.table()
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
}
##
reference.group.medians <- function(flowstate.object.reference,name.fix=NULL){
  j.match <- j.match.parameters.to.data(flowstate.object.reference)
  cols.detector <- flowstate.object.reference$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  cols.by <- flowstate.object.reference$data[,names(.SD),.SDcols = !is.numeric]
  ##generate medians
  ref.medians <- flowstate.object.reference$data[
    ,
    j = lapply(.SD,stats::median),
    .SDcols = cols.detector,
    by = cols.by
  ]
  ##ids for controls using an internal subtrahend
  ids.internal <- ref.medians[
    i = subtract.type == 'internal',
    as.character(sample.id)
  ]
  ##subtract subtrahend from internals
  for(.sample.id in ids.internal){
    ##subtraction vector
    subtrahend <- ref.medians[
      i = (sample.id == .sample.id
           & subtract.type == 'subtrahend'),
      j = unlist(.SD),
      .SDcols = cols.detector
    ]
    ##sweep out the subtrahend
    ref.medians[
      i = sample.id == .sample.id,
      j = (cols.detector) := sweep(.SD,2,subtrahend, `-`),
      .SDcols = cols.detector
    ]
  }
  ##subtract subtrahend from externals
  for(.group.type in ref.medians[,levels(group.type)]){
    ##subtraction vector
    subtrahend <- ref.medians[
      i = (group.type == .group.type
           & reference.type == 'universal negative'
           & subtract.type == 'subtrahend'),
      j = unlist(.SD),
      .SDcols = cols.detector
    ]
    ##sweep out the subtrahend
    ref.medians[
      i = (group.type == .group.type
           & subtract.type %in% c('external','subtrahend')
           & !(autofluorescence )
           & !sample.id %in% ids.internal),
      j = (cols.detector) := sweep(.SD,2,subtrahend, `-`),
      .SDcols = cols.detector
    ]
  }
  ##drop subtrahends
  ref.medians <- ref.medians[ref.medians[,.SD |> rowSums() != 0,.SDcols = cols.detector]]
  ##normalize
  ref.medians[
    ,
    j = (cols.detector) := {
      j <- .SD/max(.SD)
      j[j<0]<-0
      j
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ##modify 'ref.medians' for eventual use in an OLS fit as an overdetermined 'mix matrix'
  ref.medians[
    ,
    j = c('N','S') := {
      i <- sub(" \\(\\w+\\)","",as.character(.BY$sample.id))#subs out '(Beads|Cells)'
      splits <- strsplit(i," ")
      S <- sapply(splits,'[[',1)#marker
      N <- trimws(sub(S,"",i))#fluor
      N[grep("Unstained",S)] <- 'AF'
      S[grep("Unstained",S)] <- NA
      list(N,S)
    },
    by=sample.id
  ]
  ref.medians[,ord := seq(.N)]
  ref.medians[N=='AF',ord := max(ref.medians[['ord']]) + 1]
  data.table::setorder(ref.medians,ord)[,ord := NULL]#AF in last position
  if(!is.null(name.fix)){
    for(i in seq_along(name.fix)){
      ref.medians[N == names(name.fix[i]), N := name.fix[[i]]]
    }
  }
  ref.medians[is.na(S), alias := N]
  ref.medians[!is.na(S),alias := paste(S,N)]
  ##return the normalized reference control medians
  ref.medians[]
}
##
