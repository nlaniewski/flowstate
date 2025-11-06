#' @title Generate normalized medians/spectra from reference group 'spectral events'
#' @description
#' The 'spectral events' returned from [reference.group.spectral.events] are used to generate normalized `[0,1]` medians/spectra for use in spectral unmixing.
#'
#'
#' @param flowstate.object.reference A `flowstate` object as returned from [reference.group.spectral.events].
#' @param name.fix Named character vector -- default `NULL`; if defined, vector names should match to `sample.id`(s) and their respective values name replacements.
#' * e.g., `c('Cd8 BV786 (Cells)' = 'CD8 BV786')`.
#' @param syntactically.valid Logical -- default `FALSE`; if `TRUE`, spaces, dashes, and dots are removed from strings.
#'
#' @returns A [data.table][data.table::data.table] containing normalized reference control medians.
#' @export
#'
reference.group.medians <- function(
    flowstate.object.reference,
    name.fix=NULL,
    syntactically.valid=FALSE
)
{
  ##name fixes
  if(!is.null(name.fix)){
    for(i in names(name.fix)){
      data.table::set(
        x = flowstate.object.reference$data,
        i = flowstate.object.reference$data[,.I[sample.id == i]],
        j = 'sample.id',
        value = name.fix[[i]]
      )
    }
  }
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
  ##subtract population-specific 'universal negative'
  ##retain 'universal negative' as autofluorescence spectra
  for(.population in ref.medians[,levels(population)]){
    ref.medians[
      i = population == .population,
      j = (cols.detector) := lapply(.SD,function(j){
        af <- j[reference.type == 'universal negative']
        j <- j - af
        j[reference.type == 'universal negative'] <- af
        j
      }),
      .SDcols = cols.detector,
    ]
  }
  ##drop subtrahends
  if(any(ref.medians[,.SD |> rowSums() == 0,.SDcols = cols.detector])){
    ref.medians <- ref.medians[ref.medians[,.SD |> rowSums() != 0,.SDcols = cols.detector]]
  }
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
  ##modify 'ref.medians' for eventual use in an OLS fit as an overdetermined 'unmixing matrix'
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
  ref.medians[N == 'AF', N := paste0("AF.",population)]
  ref.medians[,ord := seq(.N)]
  ref.medians[grepl('AF',N), ord := ref.medians[,.N] + seq(.N)]
  data.table::setorder(ref.medians,ord)[,ord := NULL]#AF in last position
  if(syntactically.valid){
    ref.medians[!grepl('AF',N),N := tolower(gsub(" |-|\\.","",N))]
    ref.medians[,S := gsub("-","",S)]
  }
  ref.medians[is.na(S), alias := N]
  ref.medians[!is.na(S),alias := paste(S,N)]
  ##return the normalized reference control medians
  ref.medians[]
}
#' @title Plot Spectral Traces
#' @description
#' Plots a 'spectral trace' -- normalized `[0,1]` emission for all detectors.
#'
#' @param ref.medians The return of [reference.group.medians].
#'
#' @returns A list of [ggplot][ggplot2::ggplot] objects.
#' @export
#'
plot_spectral.trace <- function(ref.medians){
  cols.detector <- ref.medians[,names(.SD),.SDcols = is.numeric]
  cols.by <- ref.medians[,names(.SD),.SDcols = !is.numeric]
  ##
  ref.medians.split <- split(
    ref.medians,
    by=c('sample.id','population'),
    drop = TRUE,
    sorted = FALSE
  )
  ##
  ref.medians.melted <- lapply(ref.medians.split,data.table::melt,measure.vars = cols.detector)
  ref.medians.melted <- lapply(ref.medians.melted,droplevels)
  ##
  p.list <-lapply(ref.medians.melted,function(dt){
    ggplot2::ggplot(
      data = dt,
      mapping = ggplot2::aes(variable,value)
    ) +
      ggplot2::geom_line(linewidth=0.5,group=1) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        xintercept = dt[,levels(detector.peak)],
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
          unique(dt[['sample.id']])
        ),
        subtitle = paste(
          paste("Fluorophore:",unique(dt[['N']])),
          paste("Marker:", unique(dt[['S']])),
          sprintf("Group Type: %s     Population: %s",unique(dt[['group.type']]),unique(dt[['population']])),
          sep = "\n"
        ),
        x = "Detector",
        y = "Emission (Normalized)",
        caption = paste("Peak Detector:", unique(dt[['detector.peak']]))
      )
  })
}
parameters.raw.to.unmixed <- function(flowstate.object.raw,ref.medians){
  ##copy raw/overdetermined parameters
  parms <- data.table::copy(flowstate.object.raw$parameters)
  ##parameters -- non-overdetermined; essentially scatter and time
  parms.nondetector <- parms[TYPE != "Raw_Fluorescence"]
  ##merge ref.medians with parms to inherit keyword-value pairs;
  ##represents unmixed parameters; maintains the order as set in ref.medians
  ##detectors/parameters not featured (non-end members) are dropped during the merge
  parms.unmixed <- merge(
    ref.medians[,.(N.alias = detector.peak,N.ref = N,S)],
    parms,
    by = 'N.alias',
    sort = FALSE
  )
  ##modify parms.unmixed
  data.table::setnames(parms.unmixed,'N.alias','DETECTOR')
  data.table::setnames(parms.unmixed,'N.ref','N.alias')
  parms.unmixed[,TYPE := 'Unmixed_Fluorescence']
  parms.unmixed[,N := paste0(N.alias,'-A')]
  # parms.unmixed[,ord := seq(.N)]
  # parms.unmixed[N == 'AF-A',ord := ord + max(parms.unmixed[['ord']])]
  # data.table::setorder(parms.unmixed,ord)[,ord := NULL]
  parms.unmixed <- rbind(parms.nondetector,parms.unmixed,fill=TRUE)
  parms.unmixed[is.na(S),alias := N.alias]
  parms.unmixed[!is.na(S),alias := paste(S,N.alias)]
  if(any(!ref.medians[['alias']] %in% parms.unmixed[['alias']])){
    stop("Mismtach between derived parameter names and ref.medians names")
  }
  ##
  parms.unmixed[,j := paste0('$P',seq(.N)),env = list(j = 'par')]
  ##
  invisible(parms.unmixed)
}
flowstate.unmixed.object <- function(flowstate.object.raw,ref.medians){
  ##initialize a new flowstate.object;
  ##inherits keyword-value pairs from flowstate.object.raw;
  ##inherits non-detector data from flowstate.object.raw;
  ##derives 'unmixed' parameters from ref.medians

  ##initialize an empty flowstate
  fs.unmixed <- flowstate()
  ##copy/inherit keywords
  fs.unmixed$keywords <- data.table::copy(flowstate.object.raw$keywords)
  ##GUID is some sort of SpectroFlo hash -- specific to each source file;
  ##no longer applicable
  if('GUID' %in% names(fs.unmixed$keywords)){
    fs.unmixed$keywords[,j := NULL, env = list(j = "GUID")]
  }
  ##copy/inherit parameters; raw/overdetermined --> unmixed;
  ##new $PnN and $PnS keyword-value pairs as derived from ref.medians;
  ##if there are failures here, ref.medians will need to be 'fixed'
  fs.unmixed$parameters <- parameters.raw.to.unmixed(flowstate.object.raw,ref.medians)
  ##spillover matrix; derived from [['parameters']]
  spill.mat <- matrix(0,nrow = ref.medians[,.N],ncol = ref.medians[,.N])
  colnames(spill.mat) <- rownames(spill.mat) <- ref.medians[['alias']]
  diag(spill.mat) <- 1
  spill.mat[2,1] <- 1E-6#bug-fix for FlowJo
  spill <- data.table::as.data.table(spill.mat)
  rownames(spill) <- names(spill)
  data.table::setattr(spill,'applied',FALSE)
  fs.unmixed$spill <- spill
  ##names of detectors
  j.match <- j.match.parameters.to.data(flowstate.object.raw)
  cols.detector <- flowstate.object.raw$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  i <- names(flowstate.object.raw$data) %in% cols.detector#c(cols.detector,'sample.id')
  cols.retain <- names(flowstate.object.raw$data)[!i]
  ##assign to fs.unmixed[['data']] non-detector parameters; scatter and time
  fs.unmixed$data <- flowstate.object.raw$data[,.SD,.SDcols = cols.retain]
  ##
  invisible(fs.unmixed)
}
#' @title Unmix a `flowstate`
#' @description
#' A raw/overdetermined `flowstate` is unmixed using a corresponding `ref.medians`. Data is unmixed using [Ordinary Least Squares][stats::lsfit].
#'
#'
#' @param flowstate.object.raw A `flowstate` -- the return of [read.flowstate]; the `flowstate` must be raw/overdetermined.
#' @param ref.medians A [data.table][data.table::data.table] containing normalized `[0,1]` reference control medians -- the return of [reference.group.medians].
#' @param hash.historic Character string -- default `NULL`; for internal/reproducibility purposes, a hash string of the 'unmixing matrix' can be defined.
#'
#' @returns A `flowstate` containing unmixed data.
#' @export
#'
flowstate.unmix <- function(
    flowstate.object.raw,
    ref.medians,
    hash.historic = NULL
)
{
  ##initialize a flowstate.object -- unmixed
  fs.unmixed <- flowstate.unmixed.object(flowstate.object.raw,ref.medians)
  ##ref.medians as a matrix; overdetermined
  unmixing.mat <- as.matrix(ref.medians[,.SD,.SDcols = is.numeric])
  rownames(unmixing.mat) <- ref.medians[['alias']]
  if(!is.null(hash.historic)){
    hash.unmixing.mat <- apply(unmixing.mat,2,digest::digest) |> digest::digest()
    if(hash.unmixing.mat != hash.historic){
      stop("\nThe (digested) hash for the unmixing matrix does not match the historic value;\nresults will not be fully reproducible.")
    }
  }
  ##test
  if(!all(colnames(unmixing.mat) %in% names(flowstate.object.raw$data))){
    stop("Mismatch in detector names between ref.medians and flowstate.object.raw.")
  }else{
    cols.detector <- colnames(unmixing.mat)
  }
  ##ordinary least squares fit
  ##NB: a computational hit here...transposing [['data']]...then transposing the coefficients
  res <- stats::lsfit(
    x = t(unmixing.mat),
    y = t(flowstate.object.raw$data[,.SD,.SDcols = cols.detector]),
    intercept = FALSE
  )[['coefficients']] |> t()
  ##unmixed data; linear values
  ##add to fs.unmixed$data
  fs.unmixed$data[,(colnames(res)) := as.list(as.data.frame(res))]
  ##invisibly return the updated fs.unmixed
  invisible(fs.unmixed)
}
