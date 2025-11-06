raw.fluorescence.check<-function(fcs.file.paths){
  sapply(fcs.file.paths,function(i){
    parms <- readFCStext(i) |>
      parameters.to.data.table()
    if(!'TYPE' %in% names(parms)){
      stop("Are these .fcs files processed by flowstate? 'TYPE' is missing from [['parameters']].")
    }else{
      if(!"Raw_Fluorescence" %in% parms[['TYPE']]){
        stop("'Raw_Fluorescence' value not found in 'TYPE' keyword; is this a raw reference control?")
      }
    }
  }) |> invisible()
}
add.reference.keywords.to.data <- function(flowstate.object){
  ##column used to match 'sample.id' in [['data']] to keyword/value in [['keywords']]
  col.identifier <- j.match.keyword.to.data.sample.id(flowstate.object)
  ##add keywords for eventual use in data.table by/group operations
  flowstate.object$keywords[
    ,
    j = c('group.type','reference.type') := {
      list(
        stringr::str_extract(j,"[Cc]ells|[Bb]eads") |> tolower(),
        ifelse(grepl("Unstained",j),"universal negative","spectra")
      )
    },
    env = list(j = col.identifier)
  ]
  add.keywords.to.data(flowstate.object,c('group.type','reference.type'))
}
select.beads <- function(flowstate.object,height.cut.scatter = 0.5,plot=FALSE){
  if(!'population' %in% names(flowstate.object$data)){
    flowstate.object$data[,population := character()]
  }
  ##use all available scatter parameters
  cols.scatter <- grep("[FS]SC.*_[AH]",names(flowstate.object$data),value = T)
  ##trim dominant scatter peak in each scatter dim to isolate singlets;
  ##logical vector for selection
  select.population <- flowstate.object$data[
    i = group.type == 'beads',
    j = sapply(.SD,function(j){
      cuts <- peak.height.cut(
        j,
        which.peak = 1,
        height.cut = height.cut.scatter,
        plot = FALSE
      )
      data.table::`%between%`(j,cuts)
    }) |> rowSums() == length(cols.scatter),
    .SDcols = cols.scatter
  ]
  ##update [['data']] by reference
  data.table::set(
    x = flowstate.object$data,
    i = flowstate.object$data[,.I[group.type == 'beads'][select.population]],
    j = 'population',
    value = "beads"
  )
  if(plot){
    graphics::par(
      mfcol = c(2, 1)
    )
    on.exit(graphics::par(mfcol = c(1,1)))
    flowstate.object$data[
      i = group.type == 'beads',
      j = {
        for(i in c('FSC_A','SSC_A')){
          d1 <- stats::density(.SD[[i]])
          d2 <- stats::density(.SD[[i]][which(population %in% 'beads')])
          ylim <- range(d1$y,d2$y)
          plot(d1,ylim=ylim,main = sprintf("%s\nPopulation: %s",i,'beads'))
          graphics::lines(d2,col="red")
        }
      }
    ]
  }
  ##
  invisible(flowstate.object)
}
cellular.anchors <- function(
    flowstate.object,
    population.marker = NULL,
    top.N.percent = 1,
    height.cut.scatter = 0.5,
    plot = TRUE,
    ...
)
{
  ##test
  if('population' %in% names(flowstate.object$data)){
    message(
      sprintf(
        "The target column 'population' already exists;\nif it is to be overwritten/replaced: %s$data[,population := NULL]",
        deparse(substitute(flowstate.object)))
    )
  }
  ##test
  res <- flowstate.object$data[
    i = group.type == 'cells',
    j = levels(sample.id)[levels(sample.id) %in% population.marker]
  ]
  if(!all(res %in% population.marker)){
    stop("Defined 'population.marker(s)' not found among celluar samples.")
  }
  ##use all available scatter parameters -- Area and Height;
  ##FSC_A --> FSC_H --> SSC_A --> SSC_H --> SSCB_A* --> SSCB_H* (sequential cutting order;*if available)
  cols.scatter <- sort(grep("[FS]SC.*_[AH]",names(flowstate.object$data),value = T))
  ##alias column
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##detector columns
  cols.detector <- flowstate.object$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##find peak detector among each population/marker;
  ##index top n peak detector-expressing events;
  ##"anchor" scatter distributions to the indexed events;
  ##generate sequential cuts for each dominant scatter peak to isolate populations
  ##apply sequential cuts to all data; update [['data']] with a new column "population"
  lapply(names(population.marker),function(.population){
    if(plot){
      s <- seq(cols.scatter)
      mat <- matrix(
        c(
          s[s %%2 != 0],
          s[s %%2 == 0]
        ),
        nrow = 2,
        byrow = TRUE
      )
      graphics::layout(mat)
      graphics::par(oma = c(0, 0, 5, 0))
    }
    marker <- population.marker[[.population]]
    cuts.anchored <- flowstate.object$data[
      i = sample.id == marker,
      j = {
        detector.peak <- sapply(.SD,function(j){
          v <- stats::var(j)
          m <- mean(sort(j,decreasing = T)[1:200])
          c(v=v,m=m)
        }) |> apply(X = _,1,function(x) names(which.max(x)))
        if(length(unique(detector.peak)) != 1){
          stop("Edge case: different 'detector.peak' between variance and mean (sorted vector).\n",
               sprintf("%s: %s (variance) ; %s (mean)",marker,detector.peak[1],detector.peak[2])
          )
        }else{
          detector.peak <- unique(detector.peak)
        }
        top.N.percent <- ceiling(.N*(top.N.percent/100))
        cut.fluor <- sort(.SD[[detector.peak]],decreasing = T)[1:top.N.percent] |> min()
        ##initial logical
        i <- .SD[[detector.peak]] > cut.fluor
        ##initialize list
        cuts.scatter <- stats::setNames(
          vector(mode = 'list',length = length(cols.scatter)),
          nm = sort(cols.scatter)
        )
        ##store sequential cuts
        for(scatter in names(cuts.scatter)){
          cuts <- peak.height.cut(
            x = get(scatter)[i],
            which.peak = 1,
            height.cut = height.cut.scatter,
            plot = plot,
            main = sprintf("Scatter: %s",scatter),
            # main = sprintf("Population: %s\nSample: %s\nScatter: %s",.population,.BY$sample.id,scatter),
            sub = sprintf("Peak height cut: %s",height.cut.scatter)
          )
          ##store the cuts; list
          cuts.scatter[[scatter]] <- cuts
          ##update logical
          i[i] <- data.table::`%between%`(get(scatter)[i],cuts)
        }
        if(plot){
          graphics::mtext(sprintf("Population: %s\nAnchor (Sample): %s",.population,.BY$sample.id),
                side = 3,
                line = 1,
                outer = TRUE,
                cex = 1.5
          )
        }
        ##return the stored cuts
        cuts.scatter
      },
      .SDcols = cols.detector,
      by = c(flowstate.object$data[,names(.SD),.SDcols = is.factor])
    ][,population := .population]
    ##
    select.population <- flowstate.object$data[
      i = group.type == 'cells',
      j = {
        for(scatter in cuts.anchored[,names(.SD),.SDcols = is.numeric]){
          if(!exists("vec")){
            vec <- data.table::`%between%`(get(scatter),cuts.anchored[[scatter]])
          }else{
            vec[vec] <- data.table::`%between%`(get(scatter)[vec],cuts.anchored[[scatter]])
          }
        }
        vec
      }
    ]
    data.table::set(
      x = flowstate.object$data,
      i = flowstate.object$data[,.I[group.type == 'cells'][select.population]],
      j = 'population',
      value = .population
    )
  })
  ##factor
  flowstate.object$data[,population := factor(population)]
  ##return
  if(plot){on.exit(graphics::par(mfcol = c(1,1),mfrow = c(1,1),oma = c(0,0,0,0)))}
  invisible(flowstate.object)
}
plot_populations <- function(flowstate.object,pattern.scatter = c('A','AH'),bins = 300, sample.n = 5E4){
  cols.scatter <- sort(grep(sprintf("[FS]SC.*_[%s]$",match.arg(pattern.scatter)),names(flowstate.object$data),value = T))
  silhouette <- flowstate.object$data[
    ,
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n)}]),
    .SDcols = cols.scatter
  ]
  foreground <- flowstate.object$data[
    i = !is.na(population),
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n,replace = T)}]),
    .SDcols = cols.scatter,
    by=population
  ]
  plot.body <- ggplot2::ggplot(
    data = NULL,
    mapping = ggplot2::aes(value.x, value.y)
  )
  p <- plot.body + ggplot2::geom_hex(
    data = silhouette,
    bins = bins,
    fill = "gray"
  )
  p <- p + ggplot2::geom_hex(
    data = foreground,
    bins = bins,
    show.legend = FALSE) +
    viridis::scale_fill_viridis(
      option = "plasma",
      limits = c(0, 100),
      oob = scales::squish
    )
  p <- p + ggplot2::facet_grid(
    variable.y ~ variable.x,
    scales = "free")
  p <- p +
    ggplot2::scale_x_continuous(labels = scales::scientific) +
    ggplot2::scale_y_continuous(labels = scales::scientific)
  p <- p + ggplot2::labs(
    title = sprintf(
      "Populations (Anchored): %s",
      flowstate.object$data[,levels(population) |> paste0(collapse = "; ")]
    )
  )
  p
}
detector.peak <- function(flowstate.object,top.n=150){
  ##alias column
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##detector columns
  cols.detector <- flowstate.object$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##variance -- per-detector/per-sample/per-population;
  ##returns peak detector
  detectors.variance <- flowstate.object$data[
    ,
    j = .(detector.peak = sapply(.SD,stats::var) |> which.max() |> names()),
    .SDcols = cols.detector,
    keyby = .(sample.id,population)
  ]
  ##test returned peak detectors using sorted vectors;
  ##mean of top n expressing events
  ##helps resolve to maximally expressing, population-specific events
  detectors.mean <- flowstate.object$data[
    ,
    j = {
      detector <- detectors.variance[
        sample.id == .BY$sample.id][['detector.peak']] |> unique()
      list(
        detector.peak = detector,
        .mean = sapply(detector,function(.detector){
          sort(.SD[[.detector]],decreasing = T)[1:top.n] |> mean()
        })
      )
    },
    keyby = .(sample.id,population)
  ]
  ##a single peak detector per sample;
  ##population-specific;
  ##retain all 'Unstained' populations for use as negative controls and autofluorescence
  unstained <- detectors.mean[grepl("Unstained",sample.id)]
  i <- detectors.mean[
    i = !grepl("Unstained",sample.id),
    j = .I[which.max(.mean)],
    by=sample.id
  ][['V1']]
  detectors.mean <- rbind(detectors.mean[i],unstained)
  ##factor detectors
  detectors.mean[,detector.peak:=factor(detector.peak,levels=cols.detector)]
  ##add to [['data']]
  for(i in detectors.mean[,seq(.N)]){
    data.table::set(
      x = flowstate.object$data,
      i = flowstate.object$data[
        ,
        .I[sample.id == detectors.mean[i,as.character(sample.id)] &
             population == detectors.mean[i,as.character(population)]
        ]
      ],
      j = c('detector.peak','select.detector'),
      value = list(detectors.mean[i,detector.peak],TRUE)
    )
  }
}
select.spectral.events <- function(
    flowstate.object,
    n.spectral.events = 200,
    n.spectral.events.unstained = 1000,
    n.spectral.events.override = NULL
)
{
  ##initialize a logical
  flowstate.object$data[,select.spectral := FALSE]
  ##select for top n spectral events;
  ##sample-specific/population-specific/detector-specific;
  ##over-sample for unstained;
  ##override selection (if needed to avoid edge case)
  select.i <- flowstate.object$data[
    ,
    j = {
      detector.peak <- as.character(.BY$detector.peak)
      sample.id <- as.character(.BY$sample.id)
      n <- ifelse(grepl("Unstained",sample.id),n.spectral.events.unstained,n.spectral.events)
      if(sample.id %in% names(n.spectral.events.override)){
        n <- n.spectral.events.override[[sample.id]]
      }
      .I[order(.SD[[detector.peak]],decreasing = T)][1:n]
    },
    by = .(sample.id,detector.peak,population)
  ][['V1']]
  ##update [['data']]
  flowstate.object$data[select.i,select.spectral := TRUE]
}
#' @title A `flowstate` containing spectrally-associated bead/cellular events.
#' @description
#' This function is entirely dependent -- as of now -- on the following:
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
#' @param raw.reference.group \link[base]{file.path} -- one of:
#' \itemize{
#'   \item a reference group directory containing raw single-color/universal negative control .fcs files; all will be processed.
#'   \item a vector of files within that directory; used to process a subset of controls.
#' }
#' @param transform Logical -- default `FALSE`; if `TRUE`, parameters (detectors) having a keyword-value pair of  `TYPE/Raw_Fluorescence` will be transformed using [asinh(x/5000)][base::asinh].
#' @param quantiles Named \link[base]{list} -- default `NULL`; if defined, the list must conform to the following:
#' \itemize{
#'   \item `list('scatter' = probs,'fluor' = probs)`, where [probs][stats::quantile] defines a numeric vector of length 2.
#'   \itemize{
#'     \item e.g., `list('scatter' = c(0,0.99),'fluor' = c(0,0.999))`; as defined, 'scatter' (forward and side) and 'fluor' (fluorescence) associated events will be excluded from processing based on their respective `probs`.
#'   }
#' }
#' @param population.markers Named character vector -- default `NULL`; `population.markers` must be defined for the function to be successful (for cellular controls). The named character vector should take the following form:
#' \itemize{
#'   \item `c(population.name1 = 'sample.name1',population.name2 = 'sample.name2',...)`, where `population.name(s)` are lineage cell types and `'sample.name(s)'` are lineage markers used to stain those respective cell types.
#'   \itemize{
#'     \item e.g., `c(lymphocytes = 'CD3 BV510 (Cells)', monocytes = 'CD14 SB550 (Cells)'`.
#'   }
#' }
#' The defined `population.markers` will be used as 'cellular anchors': top expressing events (peak detector) in the named sample will be used to 'anchor' scatter distributions for proper assignment of the named population.
#' @param top.N.percent Numeric -- default `1`; defines the number of top expressing events (peak detector) used in conjunction with `population.markers` for detecting and assigning the named populations.
#' @param height.cut.scatter Numeric -- default `0.5`; defines the value at which scatter peak heights will be cut for fine tuning selected populations (as defined through `population.markers`).
#' @param n.detector.events Numeric -- default `100`; the number of maximally expressing events (sorted vector) used to auto-detect reference control-specific peak detectors.
#' @param detector.override Named character vector -- default `NULL`; if defined, the supplied detector name (value) will override the auto-detected peak detector on a reference control-specific (name) basis. See example.
#' @param n.spectral.events Numeric -- default `250`; the number of maximally expressing events (sorted vector) used to define 'spectral events'.
#' @param n.spectral.events.unstained Numeric -- default `1000`; the number of maximally expressing events (sorted vector) used to define '(universal) negative' events.
#' @param n.spectral.events.override Named numeric vector -- default `NULL`; if defined, the supplied numeric (value) will override `n.spectral.events` on a sample-specific (name) basis. See example.
#' @param plot Logical -- default `FALSE`; plots diagnostic/QC output for evaluating function performance.
#' @param plot.dir Character vector -- default `'$PROJ'`; the value associated with the keyword '$PROJ' (experiment name) will be used to construct a plot output directory.
#' @param verbose Logical -- default `FALSE`; print function-associated messages to the console.
#'
#' @returns A `flowstate` containing spectrally-associated bead/cellular events.
#' @export
#'
#' @examples
#' \dontrun{
#' ##single-color reference controls
#' raw.ref.files <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "Beads")
#'
#' ##simulate a SpectroFlo directory structure
#' raw.ref.dir <- file.path(tempdir(),"Raw","Reference Group")
#' if(!dir.exists(raw.ref.dir)) dir.create(raw.ref.dir,recursive = TRUE)
#' file.copy(
#'   from = raw.ref.files,
#'   to = file.path(raw.ref.dir,basename(raw.ref.files))
#' ) |> invisible()
#'
#' ##generate spectral events
#' ref <- reference.group.spectral.events(
#'   raw.reference.group.directory = raw.ref.dir,
#'   cluster.types = 'beads',
#'   plot.select = TRUE,
#'   verbose = TRUE
#' )
#'
#' ##a flowstate containing spectral events
#' class(ref)
#' ref[['data']]
#'
#' ##peak detectors
#' ref$data[,.SD,.SDcols = is.factor] |> unique()
#'
#' ##number of spectral events as defined using defaults
#' ref$data[,.N,by = sample.id]
#'
#' ##saved plot output
#' list.files(tempdir()) |> grep(pattern = "select.*.pdf",value = TRUE)
#'
#' ##detector.override
#' ##for the sake of example
#' ref <- reference.group.spectral.events(
#'   raw.reference.group.directory = raw.ref.dir,
#'   cluster.types = 'beads',
#'   detector.override = c("Unstained (Beads)" = "UV7")
#' )
#' ref$data[,.SD,.SDcols = is.factor] |> unique()
#'
#' ##n.spectral.events.override
#' ##for the sake of example
#' ref <- reference.group.spectral.events(
#'   raw.reference.group.directory = raw.ref.dir,
#'   cluster.types = 'beads',
#'   n.spectral.events.override = c("TNFa PE (Beads)" = 500)
#' )
#'
#' ##number of spectral events as defined using n.spectral.events.override
#' ref$data[,.N,by = sample.id]
#' }
#'
reference.group.spectral.events <- function(
    raw.reference.group,
    transform = FALSE,
    quantiles = NULL,
    population.markers = NULL,
    top.N.percent = 1,
    height.cut.scatter = 0.5,
    n.detector.events = 200,
    detector.override = NULL,
    n.spectral.events = 200,
    n.spectral.events.unstained = 1000,
    n.spectral.events.override = NULL,
    plot = FALSE,
    plot.dir = '$PROJ',
    verbose = FALSE
)
{
  if(!all(grepl(".fcs",raw.reference.group)) & length(raw.reference.group) == 1){
    ##get paths to raw .fcs files if a directory;
    ##reference group controls
    ref.paths <- list.files(raw.reference.group,full.names = T)
  }else if(all(grepl(".fcs",raw.reference.group)) & length(raw.reference.group) !=1){
    ref.paths <- raw.reference.group
  }
  ##get keyword identifier for use in defining 'sample.id' argument
  .sample.id <- cytometer.identifier(ref.paths)
  ##test for 'Raw' fluorescence
  raw.fluorescence.check(ref.paths)
  ##test 'population.markers' argument; need to stop the function here if not defined
  if(is.null(population.markers) & any(grepl("cells",ref.paths,ignore.case = T))){
    sample.ids <- check.keyword(ref.paths,keyword = .sample.id)
    sample.ids <- grep("cells",sample.ids,ignore.case = T,value = T)
    sample.ids <- paste0(sample.ids,collapse = " ; ")
    stop(
      paste(
        "For the effective processing of cellular controls, please define",
        "'population.markers' as documented in '?flowstate::reference.group.spectral.events'.",
        "Choose from the following:",
        sample.ids,
        sep = "\n"
      )
    )
  }
  ##read multiple files; concatenate
  if(verbose) message("Reading/concatenating .fcs files...")
  ref <- read.flowstate(
    ref.paths,
    colnames.type = 'N',
    sample.id = .sample.id,
    concatenate = TRUE
  ) |> suppressMessages()
  ##add keywords for eventual use in data.table by/group operations
  add.reference.keywords.to.data(ref)
  ##select against saturating events -- scatter and fluors; mark saturating events as FALSE
  select.nonsaturating(ref)
  ##subset the flowstate.object to retain only non-saturating events;
  ##a computational hit here (RAM) if the object is big due to assignment
  if(verbose) message("Subsetting to remove saturating events...")
  ref <- subset(ref,select.nonsaturating);ref$data[,select.nonsaturating := NULL]
  if(!is.null(quantiles) & 'scatter' %in% names(quantiles)){
    ##select against quantile events -- scatter; mark out of bounds events as FALSE
    select.quantile(ref,type='scatter',probs=quantiles[['scatter']])
    ##subset the flowstate.object to retain only in-bounds events;
    ##a computational hit here (RAM) if the object is big due to assignment
    if(verbose) message("Subsetting to remove out-of-bounds (quantiles) scatter events...")
    ref <- subset(ref,select.quantile);ref$data[,select.quantile := NULL]
    invisible(gc())
  }
  if(transform){
    ##transform detectors
    if(verbose) message("Transforming raw fluorescence (detectors)...")
    flowstate.transform(ref) |> suppressMessages()
  }
  if(!is.null(quantiles) & 'fluor' %in% names(quantiles)){
    ##select against quantile events -- fluors; mark out of bounds events as FALSE
    select.quantile(ref,type='fluor',probs=quantiles[['fluor']])
    ##subset the flowstate.object to retain only in-bounds events;
    ##a computational hit here (RAM) if the object is big due to assignment
    if(verbose) message("Subsetting to remove out-of-bounds (quantiles) fluorescent events...")
    ref <- subset(ref,select.quantile);ref$data[,select.quantile := NULL]
    invisible(gc())
  }
  ##population-specific scatter selection; outputs diagnostic plots
  if(plot){
    plot.dir <- ref$keywords[[plot.dir]] |> unique()
    if(length(plot.dir)==1){
      plot.dir <- file.path("./data_results",plot.dir)
      if(!dir.exists(plot.dir)) dir.create(plot.dir,recursive = T)
    }else{
      message("Could not resolve 'plot.dir' to a unique directory path;\n
              defaulting to tempdir().")
      plot.dir <- tempdir()
    }
  }
  if(plot){
    file.out <- file.path(
      plot.dir,
      sprintf("%s_%s_select_population_density_cuts.pdf",basename(plot.dir),Sys.Date())
    )
    grDevices::pdf(file.out,width = 8, height = 8)
  }
  ##cellular-based/marker-based fluor anchors to select appropriate scatter populations
  if(verbose) message("Using cellular-based markers/fluors to anchor scatter populations...")
  cellular.anchors(
    flowstate.object = ref,
    top.N.percent = top.N.percent,
    population.marker = population.markers,
    height.cut.scatter = height.cut.scatter,
    plot = plot
  )
  if(plot){
    if(verbose) message(sprintf("Saving %s",file.out))
    grDevices::dev.off()
  }
  ##plot selection results
  if(plot){
    file.out <- file.path(
      plot.dir,
      sprintf("%s_%s_select_population_NxN.pdf",basename(plot.dir),Sys.Date())
    )
    grDevices::pdf(file.out,width = 8, height = 8)
    print(plot_populations(ref,pattern.scatter = "AH",sample.n = 1E4))
    if(verbose) message(sprintf("Saving %s",file.out))
    grDevices::dev.off()
  }
  ##subset the flowstate.object to retain selected populations;
  ##a computational hit here (RAM) if the object is big due to assignment
  if(verbose) message("Subsetting to retain sample-specific scatter populations...")
  ref <- subset(ref,!is.na(population))
  invisible(gc())
  ##assign sample-specific/population-specific peak detectors
  ##use variance to detect peak detector;
  ##mean of sorted vectors within peak detector to resolve to specific population (max expressing)
  if(verbose) message("Variance + mean (sorted vectors) to determine peak detectors and max-expressing scatter populations...")
  detector.peak(ref,top.n = n.detector.events)
  ##subset the flowstate.object to retain selected detector-specific populations;
  ##a computational hit here (RAM) if the object is big due to assignment
  if(verbose) message("Subsetting to retain sample-specific max expressing (peak detector) scatter population...")
  ref <- subset(ref,select.detector);ref$data[,select.detector := NULL]
  invisible(gc())
  ##order [['data']] by 'detector.peak' factor levels
  data.table::setorder(ref$data,detector.peak)
  ##select top n expressing 'spectral events'
  if(verbose) message("Selecting top n 'spectral events'...")
  select.spectral.events(
    flowstate.object = ref,
    n.spectral.events = n.spectral.events,
    n.spectral.events.unstained = n.spectral.events.unstained,
    n.spectral.events.override = n.spectral.events.override
  )
  ##return -- invisibly -- the spectral events
  if(verbose) message("Returning 'spectral events'.")
  invisible(ref)
}
#' @title Plot Spectral Events
#' @description
#' Plot the results of [reference.group.spectral.events]; visual assessment of selected 'spectral events' is used to gauge function performance.
#'
#' @param flowstate.object.reference A `flowstate` object as returned from [reference.group.spectral.events].
#'
#' @returns A list of [ggplot][ggplot2::ggplot] objects.
#' @export
#'
plot_spectral.events <- function(flowstate.object.reference){
  ##alias column
  j.match <- j.match.parameters.to.data(flowstate.object.reference)
  ##detector columns
  cols.detector <- flowstate.object.reference$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##factor columns
  cols.by <- flowstate.object.reference$data[,names(.SD),.SDcols = is.factor]
  ##subset to detector.peak-specific values
  dt <- flowstate.object.reference$data[
    ,
    j = {
      detector.peak <- as.character(.BY$detector.peak)
      list(value = .SD[[detector.peak]],select.spectral)
    },
    by = cols.by
  ]
  ##add an index for plotting purposes
  dt[,index := seq(.N),by=cols.by]
  ##split the data.table and plot
  splits <- split(dt,by=c('sample.id','population'),sorted = F,drop = T)
  splits <- lapply(splits,droplevels)
  p <- lapply(splits,function(.dt){
    .dt
    ggplot2::ggplot(data = NULL) +
      ggplot2::geom_hex(
        data = .dt[(select.spectral)],
        mapping = ggplot2::aes(x = index, y = value),
        bins = 100,
        fill = "red"
      ) +
      ggplot2::geom_hex(
        data = .dt[!(select.spectral)],
        mapping = ggplot2::aes(x = index, y = value),
        bins = 100,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_gradient(low = 'lightgray',high = 'darkgray') +
      ggplot2::theme_light() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = sprintf("%s\nDetector (Peak): %s\nScatter Population: %s",
                        .dt[,levels(sample.id)],.dt[,levels(detector.peak)],.dt[,levels(population)]),
        caption = sprintf("%d 'spectral' events (out of %d) selected in %s for calculating per-detector (n = %d) medians.",
                          .dt[(select.spectral),.N], .dt[,.N], .dt[,levels(detector.peak)], length(cols.detector)),
        x = 'Index',
        y = .dt[,levels(detector.peak)]
      )
  })
  return(p)
}
##
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
