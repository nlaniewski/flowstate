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
    stop("Defined 'population.marker(s)' not found among celluar ('(Cells)') samples.")
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
  p <- p + ggplot2::labs(
    title = sprintf(
      "Populations (Anchored): %s",
      ref$data[,levels(population) |> paste0(collapse = "; ")]
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
plot_spectral.events <- function(flowstate.object){
  ##alias column
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##detector columns
  cols.detector <- flowstate.object$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##factor columns
  cols.by <- flowstate.object$data[,names(.SD),.SDcols = is.factor]
  ##subset to detector.peak-specific values
  dt <- flowstate.object$data[
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
#' @title A `flowstate` containing spectrally-associated bead/cellular events.
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
#' @param cluster.types Character vector; some combination of `c('beads','lymphocytes','monocytes')`.  An attempt will be made to auto-detect these scatter-based populations based on some rudimentary [clustering][stats::kmeans()].
#' @param detection.method Character vector `1L`; one of \link[stats]{var} or \link[base]{mean}.  The defined statistical method will be used to auto-detect reference control-specific peak detectors.
#' @param n.detector.events Numeric -- default `100`; the number of maximally expressing events (sorted vector) used to auto-detect reference control-specific peak detectors.
#' @param detector.override Named character vector -- default `NULL`; if defined, the supplied detector name (value) will override the auto-detected peak detector on a reference control-specific (name) basis. See example.
#' @param n.spectral.events Numeric -- default `250`; the number of maximally expressing events (sorted vector) used to define 'spectral events'.
#' @param n.spectral.events.unstained Numeric -- default `1000`; the number of maximally expressing events (sorted vector) used to define '(universal) negative' events.
#' @param n.spectral.events.override Named numeric vector -- default `NULL`; if defined, the supplied numeric (value) will override `n.spectral.events` on a sample-specific (name) basis. See example.
#' @param plot.select Logical; default `FALSE`. Plots diagnostic/QC output for evaluating function performance.
#' @param plot.select.dir Character vector; a file path for storing plot output.
#' @param verbose Logical; default `FALSE`. Print messages to the console.
#'
#' @returns A `flowstate` containing spectrally-associated bead/cellular events.
#' @export
#'
#' @examples
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
#'
reference.group.spectral.events <- function(
    raw.reference.group.directory,
    cluster.types = c('beads','lymphocytes','monocytes'),
    detection.method = c('var','mean'),
    n.detector.events = 100,
    detector.override = NULL,
    n.spectral.events = 250,
    n.spectral.events.unstained = 1000,
    n.spectral.events.override = NULL,
    plot.select = FALSE,
    plot.select.dir = tempdir(),
    verbose = FALSE
)
{
  ##paths to raw .fcs files; reference group controls
  ref.paths <- list.files(raw.reference.group.directory,full.names = T)
  ##test for SpectroFlo
  if(!check.keyword(ref.paths,keyword = 'CREATOR',value = "SpectroFlo")){
    message("This function has only been tested with SpectroFlo-generated .fcs files and may not perform correctly or at all when using other instrument/software-generated .fcs files.")
  }
  ##test for 'Raw' fluorescence
  sapply(ref.paths,raw.fluorescence.check) |> invisible()
  ##read multiple files; concatenate
  if(verbose) message("Reading/concatenating .fcs files...")
  ref <- read.flowstate(
    ref.paths,
    colnames.type = 'N',
    concatenate = TRUE
  ) |> suppressMessages()
  ##add keywords for eventual use in data.table by/group operations
  add.reference.keywords.to.data(ref)
  ##select against saturating events -- scatter and fluors; mark saturating events as FALSE
  select.nonsaturating(ref)
  ##alias column
  j.match <- j.match.parameters.to.data(ref)
  ##detector columns
  cols.detector <- ref$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  ##transform detectors
  if(verbose) message("Transforming raw fluorescence (detectors)...")
  flowstate.transform(ref) |> suppressMessages()
  ##rudimentary kmeans to capture 2 scatter populations;
  ##for PBMC, approximates a lymphocyte cluster and a high-scatter/monocyte cluster;
  ##peak-finding can also work, but need a well defined monocyte population...not always present
  cols.scatter <- grep("[FS]S",names(ref$data),value = T)
  if('cells' %in% ref$data[,levels(group.type)]){
    if(verbose) message("Kmeans to define lymphocytes and monocytes (scatter)...")
    ref$data[
      i = (select.nonsaturating) & group.type == 'cells',
      j = cluster := {
        set.seed(1337)
        js |>
          data.table::as.data.table() |>
          scale() |>
          stats::kmeans(x=_,centers = 2) |>
          _[['cluster']]
      },
      env = list(js = as.list(cols.scatter))
    ]
    labels <- c('lymphocytes','monocytes')
    cluster.lymphocytes <- ref$data[!is.na(cluster),.N,keyby = cluster][,which.max(N)]
    if(cluster.lymphocytes != 1){
      labels <- rev(labels)
    }
    ref$data[,cluster := factor(cluster,labels = labels)]
  }
  ##beads as is; add to 'cluster' factor level
  ref$data[group.type == 'beads', cluster := 'beads']
  ##drop and refactor if 'cluster.types' is defined
  ref$data[!cluster %in% cluster.types,cluster := NA]
  ref$data[,cluster := factor(cluster)]
  ##clean up the cluster assignments;
  ##add a logical -- 'select.scatter' -- for selecting expected distributions among either lymphocytes or monocytes
  cols.scatter.subset <- grep("_A",cols.scatter,value = T)
  clusters <- ref$data[,levels(cluster)]
  if(plot.select){
    file.out.name <- sprintf("%s_select_scatter.pdf", unique(ref$keywords$`$PROJ`))
    if(!dir.exists(plot.select.dir)) dir.create(plot.select.dir,recursive = T)
    grDevices::pdf(file.path(plot.select.dir,file.out.name),width = 6,height = 8)
  }
  graphics::par(
    mfcol = c(length(cols.scatter.subset), length(clusters)),
    oma = c(0, 0, 1, 0)
  )
  for(.cluster in clusters){
    ref$data[
      i = cluster == .cluster,
      j = select.scatter := {
        for(i in cols.scatter.subset){
          cuts <- peak.height.cut(
            x = if(!exists('vec')){
              .SD[[i]]
            }else{
              .SD[[i]][vec]
            },
            height.cut = 0.25,
            which.peak = 1,
            plot = ifelse(plot.select,TRUE,FALSE),
            main = sprintf("%s Select  (%s)",i,.BY$cluster),
            xlim = c(0,4194304)#generalize this? [['parameters']][['R']]
          )
          if(!exists('vec')){
            vec <- data.table::`%between%`(.SD[[i]],cuts)
          }else{
            vec[which(vec)] <- data.table::`%between%`(.SD[[i]][vec],cuts)
          }
        }
        vec[]
      },
      by = cluster
    ]
  }
  if(plot.select){
    grDevices::dev.off()
    graphics::par(mfcol = c(1, 1), oma = c(0, 0, 0, 0))
  }
  detection.method <- match.arg(detection.method)
  if(detection.method == 'var'){
    ##variance per-sample/per-detector to detect peak detector;
    ##will produce edge cases in certain circumstances (e.g., scatter-specific populations)
    if(verbose) message("Using per-detector variance measure to determine sample-specific peak detectors...")
    ref$data[
      i = (select.scatter),
      j = detector.peak := {
        sapply(.SD,stats::var) |> which.max() |> names()
      },
      .SDcols = cols.detector,
      by = .(sample.id)
    ]
  }else if(detection.method == 'mean'){
    ##mean of top n.detector.events per-sample/per-detector to detect peak detector;
    ##variance measure sometimes fails; due to heterogeneous cells?
    ##this seems to work for 'scatter selected' events; best so far for edge cases
    ##computational hit though, with having to sort each j (n detectors) for each control
    if(verbose) message("Using mean value of sorted per-detector vectors (top n events) to determine sample-specific peak detectors...")
    ref$data[
      i = (select.scatter),
      j = detector.peak := {
        sapply(.SD,function(j){
          sort(j,decreasing = T)[1:n.detector.events] |>
            mean()
        }) |> which.max() |> names()
      },
      .SDcols = cols.detector,
      by = .(sample.id)
    ]
  }
  ##update keywords with sample-specific peak detector; factor
  res <- ref$data[
    i = (select.scatter),
    j = .(detector.peak = unique(detector.peak)),
    by = .(sample.id)
  ]
  if(!is.null(detector.override)){
    res[
      i = sample.id %in% names(detector.override),
      j = detector.peak := detector.override
    ]
  }
  ref$keywords[,detector.peak := res[['detector.peak']]]
  ref$keywords[,detector.peak := factor(detector.peak,levels = cols.detector)]
  add.keywords.to.data(ref,'detector.peak')
  ##now grouping by sample, detector, and cluster;
  ##find within scatter population (beads, lymphs, monos) peak detector means
  res <- ref$data[
    i = (select.scatter),
    j = .(means = {
      detector.peak <- as.character(.BY$detector.peak)
      sort(.SD[[detector.peak]],decreasing = T)[1:n.detector.events] |>
        mean()
    }),
    keyby = .(sample.id,cluster,detector.peak)
  ]
  ##select max mean per-sample to select appropriate (max-expressing) scatter population
  select.i <- res[
    ,
    j = .I[which.max(means)],
    by = .(sample.id)
  ][['V1']]
  res <- res[select.i]
  ##update [['data']]
  for(i in res[,seq(.N)]){
    data.table::set(
      x = ref$data,
      i = ref$data[
        ,
        .I[(select.scatter)
           & sample.id == res[i,as.character(sample.id)]
           & cluster == res[i,as.character(cluster)]
        ]
      ],
      j = c('select.detector'),
      value = TRUE
    )
  }
  ##for unstained cells, retain scatter populations
  ref$data[
    i = grepl("Unstained",sample.id) & group.type == 'cells' & (select.scatter),
    j = select.detector := TRUE
  ]
  ##order [['data']] by 'detector.peak'
  data.table::setorder(ref$data,detector.peak)
  ##subset to detector-specific events
  ref <- subset(ref,(select.detector))
  ##NULL redudant columns
  ref$data[,grep("select",names(ref$data),value = T) := NULL]
  ##select for top n spectral events;
  ##sample-specific/scatter population-specific/detector-specific
  ##over-sample for unstained
  if(verbose) message("Selecting top n 'spectral events'...")
  ref$data[,select.spectral := FALSE]
  select.i <- ref$data[
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
    by = .(sample.id,detector.peak,cluster)
  ][['V1']]
  ref$data[select.i,select.spectral := TRUE]
  ##plot selected spectral events
  if(plot.select){
    if(verbose) message("Printing 'Spectral Events' to device (.pdf)...")
    file.out.name <- sprintf("%s_select_events.pdf",unique(ref$keywords$`$PROJ`))
    grDevices::pdf(file.path(plot.select.dir,file.out.name),width = 7,height = 7)
    ##
    ref$data[
      ,#i = sample.id == 'CD11b Spark UV387 (Cells)',
      j = {
        # message(.BY$sample.id)
        ##variables
        l <- length(select.spectral)
        index <- seq_len(l)
        l.select <- length(which(select.spectral))
        detector.peak <- as.character(.BY$detector.peak)
        ##plot
        p <- ggplot2::ggplot() +
          ggplot2::geom_hex(
            data = NULL,
            mapping = ggplot2::aes(
              x = index[!(select.spectral)],
              y = .SD[[detector.peak]][!(select.spectral)]),
            bins = 100,
            alpha = 0.6,
            show.legend = F
          ) +
          ggplot2::scale_fill_gradient(low = 'gray',high = 'black') +
          ggplot2::geom_hex(
            data = NULL,
            mapping = ggplot2::aes(
              x = index[(select.spectral)],
              y = .SD[[detector.peak]][(select.spectral)]),
            fill = 'red',
            bins = 100,
            show.legend = F
          ) +
          ggplot2::labs(
            title = sprintf("%s\nDetector (Peak): %s\nScatter Population: %s",
                            .BY$sample.id,detector.peak,.BY$cluster),
            caption = sprintf("%d 'spectral' events (out of %d) selected in %s for calculating per-detector (n = %d) medians.",
                              l.select, l, detector.peak, length(.SD)),
            x = 'Index',
            y = detector.peak
          )
        print(p)
      },
      .SDcols = cols.detector,
      by = c(ref$data[,names(.SD),.SDcols = is.factor])
    ]
    ##
    grDevices::dev.off()
  }
  ##subset to include only 'spectral events'
  ref <- subset(ref,select.spectral==TRUE)
  ref$data[,select.spectral := NULL]
  ##transform inverse; [['data']] to linear/raw
  flowstate.transform.inverse(ref)
  ##return the flowstate.object
  invisible(ref)
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
##
reference.group.medians <- function(flowstate.object.reference,name.fix=NULL,syntactically.valid=FALSE){
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
  ref.medians[,ord := seq(.N)]
  ref.medians[N=='AF',ord := max(ref.medians[['ord']]) + seq(.N)]
  data.table::setorder(ref.medians,ord)[,ord := NULL]#AF in last position
  if(!is.null(name.fix)){
    for(i in seq_along(name.fix)){
      ref.medians[N == names(name.fix[i]), N := name.fix[[i]]]
    }
  }
  if(syntactically.valid){
    ref.medians[,N := tolower(gsub(" |-|\\.","",N))]
    ref.medians[,S := gsub("-","",S)]
  }
  ref.medians[is.na(S), alias := N]
  ref.medians[!is.na(S),alias := paste(S,N)]
  ##return the normalized reference control medians
  ref.medians[]
}
##
