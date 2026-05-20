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
#' @param height.cut.scatter Numeric -- default `0.25`; defines the value at which scatter peak heights will be cut for fine tuning selected populations (as defined through `population.markers`).
#' @param n.detector.events Numeric -- default `100`; the number of maximally expressing events (sorted vector) used to auto-detect reference control-specific peak detectors.
#' @param detector.override Named character vector -- default `NULL`; if defined, the supplied detector name (value) will override the auto-detected peak detector on a reference control-specific (name) basis. See example.
#' @param n.spectral.events Numeric -- default `200`; the number of maximally expressing events (sorted vector) used to define 'spectral events'.
#' @param n.spectral.events.unstained Numeric -- default `1000`; the number of maximally expressing events (sorted vector) used to define '(universal) negative' events.
#' @param n.spectral.events.override Named numeric vector -- default `NULL`; if defined, the supplied numeric (value) will override `n.spectral.events` on a sample-specific (name) basis. See example.
#' @param plot Logical -- default `FALSE`; plots diagnostic/QC output for evaluating function performance.
#' @param plot.dir Character vector -- default `'$PROJ'`; the value associated with the keyword '$PROJ' (experiment name) will be used to construct a plot output directory.
#' @param return.spectral.events Logical -- default `TRUE`; the returned `flowstate` will be subset to include only 'spectral events' -- top expressing events used to calculate medians to define spectra. If `FALSE`, the returned `flowstate` will contain both 'spectral events' and non-'spectral events' -- useful for diagnostic/QC purposes.
#' @param verbose Logical -- default `FALSE`; print function-associated messages to the console.
#'
#' @returns A `flowstate` containing spectrally-associated bead/cellular events.
#' @keywords internal
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
    height.cut.scatter = 0.25,
    n.detector.events = 100,
    detector.override = NULL,
    n.spectral.events = 200,
    n.spectral.events.unstained = 1000,
    n.spectral.events.override = NULL,
    plot = FALSE,
    plot.dir = '$PROJ',
    return.spectral.events = TRUE,
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
    if(plot.dir == '$PROJ'){
      plot.dir <- ref$keywords[[plot.dir]] |> unique()
      if(length(plot.dir)==1){
        plot.dir <- file.path("./data_results",plot.dir)
      }else{
        message("Could not resolve 'plot.dir' to a unique directory path;\n
              defaulting to tempdir().")
        plot.dir <- tempdir()
      }
    }
    if(!dir.exists(plot.dir)) dir.create(plot.dir,recursive = T)
  }
  # if(plot){
  #   plot.dir <- ref$keywords[[plot.dir]] |> unique()
  #   if(length(plot.dir)==1){
  #     plot.dir <- file.path("./data_results",plot.dir)
  #     if(!dir.exists(plot.dir)) dir.create(plot.dir,recursive = T)
  #   }else{
  #     message("Could not resolve 'plot.dir' to a unique directory path;\n
  #             defaulting to tempdir().")
  #     plot.dir <- tempdir()
  #   }
  # }
  if(plot){
     file.out <- file.path(
      plot.dir,
      sprintf("%s_flowstate_population_density_cuts_%s.pdf",
              basename(plot.dir),Sys.Date())
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
      sprintf("%s_flowstate_population_NxN_%s.pdf",
              basename(plot.dir),Sys.Date())
    )
    grDevices::pdf(file.out,width = 8, height = 8)
    if(verbose) message(sprintf("Saving %s",file.out))
    print(plot_populations(ref,pattern.scatter = "AH",sample.n = 1E4))
    grDevices::dev.off()
  }
  ##subset the flowstate.object to retain selected populations;
  ##a computational hit here (RAM) if the object is big due to assignment
  if(verbose) message("Subsetting to retain sample-specific scatter populations...")
  ref <- subset(ref,!is.na(population))
  invisible(gc())
  if(plot){
    file.out <- file.path(
      plot.dir,
      sprintf("%s_flowstate_population_NxN_subset_%s.pdf",
              basename(plot.dir),Sys.Date())
    )
    grDevices::pdf(file.out,width = 8, height = 8)
    if(verbose) message(sprintf("Saving %s",file.out))
    print(plot_populations(ref,pattern.scatter = "AH",sample.n = 1E4))
    grDevices::dev.off()
  }
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
  if(plot){
    file.out <- file.path(
      plot.dir,
      sprintf("%s_flowstate_ref_spectral_events_%s.pdf",
              basename(plot.dir),Sys.Date())
    )
    flowstate.transform(ref) |> suppressMessages()
    p <- plot_spectral.events(ref)
    grDevices::pdf(file.out,width = 8,height = 6)
    if(verbose) message(sprintf("Saving %s",file.out))
    for(i in names(p)){
      print(p[[i]])
    }
    grDevices::dev.off()
  }
  if(return.spectral.events){
    ##subset
    ref <- subset(ref,select.spectral);ref$data[,select.spectral := NULL]
  }
  ##
  if('transform' %in% names(ref$parameters)){
    flowstate.transform.inverse(ref) |> invisible()
    if(all(is.na(ref$parameters$transform))){
      ref$parameters[
        ,
        c('transform','cofactor') := NULL
      ]
    }
  }
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
#' @keywords internal
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
