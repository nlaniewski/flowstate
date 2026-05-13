raw.fluorescence.check <- function(fcs.file.paths){
  invisible(sapply(fcs.file.paths, function(i){
    parms <- parameters.to.data.table(readFCStext(i))
    if(!'TYPE' %in% names(parms)){
      stop("Are these .fcs files processed by flowstate? 'TYPE' is missing from [['parameters']].")
    }else{
      if(!"Raw_Fluorescence" %in% parms[['TYPE']]){
        stop("'Raw_Fluorescence' value not found in 'TYPE' keyword; is this a raw reference control?")
      }
    }
  }))
}

## function to add metadata to [['keywords]] and [['data']]
reference.group.keywords <- function(flowstate){
  ## column identifier for [['data']] and [['keywords']]
  if(!'sample.id' %in% intersect(names(flowstate$data), names(flowstate$keywords))){
    stop("'sample.id' identifier not found.")
  }
  ## CREATOR/software
  software <- ref$keywords[, unique(CREATOR)]
  stopifnot(
    "Mixed software (keyword: CREATOR); can't reliably derive metadata" = length(software) == 1
  )
  if(grepl("SpectroFlo", software)) reference.group.keywords.spectroflo(flowstate)
  # ## add keyword metadata based on splitting sample.id;
  # ## following SpectroFlo naming convention: marker(S) fluorophore(N) (type) -- literal spaces separating
  # keywords.to.add <- c('tissue.type', 'N', 'S')
  # ## test
  # if(any(flowstate$data[, unique(sample.id)] != flowstate$keywords[, unique(sample.id)])){
  #   stop("Sample order does not match between [['data']] and [['keywords']].")
  # }
  # ## add keyword/value pairs based on gsub/regex; dependent on SpectrolFlo naming convention
  # ## updates [['keywords']]
  # flowstate$keywords[
  #   ,
  #   j = (keywords.to.add) := {
  #     tissue.type <- factor(gsub("^.*\\((.*?)\\).*$", "\\1", sample.id))
  #     res <- sub(" \\(.*$", "", sample.id)
  #     marker <- factor(ifelse(
  #       tissue.type %in% c("Beads", "Cells") & grepl("unstained", res, ignore.case = T),
  #       NA,
  #       sub(" +.*$", "", res)
  #     ))
  #     fluorophore <- factor(ifelse(
  #       tissue.type %in% c("Beads", "Cells") & is.na(marker),
  #       "AF",
  #       sub(".*? ", "", res)
  #     ))
  #     list(tissue.type, N = fluorophore, S = marker)
  #   }
  # ]
  ## updates [['data']]
  add.keywords.to.data(flowstate, keywords.to.add)
  ## return
  invisible(flowstate)
}

reference.group.keywords.spectroflo <- function(flowstate){
  ## add keyword metadata based on splitting sample.id;
  ## following SpectroFlo naming convention: marker(S) fluorophore(N) (type) -- literal spaces separating
  keywords.to.add <- c('tissue.type', 'N', 'S')
  ## test
  if(any(flowstate$data[, unique(sample.id)] != flowstate$keywords[, unique(sample.id)])){
    stop("Sample order does not match between [['data']] and [['keywords']].")
  }
  ## add keyword/value pairs based on gsub/regex; dependent on SpectrolFlo naming convention
  ## updates [['keywords']]
  flowstate$keywords[
    ,
    j = (keywords.to.add) := {
      tissue.type <- factor(gsub("^.*\\((.*?)\\).*$", "\\1", sample.id))
      res <- sub(" \\(.*$", "", sample.id)
      marker <- factor(ifelse(
        tissue.type %in% c("Beads", "Cells") & grepl("unstained", res, ignore.case = T),
        NA,
        sub(" +.*$", "", res)
      ))
      fluorophore <- factor(ifelse(
        tissue.type %in% c("Beads", "Cells") & is.na(marker),
        "AF",
        sub(".*? ", "", res)
      ))
      list(tissue.type, N = fluorophore, S = marker)
    }
  ]
}

#' @title Add scatter-specific populations to `[['data']]`
#' @description
#' Effective processing of single-color controls first involves defining dominant scatter-specific populations that are free of debris/aggregates and contain max-expressing events.  To efficiently select these populations in cellular controls, a defined set of cluster of differentiation (CD)/lineage surface markers are used to restrict peak-finding/cutting methods to the scatter profile of specific immune cell types.  Once these profiles/populations are defined, they are applied to all of `[['data']]`.
#'
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param population.marker Named character vector -- default `NULL`; `population.marker` must be defined for the function to be successful (for cellular controls). The named character vector should take the following form:
#' \itemize{
#'   \item `c(population.name1 = 'sample.name1',population.name2 = 'sample.name2',...)`, where `population.name(s)` are CD/lineage cell types and `'sample.name(s)'` are CD/lineage markers used to stain those respective cell types.
#'   \itemize{
#'     \item e.g., `c(lymphocytes = 'CD3 BV510 (Cells)', monocytes = 'CD14 SB550 (Cells)'`.
#'   }
#' }
#' @param top.N.percent Numeric -- default `10`; defines the percentage of top expressing events (peak detector) used in conjunction with `population.marker` for detecting and assigning the named populations.
#' @param scatter.peak.cut Numeric -- default `0.25`; defines the value at which scatter peak heights will be cut for fine tuning the scatter profile of selected populations (as defined through `population.marker`).
#' @param plot Logical -- default `FALSE`; plots diagnostic/QC output for evaluating function performance.
#' @param plot.output.dir A character vector of a full path name -- default [tempdir()][base::tempdir]; any/all generated plots will be saved to this directory.
#'
#' @returns UPDATES BY REFERENCE and invisibly returns `flowstate`; a new column (factored) named `population` is added to `[['data']]`.
#' @export
#'
select_scatter.population <- function(
    flowstate,
    population.marker,
    top.N.percent = 10,
    scatter.peak.cut = 0.25,
    plot = F,
    plot.output.dir = tempdir()
)
{
  ## relevant variables
  j.match <- j.match.parameters.to.data(flowstate)
  cols.scatter <- (flowstate$parameters
                   [grep("scatter",TYPE,ignore.case = T)][[j.match]])
  cols.scatter <- sort(cols.scatter)
  cols.detector <- (flowstate$parameters
                    [TYPE == "Raw_Fluorescence"][[j.match]])
  cols.by <- flowstate$data[,names(.SD),.SDcols = is.factor]
  proj <- factor(flowstate$parameters[,levels(PROJ)])
  ## if there are any bead controls, add to population marker;
  ## any single fluor bead control
  if(flowstate$keywords[,"Beads" %in% levels(type)]){
    population.marker <- c(population.marker,beads = "Unstained (Beads)")
  }
  ## prepare plotting device to capture graphic output
  if(plot){
    file.out <- file.path(
      plot.output.dir,
      sprintf("reference_group_scatter_populations_density_%s_%s.pdf",
              proj,Sys.Date())
    )
    col.set <- length(population.marker)
    row.set <- length(cols.scatter)
    grDevices::pdf(file.out,width = 4 * col.set, height = 2*row.set)
    graphics::par(mfcol = c(row.set, col.set))
    on.exit(grDevices::dev.off())
  }
  ## sample/type-specific (population.marker for Cells, Unstained for Beads);
  ## use the mean of top-expressing events (sorted vectors) to find peak detector;
  ## using peak detector, return the index of top.N expressing events (ordered vector)
  ## sequentially cut dominant scatter peaks; updates i
  ## return the ranges to define scatter bounds
  scatter.bounds <- flowstate$data[
    i = sample.id %in% population.marker,
    j = {
      if(.BY$type == "Cells"){
        means.detector <- sapply(cols.detector,function(j){
          mean(sort(.SD[[j]],decreasing = T)[1:100])
        })
        detector <- names(which.max(means.detector))
        top.N <- ceiling(.N*(top.N.percent/100))
        i <- order(.SD[[detector]],decreasing = T)[1:top.N]
      }else if(.BY$type == "Beads"){
        i <- sample(.N,ceiling(.N*(top.N.percent/100)))
      }
      ##
      .population <- names(which(population.marker == .BY$sample.id))
      ##
      for(scatter in cols.scatter){
        bounds <- peak.height.cut(
          x = .SD[[scatter]][i],
          # which.peak = 1,
          height.cut = scatter.peak.cut,
          plot = plot,
          main = paste0(scatter,"\n",.population)
        )
        i <- i[data.table::`%between%`(.SD[[scatter]][i],bounds)]
      }
      c(population = .population,
        sapply(cols.scatter,function(j){
          range(.SD[[j]][i])
        },simplify = F)
      )
    },
    .SDcols = c(cols.detector,cols.scatter),
    by = cols.by
  ][,population := factor(population)]
  ## apply scatter bounds to all data, type-specific
  for(.population in scatter.bounds[,levels(population)]){
    .type <- ifelse(.population == 'beads',"Beads","Cells")
    i <- flowstate$data[
      i = type == .type,
      j = {
        rowSums(sapply(cols.scatter,function(j){
          data.table::`%between%`(
            .SD[[j]],
            scatter.bounds[population == .population][[j]]
          )
        })) == length(cols.scatter)
      },
      by = cols.by
    ][['V1']]
    data.table::set(
      x = flowstate$data,
      i = flowstate$data[,.I[type == .type][i]],
      j = 'population',
      value = .population
    )
  }
  ## factor population
  flowstate$data[,population := factor(population)]
  ## return
  invisible(flowstate)
}

.reference.group.spectral.events <- function(
    flowstate,
    top.n = 200,
    top.n.override = NULL
)
{
  ## alias for selecting columns/variables; matches to N ($PnN)
  alias <- merge(
    x = flowstate$parameters[, j = .(N, TYPE)],
    y = alias_dt(flowstate),
    sort = F
  )
  ## detectors for use in .SDcols
  cols.detector <- alias[TYPE == "Raw_Fluorescence"][['alias']]
  ## factors in [['data']] for use in by
  cols.by <- flowstate$data[, names(.SD), .SDcols = is.factor]
  ## additional metadata
  cols.mdat <- names(flowstate$keywords)[names(flowstate$keywords) %in% c('$CYT', '$CYTSN', '$PROJ')]
  mdat <- flowstate$keywords[, unique(.SD), .SDcols = cols.mdat]
  ## unstained/AF 0.995 quantile;
  ## used in subsequent step to help resolve dim/AF-impacted fluors through subtraction
  q.af <- flowstate$data[
    i = N == "AF",
    j = {
      q <- lapply(.SD, stats::quantile, probs = 0.995)
      detector <- names(which.max(q))
      c(detector = detector, q)
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## use the mean of top-expressing events (sorted vectors) to find peak detector;
  ## subtract 'q.af' to resolve dim/AF-impacted fluors (otherwise obscured by AF);
  ## index top expressing events in peak detector (from ordered vector) -- 'spectral events';
  ## return the spectral events
  spectral.events <- flowstate$data[
    ,#i = type == "Beads",
    j = {
      means.detector <- sapply(.SD, function(j){
        mean(sort(j, decreasing = T)[1:100])
      })
      if(.BY$N != "AF"){
        means.detector <- means.detector - q.af[
          i = population == .BY$population,
          j = .SD,
          .SDcols = is.numeric
        ]
      }
      detector <- names(which.max(means.detector))
      if(.BY$N == "AF"){
        .top.n <- ifelse(.N >= 1000, 1000, .N)
      }else if(.BY$sample.id %in% names(top.n.override)){
        .top.n <- top.n.override[[as.character(.BY$sample.id)]]
      }else{
        .top.n <- top.n
      }
      i.top <- order(.SD[[detector]], decreasing = T)[1:.top.n]
      detector.mean <- mean(.SD[[detector]][i.top])
      ##
      c(detector = detector, detector.mean = detector.mean, .SD[i.top])
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## factor detector column; levels equal to order in [['parameters']]
  ## order by detector, N, S, and type
  spectral.events[, detector := factor(detector, levels = cols.detector)]
  data.table::setorder(spectral.events, detector, N, S, tissue.type)
  ## laser column
  spectral.events[, j = laser := gsub("\\d", "", detector)]
  spectral.events[, j = laser := factor(laser, levels = unique(laser))]
  ## reorder
  data.table::setcolorder(
    x = spectral.events,
    neworder = spectral.events[, names(.SD), .SDcols = is.factor]
  )
  ## subset to retain representative max expressing population-specific spectral events;
  ## Cells only
  i.drop <- spectral.events[
    i = N != "AF" & tissue.type == "Cells",
    j = .I[detector.mean != max(detector.mean)],
    by = sample.id
  ][['V1']]
  if(length(i.drop)!=0){
    spectral.events <- spectral.events[-i.drop]
  }
  spectral.events[, detector.mean := NULL][]
  ## add instrument metadata
  data.table::setattr(spectral.events, "mdat", mdat)
  ##
  invisible(spectral.events)
}

spectral.events.select.events <- function(spectral.events){
  ## variables
  cols.detector <- spectral.events[, names(.SD), .SDcols = is.numeric]
  cols.by <- spectral.events[, names(.SD), .SDcols = is.factor]
  ## initialize a logical
  spectral.events[N == 'AF', select.events := !vector(length = .N)]
  ##
  spectral.events[
    i = N != 'AF',
    j = select.events := {
      ## among the sorted (by peak detector) spectral events:
      ## test cosine similarity of the brightest event against the bottom ten events
      res <- sapply(seq(.N-10, .N), function(i){
        cosine.similarity(unlist(.SD[1]), unlist(.SD[i]))
      })
      ## if the mean result is less than 0.985, test by chunks of ten;
      ## set a new index and return a logical for selecting events
      if(mean(res) < 0.985){
        x <- seq(.N)
        chunks <- split(x, ceiling(seq_along(x) / 5))
        m1 <- sapply(.SD[chunks[[1]]], stats::median)
        res <- TRUE
        i <- 1
        while(res){
          i <- i + 1
          if(i > length(chunks)) break
          m2 <- sapply(.SD[chunks[[i]]], stats::median)
          res <- cosine.similarity(m1, m2) > 0.985
        }
        i <- max(chunks[[i-1]])
        select.events <- vector(length = .N)
        select.events[1:i] <- TRUE
        select.events
      }else{
        select.events <- !vector(length = .N)
      }
      ##
      select.events
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ##
  invisible(spectral.events)
}

.reference.group.spectra <- function(spectral.events){
  ##
  cols.by <- spectral.events[, names(.SD), .SDcols = is.factor]
  cols.detector <- spectral.events[, names(.SD), .SDcols = is.numeric]
  mdat <- attr(spectral.events, "mdat")
  ## derive medians
  medians.reference <- spectral.events[
    ,
    j = lapply(.SD, stats::median),
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## subtract unstained/AF; tissue.type/population-specific
  medians.reference <- medians.reference[
    ,
    j = (cols.detector) := lapply(.SD, function(j){
      af <- j[N == 'AF']
      j <- j - af
      j[N == 'AF'] <- af
      j
    }),
    .SDcols = cols.detector,
    by = population
  ]
  ## normalize [0,1] medians -- spectra
  spectra <- data.table::copy(medians.reference)[
    ,
    j = (cols.detector) := {
      j <- .SD / max(.SD)
      j[j < 0] <- 0
      j
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## AF in last position
  spectra[,ord := seq(.N)]
  spectra[grepl('AF', N), ord := spectra[, .N] + seq(.N)]
  data.table::setorder(spectra, ord)[, ord := NULL]
  ## alias column used during unmixing/merging; depending on naming convention
  spectra[, alias := trimws(paste(ifelse(is.na(S), "", as.character(S)), N))]
  spectra[, alias := factor(alias, levels = unique(alias))]
  ##
  data.table::setattr(spectra, 'mdat', mdat)
  ## return
  invisible(spectra)
}

#' Title
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
#' @param name.fix Named character vector -- default `NULL`; if defined, vector names should match to `sample.id`(s) and their respective values name replacements.
#' * e.g., `c('cd8alpha BV786 (Cells)' = 'CD8a BV786')`.
#' @param population.marker Named character vector -- default `NULL`; `population.marker` must be defined for the function to be successful (for cellular controls). The named character vector should take the following form:
#' \itemize{
#'   \item `c(population.name1 = 'sample.name1',population.name2 = 'sample.name2',...)`, where `population.name(s)` are CD/lineage cell types and `'sample.name(s)'` are CD/lineage markers used to stain those respective cell types.
#'   \itemize{
#'     \item e.g., `c(lymphocytes = 'CD3 BV510 (Cells)', monocytes = 'CD14 SB550 (Cells)'`.
#'   }
#' }
#' @param top.N.percent Numeric -- default `10`; defines the percentage of top expressing events (peak detector) used in conjunction with `population.marker` for detecting and assigning the named populations.
#' @param scatter.peak.cut Numeric -- default `0.25`; defines the value at which scatter peak heights will be cut for fine tuning the scatter profile of selected populations (as defined through `population.marker`).
#' @param top.n Numeric -- default `200`; the number of maximally expressing events (sorted vector) used to define 'spectral events' -- events used to define detector medians.
#' @param top.n.override Named numeric vector -- default `NULL`; if defined, the supplied numeric (value) will override `top.n` on a sample-specific (name) basis. See example.
#' @param output.dirs Character vector/an output directory -- default `derived`; the default will derive output directories based on the location of `raw.reference.group` file(s).
#' @param plot Logical -- default `FALSE`; plots diagnostic/QC output for evaluating function performance.
#'
#' @returns A `flowstate` containing spectrally-associated bead/cellular events.
#' @export
#'
#' @examples
#' # placeholder
reference.group.spectra <- function(
    raw.reference.group,
    name.fix = NULL,
    population.marker = NULL,
    top.N.percent = 10,
    scatter.peak.cut = 0.25,
    top.n = 200,
    top.n.override = NULL,
    output.dirs = "derived",
    plot = F
)
{
  ## alias for raw.reference.group
  dir.input <- raw.reference.group
  ## prepare/create output directories
  if(output.dirs == "derived" & all(grepl("data_source", dir.input))){
    res <- strsplit(dir.input, "/")[[1]]
    exp.root <- grep("data_source", res) + 1
    exp.name <- res[exp.root]
    dir.exp <- Reduce(file.path, res[1:exp.root])
    dirs.output <- sapply(c("modified", "results"), function(i){
      sub("source", i, dir.exp)
    })
    invisible(sapply(dirs.output, function(i){
      if(!dir.exists(i)) dir.create(i, recursive = T)
    }))
  }
  ## accepts a directory (containing .fcs files) or accepts .fcs file paths
  if(!all(grepl(".fcs", dir.input)) & length(dir.input) == 1){
    ## get paths to raw .fcs files if a directory;
    ## reference group controls
    ref.paths <- list.files(dir.input, full.names = T)
  }else if(all(grepl(".fcs", dir.input)) & length(dir.input) != 1){
    ref.paths <- dir.input
  }
  ## test for the presence of keyword-value pairs: PnTYPE/Raw_Fluorescence
  raw.fluorescence.check(ref.paths)
  ## get keyword identifier for use in defining 'sample.id' argument
  .sample.id <- cytometer.identifier(ref.paths)
  ## test 'population.marker' argument; need to stop the function here if not defined
  if(is.null(population.marker) & any(grepl("cells", ref.paths, ignore.case = T))){
    sample.ids <- check.keyword(ref.paths, keyword = .sample.id)
    sample.ids <- grep("cells", sample.ids, ignore.case = T, value = T)
    sample.ids <- paste0(sample.ids, collapse = " ; ")
    stop(
      paste(
        "For the effective processing of cellular controls, please define",
        "'population.marker' as documented in '?flowstate::reference.group.spectra'.",
        "Choose from the following:",
        sample.ids,
        sep = "\n"
      )
    )
  }
  ## read and concatenate all raw reference controls;
  ## cellular-based single stained controls (Cells);
  ## and/or bead-based single stained controls (Beads)
  ref <- read.flowstate(
    fcs.file.paths = ref.paths,
    colnames.type = 'N',
    sample.id = .sample.id,
    concatenate = T
  )
  ## update based on name.fix
  if(!is.null(name.fix)){
    for(i in names(name.fix)){
      data.table::set(
        x = ref$keywords,
        i = ref$keywords[, .I[sample.id == i]],
        j = 'sample.id',
        value = name.fix[[i]]
      )
    }
    add.keywords.to.data(ref, 'sample.id')
  }
  ## add metadata to [['keywords']] and [['data']] -- reference group-specific
  reference.group.keywords(ref)
  ## add logical for selecting against saturating events; subset
  select.nonsaturating(ref)
  ref <- subset(ref, select.nonsaturating)
  ref$data[, select.nonsaturating := NULL]
  ## add population (factor) to [['data']]
  select.scatter.population(
    flowstate = ref,
    population.marker = population.marker,
    top.N.percent = top.N.percent,
    scatter.peak.cut = scatter.peak.cut,
    plot = plot,
    plot.output.dir = dirs.output[['results']]
  )
  if(plot){
    p <- plot_populations(ref,pattern.scatter = "AH")
    for(type in names(p)){
      grDevices::png(
        filename = file.path(
          dirs.output[['results']],
          sprintf("reference_group_scatter_populations_NxN_%s_%s_%s.png",
                  exp.name,type,Sys.Date())
        ),
        width = 8, height = 8, units = "in", res = 600
      )
      print(p[[type]])
      grDevices::dev.off()
    }
  }
  ## subset to retain only selected population(s)
  ref <- subset(ref,!is.na(population))
  ## reference group spectral events
  spectral.events <- .reference.group.spectral.events(
    flowstate = ref,
    top.n = top.n,
    top.n.override = top.n.override
  )
  ## generate and plot spectral ribbons
  ## intermediate plots for diagnostics/QC
  if(plot){
    p.ribbons <- plot_spectral.ribbons(spectral.events)
    file.out <- file.path(
      dirs.output[['results']],
      sprintf("reference_group_spectral_ribbons_%s_%s.pdf",
              exp.name,Sys.Date())
    )
    grDevices::pdf(file.out,width = 10, height = 6)
    for(i in p.ribbons[,seq(.N)]){
      suppressWarnings(print(p.ribbons[,ribbon[[i]]]))
    }
    grDevices::dev.off()
  }
  ## reference group spectra
  spectra <- .reference.group.spectra(spectral.events)
  ## save the output
  file.out <- file.path(
    dirs.output[['modified']],
    sprintf("reference_group_spectra_%s_%s.rds",exp.name,Sys.Date())
  )
  saveRDS(spectra,file.out)
  data.table::fwrite(spectra,sub(".rds",".csv",file.out))
  ##
  if(plot){
    ## spectral traces plot list
    p.traces <- plot_spectral.trace(spectra)
    ## save the plots
    file.out <- file.path(
      dirs.output[['results']],
      sprintf(
        "reference_group_spectral_traces_%s_%s.pdf",exp.name,Sys.Date())
    )
    grDevices::pdf(file.out,width = 10,height = 6)
    for(i in p.traces[,seq(.N)]){
      print(p.traces[,trace[[i]]])
    }
    grDevices::dev.off()
  }
  ##
  invisible(spectra)
}
##
plot_populations <- function(flowstate,pattern.scatter = c('A','AH'),bins = 300, sample.n = 5E4){
  cols.scatter <- sort(
    grep(
      sprintf("[FS]SC.*_[%s]$",match.arg(pattern.scatter)),
      names(flowstate$data),
      value = T
    )
  )
  silhouette <- flowstate$data[
    i = is.na(population),
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n)}]),
    .SDcols = cols.scatter,
    by = type
  ]
  foreground <- flowstate$data[
    i = !is.na(population),
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n,replace = T)}]),
    .SDcols = cols.scatter,
    by=.(type,population)
  ]
  plot.body <- ggplot2::ggplot(
    data = NULL,
    mapping = ggplot2::aes(value.x, value.y)
  )
  ##
  p.list <- sapply(flowstate$data[,levels(type)],function(.type){
    p <- plot.body + ggplot2::geom_hex(
      data = silhouette[type == .type],
      bins = bins,
      fill = "darkgray"
    )
    p <- p + ggplot2::geom_hex(
      data = foreground[type == .type],
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
        "Scatter-bound Populations: %s",
        foreground[type == .type,paste0(unique(population),collapse = "; ")]
      )
    )
    p
  },simplify = F)
  ##
  return(p.list)
}

plot_spectral.ribbons <- function(spectral.events){
  ##
  cols.by <- spectral.events[, names(.SD) ,.SDcols = is.factor]
  cols.detector <- spectral.events[, names(.SD), .SDcols = is.numeric]
  mdat <- attr(spectral.events, "mdat")
  mdat <- paste(paste0(names(mdat), ":"), mdat, collapse = " ; ")
  ##
  lims <- spectral.events[, j = range(sapply(.SD, range)), .SDcols = cols.detector]
  lims[1] <- floor(lims[1]); lims[2] <- ceiling(lims[2])
  lims <- asinh(lims/5000)
  p.ribbons <- spectral.events[
    ,
    j = .(ribbon = {
      ribbon <- data.table::melt(
        data = asinh(.SD/5000),
        measure.vars = cols.detector
      )
      p <- ggplot2::ggplot(
        data = ribbon,
        mapping = ggplot2::aes(variable, value)
      ) +
        ggplot2::geom_bin_2d(
          bins = list(length(cols.detector), .N),
          boundary = 0.5,
          show.legend = F,
          na.rm = T
        ) +
        viridis::scale_fill_viridis(
          option = "viridis",
          limits = c(0,.N/10),
          oob = scales::squish
        ) +
        ggplot2::ylim(lims) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          )
        ) +
        ggplot2::labs(
          title = paste0(
            "Spectral Signature: Ribbons",
            paste0(rep(" ", 10), collapse = ""),
            .BY$sample.id
          ),
          subtitle = paste(
            paste("Fluorophore:", .BY$N),
            paste("Marker:", .BY$S),
            sprintf("Type: %s     Population: %s", .BY$tissue.type, .BY$population),
            sep = "\n"
          ),
          caption = paste(
            paste("Peak Detector:", .BY$detector),
            sprintf("N = %d binned 'spectral events'", .N),
            mdat,
            sep = "\n"
          ),
          x = "Detector",
          y = "Expression"
        )
      list(p)
    }),
    .SDcols = cols.detector,
    by = cols.by
  ]
}

plot_spectral.trace <- function(
    spectra,
    plot.type = c("ggplot", "plotly")
)
{
  ##
  cols.by <- spectra[, names(.SD) ,.SDcols = is.factor]
  cols.detector <- spectra[, names(.SD), .SDcols = is.numeric]
  mdat <- attr(spectra, "mdat")
  mdat <- paste(paste0(names(mdat), ":"), mdat, collapse = " ; ")
  ##
  pt <- match.arg(plot.type)
  ##
  switch(
    match.arg(plot.type),
    ggplot = {
      spectral.traces <- spectra[
        ,
        j = .(trace = {
          trace <- data.table::melt(
            data = .SD,
            measure.vars = cols.detector
          )
          p <- ggplot2::ggplot(
            data = trace,
            mapping = ggplot2::aes(variable, value)
          ) +
            ggplot2::geom_line(linewidth=0.5, group=1) +
            # ggplot2::geom_point(size = 1) +
            ggplot2::geom_vline(
              xintercept = as.character(.BY$detector),
              linetype = 'dashed'
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(
                color = "black",
                angle = 90,
                vjust = 0.5,
                hjust = 1
              )
            ) +
            ggplot2::labs(
              title = paste0(
                "Spectral Signature: Trace",
                paste0(rep(" ",10), collapse = ""),
                .BY$sample.id
              ),
              subtitle = paste(
                paste("Fluorophore:", .BY$N),
                paste("Marker:", .BY$S),
                sprintf("Type: %s     Population: %s", .BY$tissue.type, .BY$population),
                sep = "\n"
              ),
              x = "Detector",
              y = "Emission (Normalized [0,1])",
              caption = paste(
                paste("Peak Detector:", .BY$detector),
                mdat,
                sep = "\n"
              )
            )
          list(p)
        }),
        .SDcols = cols.detector,
        by = cols.by
      ]
    },
    plotly = {
      lapply(split(spectra[N != "AF"], by = 'laser'), function(laser){
        traces <- data.table::melt(
          data = laser,
          measure.vars = cols.detector,
          variable.name = "Detector"
        )
        n <- traces[, length(unique(alias))]
        .laser <- laser[, as.character(unique(laser))]
        ##
        p <- plotly::plot_ly(
          data = droplevels(traces),
          x = ~Detector,
          y = ~value,
          color = ~alias,
          colors = grDevices::hcl.colors(n = n, palette = "Purple-Green"),
          type = "scatter",
          mode = "lines"
        )
        ##
        p <- plotly::layout(p, xaxis = list(type = 'category'))
        ##
        p <- plotly::layout(
          p,
          title = sprintf("Laser: %s", .laser),
          xaxis = list(tickmode = 'linear', dtick = 1, tickangle = 270),
          yaxis = list(title = "Emission (Normalized [0,1])")
        )
        ##
        return(p)
      })
    }
  )
}
