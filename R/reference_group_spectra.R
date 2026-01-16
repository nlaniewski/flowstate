## function to add metadata to [['keywords]] and [['data']]
reference.group.keywords <- function(flowstate){
  ## column identifier for [['data']] and [['keywords']]
  k.match <- j.match.keyword.to.data.sample.id(flowstate)
  ## add keyword metadata based on splitting sample.id;
  ## following SpectroFlo naming convention: marker(S) fluorophore(N) (type) -- literal spaces separating
  keywords.to.add <- c('type','N','S')
  flowstate$keywords[
    ,
    j = sample.id := factor(k),
    env = list(k = k.match)
  ]
  ## test
  if(any(flowstate$data[,unique(sample.id)]!=flowstate$keywords[,unique(sample.id)])){
    stop("Sample order does not match between [['data']] and [['keywords']].")
  }
  ## add keyword/value pairs based on gsub/regex; depending on SpectrolFlo naming convention
  ## updates [['keywords']]
  flowstate$keywords[
    ,
    j = (keywords.to.add) := {
      type <- factor(gsub("^.*\\((.*?)\\).*$", "\\1", sample.id))
      res <- sub(" \\(.*$","",sample.id)
      marker <- factor(ifelse(type == "Cells" & res == "Unstained",NA,sub(" +.*$","",res)))
      fluorophore <- factor(ifelse(type == "Cells" & is.na(marker),"AF",sub(".*? ", "",res)))
      list(type,N = fluorophore,S=marker)
    }
  ]
  ## updates [['data']]
  totals <- flowstate$data[, .N, by = sample.id]
  for(j in keywords.to.add){
    data.table::set(
      x = flowstate$data,
      i = NULL,
      j = j,
      value = rep(flowstate$keywords[[j]],totals[["N"]])
    )
  }
  ## return
  invisible(flowstate)
}
## function to add population (factor) to [['data']] for selecting scatter-specific populations;
## scatter bounds for populations are defined by using 'landmark' markers
select.scatter.population <- function(
    flowstate,
    population.marker,
    scatter.cut = list(
      FSC = c("A" = 0.25,"H" = 0.05,"W" = 0.05),
      SSC = c("A" = 0.25,"H" = 0.05,"W" = 0.05)
    ),
    plot = F
)
{
  ## scatter.cut as a named vector
  scatter.cut <- unlist(lapply(names(scatter.cut),function(i){
    stats::setNames(scatter.cut[[i]],nm = paste(i,names(scatter.cut[[i]]),sep = "_"))
  }))
  ## drop any pulses not in [['data']]
  scatter.cut <- scatter.cut[names(scatter.cut) %in% names(flowstate$data)]
  ## relevant variables
  j.match <- j.match.parameters.to.data(flowstate)
  cols.detector <- (flowstate$parameters
                    [TYPE == "Raw_Fluorescence"][[j.match]])
  cols.by <- flowstate$data[,names(.SD),.SDcols = is.factor]
  ## sample-specific (population.marker);
  ## use the mean of top-expressing events (sorted vectors) to find peak detector;
  ## using peak detector, return the index of top.N expressing events (ordered vector)
  i.population <- flowstate$data[
    i = sample.id %in% population.marker,
    j = .(i = {
      means.detector <- sapply(.SD,function(j){
        mean(sort(j,decreasing = T)[1:100])
      })
      detector <- names(which.max(means.detector))
      top.N <- ceiling(.N*(1/100))
      i <- .I[order(.SD[[detector]],decreasing = T)[1:top.N]]
    }),
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## sample-specific (population.maker);
  ## indexed using 'i.population'
  ## return scatter, population name, and a 'select.scatter' logical
  scatter.subset <- flowstate$data[
    i = i.population[['i']],
    j = {
      c(population = names(which(population.marker == .BY$sample.id)),
        .SD,
        select.scatter = TRUE
      )
    },
    .SDcols = names(scatter.cut),
    by = cols.by
  ][,population := factor(population)]
  ## sequentially cut scatter peaks;
  ## update 'select.scatter' logical
  for(j in names(scatter.cut)){
    scatter.subset[
      i = (select.scatter),
      j = select.scatter := data.table::`%between%`(
        scatter,
        peak.height.cut(
          scatter,
          height.cut = scatter.cut[[j]],
          plot=plot,main=paste0(j,"\n",.BY$population))
      ),
      by = c(cols.by,'population'),
      env = list(scatter = j)
    ]
  }
  ## return the bounds
  scatter.bounds <- scatter.subset[
    i = (select.scatter),
    j = lapply(.SD,range),
    .SDcols = sort(names(scatter.cut)),
    by = c(cols.by,'population')
  ]
  ## loop through scatter bounds to build population-specific selection
  ## updates by reference -- adds 'population' to [['data']]
  for(.population in scatter.bounds[,levels(population)]){
    for(scatter in scatter.bounds[,names(.SD),.SDcols = is.numeric]){
      if(!exists("vec")){
        vec <- data.table::`%between%`(
          flowstate$data[[scatter]],
          scatter.bounds[population == .population][[scatter]]
        )
      }else{
        vec[vec] <- data.table::`%between%`(
          flowstate$data[[scatter]][vec],
          scatter.bounds[population == .population][[scatter]]
        )
      }
    }
    flowstate$data[vec,population := factor(.population)]
    rm(vec)
  }
  ## return
  invisible(flowstate)
}
## function to calculate reference group spectral medians -- population-specific;
## mean of top expressing events to determine peak detector;
## indexes top expressing events (sorted vector) in peak detector -- 'spectral events';
## uses indices to calculate medians for all detectors;
## matched (population-specific) AF is removed and medians normalized -- 'spectra'
.reference.group.spectra <- function(flowstate)
{
  ## relevant variables
  j.match <- j.match.parameters.to.data(flowstate)
  cols.detector <- (flowstate$parameters
                    [TYPE == "Raw_Fluorescence"][[j.match]])
  proj <- factor(flowstate$parameters[,levels(PROJ)])
  ## unstained/AF means;
  ## used in subsequent step to help resolve dim/AF-impacted fluors through subtraction
  means.af <- flowstate$data[
    i = N == "AF",
    j = .(mean = sapply(.SD,function(j){
      mean(sort(j,decreasing = T)[1:100])
    })),
    .SDcols = cols.detector,
    by = c(flowstate$data[,names(.SD),.SDcols = is.factor])
  ]
  ## use the mean of top-expressing events (sorted vectors) to find peak detector;
  ## for unstained/AF -- include all events to calculate mean;
  ## subtract 'means.af' to resolve dim/AF-impacted fluors (otherwise obscured by AF);
  ## index top expressing events in peak detector (from ordered vector);
  ## medians calculated using the index -- 'spectral events'
  ## include mean of peak detector per sample/per population to select representative population-specific spectra
  medians.reference <- flowstate$data[
    ,
    j = {
      if(.BY$N != "AF"){
        means.detector <- sapply(.SD,function(j){
          mean(sort(j,decreasing = T)[1:100])
        })
        means.detector <- means.detector - means.af[population == .BY$population,mean]
      }else{
        means.detector <- sapply(.SD,stats::median)
      }
      mean.detector <- max(means.detector)
      detector <- names(which.max(means.detector))
      top.n <- order(.SD[[detector]],decreasing = T)[1:100]
      medians.detector <- lapply(.SD[top.n],function(j){
        stats::median(j)
      })
      # c(detector = detector,medians.detector)
      c(mean.detector = mean.detector,detector = detector,medians.detector)
    },
    .SDcols = cols.detector,
    by = c(flowstate$data[,names(.SD),.SDcols = is.factor])
  ]
  ## factor 'detector' column; levels equal to order in [['parameters']]
  ## order by 'detector'
  medians.reference[,detector := factor(detector,levels = cols.detector)]
  data.table::setorder(medians.reference,detector)
  ## subset the medians to retain representative max expressing population-specific spectra
  i.drop <- medians.reference[
    i = N != "AF",
    j = .I[which.min(mean.detector)],
    by = sample.id
  ][['V1']]
  medians.reference <- medians.reference[-i.drop][,mean.detector := NULL]
  ## subtract unstained/AF; population-specific
  medians.reference <- medians.reference[
    ,
    j = (cols.detector) := lapply(.SD,function(j){
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
      j <- .SD/max(.SD)
      j[j<0]<-0
      j
    },
    .SDcols = cols.detector,
    by =c(medians.reference[,names(.SD),.SDcols = is.factor])
  ]
  ## AF in last position
  spectra[,ord := seq(.N)]
  spectra[grepl('AF',N), ord := spectra[,.N] + seq(.N)]
  data.table::setorder(spectra,ord)[,ord := NULL]
  ## alias column used during unmixing/merging; depending on naming convention
  spectra[,alias := trimws(paste(ifelse(is.na(S),"",as.character(S)),N))]
  ## batch/project identifier
  spectra[,PROJ := proj]
  ## return
  invisible(spectra)
}
reference.group.spectra <- function(
    raw.reference.group.directory,
    population.marker,
    output.dirs = "derived",
    plot = F
)
{
  ## alias for raw.reference.group.directory
  dir.input <- raw.reference.group.directory
  ## prepare/create output directories
  if(output.dirs == "derived" & grepl("data_source",dir.input)){
    res <- strsplit(dir.input,"/")[[1]]
    exp.root <- grep("data_source",res)+1
    exp.name <- res[exp.root]
    dir.exp <- Reduce(file.path,res[1:exp.root])
    dirs.output <- sapply(c("modified","results"),function(i){
      sub("source",i,dir.exp)
    })
    sapply(dirs.output,function(i){
      if(!dir.exists(i)) dir.create(i,recursive = T)
    }) |> invisible()
  }
  ## accepts a directory (containing .fcs files) or accepts .fcs file paths
  if(!all(grepl(".fcs",dir.input)) & length(dir.input) == 1){
    ## get paths to raw .fcs files if a directory;
    ## reference group controls
    ref.paths <- list.files(dir.input,full.names = T)
  }else if(all(grepl(".fcs",dir.input)) & length(dir.input) !=1){
    ref.paths <- dir.input
  }
  ## test for the presence of keyword-value pairs: PnTYPE/Raw_Fluorescence
  raw.fluorescence.check(ref.paths)
  ## get keyword identifier for use in defining 'sample.id' argument
  .sample.id <- cytometer.identifier(ref.paths)
  ## test 'population.marker' argument; need to stop the function here if not defined
  if(any(grepl("cells",ref.paths,ignore.case = T))){
    sample.ids <- check.keyword(ref.paths,keyword = .sample.id)
    sample.ids <- grep("cells",sample.ids,ignore.case = T,value = T)
    sample.ids <- paste0(sample.ids,collapse = " ; ")
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
  ## cellular-based (PBMC) single stained controls (Cells);
  ## bead-based (binds Ig) single stained controls (Beads)
  ref <- read.flowstate(
    fcs.file.paths = ref.paths,
    colnames.type = 'N',
    sample.id = sample.id,
    concatenate = T
  )
  ## add metadata to [['keywords]] and [['data']] -- reference group-specific
  reference.group.keywords(ref)
  ## add logical for selecting against saturating events; subset
  select.nonsaturating(ref)
  ref <- subset(ref,select.nonsaturating)
  ref$data[,select.nonsaturating := NULL]
  ## add population (factor) to [['data']]
  select.scatter.population(
    flowstate = ref,
    population.marker = population.marker
  )
  if(plot){
    p <- plot_populations(ref)
    grDevices::png(
      filename = sprintf("reference_group_scatter_populations_%s_%s.png",
                         exp.name,Sys.Date()),
      width = 8, height = 8, units = "in", res = 600
    )
    print(p)
    grDevices::dev.off()
  }
  ## subset to retain only selected population(s)
  ref <- subset(ref,!is.na(population))
  ## reference group spectra;
  ## derived from the median of 'spectral events' (top expressing sorted vectors)
  spectra <- .reference.group.spectra(ref)
  ## save the output
  file.out <- file.path(
    dirs.output[['modified']],
    sprintf("reference_group_spectral_medians_%s_%s.rds",exp.name,Sys.Date())
  )
  saveRDS(spectra,file.out)
  data.table::fwrite(spectra,sub(".rds",".csv",file.out))
}
##
