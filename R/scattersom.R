scattersom <- function(flowstate, population, bandwidth.adjust, plot.contour = F, seed.val = 1337){
  ## NULL out if derivative columns exist
  if('population' %in% names(flowstate$data)) flowstate$data[, population := NULL]
  if('select.population' %in% names(flowstate$data)) flowstate$data[, select.population := NULL]
  detector.names <- flowstate$parameters[TYPE == 'Raw_Fluorescence', N]
  cols.detector <- names(flowstate$data)[flowstate$data[, sapply(.SD, attr, which = "N") %in% detector.names]]
  ## defined through function argument:
  ## population label and the sample to use to define and select a landmark/scatter profile

  ## subset [['data']] (concatenate), indexing only those events associated with the single sample as defined in 'population'
  ref.sub <- flowstate$data[sample.id == population]
  ## data-drive method to find peak detector
  detector <- ref.sub[
    ,
    j = {
      means.detector <- sapply(cols.detector, function(j){
        mean(sort(.SD[[j]], decreasing = T)[1:100])
      })
      detector <- names(which.max(means.detector))
    }
  ]
  ## transform that single detector
  ref.sub[, (detector) := asinh(.SD[[detector]]/1000)]
  ## dimensions for FlowSOM: all scatter pulses + the transformed peak detector
  cols.scatter <- sort(grep("[FS]SC", names(ref.sub), value = T))
  dims <- c(cols.scatter, detector)
  ## a scaling function to scale [['data']];
  ## a mix of linear scatter + transformed fluor so needs scaling
  scale.func <- function(x) {
    (x - mean(x))/stats::sd(x)
  }
  ## FlowSOM on all scatter pulses + the transformed peak detector;
  ## use a default 10 x 10 grid
  ## map the result to get per-event node assignments
  fsom <- flowsomlite::SOM(
    data = as.matrix(ref.sub[, lapply(.SD, scale.func), .SDcols = dims]),
    xdim = 10,
    ydim = 10,
    map = TRUE,
    seed.val = seed.val
  )
  ## UPDATE BY REFERENCE -- add map to [['data']]
  ref.sub[
    ,
    j = c('node', 'node.dist') := as.list(as.data.frame(fsom$mapping))
  ]
  ## cluster
  d <- stats::dist(fsom$codes, method = "euclidean")
  hc <- stats::hclust(d, method = "average")
  cl <- stats::cutree(hc, k = 15)
  ##  UPDATE BY REFERENCE -- add/map cluster result to [['data']]
  ref.sub[, cluster := cl[node]]
  ## data-driven detection of optimal/representative cluster:
  ## max mean and max count for peak detector
  ## cluster-of-interest -- use data.table chains
  coi <-(
    ref.sub[,.(m = mean(.SD[[detector]]), .N), by = cluster]
    [order(-m)]
    [m > stats::median(m), as.numeric(cluster[which.max(N)])]
  )
  ## subset based on coi -- constrains the local scatter space to make contour fitting more robust;
  ## trying to eliminate user-defined parameters so this can be automated
  ## return the derived contour bounds
  bounds <- contour_pulses(ref.sub[cluster == coi], bandwidth.adjust = bandwidth.adjust, plot = plot.contour)
  ## apply the bounds to the entire concatenate [['data']]
  for(i in seq(bounds)){
    pulse.pair <- attr(bounds[[i]], which = 'pulse.pair')
    if(i == 1){
      flowstate$data[
        ,
        j = select.population := {
          as.logical(sp::point.in.polygon(ji, jj, bounds[[i]]$x, bounds[[i]]$y))
        },
        env = list(ji = pulse.pair[1], jj = pulse.pair[2])
      ]
    }else{
      flowstate$data[
        i = (select.population),
        j = select.population := {
          as.logical(sp::point.in.polygon(ji, jj, bounds[[i]]$x, bounds[[i]]$y))
        },
        env = list(ji = pulse.pair[1], jj = pulse.pair[2])
      ]
    }
  }
  ## apply the population-specific label
  flowstate$data[
    i = (select.population),
    population := factor(names(population))
  ]
  ## NULL the place-holder logical
  flowstate$data[, select.population := NULL][]
  ## return
  invisible(flowstate)
}
##
