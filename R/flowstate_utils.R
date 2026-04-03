fcs.text.primary.required.keywords <-
  c(
    '$BEGINANALYSIS',
    '$BEGINDATA',
    '$BEGINSTEXT',
    '$BYTEORD',
    '$DATATYPE',
    '$ENDANALYSIS',
    '$ENDDATA',
    '$ENDSTEXT',
    '$MODE',
    '$NEXTDATA'
  )

flowstate.parameter.keywords <-
  c(
    'B',
    'DETECTOR',
    'DISPLAY',
    'E',
    'N',
    'R',
    'S',
    'TYPE',
    'V'
  )

flowstate.transform.types <-
  c(
    "Raw_Fluorescence",
    "Unmixed_Fluorescence",
    "Ion_Count"
  )
cytometer.identifier.types <- c(
  "Aurora" = "TUBENAME",
  "ID7000" = "$CELLS",
  "DVS|FLUIDIGM|CYTOF" = "$FIL"
)

cytometer.identifier <- function(fcs.file.paths){
  cyt <- check.keyword(fcs.file.paths,keyword = '$CYT')
  i <- which(sapply(names(cytometer.identifier.types),function(i){grepl(i,cyt)}))
  if(length(i) != 1){
    message("Cytometer keyword 'identifier' could not be resolved;\nis the cytometer indexed in 'flowstate:::cytometer.identifer.types()'?")
    id <- '$FIL'
  }else{
    id <- cytometer.identifier.types[[i]]
  }
}

alias_dt <- function(flowstate){
  alias.vec <- flowstate$data[, unlist(sapply(.SD, attr, which = 'N'))]
  data.table::data.table(
    N = alias.vec,
    alias = names(alias.vec)
  )
}

j.match.keyword.to.data.sample.id <- function(flowstate.object){
  flowstate.object$keywords[
    ,
    sapply(.SD,function(j){
      sum(j %in% flowstate.object$data[,levels(sample.id)],na.rm = TRUE)
    })
  ] |> which.max() |> names()
}

#' @title Add `[['keyword']]` values to `[['data']]`
#' @description
#' Used for efficient addition of factored keyword values to `[['data']]`; updates by reference.
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param keywords.to.add Character vector; keyword names in `[['keywords']]` whose factored values will be added to `[['data']]`.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'   \item `[['data']]`; factored `keywords.to.add` are added as additional columns.
#'   \item `[['keywords']]`; `keywords.to.add` are converted to class `factor`.
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstates; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #create a new set of keywords/values
#' fs$keywords[
#' ,
#' j = c('block.id', 'block.aliquot') := data.table::tstrsplit(
#' sample.id, "_", keep = 4:5)
#' ]
#' #new columns; character class
#' fs$keywords[, .(block.id, block.aliquot)]
#'
#' #add the factored keyword values to fs[['data']]
#' add.keywords.to.data(fs, c('block.id', 'block.aliquot'))
#'
#' fs$data[, .(block.id, block.aliquot)]
#'
#' #add additional keywords to [['data']]
#' add.keywords.to.data(fs, c('$DATE', '$PROJ'))
#'
#' fs$data[, .(`$DATE`, `$PROJ`)]
#' fs$keywords[, .(`$DATE`, `$PROJ`)]
#'
add.keywords.to.data <- function(flowstate, keywords.to.add){
  ## factor keywords.to.add
  for(j in keywords.to.add){
    data.table::set(
      x = flowstate$keywords,
      j = j,
      value = factor(flowstate$keywords[[j]])
    )
  }
  ## keywords.to.add to [['data']]; index by sample.id (levels)
  for(i in flowstate$keywords[, levels(sample.id)]){
    data.table::set(
      x = flowstate$data,
      i = flowstate$data[, .I[sample.id == i]],
      j = keywords.to.add,
      value = flowstate$keywords[
        i = sample.id == i,
        j = ..keywords.to.add
      ]
    )
  }
  ##
  invisible(flowstate)
}
##
# mass.cytometry.detect <- function(flowstate.object,sample.val = 1E3,threshold.val = 0.25){
#   zeroes <- flowstate.object$data[
#     i = sample(.N,sample.val),
#     j = sapply(.SD,function(j){length(which(j==0))})
#   ]
#   res <- length(which(zeroes>0))/length(zeroes)
#   isTRUE(res>threshold.val)
# }

check.keyword <- function(fcs.file.paths, keyword = NULL, value = NULL){
  kw <- sapply(fcs.file.paths, function(i){
    readFCStext(i)[[keyword]]
  })
  if(!is.null(value)){
    all(grepl(value, kw))
  }else{
    return(unique(kw))
  }
}

#' @title Create a column (logical) for selecting against saturating events
#' @description
#' A new column (logical) named `select.nonsaturating` is created in `[['data']]` and any/all events that are at detector limits -- as defined by the '$PnR' value in `[['parameters']]` -- will be marked as `FALSE` with non-saturating events marked as `TRUE`.
#'
#' Detection of saturating events requires linear (non-transformed) values and as such should be performed before [flowstate.transform].
#'
#' @param flowstate A flowstate as returned from [read.flowstate].
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; adds a column (logical) named `select.nonsaturating`
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#'
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstates; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #UPDATES BY REFERENCE -- adds a new column named 'select.nonsaturating'
#' select_nonsaturating(fs)
#' fs$data[, .N, by = select.nonsaturating]
#'
#' #visualize
#' plot(fs,FSC_A,SSC_A) + ggplot2::facet_wrap(~select.nonsaturating)
#'
#' #subset to retain only non-saturating events
#' fs <- subset(fs, select.nonsaturating)
#' fs$data[, .N, by = select.nonsaturating]
#'
#' #after the subset, the column is now redundant
#' all(fs$data[['select.nonsaturating']])
#'
#' #NULL it out
#' fs$data[,'select.nonsaturating' := NULL]
#'
select_nonsaturating <- function(flowstate){
  ## scatter and raw/unmixed fluorescence (linear values);
  ## '$PnR/n1' value to get detector upper limit
  type.string <- c(
    sprintf("%s_Scatter",c("Forward","Side")),
    sprintf("%s_Fluorescence",c("Raw","Unmixed"))
  )
  ##
  ranges <- flowstate$parameters[
    i = TYPE %in% type.string,
    j = .(N, R)
  ][, R := as.numeric(R)]
  ##
  ranges <- merge(ranges, alias_dt(flowstate), sort = F)
  ## initialize a logical
  flowstate$data[, select.nonsaturating := TRUE]
  ## loop through and set to FALSE any event in any j that is saturating
  for(j in ranges[['alias']]){
    data.table::set(
      x = flowstate$data,
      i = flowstate$data[, .I[.SD[[j]] >= ranges[alias == j, R]]],
      j = 'select.nonsaturating',
      value = FALSE
    )
  }
  ##
  invisible(flowstate)
}

#' @title Create a column (logical) for selecting against events below/above a quantile range
#'
#' @inheritParams stats::quantile
#' @param flowstate A flowstate as returned from [read.flowstate].
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; adds a column (logical) named `select.quantile`
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#'
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstates; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #transform
#' flowstate.transform(fs)
#'
#' #UPDATES BY REFERENCE -- adds a new column named 'select.quantile'
#' select_quantile(fs, probs = c(0.0005, 0.9995))
#' fs$data[, .N, by = select.quantile]
#'
#' #visualize
#' plot(fs, CD3, CD8) + ggplot2::facet_wrap( ~ select.quantile)
#'
#' #subset to retain only events in quantile range
#' fs <- subset(fs, select.quantile)
#' fs$data[, .N, by = select.quantile]
#'
#' #after the subset, the column is now redundant
#' all(fs$data[['select.quantile']])
#'
#' #NULL it out
#' fs$data[,'select.quantile' := NULL]
#'
select_quantile <- function(flowstate, probs = c(0.0005, 0.9995)){
  ## scatter and raw/unmixed fluorescence
  type.string <- c(
    sprintf("%s_Scatter", c("Forward", "Side")),
    sprintf("%s_Fluorescence", c("Raw", "Unmixed"))
  )
  ##
  alias <- merge(
    x = flowstate$parameters[
      i = TYPE %in% type.string,
      j = .(N)
    ],
    y = alias_dt(flowstate),
    sort = F
  )[['alias']]
  ## initialize a logical
  flowstate$data[, select.quantile := TRUE]
  ## loop through and set to FALSE any event in any j that exceeds quantile probs
  for(j in alias){
    data.table::set(
      x = flowstate$data,
      i = flowstate$data[
        ,
        j = .I[!data.table::`%between%`(
          .SD[[j]],
          stats::quantile(.SD[[j]], probs = probs)
          )]
      ],
      j = 'select.quantile',
      value = FALSE
    )
  }
  ##
  invisible(flowstate)
}

cosine.similarity <- function(x,y){crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))}

kde.contour <- function(x, y, bandwidth.adjust = NULL, grid.points = 100, prob = 0.95){
  ## pre-calculate bandwidth
  h <- sapply(list(x,y), function(j) MASS::bandwidth.nrd(j))
  ## adjust bandwidth
  if(!is.null(bandwidth.adjust)){
    h <- h * bandwidth.adjust
  }
  ## kernel density estimate
  dens <- MASS::kde2d(x = x, y = y, h = h, n = grid.points)
  ## major contours
  levels <- stats::quantile(dens$z, probs = prob)
  ## major contours -- lines
  cl <- grDevices::contourLines(dens, levels = levels)
  ## major contour -- line
  cl <- cl[[which.max(sapply(cl, function(i) length(i$x)))]]
}
## SpectroFlo Reference Group-specific function
select_scatter.population.contour <- function(
    flowstate,
    population.marker,
    top.N.percent = 10,
    bandwidth.adjust = 2,
    grid.points = 100,
    prob = 0.975
){
  ## alias
  alias <- merge(
    x = flowstate$parameters[, j = .(N, TYPE)],
    y = alias_dt(flowstate),
    sort = F
  )
  ## relevant variables
  #cols.scatter <- sort(alias[grep("Scatter", TYPE)][['alias']])
  cols.scatter <- c('FSC_A', 'SSC_A')
  cols.detector <- alias[TYPE == "Raw_Fluorescence"][['alias']]
  cols.by <- flowstate$data[, names(.SD), .SDcols = is.factor]
  proj <- factor(unique(flowstate$keywords[['$PROJ']]))
  ## if there are any bead controls, add to population marker;
  ## any single fluor bead control should work...
  if(flowstate$keywords[, "Beads" %in% levels(tissue.type)]){
    ## take the first non-Unstained fluor bead control
    bead.ctrl <- flowstate$keywords[
      i = tissue.type == "Beads" & !grepl("Unstained", sample.id)
    ][1, as.character(sample.id)]
    population.marker <- c(population.marker, beads = bead.ctrl)
  }
  ##
  contour.bound <- flowstate$data[
    i = sample.id %in% population.marker,
    j = {
      ## mean of sorted vector to find peak detector
      means.detector <- sapply(cols.detector,function(j){
        mean(sort(.SD[[j]],decreasing = T)[1:100])
      })
      ## peak detector
      detector <- names(which.max(means.detector))
      ## top.N expressing events
      top.N <- ceiling(.N*(top.N.percent/100))
      ## index of ordered peak detector vector -- top expressing events
      i <- order(.SD[[detector]],decreasing = T)[1:top.N]
      ## the index allows a more discrete identification of scatter-based contour/landmark
      cl <- kde.contour(
        x = FSC_A[i],
        y = SSC_A[i],
        bandwidth.adjust = bandwidth.adjust,
        grid.points = grid.points,
        prob = prob
      )
      ## return population label(s) and major contour line (x, y)
      c(
        population =  names(which(population.marker == .BY$sample.id)),
        cl[c('x', 'y')]
      )
    },
    .SDcols = c(cols.detector,cols.scatter),
    by = cols.by
  ][,population := factor(population)][]
  ## apply tissue.type/population-specific contour bounds to all relevant data (rows)
  ## adds a population label and a logical for selection
  for(.population in contour.bound[, levels(population)]){
    .type <- ifelse(.population == 'beads', "Beads", "Cells")
    flowstate$data[
      i = tissue.type == .type,
      j = c('population', 'select.contour') := {
        cl <- sapply(c('x', 'y'), function(i){
          contour.bound[population == .population][[i]]}, simplify = F)
        res <- as.logical(sp::point.in.polygon(FSC_A, SSC_A, cl$x, cl$y))
        list(.population, res)
      }
    ][]
  }
  ##
  invisible(flowstate)
}
## SpectroFlo Reference Group-specific function
select_contour.singlets <- function(
    flowstate,
    scatter.singlets = c('FSC', 'SSC'),
    n.adjust = 0.2,
    prob = 0.90
)
{
  pulses <- paste(match.arg(scatter.singlets), c('H', 'W'), sep = "_")
  flowstate$data[
    ,
    j = select.contour := {
      .i <- sample(.N,floor(.N*n.adjust))
      cl <- kde.contour(
        x = .SD[[pulses[1]]][.i],
        y = .SD[[pulses[2]]][.i],
        prob = prob
      )
      in_contour <- as.logical(sp::point.in.polygon(.SD[[pulses[1]]], .SD[[pulses[2]]], cl$x, cl$y))
    },
    .SDcols = pulses,
    by = population
  ]
  ##
  invisible(flowstate)
}
