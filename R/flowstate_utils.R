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

parameter.alias.merged <- function(flowstate){
  ##
  alias <- attr(flowstate[['parameters']], which = 'alias', exact = TRUE)
  ##
  parameter.alias.merged <- merge(
    flowstate$parameters,
    data.table::data.table(
      N = names(alias),
      alias = alias
    ),
    sort = F
  )
  ##
  invisible(parameter.alias.merged)
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
#' select.nonsaturating(fs)
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
select.nonsaturating <- function(flowstate){
  ## scatter and raw/unmixed fluorescence (linear values);
  ## '$PnR/n1' value to get detector upper limit
  type.string <- c(
    sprintf("%s_Scatter",c("Forward","Side")),
    sprintf("%s_Fluorescence",c("Raw","Unmixed"))
  )
  ranges <- parameter.alias.merged(flowstate)[
    i = TYPE %in% type.string,
    j = .SD,
    .SDcols = c('R','alias')
  ][,R := as.numeric(R)]
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
#' select.quantile(fs, probs = c(0.0005, 0.9995))
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
select.quantile <- function(flowstate, probs = c(0.0005, 0.9995)){
  ## scatter and raw/unmixed fluorescence
  type.string <- c(
    sprintf("%s_Scatter", c("Forward", "Side")),
    sprintf("%s_Fluorescence", c("Raw", "Unmixed"))
  )
  alias <- parameter.alias.merged(flowstate)[TYPE %in% type.string][['alias']]
  ##initialize a logical
  flowstate$data[, select.quantile := TRUE]
  ##loop through and set to FALSE any event in any j that exceeds quantile probs
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
