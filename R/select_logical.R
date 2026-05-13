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

#' @title Create a column (logical) for selecting singlet (scatter pulse) events
#' @description
#' A new column (logical) named `select.singlets` is created in `[['data']]` and any/all singlet events will be marked as `TRUE` with doublet events marked as `FALSE`; singlet events are identified through scatter pulse geometry (FSC-A/FSC-H ; SSC-A/SSC-H).
#'
#' @param flowstate A flowstate as returned from [read.flowstate].
#' @param quantiles Numeric vector (length 2) -- default `c(0.85, 0.95)`; respective quantile probabilities used to select singlet events first by `FSC-A vs FSC-H`, then by `SSC-A vs SSC-H`.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; adds a column (logical) named `select.singlets`
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
#' #UPDATES BY REFERENCE -- adds a new column named 'select.singlets'
#' select_singlets(fs, quantiles = c(0.95, 0.95))
#' fs$data[, .N, by = select.singlets]
#'
#' #visualize
#' plot(fs, FSC_A, FSC_H) + ggplot2::facet_wrap( ~ select.singlets)
#' plot(fs, SSC_A, SSC_H) + ggplot2::facet_wrap( ~ select.singlets)
#'
#' #subset to retain only singlet events
#' fs <- subset(fs, select.singlets)
#' fs$data[, .N, by = select.singlets]
#'
#' #after the subset, the column is now redundant
#' all(fs$data[['select.singlets']])
#'
#' #NULL it out
#' fs$data[,'select.singlets' := NULL]
#'
select_singlets <- function(flowstate, quantiles = c(0.85, 0.95)){
  ## remove doublets -- forward scatter geometry; initial selector
  ## selects events (collinear) based on quantile [1]
  flowstate$data[
    ,
    j = select.singlets := {
      res <- FSC_A/FSC_H
      res < stats::quantile(res, probs = quantiles[1])
    }
  ]
  ## remove doublets -- side scatter geometry; indexed against the initial selector
  ## selects events (collinear) based on quantile [2]
  flowstate$data[
    i = (select.singlets),
    j = select.singlets := {
      res <- SSC_A/SSC_H
      res < stats::quantile(res, probs = quantiles[2])
    }
  ]
  ##
  invisible(flowstate)
}
