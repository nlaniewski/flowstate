spillover.update.external <- function(flowstate, spillover.flowjo.path.csv){
  ##read previously exported spillover matrix (.csv) from FlowJo
  spill.mat.external <- utils::read.csv(
    spillover.flowjo.path.csv,
    check.names = FALSE,
    row.names = NULL
  )
  ##drop column 1
  if(is.character(spill.mat.external[,1])){
    spill.mat.external <- spill.mat.external[-1]
  }
  ##split column names using the FlowJo delimiter
  cols.split <- strsplit(colnames(spill.mat.external), " :: ")
  ##which split matches the names of flowstate[['spill']]
  split.index <- which(sapply(seq(unique(sapply(cols.split, length))), function(i){
    all(sapply(cols.split, '[[', i) %in% names(flowstate$spill))
  }))
  ##update the column names
  colnames(spill.mat.external) <- sapply(cols.split, '[[', split.index)
  ##update flowstate[['spill']]
  for(j in names(flowstate$spill)){
    data.table::set(
      x = flowstate$spill,
      j = j,
      value = spill.mat.external[[j]]
    )
  }
}
#' @title Update the values of `flowstate[['spill']]`
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param i Variable (unquoted).
#' @param j Variable (unquoted).
#' @param value Numeric; correction value to be used during compensation.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['spill']]`; updates `[i, j]` with the defined correction value.
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #row index; CD4 vs CD8
#' index <- which(names(fs$spill) %in% c('CD4', 'CD8'))
#' fs$spill[index, .(CD4, CD8)]
#'
#' #update a spill value
#' spillover.update.value(fs, i = CD8, j = CD4, value = 0.03)
#'
#' fs$spill[index, .(CD4, CD8)]
#'
spillover.update.value <- function(flowstate, i, j, value){
  i.char <- deparse(substitute(i))
  j.char <- deparse(substitute(j))
  vec.i <- sapply(c(i.char, j.char), function(x) {
    which(names(flowstate$spill) %in% x)
  })
  for(i in names(vec.i)){
    l <- length(vec.i[[i]])
    if(l != 1){
      stop(
        sprintf(
          "Indexing for %s returned %s results; is %s found/unique in [['spill']]?",
          i, l, i
        )
      )
    }
  }
  ##
  data.table::set(
    x = flowstate$spill,
    i = vec.i[1],
    j = vec.i[2],
    value = value
  )
}

#' @title Compensate `flowstate[['data']]` using values stored in `flowstate[['spill']]`
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param decompensate Logical; if `TRUE`, will return data values to their original state.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; applies compensation to data using values stored in spill.
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
#' #transform
#' flowstate.transform(fs, c('CD4', 'CD8'))
#'
#' #row index; CD4 vs CD8
#' index <- which(names(fs$spill) %in% c('CD4', 'CD8'))
#'
#' #update a spill value
#' spillover.update.value(fs, CD8, CD4, 0.03)
#'
#' #over-compensated data
#' plot(fs, CD4, CD8) + ggplot2::labs(title = "Over-compensated")
#'
#' #apply compensation
#' spillover.apply(fs)
#'
#' #compensated data
#' plot(fs, CD4, CD8) + ggplot2::labs(title = "Compensated")
#'
#' #return data values to original state by decompensating
#' spillover.apply(fs, decompensate = TRUE)
#'
#' plot(fs, CD4, CD8) + ggplot2::labs(title = "Over-compensated")
#'
spillover.apply <- function(flowstate, decompensate = FALSE){
  if(attr(flowstate$spill, 'applied') & isFALSE(decompensate)){
    stop("Spillover has already been applied.")
  }
  spill.index <- sapply(c(1, 2), function(margin) {
    which(apply(flowstate$spill, margin, function(x) {
      any(x != 0 & x != 1 & x != 1E-6)
    }))
  }) |> Reduce(f = union, x = _) |> sort()
  spill.index <- names(flowstate$spill)[spill.index]
  flowstate.transform.inverse(flowstate, j = spill.index)
  spill.mat <- as.matrix(
    flowstate$spill[
      i = which(colnames(flowstate$spill) %in% spill.index),
      j = .SD,
      .SDcols = spill.index
    ]
  )
  ##
  flowstate$data[
    ,
    j = (spill.index) := as.list(
      as.data.frame(
        if(decompensate){
          as.matrix(flowstate$data[, .SD, .SDcols = spill.index]) %*% spill.mat
        }else{
          as.matrix(flowstate$data[, .SD, .SDcols = spill.index]) %*% solve(spill.mat)
        }
      )
    )
  ]
  flowstate.transform(flowstate, j = spill.index)
  ##
  if(decompensate){
    data.table::setattr(flowstate$spill, name = "applied", value = FALSE)
  }else{
    data.table::setattr(flowstate$spill, name = "applied", value = TRUE)
  }
  ##
  invisible(flowstate)
}
##
