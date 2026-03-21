#' @title Transform linear values
#' @description
#' Linear values in `[['data']]` are transformed using a defined function.  For full spectrum cytometry, transforming linear fluorescence values using \link{asinh} is preferred, with around-zero compression adjusted using a 'cofactor' of 5,000.
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param j Character vector -- default `NULL`; any/all parameters having a keyword-value pair of `TYPE/Raw|Unmixed_Fluorescence` will be transformed in `[['data']]`.  If a defined character vector: specific columns in `[['data']]` that are to be transformed.
#' @param transform.type Character vector -- default \link{asinh}; quoted function that will be used to transform `j` in `[['data']]`.
#' @param cofactor Numeric -- default 5000; if `transform.type` = \link{asinh} (default), the cofactor will be used to modify the transformation.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; transformed values
#'    \item `flowstate[['parameters']]`; an attribute -- `transformed` -- is added to record which parameters have been transformed.
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read .fcs files as a flowstate object; concatenate
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #plot and mean values of linear columns
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#'
#' #transform
#' flowstate.transform(
#'   fs,
#'   j = c('CD3', 'CD4', 'CD8'),
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['parameters']] attribute
#' #applied transform function and cofactor mapped to $PnN
#' attr(fs$parameters,'transformed')
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#'
#' #transform all/remaining
#' flowstate.transform(
#'   fs,
#'   j = NULL,
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['parameters']] attribute
#' #applied transform function and cofactor mapped to $PnN
#' attr(fs$parameters,'transformed')
#'
flowstate.transform<-function(flowstate, j = NULL, transform.type = "asinh", cofactor = 5000){
  ## merged parameters; alias column; parameters to transform
  pam <- parameter.alias.merged(flowstate)[TYPE %in% flowstate.transform.types]
  ## if j argument is defined, subset
  if(!is.null(j)){
    pam <- pam[alias %in% j]
  }
  j <- pam[['alias']]
  if(!all(j %in% names(flowstate$data))){
    stop(paste(
      paste0(j[!j %in% names(flowstate$data)], collapse = ", "),
      "not found in [['data']]")
    )
  }
  ## has transformation already been applied?
  ## check for the attribute
  if("transformed" %in% names(attributes(flowstate$parameters))){
    ## what has been transformed?
    j.transformed <- attr(flowstate$parameters, 'transformed')
    j.transformed <- pam[N %in% names(j.transformed)][['alias']]
    j <- j[!j %in% j.transformed]
    if(length(j)==0){
      return(
        message("Transformation has already been applied to j!")
      )
    }
  }
  ##
  if(flowstate$data[, length(levels(sample.id)) == 1]){
    message(paste(flowstate$data[, levels(sample.id)], "-->", "transforming..."))
  }else{
    message('flowstate --> transforming...')
  }
  trans.func <- get(transform.type)
  for(.j in j){
    data.table::set(
      x = flowstate$data,
      j = .j,
      value = trans.func(flowstate$data[[.j]] / cofactor)
    )
  }
  ##
  nm <- pam[['N']]
  val <- stats::setNames(paste0(rep(transform.type, length(nm)), ";", cofactor), nm)
  data.table::setattr(
    x = flowstate$parameters,
    name = 'transformed',
    value = val
  )
}
#' @title Transform values to linear -- inverse
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param j Character vector -- default `NULL`; any/all parameters having a stored `transformed` attribute -- from the result of [flowstate.transform] -- will be inverse transformed in `[['data']]` using [base::sinh()].  If a defined character vector: specific columns in `[['data']]` that are to be inverse transformed.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; inverse transformed values -- linear.
#'    \item `flowstate[['parameters']]`; the stored attribute -- `transformed` -- is modified, removing parameters that have been inverse transformed.
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read .fcs files as a flowstate object; concatenate
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #plot and mean values of linear columns
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' res1.linear <- fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD8')]
#' print(res1.linear)
#'
#' #transform
#' flowstate.transform(
#'   fs,
#'   j = c('CD3','CD4','CD8'),
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['parameters']] attribute
#' #applied transform function and cofactor mapped to $PnN
#' 'transformed' %in% names(attributes(fs$parameters))
#' attr(fs$parameters, 'transformed')
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' fs$data[,sapply(.SD,mean),.SDcols = c('CD3', 'CD4', 'CD8')]
#'
#' #inverse transformation; defined `j` argument
#' flowstate.transform.inverse(fs, j = c('CD3', 'CD8'))
#'
#' #updated [['parameters']] attribute
#' #returns only `PerCP-Fire 806-A (CD4)` as all other transformations have been inversed
#' 'transformed' %in% names(attributes(fs$parameters))
#' attr(fs$parameters, 'transformed')
#'
#' #plot and mean values of linear columns
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' res2.linear <- fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD8')]
#' print(res2.linear)
#'
#' #linear --> transformed --> inverse --> linear
#' res1.linear == res2.linear
#'
flowstate.transform.inverse <- function(flowstate, j = NULL){
  if(is.null(attr(flowstate$parameters, 'transformed', exact = TRUE))){
    return(message(
      sprintf(
        "attributes(%s$parameters)[['transformed']] not found;\nno inverse transformation applied.",
        deparse(substitute(flowstate))
      )
    ))
  }
  alias <- attr(flowstate$parameters,'alias')
  ##
  j.transformed <- attr(flowstate$parameters, 'transformed', exact = TRUE)
  j.inverse <- strsplit(j.transformed, ";")
  j.inverse <- lapply(j.inverse, sub, pattern = "asinh", replacement = "sinh")
  names(j.inverse) <- alias[names(j.inverse)]
  if(!is.null(j)){
    i <- names(j.inverse) %in% j
    j.inverse <- j.inverse[i]
  }
  ##
  message('flowstate --> transforming -- inverse...')
  for(.j in names(j.inverse)){
    trans.func <- get(j.inverse[[.j]][1])
    cofactor <- as.numeric(j.inverse[[.j]][2])
    ##
    data.table::set(
      x = flowstate$data,
      j = .j,
      value = trans.func(flowstate$data[[.j]]) * cofactor
    )
  }
  ##
  drop.trans <- names(alias)[alias %in% names(j.inverse)]
  j.transformed <- j.transformed[!names(j.transformed) %in% drop.trans]
  ##
  data.table::setattr(
    x = flowstate$parameters,
    name = 'transformed',
    value = if(length(j.transformed) == 0){
      NULL
    }else{
      j.transformed
    }
  )
  ##
}
