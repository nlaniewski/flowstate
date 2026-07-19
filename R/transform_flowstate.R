#' @title Transform linear values
#' @description
#' Linear values in `[['data']]` are transformed using a defined function.  For full spectrum cytometry, transforming linear fluorescence values using \link{asinh} is preferred, with around-zero compression adjusted using a 'cofactor' of 5,000.
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param j Character vector -- default `NULL`; any/all parameters having a keyword-value pair of (instrument-specific):
#' \itemize{
#'    \item `Aurora`: `TYPE %in% c("Raw_Fluorescence", "Unmixed_Fluorescence")`
#'    \item `FACSDiscover [AS]8`: `KIND %in% COLOR`
#' } will be transformed in `[['data']]`.  If a defined character vector: specific columns in `[['data']]` that are to be transformed.
#' @param transform.func Character vector -- default \link{asinh}; quoted function that will be used to transform `j` in `[['data']]`.
#' @param cofactor Numeric -- default 5000; if `transform.func` = \link{asinh} (default), the cofactor will be used to modify the transformation.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; transformed values.
#'    \item `flowstate[['data']]`; an attribute -- `transformed` -- is added to each `j`.
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
#'   transform.func = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['data']] attributes; applied transform function and cofactor
#' fs$data[, lapply(.SD, attr, which = 'transformed')]
#' fs$data[, attr(CD3, 'transformed')]
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#'
#' #transform all/remaining
#' flowstate.transform(
#'   fs,
#'   j = NULL,
#'   transform.func = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['data']] attributes; applied transform function and cofactor
#' fs$data[, lapply(.SD, attr, which = 'transformed')]
#'
flowstate.transform<-function(flowstate, j = NULL, transform.func = "asinh", cofactor = 5000){
  ## instrument-specific: parameter keyword encoding 'TYPE|KIND'
  i <- flowstate.transform.keywords %in% names(flowstate$parameters)
  ## parameters to transform -- alias
  alias <- merge(
    x = flowstate$parameters[
      i = kw %in% flowstate.transform.strings,
      j = .(N),
      env = list(kw = flowstate.transform.keywords[i])
    ],
    y = alias_dt(flowstate),
    sort = F
  )[['alias']]
  ## if j argument is defined, subset
  if(!is.null(j)){
    alias <- alias[alias %in% j]
  }
  ##
  j <- alias
  ##
  if(!all(j %in% names(flowstate$data))){
    stop(paste(
      paste0(j[!j %in% names(flowstate$data)], collapse = ", "),
      "not found in [['data']]")
    )
  }
  ## has transformation already been applied?
  ## check for the attribute
  res.transformed <- flowstate$data[, sapply(.SD, function(j)
    'transformed' %in% names(attributes(j)))]
  if(any(res.transformed)){
    ## what has been transformed?
    j.transformed <- names(which(res.transformed))
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
  # trans.func <- get(transform.func)
  for(.j in j){
    data.table::set(
      x = flowstate$data,
      j = .j,
      value = get(transform.func)(flowstate$data[[.j]] / cofactor)
    )
    data.table::setattr(
      x = flowstate$data[[.j]],
      name = "transformed",
      value = mget(c("transform.func","cofactor"))
    )
  }
  ##
  invisible(flowstate)
}
#' @title Transform values to linear -- inverse
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param j Character vector -- default `NULL`; any/all parameters having a stored `transformed` attribute -- from the result of [flowstate.transform] -- will be inverse transformed in `[['data']]` using the stored function (default [base::sinh()]).  If a defined character vector: specific columns in `[['data']]` that are to be inverse transformed.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; inverse transformed values -- linear.
#'    \item `flowstate[['data']]`; the stored attribute -- `transformed` -- is removed from each `j`.
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
#'   transform.func = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated [['data']] attributes; applied transform function and cofactor
#' fs$data[, lapply(.SD, attr, which = 'transformed')]
#' fs$data[, attr(CD3, 'transformed')]
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#' fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#'
#' #inverse transformation; defined `j` argument
#' flowstate.transform.inverse(
#'   fs,
#'   j = c('CD3', 'CD8')
#' )
#'
#' #updated [['data']] attributes; returns only `CD4` as all other transformations have been inversed
#' fs$data[, lapply(.SD, attr, which = 'transformed')]
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
  ## has transformation been applied?
  ## check for the attribute
  res.transformed <- flowstate$data[, sapply(.SD, function(j)
    'transformed' %in% names(attributes(j)))]
  ##
  if(any(res.transformed)){
    ## what has been transformed?
    j.transformed <- names(which(res.transformed))
    ## subset if j is defined
    if(!is.null(j)){
      j.transformed <- j.transformed[j.transformed %in% j]
    }
  }else{
    return(message(
      sprintf(
        "'transformed' attribute not found in %s[['data']];\nno inverse transformation applied.",
        deparse(substitute(flowstate))
      )
    ))
  }
  ##
  message('flowstate --> transforming -- inverse...')
  ##...
  for(j in j.transformed){
    ## retrieve the transform attribute
    transform.attr <- attr(flowstate$data[[j]], which = "transformed")
    ## apply the inverse or error out
    if(transform.attr$transform.func == "asinh"){
      transform.func <- get("sinh")
      cofactor <- transform.attr$cofactor
      data.table::set(
        x = flowstate$data,
        j = j,
        value = transform.func(flowstate$data[[j]]) * cofactor
      )
      data.table::setattr(
        x = flowstate$data[[j]],
        name = "transformed",
        value = NULL
      )
    }else{
      stop(
        sprintf(
          "Currently only supports 'sinh' (the inverse of 'asinh'); '%s' is stored in the 'transform' attribute for '%s'",
          transform.attr$transform.func,
          j
        )
      )
    }
  }
  ##
  invisible(flowstate)
}
