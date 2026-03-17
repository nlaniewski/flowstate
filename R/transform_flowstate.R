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
#' @title Transform `flowstate[['data']]` -- Inverse
#'
#' @param flowstate.object A flowstate object as returned from [read.flowstate].
#' @param .j Character vector -- default `NULL`; any/all parameters having a keyword-value pair of `transform/asinh` will be inverse transformed in `[['data']]` using [base::sinh()].  If a character vector: specific columns in `[['data']]` that are to be inverse transformed.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; inverse transformed values -- linear
#'    \item `flowstate[['parameters']]`; modifies two columns ('transform' and 'cofactor') -- sets to `NA`
#' }
#'
#' Invisibly returns the `flowstate.object`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read .fcs files as a flowstate object; concatenate
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   concatenate = TRUE
#' )
#'
#' #plot and mean values of linear columns
#' plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#' res1.linear <- fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
#' print(res1.linear)
#'
#' #transform
#' flowstate.transform(
#'   fs,
#'   .j = c('CD3','CD8'),
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#' #updated parameters
#' fs$parameters[!is.na(transform)]
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#' fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
#'
#' #inverse transformation
#' flowstate.transform.inverse(fs)
#'
#' #updated parameters; transform and cofactor set to NA
#' fs$parameters[S %in% c('CD3','CD8')]
#'
#' #plot and mean values of linear columns
#' plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#' res2.linear <- fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
#' print(res2.linear)
#'
#' #linear --> transformed --> inverse --> linear
#' res1.linear == res2.linear
#'
flowstate.transform.inverse<-function(flowstate.object,.j=NULL){
  ##name of [['parameters']] column that matches [['data']]
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##create a subset of [['parameters']]; needs 'transform' and 'cofactor' columns
  parameters.transformed <- flowstate.object$parameters[
    i = !is.na(transform),
    .SD,
    .SDcols = c(j.match,'transform','cofactor')
  ]
  data.table::setnames(
    parameters.transformed,
    old = j.match,
    new = 'alias'
  )
  ##
  if(!is.null(.j)){
    parameters.transformed<-parameters.transformed[alias %in% .j]
  }
  ##
  for(j in parameters.transformed[['alias']]){
    ##inverse transform [['data']]
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = if(parameters.transformed[alias == j][['transform']] %in% 'asinh'){
        sinh(flowstate.object$data[[j]])*parameters.transformed[alias == j][['cofactor']]
      }else{
        stop("In development: no inverse function defined for this transform type.")
      }
    )
    ##update [['parameters']]
    data.table::set(
      x = flowstate.object$parameters,
      i = which(flowstate.object$parameters[[j.match]] %in% j),
      j = c('transform','cofactor'),
      value = list(NA,NA)
    )
  }
}
