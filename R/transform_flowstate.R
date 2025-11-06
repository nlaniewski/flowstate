#' @title Transform `flowstate[['data']]`
#'
#' @param flowstate.object A flowstate object as returned from [read.flowstate].
#' @param .j  Character vector -- default `NULL`; any/all parameters having a keyword-value pair of `TYPE/Raw|Unmixed_Fluorescence` will be transformed in `[['data']]`.  If a character vector: specific columns in `[['data']]` that are to be transformed.
#' @param transform.type Character string -- default \link{asinh}; quoted function that will be used to transform `.j` in `[['data']]`.
#' @param cofactor Numeric -- default 5000; if `transform.type` = \link{asinh} (default), the cofactor will be used to modify the transformation.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'    \item `flowstate[['data']]`; transformed values
#'    \item `flowstate[['parameters']]`; adds two columns ('transform' and 'cofactor')
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
#' fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD4','CD8')]
#'
#' #transform
#' flowstate.transform(
#'   fs,
#'   .j = c('CD3','CD4','CD8'),
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#' #updated parameters
#' fs$parameters[!is.na(transform)]
#'
#' #plot and mean values of transformed columns from updated fs[['data']]
#' plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#' fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD4','CD8')]
#'
#' #transform all
#' flowstate.transform(
#'   fs,
#'   .j = NULL,
#'   transform.type = "asinh",
#'   cofactor = 5000
#' )
#'
#' #updated parameters
#' fs$parameters[!is.na(transform)]
#'
flowstate.transform<-function(flowstate.object,.j = NULL,transform.type="asinh",cofactor=5000){
  ##name of [['parameters']] column that matches [['data']]
  j.match <- j.match.parameters.to.data(flowstate.object)
  if(is.null(.j)){
    .j <- flowstate.object$parameters[TYPE %in% flowstate.transform.types][[j.match]]
  }
  if(!all(.j %in% names(flowstate.object$data))){
    stop(paste(
      paste0(.j[!.j %in% names(flowstate.object$data)],collapse = ", "),
      "not found in [['data']]")
    )
  }
  ##
  if('transform' %in% names(flowstate.object$parameters)){
    if(any(.j %in% flowstate.object$parameters[!is.na(transform)][[j.match]])){
      .j <- .j[!.j %in% flowstate.object$parameters[!is.na(transform)][[j.match]]]
      # stop("Transformation has already been applied to .j!")
      if(length(.j)==0){
        return(
          message("Transformation has already been applied to .j!")
        )
      }
    }
  }
  ##
  if(flowstate.object$data[,length(levels(sample.id))==1]){
    message(paste(flowstate.object$data[,levels(sample.id)],"-->","transforming..."))
  }else{
    message('flowstate.object --> transforming...')
  }
  trans.func<-get(transform.type)
  for(j in .j){
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = trans.func(flowstate.object$data[[j]]/cofactor)
    )
    ##
    data.table::set(
      x = flowstate.object$parameters,
      i = which(flowstate.object$parameters[[j.match]] %in% j),
      j = c('transform','cofactor'),
      value = list(transform.type,cofactor)
    )
  }
  ##
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
