#' @title Update the values of a flowstate spill data.table
#'
#' @param flowstate.object the return of [read.flowstate].
#' @param i Variable (unquoted).
#' @param j Variable (unquoted).
#' @param value Numeric; correction value to be used during compensation.
#'
#' @returns Updated `flowstate$spill`; !!!updates by reference -- no assignment.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = ".fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   cofactor = 5000,
#'   concatenate = TRUE
#' )
#'
#' #row index; CD4 vs CD8
#' index<-which(names(fs$spill) %in% c('CD4','CD8'))
#' fs$spill[index,.(CD4,CD8)]
#'
#' #update a spill value
#' spillover.update.value(fs,CD8,CD4,0.03)
#'
#' fs$spill[index,.(CD4,CD8)]
spillover.update.value<-function(flowstate.object,i,j,value){
  i.char<-deparse(substitute(i));j.char<-deparse(substitute(j))
  vec.i<-sapply(c(i.char,j.char),function(x){which(names(flowstate.object$spill) %in% x)})
  for(i in names(vec.i)){
    l<-length(vec.i[[i]])
    if(l!=1){
      stop(
        sprintf(
          "Indexing for %s returned %s results; is %s found/unique in [['spill']]?",
          i,l,i
        )
      )
    }
  }
  ##
  data.table::set(
    x = flowstate.object$spill,
    i = vec.i[1],
    j = vec.i[2],
    value = value
  )
}

#' @title Compensate data using values stored in spill
#'
#' @param flowstate.object the return of [read.flowstate].
#' @param decompensate Logical; if `TRUE`, will return data values to their original state.
#'
#' @returns Updated `flowstate$data`; !!!updates by reference -- no assignment.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = ".fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   cofactor = 5000,
#'   concatenate = TRUE
#' )
#'
#' #row index; CD4 vs CD8
#' index<-which(names(fs$spill) %in% c('CD4','CD8'))
#'
#' #update a spill value
#' spillover.update.value(fs,CD8,CD4,0.03)
#'
#' #over-compensated data
#' plot(fs,CD4,CD8) + ggplot2::labs(title = "Over-compensated")
#'
#' #apply compensation
#' spillover.apply(fs)
#'
#' #compensated data
#' plot(fs,CD4,CD8) + ggplot2::labs(title = "Compensated")
#'
#' #return data values to original state by decompensating
#' spillover.apply(fs,decompensate=TRUE)
#'
#' plot(fs,CD4,CD8) + ggplot2::labs(title = "Over-compensated")
spillover.apply<-function(flowstate.object,decompensate=FALSE){
  if(attr(flowstate.object$spill,'applied') & isFALSE(decompensate)){
    stop("Spillover has already been applied.")
  }
  cols.spill<-sapply(c(1,2),function(margin){
    which(apply(flowstate.object$spill,margin,function(x){any(x!=0&x!=1&x!=1E-6)}))
  })
  cols.spill<-names(flowstate.object$spill)[cols.spill]
  ##transform lookup; in development
  # col.match<-fs$parameters[,sapply(.SD,function(j){all(names(cols.spill) %in% j)})]
  # col.match<-names(which(col.match))
  # if(length(col.match)>1) col.match<-col.match[1]
  # fs$parameters[S %in% names(cols.spill)]
  ##need to reverse the transformation;
  ##need raw values when applying spillover
  for(j in cols.spill){
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = sinh(flowstate.object$data[[j]])*5000
    )
  }
  ##
  spill.mat<-as.matrix(
    flowstate.object$spill[
      i = rownames(flowstate.object$spill) %in% cols.spill,
      j = .SD,
      .SDcols = cols.spill
    ]
  )
  ##
  flowstate.object$data[
    ,
    (cols.spill) := as.list(
      as.data.frame(
        if(decompensate){
          as.matrix(flowstate.object$data[,.SD,.SDcols = cols.spill]) %*% spill.mat
        }else{
          as.matrix(flowstate.object$data[,.SD,.SDcols = cols.spill]) %*% solve(spill.mat)
        }
      )
    )
  ]
  ##
  for(j in cols.spill){
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = asinh(flowstate.object$data[[j]]/5000)
    )
  }
  ##
  if(decompensate){
    data.table::setattr(flowstate.object$spill,name = "applied",value = FALSE)
  }else{
    data.table::setattr(flowstate.object$spill,name = "applied",value = TRUE)
  }
}
##
