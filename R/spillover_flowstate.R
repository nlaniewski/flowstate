spillover.update.external<-function(flowstate.object,spillover.flowjo.path.csv){
  ##read previously exported spillover matrix (.csv) from FlowJo
  spill.mat.external<-utils::read.csv(
    spillover.flowjo.path.csv,
    check.names = FALSE,
    row.names = NULL
  )
  ##drop column 1
  if(is.character(spill.mat.external[,1])){
    spill.mat.external<-spill.mat.external[-1]
  }
  ##split column names using the FlowJo delimiter
  cols.split<-strsplit(colnames(spill.mat.external)," :: ")
  ##which split matches the names of flowstate.object[['spill']]
  split.index<-which(sapply(seq(unique(sapply(cols.split,length))),function(i){
    all(sapply(cols.split,'[[',i) %in% names(flowstate.object$spill))
  }))
  ##update the column names
  colnames(spill.mat.external)<-sapply(cols.split,'[[',split.index)
  ##update flowstate.object[['spill']]
  for(j in names(flowstate.object$spill)){
    data.table::set(
      x = flowstate.object$spill,
      j = j,
      value = spill.mat.external[[j]]
    )
  }
}
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
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
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
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   concatenate = TRUE
#' )
#'
#' #transform
#' flowstate.transform(fs,c('CD4','CD8'))
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
  spill.index<-sapply(c(1,2),function(margin){
    which(apply(flowstate.object$spill,margin,function(x){any(x!=0&x!=1&x!=1E-6)}))
  }) |> Reduce(f = union, x = _) |> sort()
  spill.index<-names(flowstate.object$spill)[spill.index]
  ##are any js in spill.index already transformed?
  if('transform' %in% names(flowstate.object$parameters)){
    j.match <- j.match.parameters.to.data(flowstate.object)
    parameters.subset<-flowstate.object$parameters[
      i = flowstate.object$parameters[[j.match]] %in% spill.index & !is.na(transform)
    ]
    ##need to reverse the transformation;
    ##need raw values when applying spillover
    flowstate.transform.inverse(flowstate.object,.j=parameters.subset[[j.match]])
  }
  ##
  spill.mat<-as.matrix(
    flowstate.object$spill[
      i = rownames(flowstate.object$spill) %in% spill.index,
      j = .SD,
      .SDcols = spill.index
    ]
  )
  ##
  flowstate.object$data[
    ,
    (spill.index) := as.list(
      as.data.frame(
        if(decompensate){
          as.matrix(flowstate.object$data[,.SD,.SDcols = spill.index]) %*% spill.mat
        }else{
          as.matrix(flowstate.object$data[,.SD,.SDcols = spill.index]) %*% solve(spill.mat)
        }
      )
    )
  ]
  if('transform' %in% names(flowstate.object$parameters)){
    ##reapply the transformation
    flowstate.transform(flowstate.object,.j=parameters.subset[[j.match]])
  }
  ##
  if(decompensate){
    data.table::setattr(flowstate.object$spill,name = "applied",value = FALSE)
  }else{
    data.table::setattr(flowstate.object$spill,name = "applied",value = TRUE)
  }
}
##
