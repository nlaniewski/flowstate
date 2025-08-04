fcs.text.primary.required.keywords<-
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

flowstate.parameter.keywords<-
  c(
    'B',
    'DISPLAY',
    'E',
    'N',
    'R',
    'S',
    'TYPE',
    'V'
  )

flowstate.transform.inverse<-function(flowstate.object){
  if(!any(c('transform','cofactor') %in% names(flowstate.object$parameters))){
    message(
      paste(
        "An inverse transformation of a flowstate.object requires the following columns in [['parameters']]: 'transform' and/or 'cofactor';",
        "these columns are not found -- the inverse transformation may have already been performed.",
        sep = "\n"
      )
    )
    return()
  }
  parms<-data.table::copy(flowstate.object$parameters)[!is.na(transform)]
  col.match<-parms[,sapply(.SD,function(j){all(j %in% names(flowstate.object$data))})]
  col.match<-names(which(col.match))[1]
  parms<-parms[grep("Raw|Unmixed_Fluorescence",TYPE),.SD,.SDcols = c(col.match,'transform','cofactor')]
  data.table::setnames(parms,old = col.match, new = 'alias')
  for(j in parms[!is.na(transform),alias]){
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = if(parms[alias == j][['transform']] %in% 'asinh'){
        sinh(flowstate.object$data[[j]])*parms[alias == j][['cofactor']]
      }else{
        stop("In development: no inverse function defined for this transform type.")
      }
    )
  }
  flowstate.object$parameters[,c('transform','cofactor') := NULL]
}

flowstate.transform<-function(flowstate.object,.j,transform.type="asinh",cofactor=5000){
  if(!all(.j %in% names(flowstate.object$data))){
    stop(".j not found in [['data']]")
  }
  ##
  col.match<-flowstate.object$parameters[,sapply(.SD,function(j){all(.j %in% j)})]
  col.match<-names(which(col.match))
  if(length(col.match)>1) col.match<-col.match[1]
  ##
  if('transform' %in% names(flowstate.object$parameters)){
    if(any(.j %in% flowstate.object$parameters[!is.na(transform)][[col.match]])){
      stop("Transformation has already been applied to .j!")
    }
  }
  ##
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
      i = which(flowstate.object$parameters[[col.match]] %in% j),
      j = c('transform','cofactor'),
      value = list(transform.type,cofactor)
    )
  }
  ##
}

add.keywords.to.data <- function(flowstate.object,keywords){
  ##
  res <- names(which(flowstate.object$keywords[
    ,
    sapply(.SD,function(j){anyNA(suppressWarnings(as.numeric(j)))}),
    .SDcols = keywords
  ]))
  if(length(res)>0){
    stop(
      paste(
        "Adding the following keywords to [['data']] will be problematic:",
        paste0(res,collapse = ", "),
        "Eventual conversion to numeric is expected.",
        sep = "\n"
      )
    )
  }
  ##
  ids<-flowstate.object$data[,levels(sample.id)]
  col.id <- names(which(flowstate.object$keywords[,sapply(.SD,function(j){all(j == ids)})]))
  if(length(col.id) != 1){
    stop("Expect a single identifier column ('sample.id') in [['data']] to match to a single identifier column in [['keywords']]")
  }
  ##
  totals<-flowstate.object$data[,.N,by=sample.id]
  for(j in keywords){
    data.table::set(
      x = flowstate.object$data,
      i = NULL,
      j = j,
      value = rep(flowstate.object$keywords[[j]],totals[['N']])
    )
  }
}
