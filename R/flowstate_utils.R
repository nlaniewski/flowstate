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

flowstate.transform.inverse<-function(flowstate.object,.j){
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
  j.match<-parms[
    ,
    names(which(sapply(.SD,function(j){all(j %in% names(flowstate.object$data))})))[1]
  ]
  parms<-parms[
    grep("Raw|Unmixed_Fluorescence",TYPE),
    .SD,
    .SDcols = c(j.match,'transform','cofactor')
  ]
  data.table::setnames(parms,old = j.match, new = 'alias')
  if(!is.null(.j)){
    parms<-parms[alias %in% .j]
  }
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
    data.table::set(
      x = flowstate.object$parameters,
      i = which(flowstate.object$parameters[[j.match]] %in% j),
      j = c('transform','cofactor'),
      value = list(NA,NA)
    )
  }
}

flowstate.transform<-function(flowstate.object,.j,transform.type="asinh",cofactor=5000){
  if(!all(.j %in% names(flowstate.object$data))){
    stop(".j not found in [['data']]")
  }
  ##
  j.match<-flowstate.object$parameters[
    ,
    names(which(sapply(.SD,function(j){all(.j %in% j)})))[1]
  ]
  ##
  if('transform' %in% names(flowstate.object$parameters)){
    if(any(.j %in% flowstate.object$parameters[!is.na(transform)][[j.match]])){
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
      i = which(flowstate.object$parameters[[j.match]] %in% j),
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
