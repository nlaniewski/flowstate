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
  parms<-data.table::copy(flowstate.object$parameters)
  col.match<-parms[,sapply(.SD,function(j){all(names(flowstate.object$data) %in% j)})]
  col.match<-names(which(col.match))
  parms<-parms[grep("Raw|Unmixed_Fluorescence",TYPE),.SD,.SDcols = c(col.match,'transform','cofactor')]
  data.table::setnames(parms,old = col.match, new = 'alias')
  for(j in parms$alias){
    data.table::set(
      x = flowstate.object$data,
      j = j,
      value = if(parms[alias == j][['transform']] == 'asinh'){
        sinh(flowstate.object$data[[j]])*parms[alias == j][['cofactor']]
      }else{
        stop("In development: no inverse function defined for this transform type.")
      }
    )
  }
  flowstate.object$parameters[,c('transform','cofactor') := NULL]
}
