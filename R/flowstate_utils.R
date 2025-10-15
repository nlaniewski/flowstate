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
    'DETECTOR',
    'DISPLAY',
    'E',
    'N',
    'R',
    'S',
    'TYPE',
    'V'
  )

j.match.parameters.to.data <- function(flowstate.object){
  names(which.max(flowstate.object$parameters[
    ,
    sapply(.SD,function(j){sum(j %in% names(flowstate.object$data),na.rm = TRUE)})
  ]))[1]
}
j.match.keyword.to.data.sample.id <- function(flowstate.object){
  flowstate.object$keywords[
    ,
    sapply(.SD,function(j){
      sum(j %in% flowstate.object$data[,levels(sample.id)],na.rm = TRUE)
    })
  ] |> which.max() |> names()
}
#' @title Add keyword values to `flowstate[['data']]`
#' @description
#' Used for efficient addition of factored keyword values to `flowstate[['data']]`; updates by reference.
#'
#' @param flowstate.object A flowstate object as returned from [read.flowstate].
#' @param keywords Character string; keyword names in `flowstate[['keywords']]` whose factored values will be added to `flowstate[['data']]`.
#' @param type.convert Logical; default `TRUE`. For any `keywords` of class `character`, conversion to `factor` will take place before addition to `[['data']]`.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'   \item `flowstate[['data']]`; factored `keywords` are added as additional columns.
#'   \item `flowstate[['keywords']]`; if `type.convert`, `keywords` are converted to class `factor`.
#' }
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
#' #create a new set of keywords/values
#' fs$keywords[
#' ,
#' j = c('block.id','block.aliquot') := data.table::tstrsplit(
#' TUBENAME,"_",type.convert = factor,keep=4:5)
#' ]
#'
#' #add the factored keyword values to fs[['data']]
#' add.keywords.to.data(fs,c('block.id','block.aliquot'))
#'
#' fs$data[,.(block.id,block.aliquot)]
#'
#' #add to keyword [['data']] using type.convert argument
#' add.keywords.to.data(fs,c('$DATE','$PROJ'),type.convert=TRUE)
#'
#' fs$data[,.(`$DATE`,`$PROJ`)]
#' fs$keywords[,.(`$DATE`,`$PROJ`)]
add.keywords.to.data <- function(flowstate.object,keywords,type.convert=TRUE){
  ##
  if(type.convert){
    cols.convert <- flowstate.object$keywords[
      ,
      names(which(sapply(.SD,is.character))),
      .SDcols = keywords
    ]
    for(j in cols.convert){
      data.table::set(
        x = flowstate.object$keywords,
        j = j,
        value = factor(flowstate.object$keywords[[j]])
      )
    }
  }
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
  col.id <- names(which(
    flowstate.object$keywords[,sapply(.SD,function(j){all(j %in% ids)})]
  ))
  if(length(col.id)==0){
    col.id <- names(which(
      flowstate.object$keywords[,sapply(.SD,function(j){all(sub(".fcs","",j) %in% ids)})]
    ))
  }
  if(length(col.id) == 0){
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
##
# mass.cytometry.detect <- function(flowstate.object,sample.val = 1E3,threshold.val = 0.25){
#   zeroes <- flowstate.object$data[
#     i = sample(.N,sample.val),
#     j = sapply(.SD,function(j){length(which(j==0))})
#   ]
#   res <- length(which(zeroes>0))/length(zeroes)
#   isTRUE(res>threshold.val)
# }
