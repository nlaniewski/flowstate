fcs.text.primary.required.keywords <-
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

fcs.text.primary.required.keywords.zero <- stats::setNames(
  nm = fcs.text.primary.required.keywords,
  rep("0", length(fcs.text.primary.required.keywords))
)

flowstate.parameter.keywords <-
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

flowstate.transform.keywords <- c(
  "KIND",#BD FACSDiscover [AS]8 -- BD Chorus
  "TYPE"#Cytek Aurora -- SpectroFlo
)

flowstate.transform.strings <- c(
  "Raw_Fluorescence",#Cytek Aurora -- SpectroFlo
  "Unmixed_Fluorescence",#Cytek Aurora -- SpectroFlo
  "Ion_Count",#DVS/Fluidigm -- CyTOF
  'COLOR'#BD FACSDiscover [AS]8 -- BD Chorus
)

cytometer.identifier.types <- c(
  "Aurora" = "TUBENAME",
  "ID7000" = "$CELLS",
  "DVS|FLUIDIGM|CYTOF" = "$FIL",
  "FACSDiscover A8" = "$SMNO",# Sample Manager Name Object (?) -- BD Chorus
  "FACSDiscover S8" = "$SMNO"# Sample Manager Name Object (?) -- BD Chorus
)

cytometer.identifier <- function(fcs.file.paths){
  cyt <- check.keyword(fcs.file.paths,keyword = '$CYT')
  i <- which(sapply(names(cytometer.identifier.types),function(i){grepl(i,cyt)}))
  if(length(i) != 1){
    message("Cytometer keyword 'identifier' could not be resolved;\nis the cytometer indexed in 'flowstate:::cytometer.identifer.types()'?")
    id <- '$FIL'
  }else{
    id <- cytometer.identifier.types[[i]]
  }
}

alias_dt <- function(flowstate){
  alias.vec <- flowstate$data[, unlist(sapply(.SD, attr, which = 'N'))]
  data.table::data.table(
    N = alias.vec,
    alias = names(alias.vec)
  )
}

j.match.keyword.to.data.sample.id <- function(flowstate.object){
  flowstate.object$keywords[
    ,
    sapply(.SD,function(j){
      sum(j %in% flowstate.object$data[,levels(sample.id)],na.rm = TRUE)
    })
  ] |> which.max() |> names()
}

#' @title Add `[['keyword']]` values to `[['data']]`
#' @description
#' Used for efficient addition of factored keyword values to `[['data']]`; updates by reference.
#'
#' @param flowstate A flowstate object as returned from [read.flowstate].
#' @param keywords.to.add Character vector; keyword names in `[['keywords']]` whose factored values will be added to `[['data']]`.
#'
#' @returns UPDATES BY REFERENCE:
#' \itemize{
#'   \item `[['data']]`; factored `keywords.to.add` are added as additional columns.
#'   \item `[['keywords']]`; `keywords.to.add` are converted to class `factor`.
#' }
#'
#' Invisibly returns `flowstate`.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstates; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #create a new set of keywords/values
#' fs$keywords[
#' ,
#' j = c('block.id', 'block.aliquot') := data.table::tstrsplit(
#' sample.id, "_", keep = 4:5)
#' ]
#' #new columns; character class
#' fs$keywords[, .(block.id, block.aliquot)]
#'
#' #add the factored keyword values to fs[['data']]
#' add.keywords.to.data(fs, c('block.id', 'block.aliquot'))
#'
#' fs$data[, .(block.id, block.aliquot)]
#'
#' #add additional keywords to [['data']]
#' add.keywords.to.data(fs, c('$DATE', '$PROJ'))
#'
#' fs$data[, .(`$DATE`, `$PROJ`)]
#' fs$keywords[, .(`$DATE`, `$PROJ`)]
#'
add.keywords.to.data <- function(flowstate, keywords.to.add){
  ## factor keywords.to.add
  for(j in keywords.to.add){
    data.table::set(
      x = flowstate$keywords,
      j = j,
      value = factor(flowstate$keywords[[j]])
    )
  }
  ## keywords.to.add to [['data']]; index by sample.id (levels)
  for(i in flowstate$keywords[, levels(sample.id)]){
    data.table::set(
      x = flowstate$data,
      i = flowstate$data[, .I[sample.id == i]],
      j = keywords.to.add,
      value = flowstate$keywords[
        i = sample.id == i,
        j = ..keywords.to.add
      ]
    )
  }
  ##
  invisible(flowstate)
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

check.keyword <- function(fcs.file.paths, keyword = NULL, value = NULL){
  kw <- sapply(fcs.file.paths, function(i){
    readFCStext(i)[[keyword]]
  })
  if(!is.null(value)){
    all(grepl(value, kw))
  }else{
    return(unique(kw))
  }
}

cosine.similarity <- function(x,y){crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))}
