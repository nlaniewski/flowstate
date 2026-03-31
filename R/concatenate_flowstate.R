concatenate.test <- function(flowstates){
  ## test to make sure 'flowstate.objects' is a list
  if(!isa(flowstates,"list")){
    stop("For concatenation: need a list of 'flowstates'.")
  }
  ## test to make sure that each list element is a 'flowstate'
  if(!all(sapply(flowstates,isa,'flowstate'))){
    stop("For concatenation: each list element of 'flowstates' must be of class 'flowstate'.")
  }
  ## does each 'flowstate' have the same number of [['data']] columns?
  ncol.data <- unique(sapply(flowstates, function(fs.obj){ncol(fs.obj[['data']])}))
  if(length(ncol.data) != 1){
    stop("For concatenation: each 'flowstate' must have the same number of data columns.")
  }
  ## does each 'flowstate' have the same [['data']] column names?
  colnames.data <- unique(lapply(flowstates,function(fs.obj){names(fs.obj[['data']])}))
  if(length(colnames.data) != 1){
    stop("For concatenation: each 'flowstate' must have the same data column names.")
  }else{
    colnames.data <- colnames.data[[1]]
  }
  ##
  message("Concatenating 'flowstates'...")
}

parameters.unique <- function(flowstates){
  ## if 'concatenate.test(flowstates)' passes then parameters will be non-unique due to possible differences in:
  ## 'R' (range -- 'Time');
  ## 'V' (volts/gain);

  ## resolve to a unique [['parameters']]
  parameters <- unique(data.table::rbindlist(lapply(flowstates, '[[', 'parameters')))
  ## duplicate names ('N') of parameters due to differing range ('R') values; 'Time'
  pars.rangefix <- parameters[, names(which(table(N) > 1))]
  ## placeholder column; logical
  parameters[, drop := FALSE]
  ## for each parameter name in 'pars.rangefix', retain the row with max 'R' value;
  ## drop the rest
  for(par in pars.rangefix) {
    i <- which(parameters[['N']] %in% par)
    pars.drop <- i[!i %in% i[parameters[i, which.max(R)]]]
    data.table::set(x = parameters,
                    i = pars.drop,
                    j = 'drop',
                    value = TRUE)
  }
  ## drop
  parameters <- parameters[drop == FALSE][, drop := NULL]
  ## reorder based on 'par'
  data.table::setorder(
    x = parameters[, ord := match(par, paste0("$P", seq(.N)))],
    "ord"
  )[,ord := NULL]
  ## return the data.table
  invisible(parameters)
}

#' @title Concatenate a list of flowstates
#'
#' @param flowstates a list; the return of [read.flowstate]
#'
#' @returns An object of class flowstate.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstates
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S"
#' )
#'
#' #a list of flowstates
#' class(fs)
#' sapply(fs, class)
#'
#' #concatenate into a single flowstate
#' fs <- concatenate.flowstate(fs)
#' class(fs)
#'
#' fs$data
#' fs$parameters
#' fs$keywords
#' fs$spill
#'
concatenate.flowstate<-function(flowstates){
  ## test if objects can be concatenated
  concatenate.test(flowstates)
  ## is there [['spill']]"
  res.spill <- any(sapply(flowstates, function(fs) 'spill' %in% names(fs)))
  ## create a flowstate (fs) S3 object ; class 'flowstate'
  fs <- flowstate(
    data = data.table::rbindlist(lapply(flowstates, '[[', 'data')),
    parameters = parameters.unique(flowstates),
    keywords = data.table::rbindlist(lapply(flowstates,'[[','keywords'), fill = TRUE),
    spill = if(res.spill){
      unique(lapply(flowstates, '[[', 'spill'))[[1]]#needs more testing
    }else{
      data.table::data.table()
    }
  )
  ## any/all concatenated .fcs files should return:
  ## [['data']], [['parameters']], and [['keywords']]
  ## may not have spill (mass cytometry)
  if(fs[['spill']][,.N]==0){fs[['spill']] <- NULL}
  ##return
  return(fs)
}
