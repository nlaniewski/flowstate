concatenate.test<-function(flowstate.objects){
  ##test to make sure 'flowstate.objects' is a list
  if(!isa(flowstate.objects,"list")){
    stop("For concatenation: need a list of 'flowstate.objects'.")
  }
  ##test to make sure that each list element is a 'flowstate.object'
  if(!all(sapply(flowstate.objects,isa,'flowstate'))){
    stop("For concatenation: each list element of 'flowstate.objects' must be of class 'flowstate'.")
  }
  ##does each 'flowstate.object' have the same number of [['data']] columns?
  ncol.data<-unique(sapply(flowstate.objects,function(fs.obj){ncol(fs.obj[['data']])}))
  if(length(ncol.data)!=1){
    stop("For concatenation: each 'flowstate.object' must have the same number of data columns.")
  }
  ##does each 'flowstate.object' have the same [['data']] column names?
  colnames.data<-unique(lapply(flowstate.objects,function(fs.obj){names(fs.obj[['data']])}))
  if(length(colnames.data)!=1){
    stop("For concatenation: each 'flowstate.object' must have the same data column names.")
  }else{
    colnames.data<-colnames.data[[1]]
  }
  message("Concatenating 'flowstate.ojects'...")
}

parameters.unique<-function(flowstate.objects){
  ##if 'concatenate.test(flowstate.objects)' passes then parameters will be non-unique due to range ('R') differences in 'Time'
  parameters<-unique(data.table::rbindlist(lapply(flowstate.objects,'[[','parameters')))
  ##the names ('N') of parameters due to differing range ('R') values; 'Time'
  pars.rangefix<-parameters[,names(which(table(N)>1))]
  ##placeholder column; logical
  parameters[,drop := FALSE]
  ##for each parameter name in 'pars.rangefix', retain the row with max 'R' value;
  ##drop the rest
  for(par in pars.rangefix){
    i<-which(parameters[['N']] %in% par)
    pars.drop<-i[!i %in% i[parameters[i,which.max(R)]]]
    data.table::set(
      x = parameters,
      i = pars.drop,
      j = 'drop',
      value = TRUE
    )
  }
  ##drop
  parameters<-parameters[drop == FALSE][,drop := NULL]
  ##reorder based on 'par'
  data.table::setorder(
    x = parameters[
      ,
      ord := match(par,paste0("$P",seq(.N)))
    ],
    "ord"
  )[,ord := NULL]
  ##return the data.table
  parameters[]
}

#' @title Concatenate a list of flowstate objects
#'
#' @param flowstate.objects a list; the return of [read.flowstate]
#'
#' @returns An object of class flowstate.
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstate objects
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   cofactor = 5000
#' )
#' #a list of flowstate objects
#' class(fs)
#' sapply(fs,class)
#'
#' #concatenate into a single flowstate object
#' fs<-concatenate.flowstate(fs)
#' class(fs)
#'
#' fs$data
#' fs$parameters
#' fs$keywords
concatenate.flowstate<-function(flowstate.objects){
  ##test if objects can be concatenated
  concatenate.test(flowstate.objects)

  ##create a flowstate (fs) S3 object ; class 'flowstate'
  fs <- flowstate(
    data = data.table::rbindlist(lapply(flowstate.objects,'[[','data')),
    parameters = parameters.unique(flowstate.objects),
    keywords = data.table::rbindlist(lapply(flowstate.objects,'[[','keywords'),fill=TRUE),
    spill = unique(lapply(flowstate.objects,'[[','spill'))[[1]],#needs more testing,
    meta = data.table::rbindlist(lapply(flowstate.objects,'[[','meta'))
  )
  # ##S3 'flowstate' object structure; NULL
  # fs<-structure(
  #   list(
  #     data = NULL,
  #     parameters = NULL,
  #     keywords = NULL,
  #     spill = NULL,
  #     meta = NULL
  #   ),
  #   class = 'flowstate'
  # )
  # ##populate the object
  # fs[['data']] <- data.table::rbindlist(lapply(flowstate.objects,'[[','data'))
  # fs[['parameters']] <- parameters.unique(flowstate.objects)
  # fs[['keywords']] <- data.table::rbindlist(lapply(flowstate.objects,'[[','keywords'))
  # fs[['spill']] <- unique(lapply(flowstate.objects,'[[','spill'))[[1]]#needs more testing
  # fs[['meta']] <- data.table::rbindlist(lapply(flowstate.objects,'[[','meta'))
  ##return
  return(fs)
}
