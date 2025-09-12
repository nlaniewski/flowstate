##S3 constructor for flowstate object
new_flowstate <- function(
    data = data.table::data.table(),
    parameters = data.table::data.table(),
    keywords = data.table::data.table(),
    spill = data.table::data.table(),
    meta = data.table::data.table(),
    conc = data.table::data.table()
){
  #validate
  sapply(ls(environment()),function(arg){
    if(!data.table::is.data.table(get(arg))){
      stop(paste(arg, "is not a data.table!"),call. = F)
    }
  })
  ##flowstate object
  structure(
    as.list(environment()),
    class = "flowstate"
  )
}

flowstate <- function(...){
  new_flowstate(...)
}

#' @title Subset `flowstate[['data']]`
#' @description
#' Subset `flowstate[['data']]` and return a new `flowstate` object.
#'
#' @param x A flowstate object as returned from [read.flowstate].
#' @param ... `subset` argument as defined by `data.table`'s \link[data.table]{subset}; logical expression indicating elements or rows to keep.
#'
#' @returns A `flowstate` object whose `[['data']]` element (a `data.table`) has been subset to contain the rows as defined through the `subset` argument. Other elements of `x` are copied.
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
#' flowstate.transform(fs,'CD3')
#'
#' #plot;visualize valley
#' fs$data[,plot(density(CD3))]
#' abline(v = 1.1)#not data-driven; just for example
#'
#' #subset
#' fs.cd3_positive <- subset(fs,CD3 > 1.1)
#'
#' #totals
#' fs$data[,.N]
#' fs.cd3_positive$data[,.N]
subset.flowstate <- function(x,...){
  ##a new (empty) flowstate object
  x.subset <- flowstate()
  ##S3 method dispatch of subset; data.table:::subset.data.table()
  ##returns a subset of [['data']] which needs assignment
  x.subset$data <- subset(x$data,...)
  ##copy other flowstate list elements (data.tables)
  for(i in names(x)[!names(x) %in% "data"]){
    x.subset[[i]] <- data.table::copy(x[[i]])
  }
  ##return
  invisible(x.subset)
}

#' @title Split a concatenated `flowstate` object into a list
#'
#' @param x A concatenated flowstate object as returned from [read.flowstate].
#' @param f Character string; a single (factored) column by which `[['data']]` will be split into individual `flowstate`s.
#' @param drop Logical; default `TRUE`. Unused factor levels in `f` are dropped -- no `flowstate` object will be returned.
#' @param ... Potential further arguments.
#'
#' @returns A list of individual `flowstate`s.
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
#' fs.split <- split(fs,'sample.id')
#'
#' length(fs.split)
#' names(fs.split)
#' lapply(fs.split,'[[','keywords')
#'
#' #remove all rows associated with sample 1; level will remain
#' fs$data <- fs$data[!sample.id %in% levels(sample.id)[1]]
#'
#' #with drop
#' fs.split <- split(fs,'sample.id',drop = TRUE)
#' length(fs.split)
#' names(fs.split)
#'
#' #without drop
#' fs.split <- split(fs,'sample.id',drop = FALSE)
#' length(fs.split)
#' names(fs.split)
#' #empty data for sample 1
#' fs.split[[1]]$data
#'
split.flowstate <- function(x,f,drop = TRUE, ...){
  ##use factor levels of 'f' argument
  names.split <- levels(x$data[[f]])
  ##is there a single column in [['keywords']] that matches the levels?
  ##as a result of using read.flowstate(...,sample.id = ...);
  ##or from adding after the fact using add.keywords.to.data()
  col.keyword <- names(which(
    x$keywords[,sapply(.SD,function(j){all(j %in% names.split)})]
  ))
  if(length(col.keyword)>1){
    col.keyword <- sapply(col.keyword,function(i){
      names.split %in% unique(x$keywords[[i]]) |> sum()
    }) |> which.max() |> names()

  }
  if(col.keyword == 'sample.id'){
    col.keyword <- NULL
  }
  ##as a list of new (empty) flowstate objects
  x.split <- stats::setNames(
    lapply(seq_along(names.split),function(x){flowstate()}),
    nm = names.split
  )
  ##populate the flowstate objects
  for(i in names(x.split)){
    ##assign the [['data']] subset
    x.split[[i]][['data']] <- subset(x[['data']],sample.id %in% i)
    ##assign [['parameters']]
    x.split[[i]][['parameters']] <- data.table::copy(x[['parameters']])
    ##assign [['keywords']] subset
    x.split[[i]][['keywords']] <- if(is.null(col.keyword)){
      subset(x[['keywords']],sample.id %in% i)
    }else{
      subset(x[['keywords']],get(col.keyword) %in% i)
    }
    ##assign [['spill']] subset
    x.split[[i]][['spill']] <- data.table::copy(x[['spill']])
  }
  ##drop empty [['data']]; due to unused factor level
  if(drop){
    i.drop <- sapply(x.split,function(i){i[['data']][,.N==0]})
  }else{
    i.drop <- NULL
  }
  ##return
  if(is.null(i.drop)){
    invisible(x.split)
  }else{
    invisible(x.split[!i.drop])
  }
}
