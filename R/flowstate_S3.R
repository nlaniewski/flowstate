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
