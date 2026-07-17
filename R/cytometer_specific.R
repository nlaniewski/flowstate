# BD FACSDiscover [AS]8 ------------------------------------------------------

facsdiscover.bdspectral.terms <- c(
  "BDSPECTRAL LAMBDA",
  "BDSPECTRAL QSPE",
  "BDSPECTRAL UNMIXING METHOD"
)

facsdiscover.json.terms <- c(
  "BDCHORUSDATARECORD",
  "BDFIRESETTINGS"
)

## parse BDSPECTRAL keywords: LAMBDA, QSPE, and UNMIXING METHOD
### construct a data.table

bdspectral.data.table <- function(keywords){
  ## using the return of 'readFCStext(...)'
  ## 'LAMBDA' and 'QSPE' keyword names
  mat <- sapply(c('LAMBDA', 'QSPE'), function(i){
    string <- keywords[[sprintf("BDSPECTRAL %s", i)]]
    string.split <- unlist(strsplit(string, ","))
    n.cols <- as.numeric(string.split[1])
    col.names <- string.split[2:(n.cols + 1)]
    vals <- as.numeric(string.split[(n.cols + 2):length(string.split)])
    stats::setNames(vals, nm = col.names)
  })
  ##
  dt <- data.table::data.table(
    N = rownames(mat),
    mat
  )
  ##
  data.table::setattr(
    x = dt,
    name = "BDSPECTRAL UNMIXING METHOD",
    value = keywords[['BDSPECTRAL UNMIXING METHOD']]
  )
  ##
  invisible(dt)
}


