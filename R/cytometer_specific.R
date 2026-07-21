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

## 'conventional' BD FACSDiscover [AS]8 .fcs file
### drop all imaging-related parameters
### drop all 'time-to-peak' parameters
### retain only 'conventional' parameters
facsdiscover.conventional.parameters <- function(flowstate){
  ##
  cyt <- flowstate$keywords[, unique(`$CYT`)]
  stopifnot(grepl("FACSDiscover [AS]8", cyt))
  ## copy parameters
  parameters <- data.table::copy(flowstate$parameters)
  ## initialize a logical
  parameters[, drop.params := FALSE]
  ##
  parameters[
    i = !FEATURE %in% c("Area", "Height", "Width"),
    j = drop.params := TRUE
  ]
  ##
  parameters[
    i = KIND == "COLOR" & MEAS != "A",
    j = drop.params := TRUE
  ]
  ##
  parameters[
    i = !(drop.params) & grepl("Img|Imaging", N),
    j = drop.params := TRUE
  ]
  ##
  parameters[
    i = N %in% c("Time", "Saturated"),
    j = drop.params := FALSE
  ]
  ## subset to drop parameters
  parameters <- (
    subset(parameters, !drop.params)
    [, drop.params := NULL]
  )
  ## needs reassignment
  return(parameters)
}

## 'conventional' BD FACSDiscover [AS]8 .fcs file
### drop all imaging-related data
### drop all 'time-to-peak' data
### retain only 'conventional' data -- greatly reduces object size
facsdiscover.conventional.data <- function(flowstate){
  ##
  cols.keep <- flowstate$parameters[, N]
  i.drop <- flowstate$data[, !sapply(.SD, attr, which = "N") %in% cols.keep]
  cols.drop <- names(flowstate$data)[i.drop]
  ##
  flowstate$data[, (cols.drop) := NULL]
  ##
  invisible(flowstate)
}
