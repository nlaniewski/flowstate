parameters_from_keywords <- function(keywords){
  ## using the return of readFCS::readTEXT()
  ### keyword-value pairs

  ## number of parameters
  par.n <- as.numeric(keywords[["$PAR"]])
  ## regex pattern to find '$P#' and 'P#'
  ### required parameters are prefixed by '$'
  ### custom/instrument-specific parameters begin with 'P'
  par.pattern <- "^\\$?P\\d+"
  parameter.names <- grep(
    pattern = par.pattern,
    x = names(keywords),
    value = TRUE
  )
  ## parameters -- named vector
  pars <- trimws(unlist(keywords[parameter.names]))
  ## split parameters vector into 'types' -- list
  par.types <- split(
    x = pars,
    f = factor(sub(par.pattern, "", names(pars)))
  )
  ## convert to numeric
  pars.numeric <- names(which(sapply(par.types, function(par) {
    !anyNA(suppressWarnings(as.numeric(par)))
  })))
  for(par in pars.numeric){
    .names <- names(par.types[[par]])
    par.types[[par]] <- stats::setNames(as.numeric(par.types[[par]]), nm = .names)
  }
  ## create a data.table
  parameters.dt <- data.table::data.table(par = paste0('$P', seq(par.n)))
  ## fill the data.table by 'par.type'
  for(j in names(par.types)){
    par.vec <- par.types[[j]]
    par.i <- as.integer(gsub("\\D+", "", names(par.vec)))
    data.table::set(
      parameters.dt,
      i = par.i,
      j = j,
      value = par.vec
    )
  }
  ## return
  invisible(parameters.dt)
}

keywords_required <- c(
  "$BEGINANALYSIS", "$BEGINDATA", "$BEGINSTEXT", "$BYTEORD",
  "$DATATYPE", "$ENDANALYSIS", "$ENDDATA", "$ENDSTEXT", "$MODE",
  "$NEXTDATA", "$PAR", "$PnB", "$PnE", "$PnN", "$PnR", "$TOT"
)

## add keywords related to originality/modification
### as defined for FCS 3.1
keywords_originality <- function(){
  list(
    '$ORIGINALITY' = "DataModified",
    '$LAST_MODIFIED' = toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS2")),
    '$LAST_MODIFIER' = sprintf("flowstate_%s", utils::packageVersion("flowstate"))
  )
}

#' @title Print method for flowstate objects
#'
#' @param x An object of class flowstate; as returned from `flowstate_object`.
#' @param ... Additional arguments passed to [data.table][data.table::print.data.table].
#'
#' @return Invisibly returns `flowstate`; prints to the console.
#' @import data.table
#' @export
print.flowstate <- function(x, ...) {
  NextMethod()
  invisible(x)
}

# as.flowstate <- function(flowstate) {
#   stopifnot(data.table::is.data.table(flowstate))
#   if (!inherits(flowstate, "flowstate")) {
#     data.table::setattr(flowstate, "class", c("flowstate", class(flowstate)))
#   }
#   flowstate
# }
#
# `[.flowstate` <- function(x, ...) {
#   res <- NextMethod()
#   if (data.table::is.data.table(res)) {
#     as.flowstate(res)
#   } else {
#     res
#   }
# }
