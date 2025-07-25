% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spillover_flowstate.R
\name{spillover.update.value}
\alias{spillover.update.value}
\title{Update the values of a flowstate spill data.table}
\usage{
spillover.update.value(flowstate.object, i, j, value)
}
\arguments{
\item{flowstate.object}{the return of \link{read.flowstate}.}

\item{i}{Variable (unquoted).}

\item{j}{Variable (unquoted).}

\item{value}{Numeric; correction value to be used during compensation.}
}
\value{
Updated \code{flowstate$spill}; !!!updates by reference -- no assignment.
}
\description{
Update the values of a flowstate spill data.table
}
\examples{
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = ".fcs")

#read all .fcs files as flowstate objects; concatenate into a single object
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type="S",
  cofactor = 5000,
  concatenate = TRUE
)

#row index; CD4 vs CD8
index<-which(names(fs$spill) \%in\% c('CD4','CD8'))
fs$spill[index,.(CD4,CD8)]

#update a spill value
spillover.update.value(fs,CD8,CD4,0.03)

fs$spill[index,.(CD4,CD8)]
}
