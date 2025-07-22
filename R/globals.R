.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
utils::globalVariables(
  c("$TOT",
    "alias",
    "j",
    "N",
    "N.alias",
    "PROJ",
    "S",
    "S.alias",
    "S_N.alias",
    "TYPE"
  )
)
