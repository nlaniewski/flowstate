.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
utils::globalVariables(
  c("$FIL", "$TOT", "alias", "barcode", "barcode.alias", "j", "N",
    "N.alias", "node_barcode", "nodes.assigned", "ord", "PROJ", "R",
    "S", "S.alias", "S_N.alias", "sample.id","total.events", "TYPE")
)
