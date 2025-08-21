.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
utils::globalVariables(
  c("$FIL", "$TOT", ".", "alias", "barcode", "barcode.alias", "detector",
    "detector.peak", "detector.select", "emission", "emission.normalized",
    "fluor", "FSC_A", "j", "marker", "N", "N.alias", "node_barcode",
    "nodes.assigned", "ord", "PROJ", "R", "S", "S.alias", "S_N.alias",
    "sample.id", "scatter.select", "SSC_A", "tissue.type", "total.events",
    "TYPE", "V1", "value", "variable")
)
