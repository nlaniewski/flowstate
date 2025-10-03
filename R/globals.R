.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(c(...),c(...)) |> sort() |> dput()
utils::globalVariables(
  c("$FIL", "$TOT", ".", "alias", "barcode", "barcode.alias", "barcode.censor",
    "combn.drop", "detector", "detector.peak", "detector.select",
    "emission", "emission.normalized", "fluor", "FSC_A", "id", "j",
    "marker", "N", "N.alias", "node_barcode", "nodes.assigned", "ord",
    "peak.values", "PROJ", "R", "S", "S.alias", "S_N.alias", "sample.id",
    "scatter.select", "SSC_A", "tissue.type", "total.events", "totals",
    "TYPE", "V1", "value", "value.x", "value.y", "variable", "variable.x",
    "variable.y")
)
