.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c("$FIL", "$TOT", ".", ".mean", "alias", "autofluorescence",
    "barcode", "barcode.alias", "barcode.censor", "cluster", "combn.drop",
    "detector", "detector.peak", "detector.select", "emission", "emission.normalized",
    "fluor", "FSC_A", "group.type", "i", "id", "index", "j", "js",
    "marker", "means", "N", "N.alias", "node_barcode", "nodes.assigned",
    "ord", "peak.values", "population", "PROJ", "R", "reference.type",
    "S", "S.alias", "S_N.alias", "sample.id", "scatter.population",
    "scatter.select", "select.detector", "select.nonsaturating",
    "select.population", "select.scatter", "select.spectral", "SSC_A",
    "SSCB_A", "subtract.type", "tissue.type", "total.events", "totals",
    "TYPE", "V1", "value", "value.x", "value.y", "variable", "variable.x",
    "variable.y")
)


