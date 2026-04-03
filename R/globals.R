.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c("$FIL", "$TOT", ".", "..keywords.to.add", ".mean", "add.reference.keywords.to.data",
    "alias", "autofluorescence", "barcode", "barcode.alias", "barcode.censor",
    "cellular.anchors", "cluster", "combn.drop", "detector", "detector.peak",
    "detector.select", "emission", "emission.normalized", "fluor",
    "FSC_A", "group.type", "i", "id", "index", "j", "j.match", "j.match.parameters.to.data",
    "js", "k", "marker", "means", "N", "N.alias", "node_barcode",
    "nodes.assigned", "ord", "peak.values", "population", "PROJ",
    "R", "reference.type", "ribbon", "S", "S.alias", "S_N.alias",
    "sample.id", "scatter.population", "scatter.select", "select.contour",
    "select.detector", "select.nonsaturating", "select.population",
    "select.quantile", "select.scatter", "select.scatter.population",
    "select.spectral", "select.spectral.events", "SSC_A", "SSCB_A",
    "subtract.type", "tissue.type", "total.events", "totals", "type",
    "TYPE", "V1", "value", "value.x", "value.y", "variable", "variable.x",
    "variable.y")
)
