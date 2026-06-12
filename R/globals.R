.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c(
    "$TOT", ".", "..keywords.to.add", "FSC_A", "FSC_H", "N", "N.alias",
    "R", "S", "S.alias", "SSC_A", "SSC_H", "TYPE", "V1", "alias",
    "barcode", "barcode.alias", "barcode.censor", "combn.drop", "detector.peak",
    "detector.select", "id", "j", "j.match", "j.match.parameters.to.data",
    "node_barcode", "nodes.assigned", "ord", "sample.id", "scatter.select",
    "select.nonsaturating", "select.population", "select.quantile",
    "select.singlets", "totals", "value", "value.x", "value.y", "variable",
    "variable.x", "variable.y"
  )
)
