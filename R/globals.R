.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
res <- utils::globalVariables(
  c(
    "$CYT", "$TOT", ".", "..keywords.to.add", "alias", "barcode",
    "barcode.alias", "barcode.censor", "combn.drop", "detector.peak",
    "detector.select", "fsc.a", "fsc.h", "FSC_A", "FSC_H", "id",
    "j", "j.match", "j.match.parameters.to.data", "kw", "N", "N.alias",
    "node_barcode", "nodes.assigned", "ord", "R", "S", "S.alias",
    "sample.id", "scatter.select", "select.nonsaturating", "select.population",
    "select.quantile", "select.singlets", "ssc.a", "ssc.h", "SSC_A",
    "SSC_H", "totals", "TYPE", "V1", "value", "value.x", "value.y",
    "variable", "variable.x", "variable.y"
  )
)
