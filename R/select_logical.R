## should probably move the other 'select_...' function here for better organization

select_singlets <- function(flowstate, quantiles = c(0.85, 0.95)){
  ## remove doublets -- forward scatter geometry; initial selector
  ## selects events (collinear) based on quantile [1]
  flowstate$data[
    ,
    j = select.singlets := {
      res <- FSC_A/FSC_H
      res < stats::quantile(res, probs = quantiles[1])
    }
  ]
  ## remove doublets -- side scatter geometry; indexed against the initial selector
  ## selects events (collinear) based on quantile [2]
  flowstate$data[
    i = (select.singlets),
    j = select.singlets := {
      res <- SSC_A/SSC_H
      res < stats::quantile(res, probs = quantiles[2])
    }
  ]
  ##
  invisible(flowstate)
}
