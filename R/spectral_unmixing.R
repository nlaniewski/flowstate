flowstate.unmixed.object <- function(flowstate.object.raw,ref.medians){
  ##initialize a new flowstate.object;
  ##inherits keyword-value pairs from flowstate.object.raw;
  ##inherits non-detector data from flowstate.object.raw;
  ##derives 'unmixed' parameters from ref.medians

  ##initialize an empty flowstate
  fs.unmixed <- flowstate()
  ##copy/inherit keywords
  fs.unmixed$keywords <- data.table::copy(flowstate.object.raw$keywords)
  ##GUID is some sort of SpectroFlo hash -- specific to each source file;
  ##no longer applicable
  if('GUID' %in% names(fs.unmixed$keywords)){
    fs.unmixed$keywords[,j := NULL, env = list(j = "GUID")]
  }
  ##copy/inherit parameters; raw/overdetermined --> unmixed;
  ##new $PnN and $PnS keyword-value pairs as derived from ref.medians;
  ##if there are failures here, ref.medians will need to be 'fixed'
  fs.unmixed$parameters <- parameters.raw.to.unmixed(flowstate.object.raw,ref.medians)
  ##spillover matrix; derived from [['parameters']]
  spill.mat <- matrix(0,nrow = ref.medians[,.N],ncol = ref.medians[,.N])
  colnames(spill.mat) <- rownames(spill.mat) <- ref.medians[['alias']]
  diag(spill.mat) <- 1
  spill.mat[2,1] <- 1E-6#bug-fix for FlowJo
  spill <- data.table::as.data.table(spill.mat)
  rownames(spill) <- names(spill)
  data.table::setattr(spill,'applied',FALSE)
  fs.unmixed$spill <- spill
  ##names of detectors
  j.match <- j.match.parameters.to.data(flowstate.object.raw)
  cols.detector <- flowstate.object.raw$parameters[TYPE == "Raw_Fluorescence"][[j.match]]
  i <- names(flowstate.object.raw$data) %in% cols.detector#c(cols.detector,'sample.id')
  cols.retain <- names(flowstate.object.raw$data)[!i]
  ##assign to fs.unmixed[['data']] non-detector parameters; scatter and time
  fs.unmixed$data <- flowstate.object.raw$data[,.SD,.SDcols = cols.retain]
  ##
  invisible(fs.unmixed)
}
