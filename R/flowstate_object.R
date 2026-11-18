flowstate_object <- function(fcs.file.path, parameter.fix = NULL){
  ##
  fcs <- readFCS::readDATA(fcs.file.path) |> data.table::as.data.table()
  ##
  keywords <- readFCS::readTEXT(fcs.file.path)
  ##
  parameters <- parameters_from_keywords(keywords)
  ##
  kw.required <- keywords[names(keywords) %in% keywords_required]
  ##
  keywords <- keywords[!names(keywords) %in% keywords_required]
  i.parameters <- grep("^\\$?P\\d+", names(keywords), ignore.case = T)
  keywords <- keywords[-i.parameters]
  keywords <- c(keywords, keywords_originality())
  ##
  fcs[, sample.id := factor(keywords[['TUBENAME']])]
  ##
  data.table::setattr(keywords, 'required', unlist(kw.required))
  data.table::setattr(fcs, 'keywords', keywords)
  ##
  if(!is.null(parameter.fix)){
    for(parm.name in names(parameter.fix)){
      parameters[N == parm.name, S := parameter.fix[[parm.name]]]
    }
  }

  if(!'TYPE' %in% names(parameters)){
    parameters[grepl("^[FS]SC", N), TYPE := "Scatter"]
    parameters[N == 'Time', TYPE := "Time"]
    parameters[!is.na(S), TYPE := "Fluorescence"]
  }
  ##
  data.table::setattr(fcs, 'parameters', parameters)
  ##
  type.transform <- parameters[
    i = grepl("fluorescence", TYPE, ignore.case = T),
    j = .(N, S)
  ] |> unlist(use.names = F) |> unique()

  for(j in type.transform[type.transform %in% names(fcs)]){
    data.table::setattr(x = fcs[[j]], name = 'transformed', value = FALSE)
  }
  ##
  # class(fcs) <- c("flowstate", class(fcs))
  ##
  invisible(fcs[])
}

names_flowstate <- function(flowstate, type = c("S", "N")){
  ## for setting names, access 'parameters' attribute
  ### N <--> S
  parameters.NS <- attr(flowstate, 'parameters')[, .(N, S)]
  ##
  name.type <- switch(
    match.arg(type),
    S = c('N', 'S'),
    N = c('S', 'N')
  )
  ##
  data.table::setnames(
    x = flowstate,
    old = parameters.NS[!is.na(S)][[name.type[1]]],
    new = parameters.NS[!is.na(S)][[name.type[2]]],
    skip_absent = T
  )
}

transform_flowstate <- function(
    flowstate,
    j = NULL,
    cofactor = 5000,
    inverse = F
)
{
  ## check transformation status, access per-column (j)'transformed' attribute
  res.transformed <- sapply(flowstate, attr, 'transformed')
  res.transformed <- res.transformed[!sapply(res.transformed, is.null)]
  ##
  if(inverse){
    for(j in names(res.transformed)){
      if(is.numeric(res.transformed[[j]])){
        data.table::set(
          x = flowstate,
          j = j,
          value = sinh(flowstate[[j]]) * res.transformed[[j]]
        )
        data.table::setattr(
          x = flowstate[[j]],
          name = 'transformed',
          value = FALSE
        )
      }
    }
    ##
    invisible(flowstate)
  }
  ##
  if(!is.null(j)){
    res.transformed <- res.transformed[intersect(names(res.transformed), j)]
  }
  ##
  for(j in names(res.transformed)){
    if(!res.transformed[[j]]){
      data.table::set(
        x = flowstate,
        j = j,
        value = asinh(flowstate[[j]]/cofactor)
      )
      data.table::setattr(
        x = flowstate[[j]],
        name = 'transformed',
        value = cofactor
      )
    }else if(res.transformed[[j]] != cofactor){
      data.table::set(
        x = flowstate,
        j = j,
        value = asinh((sinh(flowstate[[j]]) * res.transformed[[j]]) / cofactor)
      )
      data.table::setattr(
        x = flowstate[[j]],
        name = 'transformed',
        value = cofactor
      )
    }
  }
  ##
  invisible(flowstate)
}

plot_flowstate <- function(
    x,
    ...,
    bins = 300,
    limits = NULL,
    sample.n = NULL
)
{
  # dots <- as.list(match.call(expand.dots = FALSE))$...
  # .dots <- sapply(dots, deparse)

  .ncol <- ncol(x)
  ##
  if(.ncol == 1){
    p <- ggplot2::ggplot(
      data = x,
      mapping = ggplot2::aes(!!as.name(names(x)))
    )
    p <- p + ggplot2::geom_density()
  }
  if(.ncol == 2){
    p <- ggplot2::ggplot(
      data = x,
      mapping = ggplot2::aes(
        x = !!as.name(names(x)[1]),
        y = !!as.name(names(x)[2])
      )
    )
    p <- p + ggplot2::geom_hex(
      bins = bins,
      show.legend = F
    )
    p <- p +
      viridis::scale_fill_viridis(
        option = 'plasma',
        limits = limits,
        oob = scales::squish
      )
  }
  if(.ncol == 3){
    vars <- x[, sapply(.SD, class)]
    vars.numeric <- names(which(vars == "numeric"))
    vars.factor <- names(which(vars == "factor"))

    p <- ggplot2::ggplot(
      data = x,
      mapping = ggplot2::aes(
        x = !!as.name(vars.numeric[1]),
        y = !!as.name(vars.numeric[2])
      )
    )
    p <- p + ggplot2::geom_hex(
      bins = bins,
      show.legend = F
    )
    p <- p + viridis::scale_fill_viridis(
      option = 'plasma',
      limits = limits,
      oob = scales::squish
    )
    p <- p + ggplot2::facet_wrap(vars.factor)
  }
  ##
  return(p)
}
