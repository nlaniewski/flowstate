#' @title Plot flowstate data
#'
#' @param x A flowstate object as returned from [read.flowstate].
#' @param ... ... \link[ggplot2]{aes} arguments (unquoted variables); essentially x and y. If a `z` aesthetic is defined, the plot will switch to \link[ggplot2]{stat_summary_hex}, using this defined third variable as a 'color-by'.
#' @param bins \link[ggplot2]{geom_hex} argument; numeric giving the number of bins in both vertical and horizontal directions.
#' @param limits \link[ggplot2]{continuous_scale} argument.
#' @param sample.n Numeric (length 1); if defined, will randomly sample events from `[[data]]`.
#'
#' @returns A \link[ggplot2]{ggplot} object.
#' @export
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   concatenate = TRUE
#' )
#'
#' #transform
#' flowstate.transform(fs,c('CD3','CD4','CD8','Viability'))
#'
#' #plot title
#' no.fill.legend <- ggplot2::guides(fill = 'none')
#' .title1 <- paste("Batch:",fs$keywords[,unique(`$PROJ`)])
#' .title2 <- paste(
#'   "Instrument Serial#:",
#'   fs$keywords[,paste(.(unique(`$CYT`),unique(`$CYTSN`)),collapse = " ")]
#' )
#' .title <- paste(.title1,.title2,sep = "\n")
#' .title <- ggplot2::labs(title = .title)
#'
#' #plot: two variables
#' plot(fs,CD3,Viability) + no.fill.legend + .title
#' plot(fs,FSC_A,SSC_A) + no.fill.legend + .title
#'
#' #plot: two variables; facet by sample.id
#' plot(fs,CD4,CD8) + no.fill.legend + .title + ggplot2::facet_wrap(~sample.id)
#'
#' #plot: three variables; color-by
#' plot(fs,FSC_A,SSC_A,z=Viability,bins=100,limits=c(0,3)) + .title
#'
plot.flowstate <- function(x,..., bins = 200, limits = NULL, sample.n = NULL){
  p <- ggplot2::ggplot(
    data = if(!is.null(sample.n)){
      x$data[sample(.N,sample.n)]
    }else{
      x$data
    },
    mapping = ggplot2::aes(...)
  )
  if("z" %in% ...names()){
    dot.names = lapply(substitute(list(...))[-1], deparse)
    z <- dot.names[...names() == "z"]
    p <- p + ggplot2::stat_summary_hex(bins = bins) +
      viridis::scale_fill_viridis(
        option = "magma",
        limits = limits,
        oob = scales::squish
      ) +
      ggplot2::labs(fill = z)
  }
  else{
    p <- p + ggplot2::geom_hex(bins = bins) +
      viridis::scale_fill_viridis(
        option = "plasma",
        limits = limits,
        oob = scales::squish
      )
  }
  return(p)
}

plot_pairs<-function (plot.dims)
{
  pair.1 <- plot.dims[seq(1, length(plot.dims), 2)]
  pair.2 <- plot.dims[seq(2, length(plot.dims), 2)]
  if (length(pair.1) != length(pair.2)) {
    pair.2[length(pair.1)] <- pair.2[length(pair.1) - 1]
  }
  pair.list <- mapply(c, pair.1, pair.2, SIMPLIFY = F, USE.NAMES = F)
}
nxn <- function(x){
  dt.melt <- data.table::melt(
    data = x[,id := seq(.N)],
    id.vars = 'id',
    value.name = 'value.x',
    variable.name = 'variable.x'
  )
  dt.melt <- data.table::merge.data.table(
    x = dt.melt,
    y = data.table::setnames(
      data.table::copy(dt.melt),
      c("variable.x", "value.x"),
      c("variable.y", "value.y")
    ),
    allow.cartesian = TRUE
  )[,id := NULL]
  cols.plot <- dt.melt[, levels(variable.x)]
  cols.combn <- as.list(as.data.frame(utils::combn(cols.plot,2)))
  dt.melt[,combn.drop := TRUE]
  invisible(lapply(cols.combn, function(col.combn) {
    data.table::set(
      x = dt.melt,
      i = dt.melt[,.I[variable.x == col.combn[1] & variable.y == col.combn[2]]],
      j = "combn.drop",
      value = FALSE)
  }))
  dt.melt <- (
    dt.melt[combn.drop == FALSE]
    [,combn.drop := NULL]
  )
  invisible(dt.melt)
}
plot_nxn_barcodes_assigned_samples <- function(
    flowstate.object,
    censor=TRUE,
    col.pattern = "CD45BC",
    sample.n = 5E4,
    metadata.include = c("sample.id","batch","block.id","well")
)
{
  ##
  js <- grep(col.pattern,names(flowstate.object$data),value = T)
  ##
  # lims <- fs$data[
  #   i = (!barcode.censor),
  #   j = lapply(.SD, range),
  #   .SDcols = cols.barcode
  # ]
  ##background/silhouette; aggregate data
  background <- flowstate.object$data[
    i = (!barcode.censor),
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n)}]),
    .SDcols=js
  ]
  ##foreground/sample-specific data; as a keyed (barcode) data.table
  foreground <- flowstate.object$data[
    i = (!barcode.censor),
    j = nxn(.SD[{set.seed(1337);sample(.N,sample.n)}]),
    .SDcols=js,
    keyby = barcode
  ]
  ##metadata for plot
  metadata.for.plot <- flowstate.object$data[
    ,
    .(
      N = .N,
      nodes.assigned = paste0(sort(unique(node_barcode)),collapse = ","),
      barcode.alias = paste0(utils::combn(length(js),3)[, barcode], collapse = " : ")
    ),
    keyby = .(barcode,barcode.censor)
  ]
  metadata.for.plot <- metadata.for.plot[
    ,
    .(
      totals = sprintf("%d (%d)",N[(!barcode.censor)],N[(barcode.censor)]),
      nodes.assigned = unique(nodes.assigned),
      barcode.alias = unique(barcode.alias)
    ),
    by=barcode
  ]
  metadata.for.plot <- (
    flowstate.object$keywords[,.SD,.SDcols = c('barcode',metadata.include)]
    [,barcode := as.numeric(barcode)]
    [metadata.for.plot,on='barcode']
  )
  cols.convert <- metadata.for.plot[, names(.SD), .SDcols = is.factor]
  for (j in cols.convert) {
    data.table::set(
      x = metadata.for.plot,
      j = j,
      value = as.character(metadata.for.plot[[j]])
    )
  }
  ##create ggplot object
  plot.body <- ggplot2::ggplot(
    data = NULL,
    mapping = ggplot2::aes(value.x,value.y)
  )
  ##add background data
  p <- plot.body +
    ggplot2::geom_hex(
      data = background,
      bins = 100,
      fill = "gray"
    )
  ##sample-specific plots; as a list
  p.list <- sapply(foreground[,unique(barcode)],function(bc,mdat=metadata.for.plot[barcode == bc]){
    ##add foreground data
    p <- p +
      ggplot2::geom_hex(
        data = foreground[barcode == bc],
        bins = 100,
        show.legend = FALSE
      ) +
      viridis::scale_fill_viridis(
        option = "plasma",
        limits = c(0,100),
        oob = scales::squish
      )
    ##add faceting and theme
    p <- p +
      ggplot2::facet_grid(variable.y~variable.x,scales="free") +
      ggplot2::theme_void()
    ##add metadata
    p <- p +
      ggplot2::labs(
        title = mdat[,sprintf("Barcode: %d  (%s)",barcode,barcode.alias)],
        subtitle = paste0(mdat[
          ,
          j = paste(metadata.include, .SD, sep = ": ",collapse = "\n"),
          .SDcols = metadata.include
        ],'\n'),
        caption = paste(
          mdat[,sprintf("Barcode-specific Nodes: %s",nodes.assigned)],
          mdat[,sprintf("Totals: %s (# = censored events)",totals)],
          sprintf("Displayed: %d randomly sampled",sample.n),
          sep = "     "
        )
      )
    ##
    return(p)
  },simplify = F)
  ##
  return(p.list)
}
