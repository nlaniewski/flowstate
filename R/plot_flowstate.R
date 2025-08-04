#' @title Plot flowstate data
#'
#' @param x A flowstate object as returned from  [flowstate].
#' @param ... ... \link[ggplot2]{aes} arguments (unquoted variables); essentially x and y. If a `z` aesthetic is defined, the plot will switch to \link[ggplot2]{stat_summary_hex}, using this defined third variable as a 'color-by'.
#' @param bins \link[ggplot2]{geom_hex} argument; numeric giving the number of bins in both vertical and horizontal directions.
#' @param limits \link[ggplot2]{continuous_scale} argument.
#' @param sample.n Numeric (length 1); if defined, will randomly sample events from `[[data]]`.
#'
#' @returns A \link[ggplot2]{ggplot} object.
#' @export
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = ".fcs")
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   cofactor = 5000,
#'   concatenate = TRUE
#' )
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

plot_samples.barcode.assigned<-function(
    flowstate.object,
    col.pattern = "CD45BC",
    bins.N=200,
    metadata.include=c('sample.id','batch.name','block.id','well'))
{
  ##barcode columns
  barcode.dims<-grep(col.pattern,names(flowstate.object$data),value = T)
  ##plot metadata
  plot.metadata <- flowstate.object$data[
    barcode!=0 & cut == FALSE,
    j = list(
      total.events = .N,
      nodes.assigned = node_barcode |> unique() |> sort() |> paste0(collapse = ",")
    ),
    keyby = barcode
  ]
  plot.metadata[
    ,
    barcode.alias := paste0(utils::combn(length(barcode.dims),3)[,barcode],collapse = " : "),
    by=barcode
  ]
  plot.metadata <- plot.metadata[flowstate.object$meta,on='barcode']
  cols.convert<-plot.metadata[,names(.SD),.SDcols = is.factor]
  for(j in cols.convert){
    data.table::set(
      x = plot.metadata,
      j = j,
      value = as.character(plot.metadata[[j]])
    )
  }
  ##axis limits
  lims<-flowstate.object$data[
    barcode!=0 & cut==FALSE,
    lapply(.SD,range),.SDcols = barcode.dims
  ]
  ##plot pairs
  pp<-plot_pairs(barcode.dims)
  ##events to sample for background
  sample.N<-floor(1E5/flowstate.object$data[barcode!=0&cut==FALSE,length(unique(barcode))])
  ##background data
  dt.bkgd<-flowstate.object$data[,.SD[sample(.N,sample.N)],by=barcode,.SDcols = barcode.dims][,!'barcode']
  ##
  plot.bodies<-lapply(pp,function(.pp){
    ggplot2::ggplot(
      data = dt.bkgd,
      mapping = ggplot2::aes(
        x = !!as.name(.pp[1]),
        y = !!as.name(.pp[2]))
    ) +
      ggplot2::geom_hex(
        fill='gray',
        bins=bins.N
      ) +
      ggplot2::coord_cartesian(
        xlim = lims[[.pp[1]]],
        ylim = lims[[.pp[2]]]
      ) +
      ggplot2::theme_classic() +
      ggplot2::guides(fill = "none")
  })
  ##
  plotlist<-lapply(
    split(
      flowstate.object$data[
        barcode!=0 & cut == FALSE,
        .SD,
        .SDcols = c(barcode.dims,'barcode')
      ],
      by='barcode',
      sorted = T),
    function(dt){
      .barcode <- dt[,unique(barcode)]
      plot.bodies.barcode<-lapply(plot.bodies,function(.plot){
        .plot<-.plot +
          ggplot2::geom_hex(
            data = dt,
            bins = bins.N
          ) +
          viridis::scale_fill_viridis(
            option = "plasma",
            limits = NULL,
            oob = scales::squish
          )
      })
      ##
      patchwork::wrap_plots(plot.bodies.barcode,guides = 'collect') +
        patchwork::plot_annotation(
          title=paste0(
            plot.metadata[
              barcode == .barcode,
              paste("Barcode:",barcode, paste0("(",barcode.alias,")"))
            ],
            "\n",
            plot.metadata[
              barcode == .barcode,
              paste("N:", total.events)
            ],
            "\n",
            plot.metadata[
              barcode == .barcode,
              paste("Barcode-specific Nodes:", nodes.assigned)
            ],
            "\n",
            "\n",
            plot.metadata[
              barcode == .barcode,
              paste(metadata.include,.SD,sep = ": ",collapse = "\n"),
              .SDcols = metadata.include
            ]
          ),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12))
        )
    })
  ##
  return(plotlist)
}
