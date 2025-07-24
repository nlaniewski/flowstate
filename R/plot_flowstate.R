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
