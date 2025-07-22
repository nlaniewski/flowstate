#' @title Plot flowstate data
#'
#' @param x A flowstate object as returned from  [flowstate].
#' @param ... ... \link[ggplot2]{aes} arguments (unquoted); essentially x and y. If a `z` aesthetic is provided, the plot will switch to \link[ggplot2]{stat_summary_hex}, using this defined third variable as a 'color-by'.
#' @param bins \link[ggplot2]{geom_hex} argument; numeric vector giving the number of bins in both vertical and horizontal directions.
#' @param limits \link[ggplot2]{continuous_scale} argument.
#' @param sample.n Numeric (length 1); if defined, will randomly sample events from `[[data]]`.
#'
#' @returns A \link[ggplot2]{ggplot} object.
#' @export
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
