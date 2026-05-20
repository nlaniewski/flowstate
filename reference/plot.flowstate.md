# Visualize `flowstate` data

Visualize `flowstate` data

## Usage

``` r
# S3 method for class 'flowstate'
plot(x, ..., bins = 200, limits = NULL, sample.n = NULL)
```

## Arguments

- x:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- ...:

  ... [aes](https://ggplot2.tidyverse.org/reference/aes.html) arguments
  (unquoted variables); essentially x and/or y. If a `z` aesthetic is
  defined, the plot will switch to
  [stat_summary_hex](https://ggplot2.tidyverse.org/reference/stat_summary_2d.html),
  using this defined third variable as a 'color-by'.

- bins:

  [geom_hex](https://ggplot2.tidyverse.org/reference/geom_hex.html)
  argument; numeric giving the number of bins in both vertical and
  horizontal directions for bivariate and/or bivariate with 'color-by'
  plots.

- limits:

  [continuous_scale](https://ggplot2.tidyverse.org/reference/continuous_scale.html)
  argument.

- sample.n:

  Numeric (length 1); if defined, will randomly sample events from
  `[[data]]`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read all .fcs files as flowstate objects; concatenate into a single object
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type = "S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstates'...

#transform
flowstate.transform(fs, c('CD3', 'CD4', 'CD8', 'Viability'))
#> flowstate --> transforming...

#plot title
no.fill.legend <- ggplot2::guides(fill = 'none')
.title1 <- paste("Batch:", fs$keywords[, unique(`$PROJ`)])
.title2 <- paste(
  "Instrument Serial#:",
  fs$keywords[, paste(.(unique(`$CYT`), unique(`$CYTSN`)), collapse = " ")]
)
.title <- paste(.title1, .title2, sep = "\n")
.title <- ggplot2::labs(title = .title)

#plot: two variables
plot(fs,CD3,Viability) + no.fill.legend + .title
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

plot(fs,FSC_A,SSC_A) + no.fill.legend + .title
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.


#plot: two variables; facet by sample.id
plot(fs,CD4,CD8) + no.fill.legend + .title + ggplot2::facet_wrap(~sample.id)
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.


#plot: three variables; color-by
plot(fs,FSC_A,SSC_A,z=Viability,bins=100,limits=c(0,3)) + .title
#> Warning: Computation failed in `stat_summary_hex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_summary_hex()`.

```
