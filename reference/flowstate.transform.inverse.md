# Transform values to linear – inverse

Transform values to linear – inverse

## Usage

``` r
flowstate.transform.inverse(flowstate, j = NULL)
```

## Arguments

- flowstate:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- j:

  Character vector – default `NULL`; any/all parameters having a stored
  `transformed` attribute – from the result of
  [flowstate.transform](https://nlaniewski.github.io/flowstate/reference/flowstate.transform.md)
  – will be inverse transformed in `[['data']]` using the stored
  function (default
  [`base::sinh()`](https://rdrr.io/r/base/Hyperbolic.html)). If a
  defined character vector: specific columns in `[['data']]` that are to
  be inverse transformed.

## Value

UPDATES BY REFERENCE:

- `flowstate[['data']]`; inverse transformed values – linear.

- `flowstate[['data']]`; the stored attribute – `transformed` – is
  removed from each `j`.

Invisibly returns `flowstate`.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read .fcs files as a flowstate object; concatenate
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type = "S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstates'...

#plot and mean values of linear columns
plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

res1.linear <- fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD8')]
print(res1.linear)
#>      CD3      CD8 
#> 25878.14 21843.70 

#transform
flowstate.transform(
  fs,
  j = c('CD3','CD4','CD8'),
  transform.func = "asinh",
  cofactor = 5000
)
#> flowstate --> transforming...

#updated [['data']] attributes; applied transform function and cofactor
fs$data[, lapply(.SD, attr, which = 'transformed')]
#>       CD8    CD3    CD4
#>    <list> <list> <list>
#> 1:  asinh  asinh  asinh
#> 2:   5000   5000   5000
fs$data[, attr(CD3, 'transformed')]
#>    transform.func cofactor
#>            <char>    <num>
#> 1:          asinh     5000

#plot and mean values of transformed columns from updated fs[['data']]
plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#>       CD3       CD4       CD8 
#> 1.8223658 1.6040152 0.9800536 

#inverse transformation; defined `j` argument
flowstate.transform.inverse(
  fs,
  j = c('CD3', 'CD8')
)
#> flowstate --> transforming -- inverse...

#updated [['data']] attributes; returns only `CD4` as all other transformations have been inversed
fs$data[, lapply(.SD, attr, which = 'transformed')]
#>       CD4
#>    <list>
#> 1:  asinh
#> 2:   5000

#plot and mean values of linear columns
plot(fs, CD3, CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

res2.linear <- fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD8')]
print(res2.linear)
#>      CD3      CD8 
#> 25878.14 21843.70 

#linear --> transformed --> inverse --> linear
res1.linear == res2.linear
#>  CD3  CD8 
#> TRUE TRUE 
```
