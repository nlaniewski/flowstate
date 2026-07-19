# Transform linear values

Linear values in `[['data']]` are transformed using a defined function.
For full spectrum cytometry, transforming linear fluorescence values
using [asinh](https://rdrr.io/r/base/Hyperbolic.html) is preferred, with
around-zero compression adjusted using a 'cofactor' of 5,000.

## Usage

``` r
flowstate.transform(
  flowstate,
  j = NULL,
  transform.func = "asinh",
  cofactor = 5000
)
```

## Arguments

- flowstate:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- j:

  Character vector – default `NULL`; any/all parameters having a
  keyword-value pair of (instrument-specific):

  - `Aurora`: `TYPE %in% c("Raw_Fluorescence", "Unmixed_Fluorescence")`

  - `FACSDiscover [AS]8`: `KIND %in% COLOR`

  will be transformed in `[['data']]`. If a defined character vector:
  specific columns in `[['data']]` that are to be transformed.

- transform.func:

  Character vector – default
  [asinh](https://rdrr.io/r/base/Hyperbolic.html); quoted function that
  will be used to transform `j` in `[['data']]`.

- cofactor:

  Numeric – default 5000; if `transform.func` =
  [asinh](https://rdrr.io/r/base/Hyperbolic.html) (default), the
  cofactor will be used to modify the transformation.

## Value

UPDATES BY REFERENCE:

- `flowstate[['data']]`; transformed values.

- `flowstate[['data']]`; an attribute – `transformed` – is added to each
  `j`.

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

fs$data[, sapply(.SD, mean), .SDcols = c('CD3', 'CD4', 'CD8')]
#>      CD3      CD4      CD8 
#> 25878.14 19921.68 21843.70 

#transform
flowstate.transform(
  fs,
  j = c('CD3', 'CD4', 'CD8'),
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

#transform all/remaining
flowstate.transform(
  fs,
  j = NULL,
  transform.func = "asinh",
  cofactor = 5000
)
#> flowstate --> transforming...

#updated [['data']] attributes; applied transform function and cofactor
fs$data[, lapply(.SD, attr, which = 'transformed')]
#>    CD45RA CD45RO  TCRgd CD45BC1    IL2    CD8  CD197 CD45BC2 CD45BC3   CD57
#>    <list> <list> <list>  <list> <list> <list> <list>  <list>  <list> <list>
#> 1:  asinh  asinh  asinh   asinh  asinh  asinh  asinh   asinh   asinh  asinh
#> 2:   5000   5000   5000    5000   5000   5000   5000    5000    5000   5000
#>     CD193 CD45BC4  CD127   CD56  CD199  CD49a GranzymeB CD45BC5   CD95    CD3
#>    <list>  <list> <list> <list> <list> <list>    <list>  <list> <list> <list>
#> 1:  asinh   asinh  asinh  asinh  asinh  asinh     asinh   asinh  asinh  asinh
#> 2:   5000    5000   5000   5000   5000   5000      5000    5000   5000   5000
#>       CD4   CD69   TNFa  CD183   IFNg  CD103 CD45BC6  CD122  ia4b7 CD45BC7
#>    <list> <list> <list> <list> <list> <list>  <list> <list> <list>  <list>
#> 1:  asinh  asinh  asinh  asinh  asinh  asinh   asinh  asinh  asinh   asinh
#> 2:   5000   5000   5000   5000   5000   5000    5000   5000   5000    5000
#>    Viability   CD27  KLRG1
#>       <list> <list> <list>
#> 1:     asinh  asinh  asinh
#> 2:      5000   5000   5000
```
