# Transform `flowstate[['data']]` – Inverse

Transform `flowstate[['data']]` – Inverse

## Usage

``` r
flowstate.transform.inverse(flowstate.object, .j = NULL)
```

## Arguments

- flowstate.object:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- .j:

  Character vector – default `NULL`; any/all parameters having a
  keyword-value pair of `transform/asinh` will be inverse transformed in
  `[['data']]` using
  [`base::sinh()`](https://rdrr.io/r/base/Hyperbolic.html). If a
  character vector: specific columns in `[['data']]` that are to be
  inverse transformed.

## Value

UPDATES BY REFERENCE:

- `flowstate[['data']]`; inverse transformed values – linear

- `flowstate[['parameters']]`; modifies two columns ('transform' and
  'cofactor') – sets to `NA`

Invisibly returns the `flowstate.object`.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read .fcs files as a flowstate object; concatenate
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type="S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstate.ojects'...

#plot and mean values of linear columns
plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

res1.linear <- fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
print(res1.linear)
#>      CD3      CD8 
#> 25878.14 21843.70 

#transform
flowstate.transform(
  fs,
  .j = c('CD3','CD8'),
  transform.type = "asinh",
  cofactor = 5000
)
#> flowstate.object --> transforming...
#updated parameters
fs$parameters[!is.na(transform)]
#>       par      B DISPLAY      E        N       R      S                 TYPE
#>    <char> <char>  <char> <char>   <char>  <char> <char>               <char>
#> 1:   $P16     32     LOG    0,0 BUV805-A 4194304    CD8 Unmixed_Fluorescence
#> 2:   $P30     32     LOG    0,0  RB744-A 4194304    CD3 Unmixed_Fluorescence
#>         V                           PROJ N.alias S.alias  S_N.alias transform
#>    <char>                         <fctr>  <char>  <char>     <char>    <char>
#> 1:   1087 COVAIL_002_CYTOKINE_2025-02-27  BUV805     CD8 CD8_BUV805     asinh
#> 2:    253 COVAIL_002_CYTOKINE_2025-02-27   RB744     CD3  CD3_RB744     asinh
#>    cofactor
#>       <num>
#> 1:     5000
#> 2:     5000

#plot and mean values of transformed columns from updated fs[['data']]
plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
#>       CD3       CD8 
#> 1.8223658 0.9800536 

#inverse transformation
flowstate.transform.inverse(fs)

#updated parameters; transform and cofactor set to NA
fs$parameters[S %in% c('CD3','CD8')]
#>       par      B DISPLAY      E        N       R      S                 TYPE
#>    <char> <char>  <char> <char>   <char>  <char> <char>               <char>
#> 1:   $P16     32     LOG    0,0 BUV805-A 4194304    CD8 Unmixed_Fluorescence
#> 2:   $P30     32     LOG    0,0  RB744-A 4194304    CD3 Unmixed_Fluorescence
#>         V                           PROJ N.alias S.alias  S_N.alias transform
#>    <char>                         <fctr>  <char>  <char>     <char>    <char>
#> 1:   1087 COVAIL_002_CYTOKINE_2025-02-27  BUV805     CD8 CD8_BUV805      <NA>
#> 2:    253 COVAIL_002_CYTOKINE_2025-02-27   RB744     CD3  CD3_RB744      <NA>
#>    cofactor
#>       <num>
#> 1:       NA
#> 2:       NA

#plot and mean values of linear columns
plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

res2.linear <- fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD8')]
print(res2.linear)
#>      CD3      CD8 
#> 25878.14 21843.70 

#linear --> transformed --> inverse --> linear
res1.linear == res2.linear
#>  CD3  CD8 
#> TRUE TRUE 
```
