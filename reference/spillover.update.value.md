# Update the values of a flowstate spill data.table

Update the values of a flowstate spill data.table

## Usage

``` r
spillover.update.value(flowstate.object, i, j, value)
```

## Arguments

- flowstate.object:

  the return of
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- i:

  Variable (unquoted).

- j:

  Variable (unquoted).

- value:

  Numeric; correction value to be used during compensation.

## Value

Updated `flowstate$spill`; !!!updates by reference – no assignment.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read all .fcs files as flowstate objects; concatenate into a single object
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type="S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstate.ojects'...

#row index; CD4 vs CD8
index<-which(names(fs$spill) %in% c('CD4','CD8'))
fs$spill[index,.(CD4,CD8)]
#>      CD4   CD8
#>    <num> <num>
#> 1:     0     1
#> 2:     1     0

#update a spill value
spillover.update.value(fs,CD8,CD4,0.03)

fs$spill[index,.(CD4,CD8)]
#>      CD4   CD8
#>    <num> <num>
#> 1:  0.03     1
#> 2:  1.00     0
```
