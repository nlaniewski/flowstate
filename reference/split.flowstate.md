# Split a concatenated `flowstate` object into a list

Split a concatenated `flowstate` object into a list

## Usage

``` r
# S3 method for class 'flowstate'
split(x, f, drop = TRUE, ...)
```

## Arguments

- x:

  A concatenated flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- f:

  Character string; a single (factored) column by which `[['data']]`
  will be split into individual `flowstate`s.

- drop:

  Logical; default `TRUE`. Unused factor levels in `f` are dropped – no
  `flowstate` object will be returned.

- ...:

  Potential further arguments.

## Value

A list of individual `flowstate`s.

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
#> Concatenating 'flowstates'...

fs.split <- split(fs,'sample.id')

length(fs.split)
#> [1] 3
names(fs.split)
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1" "COVAIL_002_CYTOKINE_BLOCK1_2"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3"
lapply(fs.split,'[[','keywords')
#> $COVAIL_002_CYTOKINE_BLOCK1_1
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 17-JUL-2026 17:49:02.55
#>      $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>              <char>      <char>       <char> <char>
#> 1: flowstate_0.17.0 aurora user DataModified     43
#>                             $PROJ $TIMESTEP   $TOT   $VOL APPLY COMPENSATION
#>                            <char>    <char> <char> <char>             <char>
#> 1: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 326.86              FALSE
#>    CHARSET          CREATOR FSC ASF    GROUPNAME
#>     <char>           <char>  <char>       <char>
#> 1:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#>                                    GUID LASER1ASF LASER1DELAY  LASER1NAME
#>                                  <char>    <char>      <char>      <char>
#> 1: 0fffcae2-2163-49f9-98a3-8f17bb4ecc13      1.14       -41.7 YellowGreen
#>    LASER2ASF LASER2DELAY LASER2NAME LASER3ASF LASER3DELAY LASER3NAME LASER4ASF
#>       <char>      <char>     <char>    <char>      <char>     <char>    <char>
#> 1:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#>    LASER4DELAY LASER4NAME LASER5ASF LASER5DELAY LASER5NAME
#>         <char>     <char>    <char>      <char>     <char>
#> 1:      20.675        Red      1.08      41.675         UV
#>                     THRESHOLD                     TUBENAME  USERSETTINGNAME
#>                        <char>                       <char>           <char>
#> 1: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_1 *COVAIL_CYTOKINE
#>    WINDOW EXTENSION                    sample.id
#>              <char>                       <fctr>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_1
#> 
#> $COVAIL_002_CYTOKINE_BLOCK1_2
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:38:12.18 Aurora  V0299 27-Feb-2025 09:49:40.92
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 17-JUL-2026 17:49:02.63
#>      $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>              <char>      <char>       <char> <char>
#> 1: flowstate_0.17.0 aurora user DataModified     43
#>                             $PROJ $TIMESTEP   $TOT   $VOL APPLY COMPENSATION
#>                            <char>    <char> <char> <char>             <char>
#> 1: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 346.31              FALSE
#>    CHARSET          CREATOR FSC ASF    GROUPNAME
#>     <char>           <char>  <char>       <char>
#> 1:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#>                                    GUID LASER1ASF LASER1DELAY  LASER1NAME
#>                                  <char>    <char>      <char>      <char>
#> 1: 6ff1ded1-96f6-4784-8a70-3cbef2e01427      1.14       -41.7 YellowGreen
#>    LASER2ASF LASER2DELAY LASER2NAME LASER3ASF LASER3DELAY LASER3NAME LASER4ASF
#>       <char>      <char>     <char>    <char>      <char>     <char>    <char>
#> 1:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#>    LASER4DELAY LASER4NAME LASER5ASF LASER5DELAY LASER5NAME
#>         <char>     <char>    <char>      <char>     <char>
#> 1:      20.675        Red      1.08      41.675         UV
#>                     THRESHOLD                     TUBENAME  USERSETTINGNAME
#>                        <char>                       <char>           <char>
#> 1: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_2 *COVAIL_CYTOKINE
#>    WINDOW EXTENSION                    sample.id
#>              <char>                       <fctr>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_2
#> 
#> $COVAIL_002_CYTOKINE_BLOCK1_3
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:50:46.94 Aurora  V0299 27-Feb-2025 10:03:26.16
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 17-JUL-2026 17:49:02.72
#>      $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>              <char>      <char>       <char> <char>
#> 1: flowstate_0.17.0 aurora user DataModified     43
#>                             $PROJ $TIMESTEP   $TOT   $VOL APPLY COMPENSATION
#>                            <char>    <char> <char> <char>             <char>
#> 1: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 350.19              FALSE
#>    CHARSET          CREATOR FSC ASF    GROUPNAME
#>     <char>           <char>  <char>       <char>
#> 1:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#>                                    GUID LASER1ASF LASER1DELAY  LASER1NAME
#>                                  <char>    <char>      <char>      <char>
#> 1: 2c9682b4-473f-4375-a851-48f9f2fe3bf7      1.14       -41.7 YellowGreen
#>    LASER2ASF LASER2DELAY LASER2NAME LASER3ASF LASER3DELAY LASER3NAME LASER4ASF
#>       <char>      <char>     <char>    <char>      <char>     <char>    <char>
#> 1:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#>    LASER4DELAY LASER4NAME LASER5ASF LASER5DELAY LASER5NAME
#>         <char>     <char>    <char>      <char>     <char>
#> 1:      20.675        Red      1.08      41.675         UV
#>                     THRESHOLD                     TUBENAME  USERSETTINGNAME
#>                        <char>                       <char>           <char>
#> 1: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_3 *COVAIL_CYTOKINE
#>    WINDOW EXTENSION                    sample.id
#>              <char>                       <fctr>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_3
#> 

#remove all rows associated with sample 1; level will remain
fs$data <- fs$data[!sample.id %in% levels(sample.id)[1]]

#with drop
fs.split <- split(fs,'sample.id',drop = TRUE)
length(fs.split)
#> [1] 2
names(fs.split)
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_2" "COVAIL_002_CYTOKINE_BLOCK1_3"

#without drop
fs.split <- split(fs,'sample.id',drop = FALSE)
length(fs.split)
#> [1] 3
names(fs.split)
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1" "COVAIL_002_CYTOKINE_BLOCK1_2"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3"
#empty data for sample 1
fs.split[[1]]$data
#> Empty data.table (0 rows and 44 cols): Time,SSC_W,SSC_H,SSC_A,FSC_W,FSC_H...
```
