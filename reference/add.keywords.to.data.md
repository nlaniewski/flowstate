# Add `[['keyword']]` values to `[['data']]`

Used for efficient addition of factored keyword values to `[['data']]`;
updates by reference.

## Usage

``` r
add.keywords.to.data(flowstate, keywords.to.add)
```

## Arguments

- flowstate:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- keywords.to.add:

  Character vector; keyword names in `[['keywords']]` whose factored
  values will be added to `[['data']]`.

## Value

UPDATES BY REFERENCE:

- `[['data']]`; factored `keywords.to.add` are added as additional
  columns.

- `[['keywords']]`; `keywords.to.add` are converted to class `factor`.

Invisibly returns `flowstate`.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read all .fcs files as flowstates; concatenate into a single object
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type = "S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstates'...

#create a new set of keywords/values
fs$keywords[
,
j = c('block.id', 'block.aliquot') := data.table::tstrsplit(
sample.id, "_", keep = 4:5)
]
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#> 2: 09:38:12.18 Aurora  V0299 27-Feb-2025 09:49:40.92
#> 3: 09:50:46.94 Aurora  V0299 27-Feb-2025 10:03:26.16
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 12-JUN-2026 18:58:02.06
#> 2: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 12-JUN-2026 18:58:02.15
#> 3: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 12-JUN-2026 18:58:02.24
#>           $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>                   <char>      <char>       <char> <char>
#> 1: flowstate_0.16.1.9000 aurora user DataModified     43
#> 2: flowstate_0.16.1.9000 aurora user DataModified     43
#> 3: flowstate_0.16.1.9000 aurora user DataModified     43
#>                             $PROJ $TIMESTEP   $TOT   $VOL APPLY COMPENSATION
#>                            <char>    <char> <char> <char>             <char>
#> 1: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 326.86              FALSE
#> 2: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 346.31              FALSE
#> 3: COVAIL_002_CYTOKINE_2025-02-27    0.0001   2000 350.19              FALSE
#>    CHARSET          CREATOR FSC ASF    GROUPNAME
#>     <char>           <char>  <char>       <char>
#> 1:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#> 2:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#> 3:   utf-8 SpectroFlo 3.3.0    1.26 BarcodedPool
#>                                    GUID LASER1ASF LASER1DELAY  LASER1NAME
#>                                  <char>    <char>      <char>      <char>
#> 1: 0fffcae2-2163-49f9-98a3-8f17bb4ecc13      1.14       -41.7 YellowGreen
#> 2: 6ff1ded1-96f6-4784-8a70-3cbef2e01427      1.14       -41.7 YellowGreen
#> 3: 2c9682b4-473f-4375-a851-48f9f2fe3bf7      1.14       -41.7 YellowGreen
#>    LASER2ASF LASER2DELAY LASER2NAME LASER3ASF LASER3DELAY LASER3NAME LASER4ASF
#>       <char>      <char>     <char>    <char>      <char>     <char>    <char>
#> 1:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#> 2:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#> 3:      1.26     -20.375     Violet      1.18           0       Blue      1.14
#>    LASER4DELAY LASER4NAME LASER5ASF LASER5DELAY LASER5NAME
#>         <char>     <char>    <char>      <char>     <char>
#> 1:      20.675        Red      1.08      41.675         UV
#> 2:      20.675        Red      1.08      41.675         UV
#> 3:      20.675        Red      1.08      41.675         UV
#>                     THRESHOLD                     TUBENAME  USERSETTINGNAME
#>                        <char>                       <char>           <char>
#> 1: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_1 *COVAIL_CYTOKINE
#> 2: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_2 *COVAIL_CYTOKINE
#> 3: (FSC,150000)And(SSC,75000) COVAIL_002_CYTOKINE_BLOCK1_3 *COVAIL_CYTOKINE
#>    WINDOW EXTENSION                    sample.id block.id block.aliquot
#>              <char>                       <fctr>   <char>        <char>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_1   BLOCK1             1
#> 2:                3 COVAIL_002_CYTOKINE_BLOCK1_2   BLOCK1             2
#> 3:                3 COVAIL_002_CYTOKINE_BLOCK1_3   BLOCK1             3
#new columns; character class
fs$keywords[, .(block.id, block.aliquot)]
#>    block.id block.aliquot
#>      <char>        <char>
#> 1:   BLOCK1             1
#> 2:   BLOCK1             2
#> 3:   BLOCK1             3

#add the factored keyword values to fs[['data']]
add.keywords.to.data(fs, c('block.id', 'block.aliquot'))

fs$data[, .(block.id, block.aliquot)]
#>       block.id block.aliquot
#>         <fctr>        <fctr>
#>    1:   BLOCK1             1
#>    2:   BLOCK1             1
#>    3:   BLOCK1             1
#>    4:   BLOCK1             1
#>    5:   BLOCK1             1
#>   ---                       
#> 5996:   BLOCK1             3
#> 5997:   BLOCK1             3
#> 5998:   BLOCK1             3
#> 5999:   BLOCK1             3
#> 6000:   BLOCK1             3

#add additional keywords to [['data']]
add.keywords.to.data(fs, c('$DATE', '$PROJ'))

fs$data[, .(`$DATE`, `$PROJ`)]
#>             $DATE                          $PROJ
#>            <fctr>                         <fctr>
#>    1: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#>    2: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#>    3: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#>    4: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#>    5: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#>   ---                                           
#> 5996: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 5997: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 5998: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 5999: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 6000: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
fs$keywords[, .(`$DATE`, `$PROJ`)]
#>          $DATE                          $PROJ
#>         <fctr>                         <fctr>
#> 1: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 2: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
#> 3: 27-Feb-2025 COVAIL_002_CYTOKINE_2025-02-27
```
