# flowstate -- Structure

## Structure

[Flow Cytometry
Standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC2892967/ "FCS 3.1")
files are read and parsed into a `flowstate` (S3 object). A `flowstate`
is primarily comprised of the following parsed (named list) elements –
each a
[`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html):

| Name | Value |
|:---|:---|
| `data` | expression values |
| `parameters` | instrument-specific measurement parameters |
| `keywords` | instrument/sample-specific metadata |
| `spill` | instrument/sample-specific spillover matrix (if encoded in the file) |

Included with the package – for the purpose of example(s) – are a few
(subsampled) .fcs files acquired on a Cytek Aurora (5L) using SpectroFlo
software.

Read/parse them as a `flowstate` using
[`flowstate::read.flowstate()`](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

``` r

## paths to example .fcs files
fcs.files <- system.file("extdata", package = "flowstate") |> 
  list.files(full.names = T) |> 
  grep(pattern = "BLOCK.*.fcs", value = T)
fcs.files |> basename() |> print()
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1.fcs" "COVAIL_002_CYTOKINE_BLOCK1_2.fcs"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3.fcs"

## read them in and concatenate; a flowstate object
fs <- flowstate::read.flowstate(
  fcs.file.paths = fcs.files,
  colnames.type = "S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstates'...

## S3 object of class "flowstate"
class(fs)
#> [1] "flowstate"

## flowstate structure
names(fs)
#> [1] "data"       "parameters" "keywords"   "spill"
```

### `[['data']]`

Data can be accessed by name via: `fs$data` or `fs[['data']]`.

`[['data']]` contains events as rows and instrument measurements (Time,
scatter, fluorescence, etc.) as columns. By design, no transformation is
applied to the linear measurement values during reading. In addition, a
unique sample identifier (factored) is added for workflow purposes.

``` r

##
fs$data[]
```

``` r

## column/variable names
names(fs$data)
#>  [1] "Time"      "SSC_W"     "SSC_H"     "SSC_A"     "FSC_W"     "FSC_H"    
#>  [7] "FSC_A"     "SSCB_W"    "SSCB_H"    "SSCB_A"    "CD45RA"    "CD45RO"   
#> [13] "TCRgd"     "CD45BC1"   "IL2"       "CD8"       "CD197"     "CD45BC2"  
#> [19] "CD45BC3"   "CD57"      "CD193"     "CD45BC4"   "CD127"     "CD56"     
#> [25] "CD199"     "CD49a"     "GranzymeB" "CD45BC5"   "CD95"      "CD3"      
#> [31] "CD4"       "CD69"      "TNFa"      "CD183"     "IFNg"      "CD103"    
#> [37] "CD45BC6"   "CD122"     "ia4b7"     "CD45BC7"   "Viability" "CD27"     
#> [43] "KLRG1"     "sample.id"

## sample identifiers
fs$data[, levels(sample.id)]
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1" "COVAIL_002_CYTOKINE_BLOCK1_2"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3"
```

`[['data']]` also contains column/variable-specific attribute metadata –
notably, an `alias` attribute (named character vector) that serves the
purpose of linking syntactically-valid names with their original
(required to be unique) ‘\$PnN’ name. The aliases are used to rename
variables in `'[['data']]`, allowing for non-standard evaluation when
programming using `data.table`.

``` r

## example alias attribute
fs$data[, attr(CD3, which = 'alias')]
#>         N         S   N.alias   S.alias 
#> "RB744-A"     "CD3"   "RB744"     "CD3"
```

### `[['parameters']]`

Parameters can be accessed by name via: `fs$parameters` or
`fs[['parameters']]`.

`[['parameters']]` contains instrument-specific measurement parameters –
keyword-value pairs denoted by ‘\$P\|P’ (N, S, V, TYPE, etc.).

``` r

##
fs$parameters
```

`[['parameters']]` also contains attribute-level metadata – notably, a
`parameters.diff` attribute (list) that serves the purpose of preserving
sample-specific differences due to unique measurement values (e.g.,
differences in measured ‘Time’).

``` r

##
pd <- attr(fs$parameters, 'parameters.diff')
lapply(pd, '[[', 'R')
#> $COVAIL_002_CYTOKINE_BLOCK1_1.fcs
#>      Time 
#> "6557384" 
#> 
#> $COVAIL_002_CYTOKINE_BLOCK1_2.fcs
#>      Time 
#> "6888061" 
#> 
#> $COVAIL_002_CYTOKINE_BLOCK1_3.fcs
#>      Time 
#> "7593375"
```

### `[['keywords']]`

Keywords can be accessed by name via: `fs$keywords` or
`fs[['keywords']]`.

`[['keywords']]` contains instrument/sample-specific metadata. The
information contained here is leveraged for a number of different
actions, such as: platform/cytometer identification; extracting and
adding unique identifiers to `[['data']]`; adding/constructing new
keywords for the purpose of analysis/workflow; QC/diagnostics.

Following best practices as laid out by the Flow Cytometry Standard, a
series of keywords are added to mark any `flowstate` object(s) as being
modified – no longer source data. The inclusion of these keywords
(`'$LAST_MODIFIED'`, `'$LAST_MODIFIER'`, `'$ORIGINALITY'`) ensures
transparency and the indication that they are explicitly modified.

``` r

##
fs$keywords
```

``` r

## keywords to indicate/track modification
fs$keywords[, .(`$LAST_MODIFIED`,`$LAST_MODIFIER`,`$ORIGINALITY`)]
#>             $LAST_MODIFIED        $LAST_MODIFIER $ORIGINALITY
#>                     <char>                <char>       <char>
#> 1: 14-JUL-2026 20:18:53.51 flowstate_0.16.1.9002 DataModified
#> 2: 14-JUL-2026 20:18:53.60 flowstate_0.16.1.9002 DataModified
#> 3: 14-JUL-2026 20:18:53.69 flowstate_0.16.1.9002 DataModified
```

### `[['spill']]`

Spill can be accessed by name via: `fs$spill` or `fs[['spill']]`.

`[['spill']]` contains an instrument/sample-specific spill(over) matrix
– if encoded. It can be modified with correction values that are then
applied to `[['data']]` to (try to) fix unmixing/compensation errors.
Unless explicitly modified by a user, the default matrix should not
contain any correction values.

*N.B.: software manufacturers – at least SpectroFlo – write a small
value (1E-6) to the first matrix pair; I think this is to overcome a bug
with how FlowJo handles the matrix…*

``` r

##
fs$spill
```

``` r

## small value; there in the source data
fs$spill[2,1]
#>    CD45RA
#>     <num>
#> 1:  1e-06
```
