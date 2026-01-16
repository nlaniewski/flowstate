# Transform `flowstate[['data']]`

Transform `flowstate[['data']]`

## Usage

``` r
flowstate.transform(
  flowstate.object,
  .j = NULL,
  transform.type = "asinh",
  cofactor = 5000
)
```

## Arguments

- flowstate.object:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- .j:

  Character vector – default `NULL`; any/all parameters having a
  keyword-value pair of `TYPE/Raw|Unmixed_Fluorescence` will be
  transformed in `[['data']]`. If a character vector: specific columns
  in `[['data']]` that are to be transformed.

- transform.type:

  Character string – default
  [asinh](https://rdrr.io/r/base/Hyperbolic.html); quoted function that
  will be used to transform `.j` in `[['data']]`.

- cofactor:

  Numeric – default 5000; if `transform.type` =
  [asinh](https://rdrr.io/r/base/Hyperbolic.html) (default), the
  cofactor will be used to modify the transformation.

## Value

UPDATES BY REFERENCE:

- `flowstate[['data']]`; transformed values

- `flowstate[['parameters']]`; adds two columns ('transform' and
  'cofactor')

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

fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD4','CD8')]
#>      CD3      CD4      CD8 
#> 25878.14 19921.68 21843.70 

#transform
flowstate.transform(
  fs,
  .j = c('CD3','CD4','CD8'),
  transform.type = "asinh",
  cofactor = 5000
)
#> flowstate.object --> transforming...
#updated parameters
fs$parameters[!is.na(transform)]
#>       par      B DISPLAY      E                N       R      S
#>    <char> <char>  <char> <char>           <char>  <char> <char>
#> 1:   $P16     32     LOG    0,0         BUV805-A 4194304    CD8
#> 2:   $P30     32     LOG    0,0          RB744-A 4194304    CD3
#> 3:   $P31     32     LOG    0,0 PerCP-Fire 806-A 4194304    CD4
#>                    TYPE      V                           PROJ      N.alias
#>                  <char> <char>                         <fctr>       <char>
#> 1: Unmixed_Fluorescence   1087 COVAIL_002_CYTOKINE_2025-02-27       BUV805
#> 2: Unmixed_Fluorescence    253 COVAIL_002_CYTOKINE_2025-02-27        RB744
#> 3: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27 PerCPFire806
#>    S.alias        S_N.alias transform cofactor
#>     <char>           <char>    <char>    <num>
#> 1:     CD8       CD8_BUV805     asinh     5000
#> 2:     CD3        CD3_RB744     asinh     5000
#> 3:     CD4 CD4_PerCPFire806     asinh     5000

#plot and mean values of transformed columns from updated fs[['data']]
plot(fs,CD3,CD8) + ggplot2::guides(fill = 'none')
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.

fs$data[,sapply(.SD,mean),.SDcols = c('CD3','CD4','CD8')]
#>       CD3       CD4       CD8 
#> 1.8223658 1.6040152 0.9800536 

#transform all
flowstate.transform(
  fs,
  .j = NULL,
  transform.type = "asinh",
  cofactor = 5000
)
#> flowstate.object --> transforming...

#updated parameters
fs$parameters[!is.na(transform)]
#>        par      B DISPLAY      E                   N       R         S
#>     <char> <char>  <char> <char>              <char>  <char>    <char>
#>  1:   $P11     32     LOG    0,0            BUV395-A 4194304    CD45RA
#>  2:   $P12     32     LOG    0,0            BUV496-A 4194304    CD45RO
#>  3:   $P13     32     LOG    0,0            BUV563-A 4194304     TCRgd
#>  4:   $P14     32     LOG    0,0           SBUV605-A 4194304   CD45BC1
#>  5:   $P15     32     LOG    0,0            BUV737-A 4194304       IL2
#>  6:   $P16     32     LOG    0,0            BUV805-A 4194304       CD8
#>  7:   $P17     32     LOG    0,0             BV421-A 4194304     CD197
#>  8:   $P18     32     LOG    0,0      Pacific Blue-A 4194304   CD45BC2
#>  9:   $P19     32     LOG    0,0            SBV475-A 4194304   CD45BC3
#> 10:   $P20     32     LOG    0,0             BV605-A 4194304      CD57
#> 11:   $P21     32     LOG    0,0             BV650-A 4194304     CD193
#> 12:   $P22     32     LOG    0,0             BV711-A 4194304   CD45BC4
#> 13:   $P23     32     LOG    0,0             BV750-A 4194304     CD127
#> 14:   $P24     32     LOG    0,0             BV785-A 4194304      CD56
#> 15:   $P25     32     LOG    0,0  KIRAVIA Blue 520-A 4194304     CD199
#> 16:   $P26     32     LOG    0,0             RB545-A 4194304     CD49a
#> 17:   $P27     32     LOG    0,0             RB613-A 4194304 GranzymeB
#> 18:   $P28     32     LOG    0,0             PerCP-A 4194304   CD45BC5
#> 19:   $P29     32     LOG    0,0  PerCP-eFluor 710-A 4194304      CD95
#> 20:   $P30     32     LOG    0,0             RB744-A 4194304       CD3
#> 21:   $P31     32     LOG    0,0    PerCP-Fire 806-A 4194304       CD4
#> 22:   $P32     32     LOG    0,0             RB780-A 4194304      CD69
#> 23:   $P33     32     LOG    0,0                PE-A 4194304      TNFa
#> 24:   $P34     32     LOG    0,0             RY586-A 4194304     CD183
#> 25:   $P35     32     LOG    0,0     PE-Dazzle 594-A 4194304      IFNg
#> 26:   $P36     32     LOG    0,0       PE-Fire 640-A 4194304     CD103
#> 27:   $P37     32     LOG    0,0            SBY720-A 4194304   CD45BC6
#> 28:   $P38     32     LOG    0,0            PE-Cy7-A 4194304     CD122
#> 29:   $P39     32     LOG    0,0   Alexa Fluor 647-A 4194304     ia4b7
#> 30:   $P40     32     LOG    0,0     Spark NIR 685-A 4194304   CD45BC7
#> 31:   $P41     32     LOG    0,0 LIVE DEAD Scarlet-A 4194304 Viability
#> 32:   $P42     32     LOG    0,0            APC-H7-A 4194304      CD27
#> 33:   $P43     32     LOG    0,0      APC-Fire 810-A 4194304     KLRG1
#>        par      B DISPLAY      E                   N       R         S
#>     <char> <char>  <char> <char>              <char>  <char>    <char>
#>                     TYPE      V                           PROJ         N.alias
#>                   <char> <char>                         <fctr>          <char>
#>  1: Unmixed_Fluorescence    258 COVAIL_002_CYTOKINE_2025-02-27          BUV395
#>  2: Unmixed_Fluorescence    523 COVAIL_002_CYTOKINE_2025-02-27          BUV496
#>  3: Unmixed_Fluorescence    535 COVAIL_002_CYTOKINE_2025-02-27          BUV563
#>  4: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27         SBUV605
#>  5: Unmixed_Fluorescence    904 COVAIL_002_CYTOKINE_2025-02-27          BUV737
#>  6: Unmixed_Fluorescence   1087 COVAIL_002_CYTOKINE_2025-02-27          BUV805
#>  7: Unmixed_Fluorescence    301 COVAIL_002_CYTOKINE_2025-02-27           BV421
#>  8: Unmixed_Fluorescence    347 COVAIL_002_CYTOKINE_2025-02-27     PacificBlue
#>  9: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27          SBV475
#> 10: Unmixed_Fluorescence    330 COVAIL_002_CYTOKINE_2025-02-27           BV605
#> 11: Unmixed_Fluorescence    306 COVAIL_002_CYTOKINE_2025-02-27           BV650
#> 12: Unmixed_Fluorescence    254 COVAIL_002_CYTOKINE_2025-02-27           BV711
#> 13: Unmixed_Fluorescence    327 COVAIL_002_CYTOKINE_2025-02-27           BV750
#> 14: Unmixed_Fluorescence    539 COVAIL_002_CYTOKINE_2025-02-27           BV785
#> 15: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27  KIRAVIABlue520
#> 16: Unmixed_Fluorescence    639 COVAIL_002_CYTOKINE_2025-02-27           RB545
#> 17: Unmixed_Fluorescence    369 COVAIL_002_CYTOKINE_2025-02-27           RB613
#> 18: Unmixed_Fluorescence    401 COVAIL_002_CYTOKINE_2025-02-27           PerCP
#> 19: Unmixed_Fluorescence    479 COVAIL_002_CYTOKINE_2025-02-27  PerCPeFluor710
#> 20: Unmixed_Fluorescence    253 COVAIL_002_CYTOKINE_2025-02-27           RB744
#> 21: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27    PerCPFire806
#> 22: Unmixed_Fluorescence    616 COVAIL_002_CYTOKINE_2025-02-27           RB780
#> 23: Unmixed_Fluorescence    312 COVAIL_002_CYTOKINE_2025-02-27              PE
#> 24: Unmixed_Fluorescence    312 COVAIL_002_CYTOKINE_2025-02-27           RY586
#> 25: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27     PEDazzle594
#> 26: Unmixed_Fluorescence    358 COVAIL_002_CYTOKINE_2025-02-27       PEFire640
#> 27: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27          SBY720
#> 28: Unmixed_Fluorescence    296 COVAIL_002_CYTOKINE_2025-02-27           PECy7
#> 29: Unmixed_Fluorescence    159 COVAIL_002_CYTOKINE_2025-02-27   AlexaFluor647
#> 30: Unmixed_Fluorescence    159 COVAIL_002_CYTOKINE_2025-02-27     SparkNIR685
#> 31: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27 LIVEDEADScarlet
#> 32: Unmixed_Fluorescence    388 COVAIL_002_CYTOKINE_2025-02-27           APCH7
#> 33: Unmixed_Fluorescence    264 COVAIL_002_CYTOKINE_2025-02-27      APCFire810
#>                     TYPE      V                           PROJ         N.alias
#>                   <char> <char>                         <fctr>          <char>
#>       S.alias                 S_N.alias transform cofactor
#>        <char>                    <char>    <char>    <num>
#>  1:    CD45RA             CD45RA_BUV395     asinh     5000
#>  2:    CD45RO             CD45RO_BUV496     asinh     5000
#>  3:     TCRgd              TCRgd_BUV563     asinh     5000
#>  4:   CD45BC1           CD45BC1_SBUV605     asinh     5000
#>  5:       IL2                IL2_BUV737     asinh     5000
#>  6:       CD8                CD8_BUV805     asinh     5000
#>  7:     CD197               CD197_BV421     asinh     5000
#>  8:   CD45BC2       CD45BC2_PacificBlue     asinh     5000
#>  9:   CD45BC3            CD45BC3_SBV475     asinh     5000
#> 10:      CD57                CD57_BV605     asinh     5000
#> 11:     CD193               CD193_BV650     asinh     5000
#> 12:   CD45BC4             CD45BC4_BV711     asinh     5000
#> 13:     CD127               CD127_BV750     asinh     5000
#> 14:      CD56                CD56_BV785     asinh     5000
#> 15:     CD199      CD199_KIRAVIABlue520     asinh     5000
#> 16:     CD49a               CD49a_RB545     asinh     5000
#> 17: GranzymeB           GranzymeB_RB613     asinh     5000
#> 18:   CD45BC5             CD45BC5_PerCP     asinh     5000
#> 19:      CD95       CD95_PerCPeFluor710     asinh     5000
#> 20:       CD3                 CD3_RB744     asinh     5000
#> 21:       CD4          CD4_PerCPFire806     asinh     5000
#> 22:      CD69                CD69_RB780     asinh     5000
#> 23:      TNFa                   TNFa_PE     asinh     5000
#> 24:     CD183               CD183_RY586     asinh     5000
#> 25:      IFNg          IFNg_PEDazzle594     asinh     5000
#> 26:     CD103           CD103_PEFire640     asinh     5000
#> 27:   CD45BC6            CD45BC6_SBY720     asinh     5000
#> 28:     CD122               CD122_PECy7     asinh     5000
#> 29:     ia4b7       ia4b7_AlexaFluor647     asinh     5000
#> 30:   CD45BC7       CD45BC7_SparkNIR685     asinh     5000
#> 31: Viability Viability_LIVEDEADScarlet     asinh     5000
#> 32:      CD27                CD27_APCH7     asinh     5000
#> 33:     KLRG1          KLRG1_APCFire810     asinh     5000
#>       S.alias                 S_N.alias transform cofactor
#>        <char>                    <char>    <char>    <num>
```
