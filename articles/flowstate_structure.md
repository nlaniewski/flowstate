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

##paths to example .fcs files
fcs.files <- system.file("extdata", package = "flowstate") |> 
  list.files(full.names = T) |> 
  grep(pattern = "BLOCK.*.fcs", value = T)
fcs.files |> basename() |> print()
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1.fcs" "COVAIL_002_CYTOKINE_BLOCK1_2.fcs"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3.fcs"

##read them in and concatenate; a flowstate object
fs <- flowstate::read.flowstate(
  fcs.file.paths = fcs.files,
  colnames.type = "S",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstates'...

##S3 object of class "flowstate"
class(fs)
#> [1] "flowstate"

##flowstate structure
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
fs$data |> head()
#>     Time    SSC_W  SSC_H    SSC_A    FSC_W  FSC_H     FSC_A   SSCB_W SSCB_H
#>    <num>    <num>  <num>    <num>    <num>  <num>     <num>    <num>  <num>
#> 1:     0 730290.7 759933 924953.4 774044.6 997371 1286682.8 714389.4 873960
#> 2:     1 753766.8 453137 569266.0 714706.8 444274  529209.4 749070.4 411951
#> 3:     4 768520.5 599195 767489.4 729628.4 827274 1006004.3 730502.3 579392
#> 4:     4 707174.8 516909 609241.7 719956.1 819074  982828.9 700778.6 463030
#> 5:     5 743404.2 484976 600888.7 764667.1 754833  961993.2 726898.2 419073
#> 6:     8 724696.1 464423 560942.6 729720.8 794040  965712.5 703710.8 401969
#>       SSCB_A    CD45RA    CD45RO     TCRgd    CD45BC1         IL2         CD8
#>        <num>     <num>     <num>     <num>      <num>       <num>       <num>
#> 1: 1040579.6 17400.736  6296.233  3487.889  1056.3926 -2868.07104 219290.2812
#> 2:  514300.5  6299.763  1073.582  2427.801 21505.0137   666.11395    260.5484
#> 3:  705412.0 -1247.931 37782.598  2189.499   146.5772    14.37641   -377.2223
#> 4:  540802.5  3424.521  6541.122 10304.088 31037.8750 -1450.19849   1130.0894
#> 5:  507705.7 -1692.568 43001.797  5568.431 55248.4922  -638.63544    236.9170
#> 6:  471449.8 26003.760  3743.393  2530.474  -219.0977   135.40210   1106.0740
#>         CD197    CD45BC2   CD45BC3     CD57      CD193    CD45BC4     CD127
#>         <num>      <num>     <num>    <num>      <num>      <num>     <num>
#> 1:  1701.1161 25329.0879 20529.553 2525.726 1695.80396 23593.2891 1343.0120
#> 2: 18941.9355 14742.9941  4324.678 4230.007  162.88908   410.2411 -642.7087
#> 3:   699.3767  1116.3433 21791.977 2767.279  930.13574 30529.1172 2044.4512
#> 4:  1627.6073  -303.3068  3377.251 4663.360   84.06442 31405.3887 3629.7041
#> 5:  2753.4761  -910.3435  5274.306 7930.704  172.67607  1022.8995 5814.7285
#> 6:  1290.4441 23654.0664  3429.047 3077.539   74.12897   335.7986 1104.4130
#>          CD56    CD199     CD49a GranzymeB      CD45BC5       CD95         CD3
#>         <num>    <num>     <num>     <num>        <num>      <num>       <num>
#> 1:  6162.3535 1315.519 1071.0972 54200.477    -6.487877  942.84930 36690.11719
#> 2:   990.9473 1259.631 1404.4385  1982.248  -114.808891 2869.56348    94.28096
#> 3: -1573.5216 1891.648  512.3239  2113.995 15245.500000  922.26361 26607.62500
#> 4: -1813.6000 1058.180  217.9938  5147.355 19775.552734 -233.85599 51053.93750
#> 5:  1465.7168 3178.108 -758.9792  1842.105 16411.826172 4008.57715 39481.12891
#> 6:  1660.5831 1381.134  684.5980 30221.033   390.726471   62.77901   254.06075
#>           CD4       CD69      TNFa      CD183       IFNg      CD103    CD45BC6
#>         <num>      <num>     <num>      <num>      <num>      <num>      <num>
#> 1:  8299.1494  4737.0107  993.0143   473.1266 -256.25250  166.61444   821.8005
#> 2:   634.6076   912.6633  683.6307  -478.1163 -109.60402  -21.77036 12146.9199
#> 3: 52088.1602  1719.4734  816.6953  -402.7975  395.29520  -83.91235  1955.8750
#> 4:  5741.3486  6392.7744  544.5918  -881.3683 -160.21269 -994.53149  1735.5421
#> 5: 41355.3477  1742.8635 2127.6262 -3175.0767   27.80499  155.62851  -825.0389
#> 6:  6056.9761 31786.8047  939.8767  -854.2112   63.90987  636.24811  3377.4231
#>         CD122      ia4b7    CD45BC7  Viability       CD27      KLRG1
#>         <num>      <num>      <num>      <num>      <num>      <num>
#> 1:   735.8677  4871.1104  1030.7883  1429.1698 -7657.7598  5408.6689
#> 2:   713.9467  1420.3256  -409.5166 21949.6270   744.8879  3665.4756
#> 3: -1774.0292  -321.3192   478.2293   543.1359  -448.9456  8486.8936
#> 4:  -477.2218   487.7629  -924.2852  1028.7167   332.4514 14109.1025
#> 5:  -652.5281 -1965.0686 15285.1455   795.2957  -447.1280   729.1112
#> 6:  -126.8232 -1043.0923 18769.6777  1230.0640 -1377.2864  7417.7949
#>                       sample.id
#>                          <fctr>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1
#> 2: COVAIL_002_CYTOKINE_BLOCK1_1
#> 3: COVAIL_002_CYTOKINE_BLOCK1_1
#> 4: COVAIL_002_CYTOKINE_BLOCK1_1
#> 5: COVAIL_002_CYTOKINE_BLOCK1_1
#> 6: COVAIL_002_CYTOKINE_BLOCK1_1

##parameter/variable names
names(fs$data)
#>  [1] "Time"      "SSC_W"     "SSC_H"     "SSC_A"     "FSC_W"     "FSC_H"    
#>  [7] "FSC_A"     "SSCB_W"    "SSCB_H"    "SSCB_A"    "CD45RA"    "CD45RO"   
#> [13] "TCRgd"     "CD45BC1"   "IL2"       "CD8"       "CD197"     "CD45BC2"  
#> [19] "CD45BC3"   "CD57"      "CD193"     "CD45BC4"   "CD127"     "CD56"     
#> [25] "CD199"     "CD49a"     "GranzymeB" "CD45BC5"   "CD95"      "CD3"      
#> [31] "CD4"       "CD69"      "TNFa"      "CD183"     "IFNg"      "CD103"    
#> [37] "CD45BC6"   "CD122"     "ia4b7"     "CD45BC7"   "Viability" "CD27"     
#> [43] "KLRG1"     "sample.id"

##sample identifiers
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
keyword-value pairs denoted by ‘\$P’ (N, S, V, TYPE, etc.).

``` r

##
fs$parameters
#>        par      B DISPLAY      E                   N       R         S
#>     <char> <char>  <char> <char>              <char>  <char>    <char>
#>  1:    $P1     32     LOG    0,0                Time 7593375      <NA>
#>  2:    $P2     32     LIN    0,0               SSC-W 4194304      <NA>
#>  3:    $P3     32     LIN    0,0               SSC-H 4194304      <NA>
#>  4:    $P4     32     LIN    0,0               SSC-A 4194304      <NA>
#>  5:    $P5     32     LIN    0,0               FSC-W 4194304      <NA>
#>  6:    $P6     32     LIN    0,0               FSC-H 4194304      <NA>
#>  7:    $P7     32     LIN    0,0               FSC-A 4194304      <NA>
#>  8:    $P8     32     LIN    0,0             SSC-B-W 4194304      <NA>
#>  9:    $P9     32     LIN    0,0             SSC-B-H 4194304      <NA>
#> 10:   $P10     32     LIN    0,0             SSC-B-A 4194304      <NA>
#> 11:   $P11     32     LOG    0,0            BUV395-A 4194304    CD45RA
#> 12:   $P12     32     LOG    0,0            BUV496-A 4194304    CD45RO
#> 13:   $P13     32     LOG    0,0            BUV563-A 4194304     TCRgd
#> 14:   $P14     32     LOG    0,0           SBUV605-A 4194304   CD45BC1
#> 15:   $P15     32     LOG    0,0            BUV737-A 4194304       IL2
#> 16:   $P16     32     LOG    0,0            BUV805-A 4194304       CD8
#> 17:   $P17     32     LOG    0,0             BV421-A 4194304     CD197
#> 18:   $P18     32     LOG    0,0      Pacific Blue-A 4194304   CD45BC2
#> 19:   $P19     32     LOG    0,0            SBV475-A 4194304   CD45BC3
#> 20:   $P20     32     LOG    0,0             BV605-A 4194304      CD57
#> 21:   $P21     32     LOG    0,0             BV650-A 4194304     CD193
#> 22:   $P22     32     LOG    0,0             BV711-A 4194304   CD45BC4
#> 23:   $P23     32     LOG    0,0             BV750-A 4194304     CD127
#> 24:   $P24     32     LOG    0,0             BV785-A 4194304      CD56
#> 25:   $P25     32     LOG    0,0  KIRAVIA Blue 520-A 4194304     CD199
#> 26:   $P26     32     LOG    0,0             RB545-A 4194304     CD49a
#> 27:   $P27     32     LOG    0,0             RB613-A 4194304 GranzymeB
#> 28:   $P28     32     LOG    0,0             PerCP-A 4194304   CD45BC5
#> 29:   $P29     32     LOG    0,0  PerCP-eFluor 710-A 4194304      CD95
#> 30:   $P30     32     LOG    0,0             RB744-A 4194304       CD3
#> 31:   $P31     32     LOG    0,0    PerCP-Fire 806-A 4194304       CD4
#> 32:   $P32     32     LOG    0,0             RB780-A 4194304      CD69
#> 33:   $P33     32     LOG    0,0                PE-A 4194304      TNFa
#> 34:   $P34     32     LOG    0,0             RY586-A 4194304     CD183
#> 35:   $P35     32     LOG    0,0     PE-Dazzle 594-A 4194304      IFNg
#> 36:   $P36     32     LOG    0,0       PE-Fire 640-A 4194304     CD103
#> 37:   $P37     32     LOG    0,0            SBY720-A 4194304   CD45BC6
#> 38:   $P38     32     LOG    0,0            PE-Cy7-A 4194304     CD122
#> 39:   $P39     32     LOG    0,0   Alexa Fluor 647-A 4194304     ia4b7
#> 40:   $P40     32     LOG    0,0     Spark NIR 685-A 4194304   CD45BC7
#> 41:   $P41     32     LOG    0,0 LIVE DEAD Scarlet-A 4194304 Viability
#> 42:   $P42     32     LOG    0,0            APC-H7-A 4194304      CD27
#> 43:   $P43     32     LOG    0,0      APC-Fire 810-A 4194304     KLRG1
#>        par      B DISPLAY      E                   N       R         S
#>     <char> <char>  <char> <char>              <char>  <char>    <char>
#>                     TYPE      V
#>                   <char> <char>
#>  1:                 Time   <NA>
#>  2:         Side_Scatter    220
#>  3:         Side_Scatter    220
#>  4:         Side_Scatter    220
#>  5:      Forward_Scatter     40
#>  6:      Forward_Scatter     40
#>  7:      Forward_Scatter     40
#>  8:         Side_Scatter    220
#>  9:         Side_Scatter    220
#> 10:         Side_Scatter    220
#> 11: Unmixed_Fluorescence    258
#> 12: Unmixed_Fluorescence    523
#> 13: Unmixed_Fluorescence    535
#> 14: Unmixed_Fluorescence      0
#> 15: Unmixed_Fluorescence    904
#> 16: Unmixed_Fluorescence   1087
#> 17: Unmixed_Fluorescence    301
#> 18: Unmixed_Fluorescence    347
#> 19: Unmixed_Fluorescence      0
#> 20: Unmixed_Fluorescence    330
#> 21: Unmixed_Fluorescence    306
#> 22: Unmixed_Fluorescence    254
#> 23: Unmixed_Fluorescence    327
#> 24: Unmixed_Fluorescence    539
#> 25: Unmixed_Fluorescence      0
#> 26: Unmixed_Fluorescence    639
#> 27: Unmixed_Fluorescence    369
#> 28: Unmixed_Fluorescence    401
#> 29: Unmixed_Fluorescence    479
#> 30: Unmixed_Fluorescence    253
#> 31: Unmixed_Fluorescence      0
#> 32: Unmixed_Fluorescence    616
#> 33: Unmixed_Fluorescence    312
#> 34: Unmixed_Fluorescence    312
#> 35: Unmixed_Fluorescence      0
#> 36: Unmixed_Fluorescence    358
#> 37: Unmixed_Fluorescence      0
#> 38: Unmixed_Fluorescence    296
#> 39: Unmixed_Fluorescence    159
#> 40: Unmixed_Fluorescence    159
#> 41: Unmixed_Fluorescence      0
#> 42: Unmixed_Fluorescence    388
#> 43: Unmixed_Fluorescence    264
#>                     TYPE      V
#>                   <char> <char>

##parameters: N, S, and TYPE
fs$parameters[, .(N, S, TYPE)]
#>                       N         S                 TYPE
#>                  <char>    <char>               <char>
#>  1:                Time      <NA>                 Time
#>  2:               SSC-W      <NA>         Side_Scatter
#>  3:               SSC-H      <NA>         Side_Scatter
#>  4:               SSC-A      <NA>         Side_Scatter
#>  5:               FSC-W      <NA>      Forward_Scatter
#>  6:               FSC-H      <NA>      Forward_Scatter
#>  7:               FSC-A      <NA>      Forward_Scatter
#>  8:             SSC-B-W      <NA>         Side_Scatter
#>  9:             SSC-B-H      <NA>         Side_Scatter
#> 10:             SSC-B-A      <NA>         Side_Scatter
#> 11:            BUV395-A    CD45RA Unmixed_Fluorescence
#> 12:            BUV496-A    CD45RO Unmixed_Fluorescence
#> 13:            BUV563-A     TCRgd Unmixed_Fluorescence
#> 14:           SBUV605-A   CD45BC1 Unmixed_Fluorescence
#> 15:            BUV737-A       IL2 Unmixed_Fluorescence
#> 16:            BUV805-A       CD8 Unmixed_Fluorescence
#> 17:             BV421-A     CD197 Unmixed_Fluorescence
#> 18:      Pacific Blue-A   CD45BC2 Unmixed_Fluorescence
#> 19:            SBV475-A   CD45BC3 Unmixed_Fluorescence
#> 20:             BV605-A      CD57 Unmixed_Fluorescence
#> 21:             BV650-A     CD193 Unmixed_Fluorescence
#> 22:             BV711-A   CD45BC4 Unmixed_Fluorescence
#> 23:             BV750-A     CD127 Unmixed_Fluorescence
#> 24:             BV785-A      CD56 Unmixed_Fluorescence
#> 25:  KIRAVIA Blue 520-A     CD199 Unmixed_Fluorescence
#> 26:             RB545-A     CD49a Unmixed_Fluorescence
#> 27:             RB613-A GranzymeB Unmixed_Fluorescence
#> 28:             PerCP-A   CD45BC5 Unmixed_Fluorescence
#> 29:  PerCP-eFluor 710-A      CD95 Unmixed_Fluorescence
#> 30:             RB744-A       CD3 Unmixed_Fluorescence
#> 31:    PerCP-Fire 806-A       CD4 Unmixed_Fluorescence
#> 32:             RB780-A      CD69 Unmixed_Fluorescence
#> 33:                PE-A      TNFa Unmixed_Fluorescence
#> 34:             RY586-A     CD183 Unmixed_Fluorescence
#> 35:     PE-Dazzle 594-A      IFNg Unmixed_Fluorescence
#> 36:       PE-Fire 640-A     CD103 Unmixed_Fluorescence
#> 37:            SBY720-A   CD45BC6 Unmixed_Fluorescence
#> 38:            PE-Cy7-A     CD122 Unmixed_Fluorescence
#> 39:   Alexa Fluor 647-A     ia4b7 Unmixed_Fluorescence
#> 40:     Spark NIR 685-A   CD45BC7 Unmixed_Fluorescence
#> 41: LIVE DEAD Scarlet-A Viability Unmixed_Fluorescence
#> 42:            APC-H7-A      CD27 Unmixed_Fluorescence
#> 43:      APC-Fire 810-A     KLRG1 Unmixed_Fluorescence
#>                       N         S                 TYPE
#>                  <char>    <char>               <char>
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
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#> 2: 09:38:12.18 Aurora  V0299 27-Feb-2025 09:49:40.92
#> 3: 09:50:46.94 Aurora  V0299 27-Feb-2025 10:03:26.16
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 12-JUN-2026 19:21:45.63
#> 2: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 12-JUN-2026 19:21:45.73
#> 3: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 12-JUN-2026 19:21:45.81
#>           $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>                   <char>      <char>       <char> <char>
#> 1: flowstate_0.16.1.9001 aurora user DataModified     43
#> 2: flowstate_0.16.1.9001 aurora user DataModified     43
#> 3: flowstate_0.16.1.9001 aurora user DataModified     43
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
#>    WINDOW EXTENSION                    sample.id
#>              <char>                       <fctr>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_1
#> 2:                3 COVAIL_002_CYTOKINE_BLOCK1_2
#> 3:                3 COVAIL_002_CYTOKINE_BLOCK1_3

##some identifiers
fs$keywords[,.(`$CYT`,TUBENAME)]
#>      $CYT                     TUBENAME
#>    <char>                       <char>
#> 1: Aurora COVAIL_002_CYTOKINE_BLOCK1_1
#> 2: Aurora COVAIL_002_CYTOKINE_BLOCK1_2
#> 3: Aurora COVAIL_002_CYTOKINE_BLOCK1_3

##keywords to indicate/track modification
fs$keywords[,.(`$LAST_MODIFIED`,`$LAST_MODIFIER`,`$ORIGINALITY`)]
#>             $LAST_MODIFIED        $LAST_MODIFIER $ORIGINALITY
#>                     <char>                <char>       <char>
#> 1: 12-JUN-2026 19:21:45.63 flowstate_0.16.1.9001 DataModified
#> 2: 12-JUN-2026 19:21:45.73 flowstate_0.16.1.9001 DataModified
#> 3: 12-JUN-2026 19:21:45.81 flowstate_0.16.1.9001 DataModified
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
fs$spill |> head()
#>    CD45RA CD45RO TCRgd CD45BC1   IL2   CD8 CD197 CD45BC2 CD45BC3  CD57 CD193
#>     <num>  <num> <num>   <num> <num> <num> <num>   <num>   <num> <num> <num>
#> 1:  1e+00      0     0       0     0     0     0       0       0     0     0
#> 2:  1e-06      1     0       0     0     0     0       0       0     0     0
#> 3:  0e+00      0     1       0     0     0     0       0       0     0     0
#> 4:  0e+00      0     0       1     0     0     0       0       0     0     0
#> 5:  0e+00      0     0       0     1     0     0       0       0     0     0
#> 6:  0e+00      0     0       0     0     1     0       0       0     0     0
#>    CD45BC4 CD127  CD56 CD199 CD49a GranzymeB CD45BC5  CD95   CD3   CD4  CD69
#>      <num> <num> <num> <num> <num>     <num>   <num> <num> <num> <num> <num>
#> 1:       0     0     0     0     0         0       0     0     0     0     0
#> 2:       0     0     0     0     0         0       0     0     0     0     0
#> 3:       0     0     0     0     0         0       0     0     0     0     0
#> 4:       0     0     0     0     0         0       0     0     0     0     0
#> 5:       0     0     0     0     0         0       0     0     0     0     0
#> 6:       0     0     0     0     0         0       0     0     0     0     0
#>     TNFa CD183  IFNg CD103 CD45BC6 CD122 ia4b7 CD45BC7 Viability  CD27 KLRG1
#>    <num> <num> <num> <num>   <num> <num> <num>   <num>     <num> <num> <num>
#> 1:     0     0     0     0       0     0     0       0         0     0     0
#> 2:     0     0     0     0       0     0     0       0         0     0     0
#> 3:     0     0     0     0       0     0     0       0         0     0     0
#> 4:     0     0     0     0       0     0     0       0         0     0     0
#> 5:     0     0     0     0       0     0     0       0         0     0     0
#> 6:     0     0     0     0       0     0     0       0         0     0     0

##small value; there in the source data
fs$spill[2,1]
#>    CD45RA
#>     <num>
#> 1:  1e-06
```
