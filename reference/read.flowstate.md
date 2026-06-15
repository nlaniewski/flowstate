# `flowstate`: read, parse, and store .fcs data

[Flow Cytometry
Standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC2892967/) files are
read, parsed and stored as objects (S3) of
[class](https://rdrr.io/r/base/class.html) `'flowstate'`.

The individual segments (`TEXT` and `DATA`) of the FCS file are parsed
and stored as follows:

- `[['data']]` – a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing raw/linear measurement values (scatter, MFI, Time, etc.).

- `[['parameters']]` – a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing instrument-specific parameters ('\$PnN','\$PnS', etc.).

- `[['keywords']]` – a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing instrument/sample-specific keyword-value pairs (metadata).

- `[['spill']]` – a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing (if present) spillover values.

## Usage

``` r
read.flowstate(
  fcs.file.paths,
  colnames.type = c("N", "S"),
  sample.id = NULL,
  concatenate = FALSE
)
```

## Arguments

- fcs.file.paths:

  Character string; path(s) returned from
  `list.files(...,full.names=T,pattern=".fcs")`.

- colnames.type:

  Character string; one of:

  - `"N"` – `[['data']]` columns are named by using only their
    respective \$PN (name) keyword value.

  - `"S"` – `[['data']]` columns are named by using only their
    respective \$PS (stain/user-defined) keyword value.

- sample.id:

  Keyword name – default `NULL`; based on cytometer platform, the
  keyword name for `sample.id` will be automatically set and the
  respective keyword values will be added to `[['data']]` as a factored
  sample identifier. One of:

  - Aurora (Cytek) – `sample.id` = `'TUBENAME'`

  - ID7000 (Sony) – `sample.id` = `'$CELLS'`

  - Unspecified – `sample.id` = `'$FIL'`

- concatenate:

  Logical – default `FALSE`; if `TRUE`, the list of flowstates will be
  combined into a single flowstate.

## Value

For a single file: an object of class `flowstate`; for multiple files: a
named list of `flowstates`; for concatenated files: an object of class
`flowstate`

## Details

Access the individual named list elements as follows (assuming
`flowstate` object is named `fs`):

- `fs$data` or `fs[['data']]`

- `fs$parameters` or `fs[['parameters']]`

- `fs$keywords` or `fs[['keywords']]`

- `fs$spill` or `fs[['spill']]`

Other `flowstate` functions operate on the entire object (`fs`) and
access specific list elements as needed.

## References

Directly/heavily-inspired by:

flowCore:

Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J,
Jiang M, Finak G (2024). *flowCore: flowCore: Basic structures for flow
cytometry data*. doi:10.18129/B9.bioc.flowCore
<https://doi.org/10.18129/B9.bioc.flowCore>, R package version 2.22.0,
<https://bioconductor.org/packages/flowCore>.

[data.table](https://rdrr.io/pkg/data.table/man/data.table.html):

Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T,
Schwendinger B, Krylov I (2025). data.table: Extension of 'data.frame'.
R package version 1.17.99, https://r-datatable.com.

## See also

[select_nonsaturating](https://nlaniewski.github.io/flowstate/reference/select_nonsaturating.md)
;
[flowstate.transform](https://nlaniewski.github.io/flowstate/reference/flowstate.transform.md)

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read a single .fcs file as a flowstate
fs <- read.flowstate(
  fcs.file.paths[1],
  colnames.type = "S"
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
class(fs)
#> [1] "flowstate"
names(fs)
#> [1] "data"       "parameters" "keywords"   "spill"     

#.fcs DATA segment as a data.table
  fs$data
#>        Time    SSC_W  SSC_H    SSC_A    FSC_W   FSC_H     FSC_A   SSCB_W SSCB_H
#>       <num>    <num>  <num>    <num>    <num>   <num>     <num>    <num>  <num>
#>    1:     0 730290.7 759933 924953.4 774044.6  997371 1286682.8 714389.4 873960
#>    2:     1 753766.8 453137 569266.0 714706.8  444274  529209.4 749070.4 411951
#>    3:     4 768520.5 599195 767489.4 729628.4  827274 1006004.3 730502.3 579392
#>    4:     4 707174.8 516909 609241.7 719956.1  819074  982828.9 700778.6 463030
#>    5:     5 743404.2 484976 600888.7 764667.1  754833  961993.2 726898.2 419073
#>   ---                                                                          
#> 1996:  2436 665170.4 766756 850038.9 735662.5  905841 1110655.4 651364.2 621187
#> 1997:  2443 755763.9 435045 547985.5 718318.1  792184  948400.2 725049.6 377339
#> 1998:  2444 751025.3 707866 886042.2 756734.6 1088113 1372354.5 748516.9 675462
#> 1999:  2448 764364.1 551760 702909.2 734433.4  471781  577486.2 760942.9 380466
#> 2000:  2448 744009.3 636641 789444.8 744564.4  825535 1024440.0 725390.6 574194
#>          SSCB_A     CD45RA    CD45RO     TCRgd     CD45BC1         IL2
#>           <num>      <num>     <num>     <num>       <num>       <num>
#>    1: 1040579.6 17400.7363  6296.233  3487.889  1056.39258 -2868.07104
#>    2:  514300.5  6299.7627  1073.582  2427.801 21505.01367   666.11395
#>    3:  705412.0 -1247.9305 37782.598  2189.499   146.57719    14.37641
#>    4:  540802.5  3424.5208  6541.122 10304.088 31037.87500 -1450.19849
#>    5:  507705.7 -1692.5685 43001.797  5568.431 55248.49219  -638.63544
#>   ---                                                                 
#> 1996:  674364.9  -455.5892 18938.584  2339.671 30702.53320  -650.51361
#> 1997:  455982.5 45408.0547  2361.198  3722.122 57515.62891   376.72684
#> 1998:  842657.9 16636.8730  8177.999  3669.374   -37.85068   -56.72515
#> 1999:  482521.5   355.6579  3132.339  3150.341 89416.53906  -789.44360
#> 2000:  694191.5 16607.0723  4449.304  2556.906 15862.92090  -348.67285
#>               CD8      CD197    CD45BC2   CD45BC3      CD57       CD193
#>             <num>      <num>      <num>     <num>     <num>       <num>
#>    1: 219290.2812  1701.1161 25329.0879 20529.553  2525.726  1695.80396
#>    2:    260.5484 18941.9355 14742.9941  4324.678  4230.007   162.88908
#>    3:   -377.2223   699.3767  1116.3433 21791.977  2767.279   930.13574
#>    4:   1130.0894  1627.6073  -303.3068  3377.251  4663.360    84.06442
#>    5:    236.9170  2753.4761  -910.3435  5274.306  7930.704   172.67607
#>   ---                                                                  
#> 1996:  59250.2812  1092.6920 15200.3584  3081.160  3450.271  -315.22495
#> 1997:    659.7032  4738.8994   405.5616  3159.196  4285.967   470.22446
#> 1998:  41272.5664  2141.1841 21296.8184  7120.216 32015.932  -755.86090
#> 1999:   -146.6094   571.7484  1858.9985 44696.523  9054.923 -2254.48193
#> 2000:  12785.9922  1001.4847 46748.2773  3473.425  9044.300   961.13348
#>          CD45BC4      CD127       CD56    CD199     CD49a GranzymeB
#>            <num>      <num>      <num>    <num>     <num>     <num>
#>    1: 23593.2891 1343.01196  6162.3535 1315.519 1071.0972 54200.477
#>    2:   410.2411 -642.70868   990.9473 1259.631 1404.4385  1982.248
#>    3: 30529.1172 2044.45117 -1573.5216 1891.648  512.3239  2113.995
#>    4: 31405.3887 3629.70410 -1813.6000 1058.180  217.9938  5147.355
#>    5:  1022.8995 5814.72852  1465.7168 3178.108 -758.9792  1842.105
#>   ---                                                              
#> 1996:  1064.3235 2493.20142 -1785.5209 1296.361  464.7216  1550.618
#> 1997: 49806.3789 -846.77362  3121.0767 1310.188  538.4861  1325.202
#> 1998:  2133.5469  -69.45152  6356.4126 1662.545 1216.8611 52676.543
#> 1999:   840.9261  888.50812  1614.4645 3131.573 -186.4928  2652.448
#> 2000:  1132.7102 -408.14355  3265.9004 1530.620  893.4030 33198.285
#>            CD45BC5       CD95         CD3        CD4       CD69      TNFa
#>              <num>      <num>       <num>      <num>      <num>     <num>
#>    1:    -6.487877  942.84930 36690.11719  8299.1494  4737.0107  993.0143
#>    2:  -114.808891 2869.56348    94.28096   634.6076   912.6633  683.6307
#>    3: 15245.500000  922.26361 26607.62500 52088.1602  1719.4734  816.6953
#>    4: 19775.552734 -233.85599 51053.93750  5741.3486  6392.7744  544.5918
#>    5: 16411.826172 4008.57715 39481.12891 41355.3477  1742.8635 2127.6262
#>   ---                                                                    
#> 1996: 10101.198242  208.18390 32280.44922  8475.2119 15024.1406  274.7441
#> 1997:   107.310089  348.53680    91.05525  2211.5044  7191.9580 -452.7868
#> 1998:    58.550423 -371.25952 22753.05859  4531.7686  4060.4192  153.8195
#> 1999:   964.303589  774.55017 27600.72266   450.5056  3263.5374 5513.1948
#> 2000: 24346.623047   32.25357 19126.30273  6096.0869 25395.6855  870.2074
#>            CD183       IFNg      CD103    CD45BC6       CD122      ia4b7
#>            <num>      <num>      <num>      <num>       <num>      <num>
#>    1:   473.1266 -256.25250  166.61444   821.8005   735.86768  4871.1104
#>    2:  -478.1163 -109.60402  -21.77036 12146.9199   713.94666  1420.3256
#>    3:  -402.7975  395.29520  -83.91235  1955.8750 -1774.02917  -321.3192
#>    4:  -881.3683 -160.21269 -994.53149  1735.5421  -477.22177   487.7629
#>    5: -3175.0767   27.80499  155.62851  -825.0389  -652.52808 -1965.0686
#>   ---                                                                   
#> 1996:  -248.7266  431.12689  216.10272   735.4802 -1338.54907  -862.8319
#> 1997:  -493.0250  352.03424  -72.76740   602.4594  -527.90393  -588.9193
#> 1998:  -405.9241 -639.24402 -151.01843  4389.2275   149.43155   516.9016
#> 1999: -1634.7971 -859.71637 -383.70520  -329.5938   516.64771  -112.7745
#> 2000:  -791.4459 -180.43523 -951.55237   821.8984   -47.60364  -353.7823
#>          CD45BC7   Viability        CD27      KLRG1
#>            <num>       <num>       <num>      <num>
#>    1:  1030.7883  1429.16980 -7657.75977  5408.6689
#>    2:  -409.5166 21949.62695   744.88788  3665.4756
#>    3:   478.2293   543.13586  -448.94562  8486.8936
#>    4:  -924.2852  1028.71667   332.45135 14109.1025
#>    5: 15285.1455   795.29565  -447.12799   729.1112
#>   ---                                              
#> 1996:  1768.1182   188.79082 -1291.30835 10857.7451
#> 1997: 27534.8418   863.10486   671.50604  1228.5106
#> 1998:  1421.7739    73.27693 -3152.55737 11221.6455
#> 1999: 10256.8574 12299.75000 -2231.92310  8055.1738
#> 2000:   191.6911   850.90717   -64.23028   130.7855
#>                          sample.id
#>                             <fctr>
#>    1: COVAIL_002_CYTOKINE_BLOCK1_1
#>    2: COVAIL_002_CYTOKINE_BLOCK1_1
#>    3: COVAIL_002_CYTOKINE_BLOCK1_1
#>    4: COVAIL_002_CYTOKINE_BLOCK1_1
#>    5: COVAIL_002_CYTOKINE_BLOCK1_1
#>   ---                             
#> 1996: COVAIL_002_CYTOKINE_BLOCK1_1
#> 1997: COVAIL_002_CYTOKINE_BLOCK1_1
#> 1998: COVAIL_002_CYTOKINE_BLOCK1_1
#> 1999: COVAIL_002_CYTOKINE_BLOCK1_1
#> 2000: COVAIL_002_CYTOKINE_BLOCK1_1
#.fcs TEXT segment parsed and stored as three elements (data.tables):
  fs$parameters #instrument-specific measurement parameters
#>        par      B DISPLAY      E                   N       R         S
#>     <char> <char>  <char> <char>              <char>  <char>    <char>
#>  1:    $P1     32     LOG    0,0                Time 6557384      <NA>
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
  fs$keywords #instrument/sample-specific metadata
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 15-JUN-2026 22:35:54.06
#>           $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>                   <char>      <char>       <char> <char>
#> 1: flowstate_0.16.1.9001 aurora user DataModified     43
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
  fs$spill #instrument/sample-specific spillover
#>     CD45RA CD45RO TCRgd CD45BC1   IL2   CD8 CD197 CD45BC2 CD45BC3  CD57 CD193
#>      <num>  <num> <num>   <num> <num> <num> <num>   <num>   <num> <num> <num>
#>  1:  1e+00      0     0       0     0     0     0       0       0     0     0
#>  2:  1e-06      1     0       0     0     0     0       0       0     0     0
#>  3:  0e+00      0     1       0     0     0     0       0       0     0     0
#>  4:  0e+00      0     0       1     0     0     0       0       0     0     0
#>  5:  0e+00      0     0       0     1     0     0       0       0     0     0
#>  6:  0e+00      0     0       0     0     1     0       0       0     0     0
#>  7:  0e+00      0     0       0     0     0     1       0       0     0     0
#>  8:  0e+00      0     0       0     0     0     0       1       0     0     0
#>  9:  0e+00      0     0       0     0     0     0       0       1     0     0
#> 10:  0e+00      0     0       0     0     0     0       0       0     1     0
#> 11:  0e+00      0     0       0     0     0     0       0       0     0     1
#> 12:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 13:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 14:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 15:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 16:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 17:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 18:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 19:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 20:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 21:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 22:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 23:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 24:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 25:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 26:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 27:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 28:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 29:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 30:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 31:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 32:  0e+00      0     0       0     0     0     0       0       0     0     0
#> 33:  0e+00      0     0       0     0     0     0       0       0     0     0
#>     CD45RA CD45RO TCRgd CD45BC1   IL2   CD8 CD197 CD45BC2 CD45BC3  CD57 CD193
#>      <num>  <num> <num>   <num> <num> <num> <num>   <num>   <num> <num> <num>
#>     CD45BC4 CD127  CD56 CD199 CD49a GranzymeB CD45BC5  CD95   CD3   CD4  CD69
#>       <num> <num> <num> <num> <num>     <num>   <num> <num> <num> <num> <num>
#>  1:       0     0     0     0     0         0       0     0     0     0     0
#>  2:       0     0     0     0     0         0       0     0     0     0     0
#>  3:       0     0     0     0     0         0       0     0     0     0     0
#>  4:       0     0     0     0     0         0       0     0     0     0     0
#>  5:       0     0     0     0     0         0       0     0     0     0     0
#>  6:       0     0     0     0     0         0       0     0     0     0     0
#>  7:       0     0     0     0     0         0       0     0     0     0     0
#>  8:       0     0     0     0     0         0       0     0     0     0     0
#>  9:       0     0     0     0     0         0       0     0     0     0     0
#> 10:       0     0     0     0     0         0       0     0     0     0     0
#> 11:       0     0     0     0     0         0       0     0     0     0     0
#> 12:       1     0     0     0     0         0       0     0     0     0     0
#> 13:       0     1     0     0     0         0       0     0     0     0     0
#> 14:       0     0     1     0     0         0       0     0     0     0     0
#> 15:       0     0     0     1     0         0       0     0     0     0     0
#> 16:       0     0     0     0     1         0       0     0     0     0     0
#> 17:       0     0     0     0     0         1       0     0     0     0     0
#> 18:       0     0     0     0     0         0       1     0     0     0     0
#> 19:       0     0     0     0     0         0       0     1     0     0     0
#> 20:       0     0     0     0     0         0       0     0     1     0     0
#> 21:       0     0     0     0     0         0       0     0     0     1     0
#> 22:       0     0     0     0     0         0       0     0     0     0     1
#> 23:       0     0     0     0     0         0       0     0     0     0     0
#> 24:       0     0     0     0     0         0       0     0     0     0     0
#> 25:       0     0     0     0     0         0       0     0     0     0     0
#> 26:       0     0     0     0     0         0       0     0     0     0     0
#> 27:       0     0     0     0     0         0       0     0     0     0     0
#> 28:       0     0     0     0     0         0       0     0     0     0     0
#> 29:       0     0     0     0     0         0       0     0     0     0     0
#> 30:       0     0     0     0     0         0       0     0     0     0     0
#> 31:       0     0     0     0     0         0       0     0     0     0     0
#> 32:       0     0     0     0     0         0       0     0     0     0     0
#> 33:       0     0     0     0     0         0       0     0     0     0     0
#>     CD45BC4 CD127  CD56 CD199 CD49a GranzymeB CD45BC5  CD95   CD3   CD4  CD69
#>       <num> <num> <num> <num> <num>     <num>   <num> <num> <num> <num> <num>
#>      TNFa CD183  IFNg CD103 CD45BC6 CD122 ia4b7 CD45BC7 Viability  CD27 KLRG1
#>     <num> <num> <num> <num>   <num> <num> <num>   <num>     <num> <num> <num>
#>  1:     0     0     0     0       0     0     0       0         0     0     0
#>  2:     0     0     0     0       0     0     0       0         0     0     0
#>  3:     0     0     0     0       0     0     0       0         0     0     0
#>  4:     0     0     0     0       0     0     0       0         0     0     0
#>  5:     0     0     0     0       0     0     0       0         0     0     0
#>  6:     0     0     0     0       0     0     0       0         0     0     0
#>  7:     0     0     0     0       0     0     0       0         0     0     0
#>  8:     0     0     0     0       0     0     0       0         0     0     0
#>  9:     0     0     0     0       0     0     0       0         0     0     0
#> 10:     0     0     0     0       0     0     0       0         0     0     0
#> 11:     0     0     0     0       0     0     0       0         0     0     0
#> 12:     0     0     0     0       0     0     0       0         0     0     0
#> 13:     0     0     0     0       0     0     0       0         0     0     0
#> 14:     0     0     0     0       0     0     0       0         0     0     0
#> 15:     0     0     0     0       0     0     0       0         0     0     0
#> 16:     0     0     0     0       0     0     0       0         0     0     0
#> 17:     0     0     0     0       0     0     0       0         0     0     0
#> 18:     0     0     0     0       0     0     0       0         0     0     0
#> 19:     0     0     0     0       0     0     0       0         0     0     0
#> 20:     0     0     0     0       0     0     0       0         0     0     0
#> 21:     0     0     0     0       0     0     0       0         0     0     0
#> 22:     0     0     0     0       0     0     0       0         0     0     0
#> 23:     1     0     0     0       0     0     0       0         0     0     0
#> 24:     0     1     0     0       0     0     0       0         0     0     0
#> 25:     0     0     1     0       0     0     0       0         0     0     0
#> 26:     0     0     0     1       0     0     0       0         0     0     0
#> 27:     0     0     0     0       1     0     0       0         0     0     0
#> 28:     0     0     0     0       0     1     0       0         0     0     0
#> 29:     0     0     0     0       0     0     1       0         0     0     0
#> 30:     0     0     0     0       0     0     0       1         0     0     0
#> 31:     0     0     0     0       0     0     0       0         1     0     0
#> 32:     0     0     0     0       0     0     0       0         0     1     0
#> 33:     0     0     0     0       0     0     0       0         0     0     1
#>      TNFa CD183  IFNg CD103 CD45BC6 CD122 ia4b7 CD45BC7 Viability  CD27 KLRG1
#>     <num> <num> <num> <num>   <num> <num> <num>   <num>     <num> <num> <num>

#read all .fcs files as a named list containing individual flowstates
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type = "S",
  concatenate = FALSE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
class(fs);class(fs[[1]])
#> [1] "list"
#> [1] "flowstate"
names(fs);names(fs[[1]])
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1.fcs" "COVAIL_002_CYTOKINE_BLOCK1_2.fcs"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3.fcs"
#> [1] "data"       "parameters" "keywords"   "spill"     

fs[[1]]$keywords
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 15-JUN-2026 22:35:54.16
#>           $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>                   <char>      <char>       <char> <char>
#> 1: flowstate_0.16.1.9001 aurora user DataModified     43
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
fs[[1]]$data[, levels(sample.id)]
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1"

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
class(fs)
#> [1] "flowstate"
names(fs)
#> [1] "data"       "parameters" "keywords"   "spill"     
fs$keywords
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#> 2: 09:38:12.18 Aurora  V0299 27-Feb-2025 09:49:40.92
#> 3: 09:50:46.94 Aurora  V0299 27-Feb-2025 10:03:26.16
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 15-JUN-2026 22:35:54.42
#> 2: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 15-JUN-2026 22:35:54.51
#> 3: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 15-JUN-2026 22:35:54.59
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
fs$data[, levels(sample.id)]
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1" "COVAIL_002_CYTOKINE_BLOCK1_2"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3"
```
