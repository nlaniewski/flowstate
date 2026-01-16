# flowstate -- Structure

## Structure

[Flow Cytometry
Standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC2892967/ "FCS 3.1")
files are read and parsed into a `flowstate` (S3 object). A `flowstate`
is primarily comprised of the following parsed (named list) elements –
each a
[`data.table::data.table()`](https://rdatatable.gitlab.io/data.table/reference/data.table.html):

| Name         | Value                                                                |
|:-------------|:---------------------------------------------------------------------|
| `data`       | expression values                                                    |
| `parameters` | instrument-specific measurement parameters                           |
| `keywords`   | instrument/sample-specific metadata                                  |
| `spill`      | instrument/sample-specific spillover matrix (if encoded in the file) |

Included with the package – for the purpose of example(s) – are a few
(subsampled) .fcs files acquired on a Cytek Aurora (5L) using SpectroFlo
software.

Read/parse them as a `flowstate` using
[`flowstate::read.flowstate()`](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

``` r
##paths to example .fcs files
fcs.files <- system.file("extdata", package = "flowstate") |> 
  list.files(full.names = T) |>
  grep(pattern = "BLOCK.*.fcs",value = T)
fcs.files |> basename() |> print()
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1.fcs" "COVAIL_002_CYTOKINE_BLOCK1_2.fcs"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3.fcs"

##read them in and concatenate; a flowstate object
fs <- flowstate::read.flowstate(
  fcs.file.paths = fcs.files,
  colnames.type = "S_N",
  concatenate = TRUE
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate
#> Concatenating 'flowstate.ojects'...

##S3 object of class "flowstate"
class(fs)
#> [1] "flowstate"
##flowstate structure
names(fs)
#> [1] "data"       "parameters" "keywords"   "spill"
```

### `[['data']]`

Data can be accessed by name via: `fs$data` or `fs[['data']]`.

`[['data']]` contains events as rows and scatter, time, and fluorescent
measurements as columns. By design, no transformation is applied to the
linear measurement values during reading. In addition, a unique sample
identifier (factored) is added for workflow purposes.

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
#>       SSCB_A CD45RA_BUV395 CD45RO_BUV496 TCRgd_BUV563 CD45BC1_SBUV605
#>        <num>         <num>         <num>        <num>           <num>
#> 1: 1040579.6     17400.736      6296.233     3487.889       1056.3926
#> 2:  514300.5      6299.763      1073.582     2427.801      21505.0137
#> 3:  705412.0     -1247.931     37782.598     2189.499        146.5772
#> 4:  540802.5      3424.521      6541.122    10304.088      31037.8750
#> 5:  507705.7     -1692.568     43001.797     5568.431      55248.4922
#> 6:  471449.8     26003.760      3743.393     2530.474       -219.0977
#>     IL2_BUV737  CD8_BUV805 CD197_BV421 CD45BC2_PacificBlue CD45BC3_SBV475
#>          <num>       <num>       <num>               <num>          <num>
#> 1: -2868.07104 219290.2812   1701.1161          25329.0879      20529.553
#> 2:   666.11395    260.5484  18941.9355          14742.9941       4324.678
#> 3:    14.37641   -377.2223    699.3767           1116.3433      21791.977
#> 4: -1450.19849   1130.0894   1627.6073           -303.3068       3377.251
#> 5:  -638.63544    236.9170   2753.4761           -910.3435       5274.306
#> 6:   135.40210   1106.0740   1290.4441          23654.0664       3429.047
#>    CD57_BV605 CD193_BV650 CD45BC4_BV711 CD127_BV750 CD56_BV785
#>         <num>       <num>         <num>       <num>      <num>
#> 1:   2525.726  1695.80396    23593.2891   1343.0120  6162.3535
#> 2:   4230.007   162.88908      410.2411   -642.7087   990.9473
#> 3:   2767.279   930.13574    30529.1172   2044.4512 -1573.5216
#> 4:   4663.360    84.06442    31405.3887   3629.7041 -1813.6000
#> 5:   7930.704   172.67607     1022.8995   5814.7285  1465.7168
#> 6:   3077.539    74.12897      335.7986   1104.4130  1660.5831
#>    CD199_KIRAVIABlue520 CD49a_RB545 GranzymeB_RB613 CD45BC5_PerCP
#>                   <num>       <num>           <num>         <num>
#> 1:             1315.519   1071.0972       54200.477     -6.487877
#> 2:             1259.631   1404.4385        1982.248   -114.808891
#> 3:             1891.648    512.3239        2113.995  15245.500000
#> 4:             1058.180    217.9938        5147.355  19775.552734
#> 5:             3178.108   -758.9792        1842.105  16411.826172
#> 6:             1381.134    684.5980       30221.033    390.726471
#>    CD95_PerCPeFluor710   CD3_RB744 CD4_PerCPFire806 CD69_RB780   TNFa_PE
#>                  <num>       <num>            <num>      <num>     <num>
#> 1:           942.84930 36690.11719        8299.1494  4737.0107  993.0143
#> 2:          2869.56348    94.28096         634.6076   912.6633  683.6307
#> 3:           922.26361 26607.62500       52088.1602  1719.4734  816.6953
#> 4:          -233.85599 51053.93750        5741.3486  6392.7744  544.5918
#> 5:          4008.57715 39481.12891       41355.3477  1742.8635 2127.6262
#> 6:            62.77901   254.06075        6056.9761 31786.8047  939.8767
#>    CD183_RY586 IFNg_PEDazzle594 CD103_PEFire640 CD45BC6_SBY720 CD122_PECy7
#>          <num>            <num>           <num>          <num>       <num>
#> 1:    473.1266       -256.25250       166.61444       821.8005    735.8677
#> 2:   -478.1163       -109.60402       -21.77036     12146.9199    713.9467
#> 3:   -402.7975        395.29520       -83.91235      1955.8750  -1774.0292
#> 4:   -881.3683       -160.21269      -994.53149      1735.5421   -477.2218
#> 5:  -3175.0767         27.80499       155.62851      -825.0389   -652.5281
#> 6:   -854.2112         63.90987       636.24811      3377.4231   -126.8232
#>    ia4b7_AlexaFluor647 CD45BC7_SparkNIR685 Viability_LIVEDEADScarlet CD27_APCH7
#>                  <num>               <num>                     <num>      <num>
#> 1:           4871.1104           1030.7883                 1429.1698 -7657.7598
#> 2:           1420.3256           -409.5166                21949.6270   744.8879
#> 3:           -321.3192            478.2293                  543.1359  -448.9456
#> 4:            487.7629           -924.2852                 1028.7167   332.4514
#> 5:          -1965.0686          15285.1455                  795.2957  -447.1280
#> 6:          -1043.0923          18769.6777                 1230.0640 -1377.2864
#>    KLRG1_APCFire810                    sample.id
#>               <num>                       <fctr>
#> 1:        5408.6689 COVAIL_002_CYTOKINE_BLOCK1_1
#> 2:        3665.4756 COVAIL_002_CYTOKINE_BLOCK1_1
#> 3:        8486.8936 COVAIL_002_CYTOKINE_BLOCK1_1
#> 4:       14109.1025 COVAIL_002_CYTOKINE_BLOCK1_1
#> 5:         729.1112 COVAIL_002_CYTOKINE_BLOCK1_1
#> 6:        7417.7949 COVAIL_002_CYTOKINE_BLOCK1_1

##variable names
names(fs$data)
#>  [1] "Time"                      "SSC_W"                    
#>  [3] "SSC_H"                     "SSC_A"                    
#>  [5] "FSC_W"                     "FSC_H"                    
#>  [7] "FSC_A"                     "SSCB_W"                   
#>  [9] "SSCB_H"                    "SSCB_A"                   
#> [11] "CD45RA_BUV395"             "CD45RO_BUV496"            
#> [13] "TCRgd_BUV563"              "CD45BC1_SBUV605"          
#> [15] "IL2_BUV737"                "CD8_BUV805"               
#> [17] "CD197_BV421"               "CD45BC2_PacificBlue"      
#> [19] "CD45BC3_SBV475"            "CD57_BV605"               
#> [21] "CD193_BV650"               "CD45BC4_BV711"            
#> [23] "CD127_BV750"               "CD56_BV785"               
#> [25] "CD199_KIRAVIABlue520"      "CD49a_RB545"              
#> [27] "GranzymeB_RB613"           "CD45BC5_PerCP"            
#> [29] "CD95_PerCPeFluor710"       "CD3_RB744"                
#> [31] "CD4_PerCPFire806"          "CD69_RB780"               
#> [33] "TNFa_PE"                   "CD183_RY586"              
#> [35] "IFNg_PEDazzle594"          "CD103_PEFire640"          
#> [37] "CD45BC6_SBY720"            "CD122_PECy7"              
#> [39] "ia4b7_AlexaFluor647"       "CD45BC7_SparkNIR685"      
#> [41] "Viability_LIVEDEADScarlet" "CD27_APCH7"               
#> [43] "KLRG1_APCFire810"          "sample.id"

##sample identifiers
fs$data[,levels(sample.id)]
#> [1] "COVAIL_002_CYTOKINE_BLOCK1_1" "COVAIL_002_CYTOKINE_BLOCK1_2"
#> [3] "COVAIL_002_CYTOKINE_BLOCK1_3"
```

### `[['parameters']]`

Parameters can be accessed by name via: `fs$parameters` or
`fs[['parameters']]`.

`[['parameters']]` contains instrument-specific measurement…parameters.

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
#>                     TYPE      V                           PROJ         N.alias
#>                   <char> <char>                         <fctr>          <char>
#>  1:                 Time   <NA> COVAIL_002_CYTOKINE_2025-02-27            Time
#>  2:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27           SSC_W
#>  3:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27           SSC_H
#>  4:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27           SSC_A
#>  5:      Forward_Scatter     40 COVAIL_002_CYTOKINE_2025-02-27           FSC_W
#>  6:      Forward_Scatter     40 COVAIL_002_CYTOKINE_2025-02-27           FSC_H
#>  7:      Forward_Scatter     40 COVAIL_002_CYTOKINE_2025-02-27           FSC_A
#>  8:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27          SSCB_W
#>  9:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27          SSCB_H
#> 10:         Side_Scatter    220 COVAIL_002_CYTOKINE_2025-02-27          SSCB_A
#> 11: Unmixed_Fluorescence    258 COVAIL_002_CYTOKINE_2025-02-27          BUV395
#> 12: Unmixed_Fluorescence    523 COVAIL_002_CYTOKINE_2025-02-27          BUV496
#> 13: Unmixed_Fluorescence    535 COVAIL_002_CYTOKINE_2025-02-27          BUV563
#> 14: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27         SBUV605
#> 15: Unmixed_Fluorescence    904 COVAIL_002_CYTOKINE_2025-02-27          BUV737
#> 16: Unmixed_Fluorescence   1087 COVAIL_002_CYTOKINE_2025-02-27          BUV805
#> 17: Unmixed_Fluorescence    301 COVAIL_002_CYTOKINE_2025-02-27           BV421
#> 18: Unmixed_Fluorescence    347 COVAIL_002_CYTOKINE_2025-02-27     PacificBlue
#> 19: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27          SBV475
#> 20: Unmixed_Fluorescence    330 COVAIL_002_CYTOKINE_2025-02-27           BV605
#> 21: Unmixed_Fluorescence    306 COVAIL_002_CYTOKINE_2025-02-27           BV650
#> 22: Unmixed_Fluorescence    254 COVAIL_002_CYTOKINE_2025-02-27           BV711
#> 23: Unmixed_Fluorescence    327 COVAIL_002_CYTOKINE_2025-02-27           BV750
#> 24: Unmixed_Fluorescence    539 COVAIL_002_CYTOKINE_2025-02-27           BV785
#> 25: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27  KIRAVIABlue520
#> 26: Unmixed_Fluorescence    639 COVAIL_002_CYTOKINE_2025-02-27           RB545
#> 27: Unmixed_Fluorescence    369 COVAIL_002_CYTOKINE_2025-02-27           RB613
#> 28: Unmixed_Fluorescence    401 COVAIL_002_CYTOKINE_2025-02-27           PerCP
#> 29: Unmixed_Fluorescence    479 COVAIL_002_CYTOKINE_2025-02-27  PerCPeFluor710
#> 30: Unmixed_Fluorescence    253 COVAIL_002_CYTOKINE_2025-02-27           RB744
#> 31: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27    PerCPFire806
#> 32: Unmixed_Fluorescence    616 COVAIL_002_CYTOKINE_2025-02-27           RB780
#> 33: Unmixed_Fluorescence    312 COVAIL_002_CYTOKINE_2025-02-27              PE
#> 34: Unmixed_Fluorescence    312 COVAIL_002_CYTOKINE_2025-02-27           RY586
#> 35: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27     PEDazzle594
#> 36: Unmixed_Fluorescence    358 COVAIL_002_CYTOKINE_2025-02-27       PEFire640
#> 37: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27          SBY720
#> 38: Unmixed_Fluorescence    296 COVAIL_002_CYTOKINE_2025-02-27           PECy7
#> 39: Unmixed_Fluorescence    159 COVAIL_002_CYTOKINE_2025-02-27   AlexaFluor647
#> 40: Unmixed_Fluorescence    159 COVAIL_002_CYTOKINE_2025-02-27     SparkNIR685
#> 41: Unmixed_Fluorescence      0 COVAIL_002_CYTOKINE_2025-02-27 LIVEDEADScarlet
#> 42: Unmixed_Fluorescence    388 COVAIL_002_CYTOKINE_2025-02-27           APCH7
#> 43: Unmixed_Fluorescence    264 COVAIL_002_CYTOKINE_2025-02-27      APCFire810
#>                     TYPE      V                           PROJ         N.alias
#>                   <char> <char>                         <fctr>          <char>
#>       S.alias                 S_N.alias
#>        <char>                    <char>
#>  1:      Time                      Time
#>  2:     SSC_W                     SSC_W
#>  3:     SSC_H                     SSC_H
#>  4:     SSC_A                     SSC_A
#>  5:     FSC_W                     FSC_W
#>  6:     FSC_H                     FSC_H
#>  7:     FSC_A                     FSC_A
#>  8:    SSCB_W                    SSCB_W
#>  9:    SSCB_H                    SSCB_H
#> 10:    SSCB_A                    SSCB_A
#> 11:    CD45RA             CD45RA_BUV395
#> 12:    CD45RO             CD45RO_BUV496
#> 13:     TCRgd              TCRgd_BUV563
#> 14:   CD45BC1           CD45BC1_SBUV605
#> 15:       IL2                IL2_BUV737
#> 16:       CD8                CD8_BUV805
#> 17:     CD197               CD197_BV421
#> 18:   CD45BC2       CD45BC2_PacificBlue
#> 19:   CD45BC3            CD45BC3_SBV475
#> 20:      CD57                CD57_BV605
#> 21:     CD193               CD193_BV650
#> 22:   CD45BC4             CD45BC4_BV711
#> 23:     CD127               CD127_BV750
#> 24:      CD56                CD56_BV785
#> 25:     CD199      CD199_KIRAVIABlue520
#> 26:     CD49a               CD49a_RB545
#> 27: GranzymeB           GranzymeB_RB613
#> 28:   CD45BC5             CD45BC5_PerCP
#> 29:      CD95       CD95_PerCPeFluor710
#> 30:       CD3                 CD3_RB744
#> 31:       CD4          CD4_PerCPFire806
#> 32:      CD69                CD69_RB780
#> 33:      TNFa                   TNFa_PE
#> 34:     CD183               CD183_RY586
#> 35:      IFNg          IFNg_PEDazzle594
#> 36:     CD103           CD103_PEFire640
#> 37:   CD45BC6            CD45BC6_SBY720
#> 38:     CD122               CD122_PECy7
#> 39:     ia4b7       ia4b7_AlexaFluor647
#> 40:   CD45BC7       CD45BC7_SparkNIR685
#> 41: Viability Viability_LIVEDEADScarlet
#> 42:      CD27                CD27_APCH7
#> 43:     KLRG1          KLRG1_APCFire810
#>       S.alias                 S_N.alias
#>        <char>                    <char>

##internal flowstate functions add a few 'alias' columns;
##they serve the purpose of renaming the columns in [['data']];
##the colnames.type argument in read.flowstate uses these alias columns;
##some combination of N and S -- made syntactically valid

##original (encoded) parameters: N and S
fs$parameters[,.(N,S)]
#>                       N         S
#>                  <char>    <char>
#>  1:                Time      <NA>
#>  2:               SSC-W      <NA>
#>  3:               SSC-H      <NA>
#>  4:               SSC-A      <NA>
#>  5:               FSC-W      <NA>
#>  6:               FSC-H      <NA>
#>  7:               FSC-A      <NA>
#>  8:             SSC-B-W      <NA>
#>  9:             SSC-B-H      <NA>
#> 10:             SSC-B-A      <NA>
#> 11:            BUV395-A    CD45RA
#> 12:            BUV496-A    CD45RO
#> 13:            BUV563-A     TCRgd
#> 14:           SBUV605-A   CD45BC1
#> 15:            BUV737-A       IL2
#> 16:            BUV805-A       CD8
#> 17:             BV421-A     CD197
#> 18:      Pacific Blue-A   CD45BC2
#> 19:            SBV475-A   CD45BC3
#> 20:             BV605-A      CD57
#> 21:             BV650-A     CD193
#> 22:             BV711-A   CD45BC4
#> 23:             BV750-A     CD127
#> 24:             BV785-A      CD56
#> 25:  KIRAVIA Blue 520-A     CD199
#> 26:             RB545-A     CD49a
#> 27:             RB613-A GranzymeB
#> 28:             PerCP-A   CD45BC5
#> 29:  PerCP-eFluor 710-A      CD95
#> 30:             RB744-A       CD3
#> 31:    PerCP-Fire 806-A       CD4
#> 32:             RB780-A      CD69
#> 33:                PE-A      TNFa
#> 34:             RY586-A     CD183
#> 35:     PE-Dazzle 594-A      IFNg
#> 36:       PE-Fire 640-A     CD103
#> 37:            SBY720-A   CD45BC6
#> 38:            PE-Cy7-A     CD122
#> 39:   Alexa Fluor 647-A     ia4b7
#> 40:     Spark NIR 685-A   CD45BC7
#> 41: LIVE DEAD Scarlet-A Viability
#> 42:            APC-H7-A      CD27
#> 43:      APC-Fire 810-A     KLRG1
#>                       N         S
#>                  <char>    <char>

##alias columns used to rename [['data']] columns
fs$parameters[,.(N.alias,S.alias)]
#>             N.alias   S.alias
#>              <char>    <char>
#>  1:            Time      Time
#>  2:           SSC_W     SSC_W
#>  3:           SSC_H     SSC_H
#>  4:           SSC_A     SSC_A
#>  5:           FSC_W     FSC_W
#>  6:           FSC_H     FSC_H
#>  7:           FSC_A     FSC_A
#>  8:          SSCB_W    SSCB_W
#>  9:          SSCB_H    SSCB_H
#> 10:          SSCB_A    SSCB_A
#> 11:          BUV395    CD45RA
#> 12:          BUV496    CD45RO
#> 13:          BUV563     TCRgd
#> 14:         SBUV605   CD45BC1
#> 15:          BUV737       IL2
#> 16:          BUV805       CD8
#> 17:           BV421     CD197
#> 18:     PacificBlue   CD45BC2
#> 19:          SBV475   CD45BC3
#> 20:           BV605      CD57
#> 21:           BV650     CD193
#> 22:           BV711   CD45BC4
#> 23:           BV750     CD127
#> 24:           BV785      CD56
#> 25:  KIRAVIABlue520     CD199
#> 26:           RB545     CD49a
#> 27:           RB613 GranzymeB
#> 28:           PerCP   CD45BC5
#> 29:  PerCPeFluor710      CD95
#> 30:           RB744       CD3
#> 31:    PerCPFire806       CD4
#> 32:           RB780      CD69
#> 33:              PE      TNFa
#> 34:           RY586     CD183
#> 35:     PEDazzle594      IFNg
#> 36:       PEFire640     CD103
#> 37:          SBY720   CD45BC6
#> 38:           PECy7     CD122
#> 39:   AlexaFluor647     ia4b7
#> 40:     SparkNIR685   CD45BC7
#> 41: LIVEDEADScarlet Viability
#> 42:           APCH7      CD27
#> 43:      APCFire810     KLRG1
#>             N.alias   S.alias
#>              <char>    <char>
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
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 16-JAN-2026 18:46:36.80
#> 2: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 16-JAN-2026 18:46:36.82
#> 3: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 16-JAN-2026 18:46:36.84
#>    $LAST_MODIFIER         $OP $ORIGINALITY   $PAR
#>            <char>      <char>       <char> <char>
#> 1:                aurora user DataModified     43
#> 2:                aurora user DataModified     43
#> 3:                aurora user DataModified     43
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
#>    WINDOW EXTENSION
#>              <char>
#> 1:                3
#> 2:                3
#> 3:                3

##some identifiers
fs$keywords[,.(`$CYT`,TUBENAME)]
#>      $CYT                     TUBENAME
#>    <char>                       <char>
#> 1: Aurora COVAIL_002_CYTOKINE_BLOCK1_1
#> 2: Aurora COVAIL_002_CYTOKINE_BLOCK1_2
#> 3: Aurora COVAIL_002_CYTOKINE_BLOCK1_3

##keywords to indicate/track modification
fs$keywords[,.(`$LAST_MODIFIED`,`$LAST_MODIFIER`,`$ORIGINALITY`)]
#>             $LAST_MODIFIED $LAST_MODIFIER $ORIGINALITY
#>                     <char>         <char>       <char>
#> 1: 16-JAN-2026 18:46:36.80                DataModified
#> 2: 16-JAN-2026 18:46:36.82                DataModified
#> 3: 16-JAN-2026 18:46:36.84                DataModified
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
#>    CD45RA_BUV395 CD45RO_BUV496 TCRgd_BUV563 CD45BC1_SBUV605 IL2_BUV737
#>            <num>         <num>        <num>           <num>      <num>
#> 1:         1e+00             0            0               0          0
#> 2:         1e-06             1            0               0          0
#> 3:         0e+00             0            1               0          0
#> 4:         0e+00             0            0               1          0
#> 5:         0e+00             0            0               0          1
#> 6:         0e+00             0            0               0          0
#>    CD8_BUV805 CD197_BV421 CD45BC2_PacificBlue CD45BC3_SBV475 CD57_BV605
#>         <num>       <num>               <num>          <num>      <num>
#> 1:          0           0                   0              0          0
#> 2:          0           0                   0              0          0
#> 3:          0           0                   0              0          0
#> 4:          0           0                   0              0          0
#> 5:          0           0                   0              0          0
#> 6:          1           0                   0              0          0
#>    CD193_BV650 CD45BC4_BV711 CD127_BV750 CD56_BV785 CD199_KIRAVIABlue520
#>          <num>         <num>       <num>      <num>                <num>
#> 1:           0             0           0          0                    0
#> 2:           0             0           0          0                    0
#> 3:           0             0           0          0                    0
#> 4:           0             0           0          0                    0
#> 5:           0             0           0          0                    0
#> 6:           0             0           0          0                    0
#>    CD49a_RB545 GranzymeB_RB613 CD45BC5_PerCP CD95_PerCPeFluor710 CD3_RB744
#>          <num>           <num>         <num>               <num>     <num>
#> 1:           0               0             0                   0         0
#> 2:           0               0             0                   0         0
#> 3:           0               0             0                   0         0
#> 4:           0               0             0                   0         0
#> 5:           0               0             0                   0         0
#> 6:           0               0             0                   0         0
#>    CD4_PerCPFire806 CD69_RB780 TNFa_PE CD183_RY586 IFNg_PEDazzle594
#>               <num>      <num>   <num>       <num>            <num>
#> 1:                0          0       0           0                0
#> 2:                0          0       0           0                0
#> 3:                0          0       0           0                0
#> 4:                0          0       0           0                0
#> 5:                0          0       0           0                0
#> 6:                0          0       0           0                0
#>    CD103_PEFire640 CD45BC6_SBY720 CD122_PECy7 ia4b7_AlexaFluor647
#>              <num>          <num>       <num>               <num>
#> 1:               0              0           0                   0
#> 2:               0              0           0                   0
#> 3:               0              0           0                   0
#> 4:               0              0           0                   0
#> 5:               0              0           0                   0
#> 6:               0              0           0                   0
#>    CD45BC7_SparkNIR685 Viability_LIVEDEADScarlet CD27_APCH7 KLRG1_APCFire810
#>                  <num>                     <num>      <num>            <num>
#> 1:                   0                         0          0                0
#> 2:                   0                         0          0                0
#> 3:                   0                         0          0                0
#> 4:                   0                         0          0                0
#> 5:                   0                         0          0                0
#> 6:                   0                         0          0                0

##small value; there in the source data
fs$spill[2,1]
#>    CD45RA_BUV395
#>            <num>
#> 1:         1e-06
```
