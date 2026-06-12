# Concatenate a list of flowstates

Concatenate a list of flowstates

## Usage

``` r
concatenate.flowstate(flowstates)
```

## Arguments

- flowstates:

  a list; the return of
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md)

## Value

An object of class flowstate.

## Examples

``` r
fcs.file.paths <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")

#read all .fcs files as flowstates
fs <- read.flowstate(
  fcs.file.paths,
  colnames.type = "S"
)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_2.fcs --> flowstate
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs --> flowstate

#a list of flowstates
class(fs)
#> [1] "list"
sapply(fs, class)
#> COVAIL_002_CYTOKINE_BLOCK1_1.fcs COVAIL_002_CYTOKINE_BLOCK1_2.fcs 
#>                      "flowstate"                      "flowstate" 
#> COVAIL_002_CYTOKINE_BLOCK1_3.fcs 
#>                      "flowstate" 

#concatenate into a single flowstate
fs <- concatenate.flowstate(fs)
#> Concatenating 'flowstates'...
class(fs)
#> [1] "flowstate"

fs$data
#>        Time    SSC_W  SSC_H    SSC_A    FSC_W  FSC_H     FSC_A   SSCB_W SSCB_H
#>       <num>    <num>  <num>    <num>    <num>  <num>     <num>    <num>  <num>
#>    1:     0 730290.7 759933 924953.4 774044.6 997371 1286682.8 714389.4 873960
#>    2:     1 753766.8 453137 569266.0 714706.8 444274  529209.4 749070.4 411951
#>    3:     4 768520.5 599195 767489.4 729628.4 827274 1006004.3 730502.3 579392
#>    4:     4 707174.8 516909 609241.7 719956.1 819074  982828.9 700778.6 463030
#>    5:     5 743404.2 484976 600888.7 764667.1 754833  961993.2 726898.2 419073
#>   ---                                                                         
#> 5996:  2454 755967.8 529870 667607.8 752809.6 781192  980148.1 748563.8 386375
#> 5997:  2454 681458.6 362164 411332.9 704912.6 856891 1006722.2 675241.8 333285
#> 5998:  2455 677347.4 701711 792170.2 702094.4 810637  948572.9 668789.1 506809
#> 5999:  2455 696805.0 546994 635246.9 734241.6 960241 1175081.5 696149.9 465449
#> 6000:  2455 717597.9 481212 575527.9 713571.4 772944  919251.3 712120.8 452721
#>          SSCB_A    CD45RA    CD45RO     TCRgd    CD45BC1         IL2
#>           <num>     <num>     <num>     <num>      <num>       <num>
#>    1: 1040579.6 17400.736  6296.233  3487.889  1056.3926 -2868.07104
#>    2:  514300.5  6299.763  1073.582  2427.801 21505.0137   666.11395
#>    3:  705412.0 -1247.931 37782.598  2189.499   146.5772    14.37641
#>    4:  540802.5  3424.521  6541.122 10304.088 31037.8750 -1450.19849
#>    5:  507705.7 -1692.568 43001.797  5568.431 55248.4922  -638.63544
#>   ---                                                               
#> 5996:  482043.9 10375.593  2197.690  2926.014 16656.1309 -1116.13855
#> 5997:  375079.9 17660.221  2524.003  3396.743 40435.3008  -350.11267
#> 5998:  564913.9  3601.014  3085.314  2354.051 -1026.0671  -749.64618
#> 5999:  540037.1  1594.400 31943.914 -6318.851 33177.8359 -1760.08154
#> 6000:  537320.1  4121.758  6534.809    80.228  7164.7236 -1993.13062
#>               CD8      CD197     CD45BC2   CD45BC3      CD57      CD193
#>             <num>      <num>       <num>     <num>     <num>      <num>
#>    1: 219290.2812  1701.1161 25329.08789 20529.553 2525.7261 1695.80396
#>    2:    260.5484 18941.9355 14742.99414  4324.678 4230.0068  162.88908
#>    3:   -377.2223   699.3767  1116.34326 21791.977 2767.2786  930.13574
#>    4:   1130.0894  1627.6073  -303.30679  3377.251 4663.3604   84.06442
#>    5:    236.9170  2753.4761  -910.34351  5274.306 7930.7036  172.67607
#>   ---                                                                  
#> 5996:  64804.6055  6916.9336 32589.94531 12607.492 3992.2173  442.76263
#> 5997:   5102.4351  6736.9619  -103.96835  3328.096 3078.8726  343.55698
#> 5998:   -533.4910  1094.5952   -35.77078 11619.281 6814.1572 -339.72955
#> 5999: 109734.0469   830.3315  4554.74463 12002.997 1818.2131  570.69165
#> 6000: 122607.1562  1581.7339  1644.67358 17445.123 -637.4016  813.78949
#>          CD45BC4     CD127       CD56      CD199     CD49a  GranzymeB
#>            <num>     <num>      <num>      <num>     <num>      <num>
#>    1: 23593.2891 1343.0120  6162.3535 1315.51941 1071.0972  54200.477
#>    2:   410.2411 -642.7087   990.9473 1259.63062 1404.4385   1982.248
#>    3: 30529.1172 2044.4512 -1573.5216 1891.64819  512.3239   2113.995
#>    4: 31405.3887 3629.7041 -1813.6000 1058.17969  217.9938   5147.355
#>    5:  1022.8995 5814.7285  1465.7168 3178.10767 -758.9792   1842.105
#>   ---                                                                
#> 5996:   747.4293 1195.2950  1832.6479 1560.49304  289.7756   1733.709
#> 5997: 26242.4258 1008.8782  1341.5767  769.40173  916.8894   2277.240
#> 5998: 11073.8740 -528.7997  1454.7413  926.90973 1844.4650 163543.656
#> 5999:  1443.5271 1471.8104  2837.1677   69.61098 1223.5042   2391.115
#> 6000: 23122.0137 3150.3503 -2212.3643 1250.06580 1173.8098   1466.612
#>            CD45BC5       CD95         CD3        CD4      CD69      TNFa
#>              <num>      <num>       <num>      <num>     <num>     <num>
#>    1:    -6.487877   942.8493 36690.11719  8299.1494 4737.0107  993.0143
#>    2:  -114.808891  2869.5635    94.28096   634.6076  912.6633  683.6307
#>    3: 15245.500000   922.2636 26607.62500 52088.1602 1719.4734  816.6953
#>    4: 19775.552734  -233.8560 51053.93750  5741.3486 6392.7744  544.5918
#>    5: 16411.826172  4008.5771 39481.12891 41355.3477 1742.8635 2127.6262
#>   ---                                                                   
#> 5996:    90.256508   521.3145 48319.74219  5864.3413 1021.8483  362.6751
#> 5997:   853.720154  1287.4240  1114.44263  2816.0164 4534.4888 -297.5424
#> 5998:  -282.902435 -2481.5879  1108.75098  1231.4691  141.7332 1626.9031
#> 5999:   255.437500   178.0732 23065.99414  4533.7129 3185.1895 5428.0649
#> 6000:   436.824249   752.9613 34094.38672  7511.4536 3074.0154 2097.9915
#>            CD183        IFNg       CD103    CD45BC6       CD122      ia4b7
#>            <num>       <num>       <num>      <num>       <num>      <num>
#>    1:   473.1266  -256.25250   166.61444   821.8005   735.86768  4871.1104
#>    2:  -478.1163  -109.60402   -21.77036 12146.9199   713.94666  1420.3256
#>    3:  -402.7975   395.29520   -83.91235  1955.8750 -1774.02917  -321.3192
#>    4:  -881.3683  -160.21269  -994.53149  1735.5421  -477.22177   487.7629
#>    5: -3175.0767    27.80499   155.62851  -825.0389  -652.52808 -1965.0686
#>   ---                                                                     
#> 5996:  -519.9846   300.67377  -334.04483   332.5317   121.39590  1214.4460
#> 5997:  -286.4161   206.22299  -241.14139 19857.8887   595.34851  2625.4993
#> 5998: -1262.1187 -1322.14587 -3573.42236 10710.0586   420.53311  2082.1335
#> 5999:  -787.6249 -1464.29724   -37.02888   881.4003   112.57523   200.0300
#> 6000:  -790.8807 -1168.62866   -23.95932 14355.7578    93.75003   896.9961
#>          CD45BC7  Viability       CD27       KLRG1                    sample.id
#>            <num>      <num>      <num>       <num>                       <fctr>
#>    1:  1030.7883  1429.1698 -7657.7598  5408.66895 COVAIL_002_CYTOKINE_BLOCK1_1
#>    2:  -409.5166 21949.6270   744.8879  3665.47559 COVAIL_002_CYTOKINE_BLOCK1_1
#>    3:   478.2293   543.1359  -448.9456  8486.89355 COVAIL_002_CYTOKINE_BLOCK1_1
#>    4:  -924.2852  1028.7167   332.4514 14109.10254 COVAIL_002_CYTOKINE_BLOCK1_1
#>    5: 15285.1455   795.2957  -447.1280   729.11121 COVAIL_002_CYTOKINE_BLOCK1_1
#>   ---                                                                          
#> 5996:   448.0044   270.1838  2724.5115   -17.58312 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5997:   470.1801   755.3552  -773.1400   808.78845 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5998:  -462.3215  4186.0117  -424.1915  1049.87292 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5999: 12610.7148  1859.3572 -1463.3777  3298.12988 COVAIL_002_CYTOKINE_BLOCK1_3
#> 6000:  -454.0528  1697.2017 -2134.2859  5035.73535 COVAIL_002_CYTOKINE_BLOCK1_3
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
fs$keywords
#>          $BTIM   $CYT $CYTSN       $DATE       $ETIM
#>         <char> <char> <char>      <char>      <char>
#> 1: 09:22:51.62 Aurora  V0299 27-Feb-2025 09:33:47.32
#> 2: 09:38:12.18 Aurora  V0299 27-Feb-2025 09:49:40.92
#> 3: 09:50:46.94 Aurora  V0299 27-Feb-2025 10:03:26.16
#>                                $FIL $FLOWRATE    $INST          $LAST_MODIFIED
#>                              <char>    <char>   <char>                  <char>
#> 1: COVAIL_002_CYTOKINE_BLOCK1_1.fcs    Medium Cytekbio 12-JUN-2026 18:58:03.63
#> 2: COVAIL_002_CYTOKINE_BLOCK1_2.fcs    Medium Cytekbio 12-JUN-2026 18:58:03.71
#> 3: COVAIL_002_CYTOKINE_BLOCK1_3.fcs    Medium Cytekbio 12-JUN-2026 18:58:03.79
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
#>    WINDOW EXTENSION                    sample.id
#>              <char>                       <fctr>
#> 1:                3 COVAIL_002_CYTOKINE_BLOCK1_1
#> 2:                3 COVAIL_002_CYTOKINE_BLOCK1_2
#> 3:                3 COVAIL_002_CYTOKINE_BLOCK1_3
fs$spill
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
```
