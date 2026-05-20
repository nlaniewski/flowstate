# Create a column (logical) for selecting against saturating events

A new column (logical) named `select.nonsaturating` is created in
`[['data']]` and any/all events that are at detector limits – as defined
by the '\$PnR' value in `[['parameters']]` – will be marked as `FALSE`
with non-saturating events marked as `TRUE`.

Detection of saturating events requires linear (non-transformed) values
and as such should be performed before
[flowstate.transform](https://nlaniewski.github.io/flowstate/reference/flowstate.transform.md).

## Usage

``` r
select_nonsaturating(flowstate)
```

## Arguments

- flowstate:

  A flowstate as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

## Value

UPDATES BY REFERENCE:

- `flowstate[['data']]`; adds a column (logical) named
  `select.nonsaturating`

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

#UPDATES BY REFERENCE -- adds a new column named 'select.nonsaturating'
select_nonsaturating(fs)
fs$data[, .N, by = select.nonsaturating]
#>    select.nonsaturating     N
#>                  <lgcl> <int>
#> 1:                 TRUE  5978
#> 2:                FALSE    22

#visualize
plot(fs,FSC_A,SSC_A) + ggplot2::facet_wrap(~select.nonsaturating)
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.
#> Warning: Computation failed in `stat_binhex()`.
#> Caused by error in `compute_group()`:
#> ! The package "hexbin" is required for `stat_bin_hex()`.


#subset to retain only non-saturating events
fs <- subset(fs, select.nonsaturating)
fs$data[, .N, by = select.nonsaturating]
#>    select.nonsaturating     N
#>                  <lgcl> <int>
#> 1:                 TRUE  5978

#after the subset, the column is now redundant
all(fs$data[['select.nonsaturating']])
#> [1] TRUE

#NULL it out
fs$data[,'select.nonsaturating' := NULL]
#>        Time    SSC_W  SSC_H    SSC_A    FSC_W  FSC_H     FSC_A   SSCB_W SSCB_H
#>       <num>    <num>  <num>    <num>    <num>  <num>     <num>    <num>  <num>
#>    1:     0 730290.7 759933 924953.4 774044.6 997371 1286682.8 714389.4 873960
#>    2:     1 753766.8 453137 569266.0 714706.8 444274  529209.4 749070.4 411951
#>    3:     4 768520.5 599195 767489.4 729628.4 827274 1006004.3 730502.3 579392
#>    4:     4 707174.8 516909 609241.7 719956.1 819074  982828.9 700778.6 463030
#>    5:     5 743404.2 484976 600888.7 764667.1 754833  961993.2 726898.2 419073
#>   ---                                                                         
#> 5974:  2454 755967.8 529870 667607.8 752809.6 781192  980148.1 748563.8 386375
#> 5975:  2454 681458.6 362164 411332.9 704912.6 856891 1006722.2 675241.8 333285
#> 5976:  2455 677347.4 701711 792170.2 702094.4 810637  948572.9 668789.1 506809
#> 5977:  2455 696805.0 546994 635246.9 734241.6 960241 1175081.5 696149.9 465449
#> 5978:  2455 717597.9 481212 575527.9 713571.4 772944  919251.3 712120.8 452721
#>          SSCB_A    CD45RA    CD45RO     TCRgd    CD45BC1         IL2
#>           <num>     <num>     <num>     <num>      <num>       <num>
#>    1: 1040579.6 17400.736  6296.233  3487.889  1056.3926 -2868.07104
#>    2:  514300.5  6299.763  1073.582  2427.801 21505.0137   666.11395
#>    3:  705412.0 -1247.931 37782.598  2189.499   146.5772    14.37641
#>    4:  540802.5  3424.521  6541.122 10304.088 31037.8750 -1450.19849
#>    5:  507705.7 -1692.568 43001.797  5568.431 55248.4922  -638.63544
#>   ---                                                               
#> 5974:  482043.9 10375.593  2197.690  2926.014 16656.1309 -1116.13855
#> 5975:  375079.9 17660.221  2524.003  3396.743 40435.3008  -350.11267
#> 5976:  564913.9  3601.014  3085.314  2354.051 -1026.0671  -749.64618
#> 5977:  540037.1  1594.400 31943.914 -6318.851 33177.8359 -1760.08154
#> 5978:  537320.1  4121.758  6534.809    80.228  7164.7236 -1993.13062
#>               CD8      CD197     CD45BC2   CD45BC3      CD57      CD193
#>             <num>      <num>       <num>     <num>     <num>      <num>
#>    1: 219290.2812  1701.1161 25329.08789 20529.553 2525.7261 1695.80396
#>    2:    260.5484 18941.9355 14742.99414  4324.678 4230.0068  162.88908
#>    3:   -377.2223   699.3767  1116.34326 21791.977 2767.2786  930.13574
#>    4:   1130.0894  1627.6073  -303.30679  3377.251 4663.3604   84.06442
#>    5:    236.9170  2753.4761  -910.34351  5274.306 7930.7036  172.67607
#>   ---                                                                  
#> 5974:  64804.6055  6916.9336 32589.94531 12607.492 3992.2173  442.76263
#> 5975:   5102.4351  6736.9619  -103.96835  3328.096 3078.8726  343.55698
#> 5976:   -533.4910  1094.5952   -35.77078 11619.281 6814.1572 -339.72955
#> 5977: 109734.0469   830.3315  4554.74463 12002.997 1818.2131  570.69165
#> 5978: 122607.1562  1581.7339  1644.67358 17445.123 -637.4016  813.78949
#>          CD45BC4     CD127       CD56      CD199     CD49a  GranzymeB
#>            <num>     <num>      <num>      <num>     <num>      <num>
#>    1: 23593.2891 1343.0120  6162.3535 1315.51941 1071.0972  54200.477
#>    2:   410.2411 -642.7087   990.9473 1259.63062 1404.4385   1982.248
#>    3: 30529.1172 2044.4512 -1573.5216 1891.64819  512.3239   2113.995
#>    4: 31405.3887 3629.7041 -1813.6000 1058.17969  217.9938   5147.355
#>    5:  1022.8995 5814.7285  1465.7168 3178.10767 -758.9792   1842.105
#>   ---                                                                
#> 5974:   747.4293 1195.2950  1832.6479 1560.49304  289.7756   1733.709
#> 5975: 26242.4258 1008.8782  1341.5767  769.40173  916.8894   2277.240
#> 5976: 11073.8740 -528.7997  1454.7413  926.90973 1844.4650 163543.656
#> 5977:  1443.5271 1471.8104  2837.1677   69.61098 1223.5042   2391.115
#> 5978: 23122.0137 3150.3503 -2212.3643 1250.06580 1173.8098   1466.612
#>            CD45BC5       CD95         CD3        CD4      CD69      TNFa
#>              <num>      <num>       <num>      <num>     <num>     <num>
#>    1:    -6.487877   942.8493 36690.11719  8299.1494 4737.0107  993.0143
#>    2:  -114.808891  2869.5635    94.28096   634.6076  912.6633  683.6307
#>    3: 15245.500000   922.2636 26607.62500 52088.1602 1719.4734  816.6953
#>    4: 19775.552734  -233.8560 51053.93750  5741.3486 6392.7744  544.5918
#>    5: 16411.826172  4008.5771 39481.12891 41355.3477 1742.8635 2127.6262
#>   ---                                                                   
#> 5974:    90.256508   521.3145 48319.74219  5864.3413 1021.8483  362.6751
#> 5975:   853.720154  1287.4240  1114.44263  2816.0164 4534.4888 -297.5424
#> 5976:  -282.902435 -2481.5879  1108.75098  1231.4691  141.7332 1626.9031
#> 5977:   255.437500   178.0732 23065.99414  4533.7129 3185.1895 5428.0649
#> 5978:   436.824249   752.9613 34094.38672  7511.4536 3074.0154 2097.9915
#>            CD183        IFNg       CD103    CD45BC6       CD122      ia4b7
#>            <num>       <num>       <num>      <num>       <num>      <num>
#>    1:   473.1266  -256.25250   166.61444   821.8005   735.86768  4871.1104
#>    2:  -478.1163  -109.60402   -21.77036 12146.9199   713.94666  1420.3256
#>    3:  -402.7975   395.29520   -83.91235  1955.8750 -1774.02917  -321.3192
#>    4:  -881.3683  -160.21269  -994.53149  1735.5421  -477.22177   487.7629
#>    5: -3175.0767    27.80499   155.62851  -825.0389  -652.52808 -1965.0686
#>   ---                                                                     
#> 5974:  -519.9846   300.67377  -334.04483   332.5317   121.39590  1214.4460
#> 5975:  -286.4161   206.22299  -241.14139 19857.8887   595.34851  2625.4993
#> 5976: -1262.1187 -1322.14587 -3573.42236 10710.0586   420.53311  2082.1335
#> 5977:  -787.6249 -1464.29724   -37.02888   881.4003   112.57523   200.0300
#> 5978:  -790.8807 -1168.62866   -23.95932 14355.7578    93.75003   896.9961
#>          CD45BC7  Viability       CD27       KLRG1                    sample.id
#>            <num>      <num>      <num>       <num>                       <fctr>
#>    1:  1030.7883  1429.1698 -7657.7598  5408.66895 COVAIL_002_CYTOKINE_BLOCK1_1
#>    2:  -409.5166 21949.6270   744.8879  3665.47559 COVAIL_002_CYTOKINE_BLOCK1_1
#>    3:   478.2293   543.1359  -448.9456  8486.89355 COVAIL_002_CYTOKINE_BLOCK1_1
#>    4:  -924.2852  1028.7167   332.4514 14109.10254 COVAIL_002_CYTOKINE_BLOCK1_1
#>    5: 15285.1455   795.2957  -447.1280   729.11121 COVAIL_002_CYTOKINE_BLOCK1_1
#>   ---                                                                          
#> 5974:   448.0044   270.1838  2724.5115   -17.58312 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5975:   470.1801   755.3552  -773.1400   808.78845 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5976:  -462.3215  4186.0117  -424.1915  1049.87292 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5977: 12610.7148  1859.3572 -1463.3777  3298.12988 COVAIL_002_CYTOKINE_BLOCK1_3
#> 5978:  -454.0528  1697.2017 -2134.2859  5035.73535 COVAIL_002_CYTOKINE_BLOCK1_3
```
