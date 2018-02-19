---
layout: default
---

Which CommonMind (CMC) samples were sequenced under the BSM grant?

From Mette Peters' message from Feb 15, 2018, at 4:27 PM to Andy Chess:

>Jessica and I are working on a master table of all CMC samples. Can you please verify if these are [all samples generated through the BSMN grant on CMC individuals][link1]. Will there be additional samples?


```r
library(synapseClient)
```

```
## Loading required package: methods
```

```
## 
## TERMS OF USE NOTICE:
##     When using Synapse, remember that the terms and conditions of use require that you:
##     1) Attribute data contributors when discussing these data or results from these data.
##     2) Not discriminate, identify, or recontact individuals or groups represented by the data.
##     3) Use and contribute only data de-identified to HIPAA standards.
##     4) Redistribute data only under these same terms of use.
```

```r
synapseLogin()
```

```
## Goodbye.
```

```
## Welcome attilagk!
```


```r
cmc <- list(synapse = read.csv("Job-53521223405979528604345739.csv"))
cmc$chesslab <- read.csv("../../data/dnalib/BSM_Project_Chess.csv")
```

### Samples missing from the Synapse table


```r
dissections <- list(synapse = levels(cmc$synapse$Dissection_ID), chesslab = levels(cmc$chesslab$Institution.Dissection.ID))
(dissections.in.chesslab.not.in.synapse <- with(dissections, setdiff(chesslab, synapse)))
```

```
##  [1] ""                                 
##  [2] "CMC_PsychENCODE_MSSM_DLPFC_1150"  
##  [3] "CMC_PsychENCODE_MSSM_DLPFC_1155"  
##  [4] "CMC_PsychENCODE_MSSM_DLPFC_1206"  
##  [5] "CMC_PsychENCODE_MSSM_DLPFC_1225"  
##  [6] "CMC_PsychENCODE_MSSM_DLPFC_1247"  
##  [7] "CMC_PsychENCODE_MSSM_DLPFC_1255"  
##  [8] "CMC_PsychENCODE_MSSM_DLPFC_1268"  
##  [9] "CMC_PsychENCODE_MSSM_DLPFC_1272"  
## [10] "CMC_PsychENCODE_MSSM_DLPFC_1277"  
## [11] "CMC_PsychENCODE_MSSM_DLPFC_1300"  
## [12] "CMC_PsychENCODE_MSSM_DLPFC_1305"  
## [13] "CMC_PsychENCODE_MSSM_DLPFC_1333"  
## [14] "CMC_PsychENCODE_MSSM_DLPFC_1338"  
## [15] "CMC_PsychENCODE_MSSM_DLPFC_1346"  
## [16] "CMC_PsychENCODE_MSSM_DLPFC_1347"  
## [17] "CMC_PsychENCODE_MSSM_DLPFC_1357"  
## [18] "CMC_PsychENCODE_MSSM_DLPFC_1361"  
## [19] "CMC_PsychENCODE_MSSM_DLPFC_1375"  
## [20] "CMC_PsychENCODE_MSSM_DLPFC_1391"  
## [21] "CMC_PsychENCODE_MSSM_DLPFC_1417"  
## [22] "CMC_PsychENCODE_MSSM_DLPFC_1429"  
## [23] "CMC_PsychENCODE_PITT_DRPC10026"   
## [24] "CMC_PsychENCODE_PITT_DRPC1088"    
## [25] "CMC_PsychENCODE_PITT_DRPC1099"    
## [26] "CMC_PsychENCODE_PITT_DRPC1189"    
## [27] "CMC_PsychENCODE_PITT_DRPC1211"    
## [28] "CMC_PsychENCODE_PITT_DRPC1307"    
## [29] "CMC_PsychENCODE_PITT_DRPC1374"    
## [30] "CMC_PsychENCODE_PITT_DRPC1386"    
## [31] "CMC_PsychENCODE_PITT_DRPC1391"    
## [32] "CMC_PsychENCODE_PITT_DRPC1454"    
## [33] "CMC_PsychENCODE_PITT_DRPC1472"    
## [34] "CMC_PsychENCODE_PITT_DRPC1542"    
## [35] "CMC_PsychENCODE_PITT_DRPC1558"    
## [36] "CMC_PsychENCODE_PITT_DRPC547"     
## [37] "CMC_PsychENCODE_PITT_DRPC625"     
## [38] "CMC_PsychENCODE_PITT_DRPC700"     
## [39] "CMC_PsychENCODE_PITT_DRPC727"     
## [40] "CMC_PsychENCODE_PITT_DRPC739"     
## [41] "CMC_PsychENCODE_PITT_DRPC822"     
## [42] "CMC_PsychENCODE_PITT_DRPC917"     
## [43] "CMC_PsychENCODE_PITT_DRPC988"     
## [44] "CMC_PsychENCODE_PITT_PFC_BP_10003"
## [45] "CMC_PsychENCODE_PITT_PFC_BP_1086"
```

### Samples missing from the Chesslab table


```r
(dissections.in.synapse.not.in.chesslab <- with(dissections, setdiff(synapse, chesslab)))
```

```
##  [1] "CMC_PsychENCODE_Pitt_DRPC1454" "CMC_PsychENCODE_Pitt_DRPC1472"
##  [3] "CMC_PsychENCODE_Pitt_DRPC1542" "CMC_PsychENCODE_Pitt_DRPC1558"
##  [5] "CMC_PsychENCODE_Pitt_DRPC547"  "CMC_PsychENCODE_Pitt_DRPC625" 
##  [7] "CMC_PsychENCODE_Pitt_DRPC700"  "CMC_PsychENCODE_Pitt_DRPC727" 
##  [9] "CMC_PsychENCODE_Pitt_DRPC739"  "CMC_PsychENCODE_Pitt_DRPC917" 
## [11] "CMC_PsychENCODE_Pitt_DRPC988"
```

Based on the `Dissection_ID` let's inspect the samples in `synapse` but not in `chesslab`!  In particular, which `Brain_Bank` did they come from?


```r
subset(cmc$synapse, subset = Dissection_ID %in% dissections.in.synapse.not.in.chesslab, select =  Brain_Bank, drop = TRUE)
```

```
##  [1] PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT
## [15] PITT PITT PITT PITT PITT PITT PITT PITT PITT PITT
## Levels: MSSM PITT
```

The all came from `PITT`.

Now let's see the rest of the samples in `chesslab`.  Are they all from `MSSM`?


```r
all("MSSM" == subset(cmc$synapse, subset = ! Dissection_ID %in% dissections.in.synapse.not.in.chesslab, select =  Brain_Bank, drop = TRUE))
```

```
## [1] TRUE
```

Yes.

[link1]: https://urldefense.proofpoint.com/v2/url?u=https-3A__www.synapse.org_-23-21Synapse-3Asyn11801978_tables_query_eyJzcWwiOiJTRUxFQ1QgKiBGUk9NIHN5bjExODAxOTc4IE9SREVSIEJZIFwiVGlzc3VlXCIgQVNDLCBcIkFTU0FZXCIgQVNDLCBcIkRpc3NlY3Rpb25fSURcIiBBU0MiLCAic2VsZWN0ZWRGYWNldHMiOlt7ImNvbmNyZXRlVHlwZSI6Im9yZy5zYWdlYmlvbmV0d29ya3MucmVwby5tb2RlbC50YWJsZS5GYWNldENvbHVtblZhbHVlc1JlcXVlc3QiLCAiY29sdW1uTmFtZSI6IkZ1bmRpbmciLCAiZmFjZXRWYWx1ZXMiOlsiVTAxTUgxMDY4OTEiXX1dLCAiaW5jbHVkZUVudGl0eUV0YWciOnRydWUsICJpc0NvbnNpc3RlbnQiOnRydWUsICJvZmZzZXQiOjAsICJsaW1pdCI6MjV9&d=DwMFaQ&c=shNJtf5dKgNcPZ6Yh64b-A&r=zYUNWLTGp4iMUkBLv5yRxoKLmEITRT5pxh58JUSn5VM&m=BCAIE0gpekjgq-YM68EFSGNwnAM_GdlTLGh4XLq0USI&s=RM7G5AwlrlKDSNj_XJ6CZCb9pNrPPhdWadzZ2-ZQjM4&e=
<!-- MathJax scripts -->
<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
