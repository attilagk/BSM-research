---
layout: default
featimg: "tot-dna-histogram-1.png"
---

Validation of somatic variants called from next generation sequencing data is possible with ddPCR.  But ddPCR requires substantial DNA.  Here our DNA libraries are analyzed in terms of their measured total DNA content (in nanograms).  The result is a ranking of study individuals according to DNA content based on all combinations of brain samples and nucleus types (either NeuN+ and NeuN-) and all replicates of DNA libraries.

## Total DNA / library



Load table of DNA libraries and summarize the distribution of total DNA (in ng) across libraries for a given tissue sample type (muscle, NeuN+, NeuN-).  Libraries from muscle tend to contain a lot more DNA than those from NeuN+ and NeuN- nuclei:


```r
dnalib <- read.csv("~/projects/bsm/data/dnalib/BSM_Project_Chess.csv")
l <- lapply(list(muscle = "mu", neunpos = "np", neunneg = "nn"),
            function(x)
                with(dnalib, subset(dnalib, subset = Sample == x,
                                    select = c("Library.name", "Individual.ID", "Sample", "Dx", "Total.DNA..ng."))))
lapply(l, function(x) summary(x$Total.DNA..ng.))
```

```
## $muscle
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    4450    6712    8650   12757   12050   53550 
## 
## $neunpos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##   24.65  690.81 1126.75 1099.14 1468.62 2565.00       4 
## 
## $neunneg
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     0.0   587.8   927.5   954.7  1322.4  2265.0       4
```


```r
histogram(~ Total.DNA..ng. | Sample, data = dnalib, subset = Sample  %in% c("np", "nn"))
```

![plot of chunk tot-dna-histogram](figure/tot-dna-histogram-1.png)

Get best library replicate within an individual


```r
lb <- lapply(l, best.replicate)
indiv <- intersect(as.character(lb$neunp$Individual.ID), as.character(lb$neunn$Individual.ID))
lbest <- cbind(lb$neunneg[indiv, c(1, 4, 5) ], lb$neunpos[indiv, c(1, 4, 5)])
names(lbest) <- paste0(rep(c("neunneg.", "neunpos."), each = 3), rep(c("Library.name", "Dx", "Total.DNA..ng."), 2))
lbest$neunneg.Library.name <- sub("^[^.]+\\.", "", lbest$neunneg.Library.name)
lbest$neunpos.Library.name <- sub("^[^.]+\\.", "", lbest$neunpos.Library.name)
```

Look at the correlation between NeuN+ and NeuN- DNA content:


```r
xyplot(neunneg.Total.DNA..ng. ~ neunpos.Total.DNA..ng., data = lbest)
```

![plot of chunk neu-neg-pos-dna-correlation](figure/neu-neg-pos-dna-correlation-1.png)

Sort individuals according to the minimum DNA within the NeuN- and NeuN+ sample (using the best library replicate given each tissue type) and look at the best individuals:


```r
lbest$min.dna <- apply(as.matrix(lbest[c(3, 6)]), 1, min)
lbest.ordered <- lbest[order(lbest$min.dna, decreasing = TRUE), sel.col <- c("neunneg.Dx", "neunneg.Library.name", "neunpos.Library.name", "min.dna")]
head(lbest.ordered, n = 10L)
```

```
##              neunneg.Dx neunneg.Library.name neunpos.Library.name min.dna
## CMC_MSSM_363        SCZ       DLPFC_1144.nn2       DLPFC_1144.np1  2090.0
## CMC_MSSM_364        SCZ       DLPFC_1141.nn1       DLPFC_1141.np2  1720.0
## CMC_MSSM_281    Control       DLPFC_1358.nn2       DLPFC_1358.np1  1593.0
## CMC_PITT_123    Control          DRPC988.nn1          DRPC988.np2  1570.0
## CMC_MSSM_185    Control       DLPFC_1430.nn1       DLPFC_1430.np1  1475.0
## CMC_MSSM_347        SCZ       DLPFC_1353.nn1       DLPFC_1353.np1  1445.0
## CMC_MSSM_130    Control       DLPFC_1150.nn1       DLPFC_1150.np2  1443.5
## CMC_MSSM_142        SCZ       DLPFC_1286.nn2       DLPFC_1286.np1  1440.0
## CMC_MSSM_107    Control       DLPFC_1400.nn2       DLPFC_1400.np2  1415.0
## CMC_MSSM_352        SCZ       DLPFC_1440.nn1       DLPFC_1440.np1  1390.0
```

The "worst individuals":


```r
tail(lbest.ordered)
```

```
##              neunneg.Dx neunneg.Library.name neunpos.Library.name min.dna
## CMC_MSSM_222        SCZ       DLPFC_1441.nn1       DLPFC_1441.np1  381.50
## CMC_MSSM_297        SCZ       DLPFC_1172.nn1       DLPFC_1172.np1  343.50
## CMC_MSSM_065    Control       DLPFC_1334.nn1       DLPFC_1334.np1  332.50
## CMC_MSSM_061    Control       DLPFC_1188.nn1       DLPFC_1188.np1  188.75
## CMC_MSSM_305        SCZ       DLPFC_1160.nn1       DLPFC_1160.np1  184.50
## CMC_MSSM_133        SCZ       DLPFC_1145.nn2       DLPFC_1145.np1   38.65
```

MSSM_179 (the first individual analyzed in the pilot phase) is in the middle:

```r
lbest.ordered["CMC_MSSM_179", ]
```

```
##              neunneg.Dx neunneg.Library.name neunpos.Library.name min.dna
## CMC_MSSM_179    Control       DLPFC_1224.nn1       DLPFC_1224.np1     832
```
