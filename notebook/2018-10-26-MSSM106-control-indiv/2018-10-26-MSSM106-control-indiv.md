---
layout: default
#title: 
tags: [ strelka2, mutect2, tnseq, lofreq, somaticsniper ]
featimg: "venn-common-sample-snvs-1.png"
---

A summary of called variants by multiple callers is presented here on control individual MSSM106.

## Title




```r
library(VennDiagram)
source("../../src/vcf.R")
```

Import VCF for each caller (separately for SNVs and indels) and store only the called variants' identity by position.


```r
results <- get.calls4indiv(indiv = "MSSM_106", callers = c("lofreqSomatic", "somaticSniper", "strelka2Germline2s", "strelka2Somatic", "TNseq"))
```


```r
barchart(caller ~ set.size | var.type, data = results$ssize, xlim = c(0, ceiling(max(ssize$set.size) / 1e6) * 1e6))
```

```
## Error in limitsFromLimitlist(have.lim = have.xlim, lim = xlim, relation = x.relation, : object 'ssize' not found
```

```r
barchart(caller ~ log10(set.size) | var.type, data = results$ssize, xlim = c(0, 7))
```

<img src="figure/call-set-size-1.png" title="plot of chunk call-set-size" alt="plot of chunk call-set-size" width="600px" />


```r
barchart(caller ~ log10(set.size) | var.type, data = results$ssize, xlim = c(0, 7))[1]
```

<img src="figure/call-set-size-snvs-1.png" title="plot of chunk call-set-size-snvs" alt="plot of chunk call-set-size-snvs" width="600px" />


```r
my.par <- list(main.cex = 1.8, fill = trellis.par.get("superpose.line")$col[seq_along(callers)], col = "gray", cat.cex = 1.4)
```

```
## Error in eval(expr, envir, enclos): object 'callers' not found
```

```r
grid.draw(venn.diagram(results$vcf$snvs, NULL, main = "SNVs", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

```
## Error in VennDiagram::draw.quintuple.venn(area1 = length(A), area2 = length(B), : object 'my.par' not found
```
