---
layout: default
#title: 
tags: [ strelka2, mutect2, tnseq, lofreq, somaticsniper ]
featimg: "venn-snvs-1.png"
---

A summary of called variants by multiple callers is presented here on control individual MSSM106.

## Title




```r
library(VennDiagram)
source("../../src/vcf.R")
```

Import VCF for each caller (separately for SNVs and indels) and store only the called variants' identity by position.


```r
callers <- c("lofreqSomatic", "somaticSniper", "strelka2Germline2s", "strelka2Somatic", "TNseq")
unfiltered <- get.calls4indiv(indiv = "MSSM_106", callers = callers, PASS = FALSE)
PASS <- get.calls4indiv(indiv = "MSSM_106", callers = callers, PASS = TRUE)
```


```r
unfiltered$vp <- lapply(unfiltered$vcf, get.venn.partitions)
PASS$vp <- lapply(PASS$vcf, get.venn.partitions)
```


```r
ssize <- rbind(cbind(unfiltered$ssize, data.frame(filter = "unfiltered")),
               cbind(PASS$ssize, data.frame(filter = "PASS")))
my.barchart <- function(do.log = FALSE, do.groups = TRUE, ...) {
    fm <- if(do.log)
              formula(caller ~ log10(set.size) | var.type)
          else
              formula(caller ~ set.size | var.type)
    bc <- if(do.groups)
              barchart(fm, data = ssize, groups = filter)
          else
              barchart(fm, data = ssize, subset = ssize$filter == "unfiltered")
    xright <- if(do.log) 7 else ceiling(max(unfiltered$ssize$set.size) / 1e6) * 1e6
    update(bc, xlim = c(0, xright), ...)
}
my.barchart(do.log = FALSE, do.groups = TRUE, auto.key = list(points = FALSE, rectangles = TRUE))
```

<img src="figure/call-set-size-1.png" title="plot of chunk call-set-size" alt="plot of chunk call-set-size" width="600px" />

```r
my.barchart(do.log = TRUE, do.groups = TRUE, auto.key = list(points = FALSE, rectangles = TRUE))
```

<img src="figure/call-set-size-2.png" title="plot of chunk call-set-size" alt="plot of chunk call-set-size" width="600px" />

```r
#my.barchart(do.log = TRUE, do.groups = FALSE)
```


```r
my.barchart(do.log = TRUE, do.groups = FALSE, auto.key = list(points = FALSE, rectangles = TRUE))[1]
```

<img src="figure/call-set-size-snvs-1.png" title="plot of chunk call-set-size-snvs" alt="plot of chunk call-set-size-snvs" width="600px" />

```r
my.barchart(do.log = TRUE, do.groups = TRUE, auto.key = list(points = FALSE, rectangles = TRUE))[1]
```

<img src="figure/call-set-size-snvs-2.png" title="plot of chunk call-set-size-snvs" alt="plot of chunk call-set-size-snvs" width="600px" />


```r
my.par <- list(main.cex = 1.8, fill = trellis.par.get("superpose.line")$col[seq_along(callers)], col = "gray", cat.cex = 1.4)
grid.draw(venn.diagram(unfiltered$vcf$snvs, NULL, main = "SNVs", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-snvs-1.png" title="plot of chunk venn-snvs" alt="plot of chunk venn-snvs" width="600px" />

The intersection of all callers contains the following variant(s)


```r
unfiltered$vp$snvs[["..values.."]][1]
```

```
## $`1`
## [1] "X:140336646_G/T"
```

The Venn diagram for PASS filtered callsets:


```r
my.par <- list(main.cex = 1.8, fill = trellis.par.get("superpose.line")$col[seq_along(callers)], col = "gray", cat.cex = 1.4)
grid.draw(venn.diagram(PASS$vcf$snvs, NULL, main = "SNVs", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-snvs-PASS-1.png" title="plot of chunk venn-snvs-PASS" alt="plot of chunk venn-snvs-PASS" width="600px" />
