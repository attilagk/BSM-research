---
layout: default
title: Call set concordance---test data
tags: [ strelka2, mutect2, tnseq, lofreq, somaticsniper ]
featimg: "venn-common-sample-wgs-snvs-1.png"
---

Attempt to represent the concordance among call sets corresponding to various variant callers (strelka2Somatic,strelka2Germline, Tnseq/Mutect2, lofreq, somaticSniper.





```r
library(VennDiagram)
source("2018-04-08-call-set-concordance.R")
```


```r
samples <- c("common.sample", "benchmark")
names(samples) <- samples
segs <- c(paste0(c(1, 3, 10, 30, 100), "MB"), "wgs")
names(segs) <- segs
vartypes <- c("snvs", "indels")
names(vartypes) <- vartypes
callers <- c("lofreqSomatic", "somaticSniper", "strelka2Germline", "strelka2Somatic", "Tnseq")
names(callers) <- sub("Tnseq", "Tnseq.Mutect2", callers)
```


```r
# common sample
cs.vcf <- import.across.segments(samples[1], sgs = segs)
```


```r
# benchmark study
bm.vcf <- import.across.segments(samples[2], sgs = segs[c(1, 3, 5, 6)])
```


```r
part.wgs <-
    rbind(cbind(data.frame(sample = "common.sample"), partition.sizes(cs.vcf[["wgs"]]$snvs)),
          cbind(data.frame(sample = "benchmark"), partition.sizes(bm.vcf[["wgs"]]$snvs)))
```


```r
part.100MB <-
    rbind(cbind(data.frame(sample = "common.sample"), partition.sizes(cs.vcf[["100MB"]]$snvs)),
          cbind(data.frame(sample = "benchmark"), partition.sizes(bm.vcf[["100MB"]]$snvs)))
```

## Number of unfiltered calls in the whole genome


```r
ssize <- rbind(cbind(data.frame(sample = "common.sample"), set.size.length(cs.vcf)),
               cbind(data.frame(sample = "benchmark"), set.size.length(bm.vcf)))
barchart(caller ~ set.size | sample, data = ssize, subset = length.Mb == 3235, xlim = c(0, 5e6))
```

<img src="figure/call-set-size-1.png" title="plot of chunk call-set-size" alt="plot of chunk call-set-size" width="600px" />

```r
barchart(caller ~ log10(set.size) | sample, data = ssize, subset = length.Mb == 3235, xlim = c(0, 7))
```

<img src="figure/call-set-size-2.png" title="plot of chunk call-set-size" alt="plot of chunk call-set-size" width="600px" />

## Concordance of unfiltered calls

### SNVs in the whole genome


```r
my.par <- list(main.cex = 1.8, fill = trellis.par.get("superpose.line")$col[seq_along(callers)], col = "gray", cat.cex = 1.4)
grid.draw(venn.diagram(cs.vcf[["wgs"]]$snvs, NULL, main = "common sample", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-common-sample-wgs-snvs-1.png" title="plot of chunk venn-common-sample-wgs-snvs" alt="plot of chunk venn-common-sample-wgs-snvs" width="600px" />


```r
grid.draw(venn.diagram(bm.vcf[["wgs"]]$snvs, NULL, main = "benchmark", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-benchmark-wgs-snvs-1.png" title="plot of chunk venn-benchmark-wgs-snvs" alt="plot of chunk venn-benchmark-wgs-snvs" width="600px" />


```r
xyplot(callers.in.partition ~ log10(calls.in.partition) | sample, data = part.wgs, pch = "|", col = "red", grid = TRUE, cex = 2)
```

<img src="figure/part-sizes-wgs-1.png" title="plot of chunk part-sizes-wgs" alt="plot of chunk part-sizes-wgs" width="600px" />

### Indels in the whole genome


```r
grid.draw(venn.diagram(cs.vcf[["wgs"]]$indels, NULL, main = "common sample", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-common-sample-wgs-indels-1.png" title="plot of chunk venn-common-sample-wgs-indels" alt="plot of chunk venn-common-sample-wgs-indels" width="600px" />


```r
grid.draw(venn.diagram(bm.vcf[["wgs"]]$indels, NULL, main = "benchmark", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-benchmark-wgs-indels-1.png" title="plot of chunk venn-benchmark-wgs-indels" alt="plot of chunk venn-benchmark-wgs-indels" width="600px" />

## Various genomic segments


```r
xyplot(log10(set.size) ~ log10(length.Mb) | sample, data = ssize, groups = caller, type = "b", lty = 2, auto.key = TRUE, grid = TRUE)
```

<img src="figure/set-size-seg-length-1.png" title="plot of chunk set-size-seg-length" alt="plot of chunk set-size-seg-length" width="600px" />

### SNVs in a 100Mb segment


```r
grid.draw(venn.diagram(cs.vcf[["100MB"]]$snvs, NULL, main = "common sample", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-common-sample-100MB-snvs-1.png" title="plot of chunk venn-common-sample-100MB-snvs" alt="plot of chunk venn-common-sample-100MB-snvs" width="600px" />


```r
grid.draw(venn.diagram(bm.vcf[["100MB"]]$snvs, NULL, main = "benchmark", main.cex = my.par$main.cex, fill = my.par$fill, col = my.par$col, cat.cex = my.par$cat.cex))
```

<img src="figure/venn-benchmark-100MB-snvs-1.png" title="plot of chunk venn-benchmark-100MB-snvs" alt="plot of chunk venn-benchmark-100MB-snvs" width="600px" />


```r
xyplot(callers.in.partition ~ log10(calls.in.partition) | sample, data = part.100MB, pch = "|", col = "red", grid = TRUE, cex = 2)
```

<img src="figure/part-sizes-100MB-1.png" title="plot of chunk part-sizes-100MB" alt="plot of chunk part-sizes-100MB" width="600px" />

