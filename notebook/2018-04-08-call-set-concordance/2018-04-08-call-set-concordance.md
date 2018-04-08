---
layout: default
title: Call set concordance---test data
tags: [ strelka2, mutect2, tnseq, lofreq, somaticsniper ]
#featimg: "runtime-length-2.png"
---

Attempt to represent the concordance among call sets corresponding to various variant callers (strelka2Somatic,strelka2Germline, Tnseq/Mutect2, lofreq, somaticSniper.




```r
library(VariantAnnotation)
library(VennDiagram)
```


```r
callers <- c("lofreqSomatic", "somaticSniper", "strelka2Germline", "strelka2Somatic", "Tnseq")
names(callers) <- sub("Tnseq", "Tnseq.Mutect2", callers)
set.types <- paste(rep(c("common-sample", "benchmark"), each = 2), c("snvs", "indels"), sep = "-")
names(set.types) <- gsub("-", ".", set.types)
fl <- lapply(set.types, paste0, "/", callers, ".vcf.gz")
fl <- lapply(fl, `names<-`, names(callers))
vcf <- lapply(fl[1:2], lapply, readVcf, "hg19")
```


```r
grid.draw(venn.diagram(lapply(vcf$common.sample.snvs, function(x) row.names(x)), NULL))
```

<img src="figure/venn-common-sample-snvs-1.png" title="plot of chunk venn-common-sample-snvs" alt="plot of chunk venn-common-sample-snvs" width="600px" />


```r
grid.draw(venn.diagram(lapply(vcf$common.sample.indels, function(x) row.names(x)), NULL))
```

<img src="figure/venn-common-sample-indels-1.png" title="plot of chunk venn-common-sample-indels" alt="plot of chunk venn-common-sample-indels" width="600px" />
