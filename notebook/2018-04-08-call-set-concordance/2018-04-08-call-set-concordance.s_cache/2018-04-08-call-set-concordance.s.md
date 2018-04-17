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





























