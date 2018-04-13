---
layout: default
title: Test runs for the benchmark study
tags: [ strelka2, mutect2, tnseq ]
featimg: "runtimes-dotplot-1.png"
---

These testruns of variant callers involves a pair of BAM files for the benchmark study.  The pair is composed of Mix1A and Mix3A.




Parse log files and collect runtimes


```bash
cd ../../results/2018-03-20-benchmark-testdata/
./runtime-length-ncores 2>/dev/null 1>runtimes.csv
```

Import runtimes along with chromosomal segment length and other information


```r
runt <- read.csv("~/projects/bsm/results/2018-03-20-benchmark-testdata/runtimes.csv")
runt$length.MB <- as.integer(sub("([[:digit:]]+)MB", "\\1", as.character(runt$length)))
```

```
## Warning: NAs introduced by coercion
```

```r
# the whole genome is 3235 MB long according to https://en.wikipedia.org/wiki/Human_genome
runt$length.MB <- ifelse(is.na(runt$length.MB), 3235, runt$length.MB)
runt <- runt[with(runt, order(numcores, length.MB)), ]
levels(runt$caller) <- sub("Tnseq", "Tnseq.Mutect2", levels(runt$caller))
```

Below are the runtimes.  Note that all callers were run on 32 cores except for somaticSniper, which has only single threaded mode.


```r
dotplot(caller ~ runtime / 3600 | length, data = runt, scales = list(x = list(relation = "free")), grid = TRUE, xlab = "runtime, hours", layout = c(1, 3))
```

<img src="figure/runtimes-dotplot-1.png" title="plot of chunk runtimes-dotplot" alt="plot of chunk runtimes-dotplot" width="600px" />

Tnseq is much slower than strelka2Somatic or strelka2Germline on both the small and the WGS data set.  At the shallower coverage of the the difference is less striking although Tnseq is still slower; see [2018-02-22-ref-tissue-proj-testdata]({{ site.baseurl }}{% post_url 2018-02-22-ref-tissue-proj-testdata %}).
