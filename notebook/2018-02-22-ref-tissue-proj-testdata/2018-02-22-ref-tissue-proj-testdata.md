---
layout: default
title: Test runs on the common sample BAMs
tags: [ strelka2, mutect2, tnseq ]
featimg: "runtime-length-2.png"
---

Test runs of variant callers were performed using as input data the BAM files from the common sample from the reference tissue project.  Runtime were measured as a function of the length and number of CPU cores.

The analysis had the following major steps
1. take the BAMs for the paired common sample: `Common_7_NeuN_DO16090243-final-all.bam` and `Fibro_Common_7_DO16090243-final-all.bam` generated in Peter Park's group
1. excise segments of various lengths (between 1 and 100 MB) starting at chr1:50,000,000 using `do-aln-subregion` 
1. run callers e.g. using `do-mystrelka2Somatic` and store runtime results
1. parse runtime results using `runtime-length-ncores`



Parse the runtime files:


```bash
cd ../../results/2018-02-22-ref-tissue-proj-testdata/
./runtime-length-ncores 2>/dev/null 1>runtimes.csv
```

Import and prepare data


```r
runt <- read.csv("~/projects/bsm/results/2018-02-22-ref-tissue-proj-testdata/runtimes.csv")
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
# somatic sniper has no parallel mode
for(l in levels(runt$length))
    runt[with(runt, caller == "somaticSniper" & length == l & numcores != 1), "runtime"] <- subset(runt, caller == "somaticSniper" & length == l & numcores == 1, runtime, drop = TRUE)
```

The following plots show how runtime scales with the size of the data.












