---
layout: default
tags: [BAM, alignment, genome]
featimg: "coverage-wgs-1.png"
---

Coverage and other statistics on the aligned reads based on the CEPH mixture samples.


```
## Loading required package: Matrix
```

```
## Loading required package: RColorBrewer
```

The depth was obtained with `samtools depth` based on alignments such as `Mix1A.bam`.  Every $1000$th nucleotide position in the genome was retained and written in files such as `Mix1A.bam.1000.depth`.  The following expressions filter this by 10 to arrive at every $10,000$th position being retained.


```r
get.depth <- function(fname = "../../results/2018-03-12-ceph-alignment-stats/Mix1A.bam.1000.depth", filter.by = 10) {
    df <- read.delim(fname, header = FALSE, col.names = c("contig", "pos", "depth"))
    df$contig <- factor(df$contig)
    return(df[seq(from = 1, to = nrow(df), by = filter.by), ])
}
samples <- c("Mix1A", "Mix1B")
fnames <- paste0("../../results/2018-03-12-ceph-alignment-stats/", samples, ".bam.1000.depth")
names(fnames) <- samples
depths <- lapply(fnames, get.depth)
depths <- do.call(rbind, lapply(names(depths), function(x) cbind(data.frame(sample = x), depths[[x]])))
depths$contig <- factor(depths$contig, ordered = TRUE, levels = c(1:22, "X", "Y", "MT", levels(depths$contig)[c(23:82, 84)]))
```

The following empirical densities of depth are based on all contigs (chromosomes) and hence show overall coverage.


```r
densityplot(~ depth, data = depths, groups = sample, subset = depth < 300, plot.points = FALSE, ylab = "density", par.settings = list(superpose.line = list(col = c("magenta", "darkgreen"))), auto.key = TRUE)
```

<img src="figure/coverage-wgs-1.png" title="plot of chunk coverage-wgs" alt="plot of chunk coverage-wgs" width="700px" />

Here coverage is shown separately for contigs (only chromosomal contigs are displayed excluding the mitocondrium, the decoy sequence and other contigs).


```r
densityplot(~ depth | contig, data = depths, groups = sample, subset = depth < 300 & contig %in% c(as.character(1:22), "X", "Y"), plot.points = FALSE, ylab = "density", par.settings = list(superpose.line = list(col = c("magenta", "darkgreen"))), auto.key = TRUE, ylim = c(0, 0.03))
```

<img src="figure/coverage-contigs-1.png" title="plot of chunk coverage-contigs" alt="plot of chunk coverage-contigs" width="700px" />

The coverage of the Y chromosome is expanded below:


```r
densityplot(~ depth | contig, data = depths, groups = sample, subset = contig %in% c("Y"), plot.points = FALSE, ylab = "density", par.settings = list(superpose.line = list(col = c("magenta", "darkgreen"))), auto.key = TRUE, xlim = c(0, 35))
```

<img src="figure/coverage-Y-1.png" title="plot of chunk coverage-Y" alt="plot of chunk coverage-Y" width="700px" />

The coverage of the mitochondrium is shown here:


```r
densityplot(~ depth | contig, data = depths, groups = sample, subset = contig %in% c("MT"), plot.points = FALSE, ylab = "density", par.settings = list(superpose.line = list(col = c("magenta", "darkgreen"))), auto.key = TRUE)
```

<img src="figure/coverage-MT-1.png" title="plot of chunk coverage-MT" alt="plot of chunk coverage-MT" width="700px" />

The depth is roughly evenly distributed along and across chromosomes, except---as expected---at pericentromeric repeat regions.


```r
useOuterStrips(xyplot(depth ~ pos | sample * contig, data = depths, subset = contig %in% c(as.character(1:22), "X", "Y"), ylim = c(0, 300), pch = "."))
```

<img src="figure/depth-contigs-1.png" title="plot of chunk depth-contigs" alt="plot of chunk depth-contigs" width="700px" />
