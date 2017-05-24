


```r
fai <- read.delim("Homo_sapiens.GRCh37.dna.fa.fai", header = FALSE, stringsAsFactors = FALSE)[1:2]
names(fai) <- c("contig", "length")
rdepth <- lapply(tissues <- c("NeuN_pl", "NeuN_mn", "muscle"),
                 function(x) {
                     read.delim("../../data/MSSM_179/aln-stats/MSSM179_muscle.bam.depth.1000", header = FALSE)
                     df <- read.delim(paste0("../../data/MSSM_179/aln-stats/MSSM179_", x , ".bam.depth.1000"),
                                      header = FALSE, col.names = c("contig", "pos", "depth"))
                     df$tissue <- factor(x, levels = tissues, ordered = TRUE)
                     return(df)
                 })
rdepth <- do.call(rbind, rdepth)
rdepth$contig <- factor(rdepth$contig, level = fai$contig, ordered = TRUE)
```

Read depth


```r
rdepth.100 <- rdepth[seq(1, nrow(rdepth), by = 100), ] # for test purposes
chromosomes <- grep("^GL", fai$contig, value = TRUE, invert = TRUE)
```

<img src="figure/read-depth-muscle-1.png" title="plot of chunk read-depth-muscle" alt="plot of chunk read-depth-muscle" width="700px" />

<img src="figure/read-depth-1.png" title="plot of chunk read-depth" alt="plot of chunk read-depth" width="700px" />
