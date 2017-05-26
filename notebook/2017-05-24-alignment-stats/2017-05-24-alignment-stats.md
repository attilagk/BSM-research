


```r
tissues <- c("NeuN_pl", "NeuN_mn", "muscle")
fai <- get.fai()
chromosomes <- grep("^GL", fai$contig, value = TRUE, invert = TRUE)
rd.genome <- get.read.depths(fai, tissues)
rd.repres <- get.read.depths(fai, tissues, suffices = "chr1_119,500,001-121,500,000")
```


```r
rd.repres <-
    do.call(rbind,
            lapply(regions <- c("2Mb", "10kb baseline", "10kb spike"),
                   function(x) {
                       r <- rd.repres
                       r$region <- factor(x, levels = regions, ordered = TRUE)
                       r
                   }))
rd.repres.100 <- rd.repres[seq(1, nrow(rd.repres), by = 100), ] # for test purposes
rd.genome.100 <- rd.genome[seq(1, nrow(rd.genome), by = 100), ] # for test purposes
fastaix <- "~/data/GRCh37/dna/Homo_sapiens.GRCh37.dna.fa.fai"
```


```r
depth.hist <-
    lapply(tissues,
           function(x) {
               df <- get.samstats(paste0("../../results/2017-05-24-alignment-stats/MSSM179_", x, ".bam.stats"))
               df <- df[3:4]
               names(df) <- c("depth", "frequency")
               df$tissue <- factor(x, levels = tissues, ordered = TRUE)
               return(df)
           })
depth.hist <- do.call(rbind, depth.hist)
xyplot(frequency ~ depth | tissue, groups = tissue, data = depth.hist, type = "s", xlim = c(-10, 300), layout = c(1, 3), grid = TRUE)
```

<img src="figure/depth-hist-1.png" title="plot of chunk depth-hist" alt="plot of chunk depth-hist" width="700px" />

<img src="figure/read-depth-genome-log-1.png" title="plot of chunk read-depth-genome-log" alt="plot of chunk read-depth-genome-log" width="700px" />

<img src="figure/read-depth-genome-1.png" title="plot of chunk read-depth-genome" alt="plot of chunk read-depth-genome" width="700px" />

<img src="figure/read-depth-repres-1.png" title="plot of chunk read-depth-repres" alt="plot of chunk read-depth-repres" width="700px" /><img src="figure/read-depth-repres-2.png" title="plot of chunk read-depth-repres" alt="plot of chunk read-depth-repres" width="700px" />

<img src="figure/read-depth-nonchromosomal-log-1.png" title="plot of chunk read-depth-nonchromosomal-log" alt="plot of chunk read-depth-nonchromosomal-log" width="700px" />

<img src="figure/read-depth-nonchromosomal-1.png" title="plot of chunk read-depth-nonchromosomal" alt="plot of chunk read-depth-nonchromosomal" width="700px" />
