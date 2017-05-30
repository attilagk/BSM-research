

## Goals

Assess concordance of call sets

1. between callers
   * compare call set by mutect2 and that by strelka for a given tissue-pair and variant type (snvs or indels)
1. between reference tissues
   * compare call set with muscle reference to that with NeuN_mn reference given a variant type
   * use call sets obtained by mutect2 $\cap$ strelka given a tissue pair and variant type

The first

## Results

### Computation and summary


```bash
./myisec.sh &&
    find results/ -name callset-sizes.tsv | xargs head
```

```
## ==> results/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 710	0000.vcf	private to	mutect2.bcf
## 80	0001.vcf	private to	strelka.bcf
## 13	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 11251	0000.vcf	private to	mutect2.bcf
## 258	0001.vcf	private to	strelka.bcf
## 334	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 693	0000.vcf	private to	mutect2.bcf
## 90	0001.vcf	private to	strelka.bcf
## 9	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 11274	0000.vcf	private to	mutect2.bcf
## 265	0001.vcf	private to	strelka.bcf
## 368	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/2_cmp-reftissues/indels/callset-sizes.tsv <==
## 7	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 3	0001.vcf	private to	muscle-NeuN_pl.bcf
## 6	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
## 
## ==> results/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## 109	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 143	0001.vcf	private to	muscle-NeuN_pl.bcf
## 225	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
```


```r
indirs <- c(paste0("results/1_isec-callers/", c("NeuN_mn-NeuN_pl", "muscle-NeuN_pl")),
       "results/2_cmp-reftissues")
indirs <- paste0(rep(indirs, each = 2), c("/indels/", "/snvs/"))
names(indirs) <- LETTERS[seq_along(indirs)]
clsets <-
    lapply(indirs, function(x) {
           df <- read.delim(paste0(x, "callset-sizes.tsv"), header = FALSE,
                            col.names = c("set.size", "file.vcf", "set.operator", "set.operand"))
           df$directory <- factor(x)
           return(df)
       })
clsets <- do.call(rbind, clsets)
```

### Comparing callers

#### NeuN_mn reference tissue, SNVs


```r
my.grid.draw(dir = "results/1_isec-callers/NeuN_mn-NeuN_pl/snvs/", cls = clsets,
             category = c("mutect2", "strelka"), col = my.cols <- c("cyan", "magenta"), fill = my.cols,
             cat.pos = 165 * c(-1, 1))
```

<img src="figure/venn-caller-neun_mn-snvs-1.png" title="plot of chunk venn-caller-neun_mn-snvs" alt="plot of chunk venn-caller-neun_mn-snvs" width="500px" />

#### NeuN_mn reference tissue, indels

<img src="figure/venn-caller-neun_mn-indels-1.png" title="plot of chunk venn-caller-neun_mn-indels" alt="plot of chunk venn-caller-neun_mn-indels" width="500px" />


#### muscle reference tissue, SNVs

<img src="figure/venn-caller-muscle-snvs-1.png" title="plot of chunk venn-caller-muscle-snvs" alt="plot of chunk venn-caller-muscle-snvs" width="500px" />

#### muscle reference tissue, indels

<img src="figure/venn-caller-muscle-indels-1.png" title="plot of chunk venn-caller-muscle-indels" alt="plot of chunk venn-caller-muscle-indels" width="500px" />

### Comparing reference tissues: SNVs

<img src="figure/venn-ref-tissue-snvs-1.png" title="plot of chunk venn-ref-tissue-snvs" alt="plot of chunk venn-ref-tissue-snvs" width="500px" />

### Comparing reference tissues: indels

<img src="figure/venn-ref-tissue-indels-1.png" title="plot of chunk venn-ref-tissue-indels" alt="plot of chunk venn-ref-tissue-indels" width="500px" />
