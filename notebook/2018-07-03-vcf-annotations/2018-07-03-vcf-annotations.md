---
layout: default
#title: 
featimg: density-1.png
---

Annotations in VCF files provide information for variant calling and filtering.  Therefore (joint) distribution of annotations across candidate variants (VCF records) is of interest.  Since each variant caller produces its own set of annotations, those annotations must be examined separately for each caller.

## Preparations



Making a list of data frames based on a list of callers.  Each data frame contains a random subsample of VCF records and a subset of VCF fields.


```bash
# SNVs
cd $HOME/projects/bsm/results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/snvs
for vcf in vcf/*.vcf.gz; do vcfinfo2tsv $vcf& done
# indels
cd $HOME/projects/bsm/results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/indels
for vcf in vcf/*.vcf.gz; do vcfinfo2tsv $vcf& done
```


```r
callers <- c("strelka2Somatic", "strelka2Germline", "lofreqSomatic", "somaticSniper", "Tnseq")
names(callers) <- callers
# helper function
import.subsample <- function(caller, ssize = 1e4, indir = "../../results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/snvs/") {
    df <- read.delim(paste0(indir, caller, ".vcf.gz.info.tsv"), na.strings = c("NA", "."))
    mycols <- sapply(df, function(y) ! all(is.na(y)))
    mysample <- sample(seq.int(nr <- nrow(df)), size = min(nr, ssize), replace = FALSE)
    df[mysample, mycols]
}
info <- list() # list of lists
info$snvs <- lapply(callers, import.subsample, ssize = 1e4, indir = "../../results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/snvs/")
info$indels <- lapply(callers, import.subsample, ssize = 1e4, indir = "../../results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/indels/")
```


```bash
cd $HOME/projects/bsm/results/2018-07-03-vcf-annotations/benchmark-mix1a-mix3a/
./parse-header-info
```

Reshaping data frames into long format.


```r
myreshape <- function(caller, infolist = info$snvs) {
    df <- infolist[[caller]]
    df <- df[sapply(df, is.numeric)]
    foo <- function(key) data.frame(Annotation = key, Value = df[[key]])
    l <- lapply(names(df), foo)
    do.call(rbind, l)
}
info.long <- list()
info.long$snvs <- lapply(callers, myreshape, info$snvs)
info.long$indels <- lapply(callers, myreshape, info$indels)
```

```
## Error in `[.default`(df, sapply(df, is.numeric)): invalid subscript type 'list'
```

## Marginal distributions

### SNVs


```r
lapply(callers, function(y) densityplot(~ Value | Annotation, data = info.long$snvs[[y]], pch = ".", main = y, scales = list(relation = "free")))
```

```
## $strelka2Somatic
```

<img src="figure/density-1.png" title="plot of chunk density" alt="plot of chunk density" width="600px" />

```
## 
## $strelka2Germline
```

<img src="figure/density-2.png" title="plot of chunk density" alt="plot of chunk density" width="600px" />

```
## 
## $lofreqSomatic
```

<img src="figure/density-3.png" title="plot of chunk density" alt="plot of chunk density" width="600px" />

```
## 
## $somaticSniper
```

<img src="figure/density-4.png" title="plot of chunk density" alt="plot of chunk density" width="600px" />

```
## 
## $Tnseq
```

<img src="figure/density-5.png" title="plot of chunk density" alt="plot of chunk density" width="600px" />

## Pairwise joint distributions

### SNVs


```r
mysplom <- function(caller, infolist = info$snvs) {
    df <- infolist[[caller]]
    df <- df[sapply(df, is.numeric)]
    splom(df, pch = '.', main = caller)
}
lapply(callers["somaticSniper" != callers], mysplom, info$snvs)
```

```
## $strelka2Somatic
```

<img src="figure/splom-1.png" title="plot of chunk splom" alt="plot of chunk splom" width="600px" />

```
## 
## $strelka2Germline
```

<img src="figure/splom-2.png" title="plot of chunk splom" alt="plot of chunk splom" width="600px" />

```
## 
## $lofreqSomatic
```

<img src="figure/splom-3.png" title="plot of chunk splom" alt="plot of chunk splom" width="600px" />

```
## 
## $Tnseq
```

<img src="figure/splom-4.png" title="plot of chunk splom" alt="plot of chunk splom" width="600px" />

