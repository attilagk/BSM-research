---
layout: default
title: Test runs on the common sample BAMs
tags: [ strelka2 ]
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
```

The following plots show how runtime scales with the size of the data.


```r
arg <- list(xlab = "length (bases)", ylab = "real time",
            labels = c("1s", "7.8s", "1m", "7.8m", "1h", "7.8h", "2.5d", "2.8w", "5m"),
            at.lin = 0:8 / 2, at.log = 60 ^ c(0:8 / 2))
xyplot(runtime / 60 ~ length.MB, data = runt, groups = caller, subset = numcores == 6 & length != "wgs", ylab = "runtime, min", xlab = "length, MB", grid = TRUE, type = "b", auto.key = list(lines = FALSE, points = TRUE), main = "Runtime at 6 cores", pch = 21, lty = 2)
```

<img src="figure/runtime-length-1.png" title="plot of chunk runtime-length" alt="plot of chunk runtime-length" width="600px" />

```r
xyplot(log(runtime, base = 60) ~ log10(length.MB * 1e6), data = runt, groups = caller, subset = numcores == 6, scales = list(y = list(at = arg$at.lin, labels = arg$labels)), ylim = c(1, 3.0), grid = TRUE, ylab = "runtime (log scale)", xlab = "log length (base)", type = "b", auto.key = list(lines = FALSE, points = TRUE), main = "Runtime at 6 cores", pch = 21, lty = 2)
```

<img src="figure/runtime-length-2.png" title="plot of chunk runtime-length" alt="plot of chunk runtime-length" width="600px" />

The next graph presents runtimes at different number of CPU cores.

<img src="figure/runtime-ncores-1.png" title="plot of chunk runtime-ncores" alt="plot of chunk runtime-ncores" width="600px" />

The number of all SNV candidates in the whole genome by strelka2Somatic caller and of those that `PASS`ed the default filters:


```bash
vcf="../../results/2018-02-22-ref-tissue-proj-testdata/WGS/6proc/strelka2Somatic/results/variants/somatic.snvs.vcf.gz"
zcat $vcf | sed -n '/^#CHROM/,$ p' | wc -l
zcat $vcf | grep -c PASS
```

```
## gzip: ../../results/2018-02-22-ref-tissue-proj-testdata/WGS/6proc/strelka2Somatic/results/variants/somatic.snvs.vcf.gz: No such file or directory
## 0
## gzip: ../../results/2018-02-22-ref-tissue-proj-testdata/WGS/6proc/strelka2Somatic/results/variants/somatic.snvs.vcf.gz: No such file or directory
## 0
```

## Caller specific notes

### Sentieon Tnseq

This caller is a reimplementation of GATK's mutect2 and like GATK has strict requirements for read groups (RG) in the input BAM files.  Our input files (`Common_7_NeuN_DO16090243-final-all.bam` and `Fibro_Common_7_DO16090243-final-all.bam`) do not fulfill the requirements and therefore it was necessary to create new input files with appropriately altered read groups.  See `correct-commons-bam` that implements the correction and produces output file `filename-correct-rg.bam` from `filename.bam`.

The following result shows that on the 1MB long data Tnseq and GATK v3.8.0 mutect2 (both with default settings) produce identical results


```bash
cd ../../results/2018-02-22-ref-tissue-proj-testdata/1MB/4proc/isec
cat README.txt
echo -e "\ncounting records in each VCF file:"
grep -c '^[^#]' 000?.vcf
```

```
## This file was produced by vcfisec.
## The command line was:	bcftools isec  -p isec mutect2/out.vcf.gz Tnseq/tnseq.vcf.gz
## 
## Using the following file names:
## isec/0000.vcf	for records private to	mutect2/out.vcf.gz
## isec/0001.vcf	for records private to	Tnseq/tnseq.vcf.gz
## isec/0002.vcf	for records from mutect2/out.vcf.gz shared by both	mutect2/out.vcf.gz Tnseq/tnseq.vcf.gz
## isec/0003.vcf	for records from Tnseq/tnseq.vcf.gz shared by both	mutect2/out.vcf.gz Tnseq/tnseq.vcf.gz
## 
## counting records in each VCF file:
## 0000.vcf:0
## 0001.vcf:0
## 0002.vcf:21
## 0003.vcf:21
```

Interestingly, GATK v3.8.0 mutect2 ran successfully both with and without the read group correction implemented by `correct-commons-bam`. Right after this I switched to GATK v4.0.2.0.
