

## Definitions

**(all) kept variants**
* these were kept by the quality filter (see `PASS` in the `FILTER` field of the vcf file)
* they include specific and non-specific variants

**sample-specific variants**
* these are specific to either of the three samples
* thus they represent a certain type of somatic mozaicism (see *Discussion* below)

## Analysis

Extract the number of kept variants:


```bash
# after running "runme.sh"
cd ../../results/common-experiment/snv/2017-01-24-alison
BN="snv-2017-01-24"
grep '^After filtering.*Sites' ${BN}.log | sed 's/.*kept \([[:digit:]]\+\) out of.*/\1/' > kept.vars
cat ${BN}.log
```

```
## 
## VCFtools - UNKNOWN
## (C) Adam Auton and Anthony Marcketta 2009
## 
## Parameters as interpreted:
## 	--gzvcf snv-2017-01-24.vcf.gz
## 	--out snv-2017-01-24
## 	--singletons
## 	--keep-filtered PASS
## 
## Using zlib version: 1.2.8
## After filtering, kept 3 out of 3 Individuals
## Outputting Singleton Locations
## After filtering, kept 4210739 out of a possible 4896434 Sites
## Run Time = 32.00 seconds
```

Read number of kept variants and also read the sites of singletons:


```r
kept.vars <- scan(file = "../../results/common-experiment/snv/2017-01-24-alison/kept.vars")
singletons <- read.delim("../../results/common-experiment/snv/2017-01-24-alison/snv-2017-01-24.singletons")
head(singletons)
```

```
##   CHROM   POS SINGLETON.DOUBLETON ALLELE            INDV
## 1     1 28863                   S      A LBRCE-pfc-1b123
## 2     1 69063                   S      C LBRCE-pfc-1b123
## 3     1 69270                   S      A  Fibro_Common_7
## 4     1 76423                   S      A   Common_7_NeuN
## 5     1 83974                   S  AAAAG  Fibro_Common_7
## 6     1 83974                   S      A  Fibro_Common_7
```

Print first the number sample-specific variants:


```r
(singleton.freq <- with(singletons, table(INDV)))
```

```
## INDV
##   Common_7_NeuN  Fibro_Common_7 LBRCE-pfc-1b123 
##           13379           81162           15013
```

Now print their fraction in all 4.210739 &times; 10<sup>6</sup> kept variants:


```r
singleton.freq / kept.vars
```

```
## INDV
##   Common_7_NeuN  Fibro_Common_7 LBRCE-pfc-1b123 
##     0.003177352     0.019275001     0.003565407
```

The same information as a barchart:


```r
barchart(singleton.freq, main = "Sample-specific variants")
```

<img src="figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="700px" />

## Discussion
