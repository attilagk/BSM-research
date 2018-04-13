---
layout: default
title: "Comparing Callsets: Callers and Reference Tissues"
featimg: "venn-caller-neun_mn-neun_pl-snvs-1.png"
---

The MuTect2 and Strelka somatic variant callers are used here to produced unfiltered and heuristically filtered call sets.  Two pairwise tissue comparisons are used: (1) NeuN+ vs NeuN- reference or (2) NeuN+ vs muscle reference.  These result in substantially concordant call sets.



## Goals

Assess concordance of call sets

1. between callers
   * compare call set by mutect2 and that by strelka for a given tissue-pair, variant type (snvs or indels), and filter setting
1. between reference tissues
   * compare call set with muscle reference to that with NeuN_mn reference given a variant type
   * use call sets obtained by mutect2 $$\cap$$ strelka given a tissue pair and variant type and filter setting

## Results

### Computation and summary

The `myisec.sh` script does all the job: various filtering, conversion, indexing and comparison using `bcftools`, as well as summarizing the results by calculating the size of call sets in files named `callset-sizes.tsv`, which are printed below.  The only filter setting tested was the *PASS* value of mutect2; otherwise unfiltered sets are reported.  Note that none of the strelka records have *PASS* filter value.


```bash
./myisec.sh $HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB
./myisec.sh $HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB PASS
find results/ -name callset-sizes.tsv | xargs head
```

```
## ==> results/mutect2-unfilt/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## nrec	ABC	tissue_pair_A	tissue_pair_B	tissue_pair_C
## 120	100	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	010	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	001	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 225	110	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 23	101	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	011	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	111	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 
## ==> results/mutect2-unfilt/2_cmp-reftissues/indels/callset-sizes.tsv <==
## nrec	ABC	tissue_pair_A	tissue_pair_B	tissue_pair_C
## 2	100	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	010	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	001	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 6	110	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 1	101	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	011	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	111	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 11274	/0000.vcf	private to	/mutect2.bcf
## 265	/0001.vcf	private to	/strelka.bcf
## 368	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 693	/0000.vcf	private to	/mutect2.bcf
## 90	/0001.vcf	private to	/strelka.bcf
## 9	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_mn/snvs/callset-sizes.tsv <==
## 3882	/0000.vcf	private to	/mutect2.bcf
## 129	/0001.vcf	private to	/strelka.bcf
## 256	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_mn/indels/callset-sizes.tsv <==
## 191	/0000.vcf	private to	/mutect2.bcf
## 33	/0001.vcf	private to	/strelka.bcf
## 7	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 11251	/0000.vcf	private to	/mutect2.bcf
## 258	/0001.vcf	private to	/strelka.bcf
## 334	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 710	/0000.vcf	private to	/mutect2.bcf
## 80	/0001.vcf	private to	/strelka.bcf
## 13	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## nrec	ABC	tissue_pair_A	tissue_pair_B	tissue_pair_C
## 8	100	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	010	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	001	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 31	110	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 3	101	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	011	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	111	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 
## ==> results/mutect2-PASS/2_cmp-reftissues/indels/callset-sizes.tsv <==
## nrec	ABC	tissue_pair_A	tissue_pair_B	tissue_pair_C
## 0	100	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	010	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	001	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	110	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	101	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	011	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 0	111	muscle-NeuN_pl	NeuN_mn-NeuN_pl	muscle-NeuN_mn
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 408	/0000.vcf	private to	/mutect2.bcf
## 591	/0001.vcf	private to	/strelka.bcf
## 42	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 17	/0000.vcf	private to	/mutect2.bcf
## 99	/0001.vcf	private to	/strelka.bcf
## 0	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_mn/snvs/callset-sizes.tsv <==
## 57	/0000.vcf	private to	/mutect2.bcf
## 302	/0001.vcf	private to	/strelka.bcf
## 83	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_mn/indels/callset-sizes.tsv <==
## 5	/0000.vcf	private to	/mutect2.bcf
## 40	/0001.vcf	private to	/strelka.bcf
## 0	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 402	/0000.vcf	private to	/mutect2.bcf
## 551	/0001.vcf	private to	/strelka.bcf
## 41	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 19	/0000.vcf	private to	/mutect2.bcf
## 92	/0001.vcf	private to	/strelka.bcf
## 1	/0002.vcf	shared by both	/mutect2.bcf /strelka.bcf
```

Import sizes of various callsets to R:


```r
indirs <- paste0("results/mutect2-", c("unfilt", "PASS"), "/")
indirs <- paste0(rep(indirs, each = 3), c(paste0("1_isec-callers/", c("muscle-NeuN_pl", "NeuN_mn-NeuN_pl", "muscle-NeuN_mn"))))
indirs <- paste0(rep(indirs, each = 2), c("/indels/", "/snvs/"))
names(indirs) <- LETTERS[seq_along(indirs)]
clsets1 <-
    lapply(indirs, function(x) {
           df <- read.delim(paste0(x, "callset-sizes.tsv"), header = FALSE,
                            col.names = c("set.size", "file.vcf", "set.operator", "set.operand"))
           df$directory <- factor(x)
           return(df)
       })
clsets1 <- do.call(rbind, clsets1)
```

### Comparing callers

#### NeuN_mn--NeuN_pl, SNVs

Without filtering mutect2 calls many more variants than strelka.  With filtering (PASS) mutect2 still calls slightly more variants than stelka without filtering.  The overlap between the two callers is substantial and filtering makes mutect2 calls even more concordant with strelka calls.


```r
my.venn3a(dir = "1_isec-callers/NeuN_mn-NeuN_pl/snvs/", cls = clsets1)
```

<img src="figure/venn-caller-neun_mn-neun_pl-snvs-1.png" title="plot of chunk venn-caller-neun_mn-neun_pl-snvs" alt="plot of chunk venn-caller-neun_mn-neun_pl-snvs" width="400px" />

#### NeuN_mn--NeuN_pl, indels

Callsets of indels share some but not all of the above properties with SNV callsets.  Importantly, the concordance seems lower for indels than for NSPs, but the much lower number of indel records leaves a greater room for certain patterns emerging by chance.

<img src="figure/venn-caller-neun_mn-neun_pl-indels-1.png" title="plot of chunk venn-caller-neun_mn-neun_pl-indels" alt="plot of chunk venn-caller-neun_mn-neun_pl-indels" width="400px" />


#### muscle--NeuN_pl, SNVs

These results are very similar to the previous ones; the muscle reference gives very slightly more candidate variants than NeuN_mn reference, which aligns with the expectation that muscle cell lineage is further away from the NeuN_pl lineage compared to NeuN_mn.  But the differences are very subtle.

<img src="figure/venn-caller-muscle-neun_pl-snvs-1.png" title="plot of chunk venn-caller-muscle-neun_pl-snvs" alt="plot of chunk venn-caller-muscle-neun_pl-snvs" width="400px" />

#### muscle--NeuN_pl, indels

<img src="figure/venn-caller-muscle-neun_pl-indels-1.png" title="plot of chunk venn-caller-muscle-neun_pl-indels" alt="plot of chunk venn-caller-muscle-neun_pl-indels" width="400px" />

#### muscle--NeuN_mn, SNVs

<img src="figure/venn-caller-muscle-neun-mn-snvs-1.png" title="plot of chunk venn-caller-muscle-neun-mn-snvs" alt="plot of chunk venn-caller-muscle-neun-mn-snvs" width="400px" />

#### muscle--NeuN_mn, indels

<img src="figure/venn-caller-muscle-neun-mn-indels-1.png" title="plot of chunk venn-caller-muscle-neun-mn-indels" alt="plot of chunk venn-caller-muscle-neun-mn-indels" width="400px" />
