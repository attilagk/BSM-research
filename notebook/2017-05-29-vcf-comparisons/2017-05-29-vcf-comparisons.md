

## Goals

Assess concordance of call sets

1. between callers
   * compare call set by mutect2 and that by strelka for a given tissue-pair, variant type (snvs or indels), and filter setting
1. between reference tissues
   * compare call set with muscle reference to that with NeuN_mn reference given a variant type
   * use call sets obtained by mutect2 $\cap$ strelka given a tissue pair and variant type and filter setting

## Results

### Computation and summary

The `myisec.sh` script does all the job: various filtering, conversion, indexing and comparison using `bcftools`, as well as summarizing the results by calculating the size of call sets in files named `callset-sizes.tsv`, which are printed below.  The only filter setting tested was the *PASS* value of mutect2; otherwise unfiltered sets are reported.  Note that none of the strelka records have *PASS* filter value.


```bash
./myisec.sh $HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB
./myisec.sh $HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB PASS
find results/ -name callset-sizes.tsv | xargs head
```

```
## ==> results/mutect2-PASS/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 19	0000.vcf	private to	mutect2.bcf
## 92	0001.vcf	private to	strelka.bcf
## 1	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 402	0000.vcf	private to	mutect2.bcf
## 551	0001.vcf	private to	strelka.bcf
## 41	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 17	0000.vcf	private to	mutect2.bcf
## 99	0001.vcf	private to	strelka.bcf
## 0	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-PASS/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 408	0000.vcf	private to	mutect2.bcf
## 591	0001.vcf	private to	strelka.bcf
## 42	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-PASS/2_cmp-reftissues/indels/callset-sizes.tsv <==
## 1	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 0	0001.vcf	private to	muscle-NeuN_pl.bcf
## 0	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
## 
## ==> results/mutect2-PASS/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## 10	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 11	0001.vcf	private to	muscle-NeuN_pl.bcf
## 31	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 710	0000.vcf	private to	mutect2.bcf
## 80	0001.vcf	private to	strelka.bcf
## 13	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 11251	0000.vcf	private to	mutect2.bcf
## 258	0001.vcf	private to	strelka.bcf
## 334	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 693	0000.vcf	private to	mutect2.bcf
## 90	0001.vcf	private to	strelka.bcf
## 9	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-unfilt/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 11274	0000.vcf	private to	mutect2.bcf
## 265	0001.vcf	private to	strelka.bcf
## 368	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/mutect2-unfilt/2_cmp-reftissues/indels/callset-sizes.tsv <==
## 7	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 3	0001.vcf	private to	muscle-NeuN_pl.bcf
## 6	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
## 
## ==> results/mutect2-unfilt/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## 109	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 143	0001.vcf	private to	muscle-NeuN_pl.bcf
## 225	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
```

Import sizes of various callsets to R:


```r
indirs <- paste0("results/mutect2-", c("unfilt", "PASS"), "/")
indirs <- paste0(rep(indirs, each = 3), c(paste0("1_isec-callers/", c("NeuN_mn-NeuN_pl", "muscle-NeuN_pl")), "2_cmp-reftissues"))
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

Without filtering mutect2 calls many more variants than strelka.  With filtering (PASS) mutect2 still calls slightly more variants than stelka without filtering.  The overlap between the two callers is substantial and filtering makes mutect2 calls even more concordant with strelka calls.


```r
my.venn3(dir = "1_isec-callers/NeuN_mn-NeuN_pl/snvs/", cls = clsets)
```

<img src="figure/venn-caller-neun_mn-snvs-1.png" title="plot of chunk venn-caller-neun_mn-snvs" alt="plot of chunk venn-caller-neun_mn-snvs" width="400px" />

#### NeuN_mn reference tissue, indels

Callsets of indels share some but not all of the above properties with SNV callsets.  Importantly, the concordance seems lower for indels than for NSPs, but the much lower number of indel records leaves a greater room for certain patterns emerging by chance.

<img src="figure/venn-caller-neun_mn-indels-1.png" title="plot of chunk venn-caller-neun_mn-indels" alt="plot of chunk venn-caller-neun_mn-indels" width="400px" />


#### muscle reference tissue, SNVs

These results are very similar to the previous ones; the muscle reference gives very slightly more candidate variants than NeuN_mn reference, which aligns with the expectation that muscle cell lineage is further away from the NeuN_pl lineage compared to NeuN_mn.  But the differences are very subtle.

<img src="figure/venn-caller-muscle-snvs-1.png" title="plot of chunk venn-caller-muscle-snvs" alt="plot of chunk venn-caller-muscle-snvs" width="400px" />

#### muscle reference tissue, indels

<img src="figure/venn-caller-muscle-indels-1.png" title="plot of chunk venn-caller-muscle-indels" alt="plot of chunk venn-caller-muscle-indels" width="400px" />

### Comparing reference tissues: unfiltered

The next Venn-Euler diagrams feature two sets corresponding to muscle and NeuN_mn reference tissue.  Both sets are **intersections** of mutect2 and strelka call sets.  This first set of diagrams was obtained **without** filtering mutect2 calls.  The overlap is relatively larger than either set difference---records that are private to either muscle or NeuN_mn reference.  This suggests that the choice of reference tissue is not critical.  Fewer records are private to NeuN_mn than to muscle.  This is again consistent with the muscle lineage being more distant to the NeuN_pl lineage.

### SNVs

<img src="figure/venn-ref-tissue-snvs-1.png" title="plot of chunk venn-ref-tissue-snvs" alt="plot of chunk venn-ref-tissue-snvs" width="400px" />

#### indels

<img src="figure/venn-ref-tissue-indels-1.png" title="plot of chunk venn-ref-tissue-indels" alt="plot of chunk venn-ref-tissue-indels" width="400px" />

### Comparing reference tissues: filtered

The next set of diagrams was obtained **with** filtered mutect2 calls.  These results show the same qualitative picture as the ones without filtering.

### SNVs

<img src="figure/venn-ref-tissue-snvs-PASS-1.png" title="plot of chunk venn-ref-tissue-snvs-PASS" alt="plot of chunk venn-ref-tissue-snvs-PASS" width="400px" />

#### indels

<img src="figure/venn-ref-tissue-indels-PASS-1.png" title="plot of chunk venn-ref-tissue-indels-PASS" alt="plot of chunk venn-ref-tissue-indels-PASS" width="400px" />
