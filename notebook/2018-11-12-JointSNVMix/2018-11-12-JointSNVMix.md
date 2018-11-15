---
layout: default
title: JointSNVMix2
tags: [ JointSNVMix2 ]
#featimg: "all-1.png"
---

Testing JointSNVMix and inspecting its output.

The `bash` code below shows how the `jsm.tsv` files have been generated.  See details in each `runme` file.


```bash
# run JointSNVMix2 
cd $HOME/projects/bsm/results/2018-11-12-JointSNVMix
for wd in chr22 1Mb 10Mb 100Mb; do
    $wd/runme & # asynchronous evaluation (each process runs on thread)
done
# ...wait for all processes to complete...
# downsample rows (genome positions): keep every 1000th
$HOME/projects/bsm/notebook/2018-11-12-JointSNVMix/sample-every chr22/jsm.tsv 1000
```


```bash
# run JointSNVMix2 
cd $HOME/projects/bsm/results/2018-11-12-JointSNVMix
head -n2 chr22/jsm.tsv.every-1000
```

```
## chrom	position	ref_base	var_base	normal_counts_a	normal_counts_b	tumour_counts_a	tumour_counts_b	p_AA_AA	p_AA_AB	p_AA_BB	p_AB_AA	p_AB_AB	p_AB_BB	p_BB_AA	p_BB_AB	p_BB_BB
## 22	16053036	A	C	96	0	117	1	1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000
```





```r
tsv <- "~/projects/bsm/results/2018-11-12-JointSNVMix/chr22/jsm.tsv.every-1000"
records <- read.delim(tsv)
# aggregate probabilities
# across genotypes for the case/tumor sample
records$p_AA_XX <- with(records, p_AA_AA + p_AA_AB + p_AA_BB)
records$p_AB_XX <- with(records, p_AB_AA + p_AB_AB + p_AB_BB)
records$p_BB_XX <- with(records, p_BB_AA + p_BB_AB + p_BB_BB)
# across genotypes for the control/normal sample
records$p_XX_AA <- with(records, p_AA_AA + p_AB_AA + p_BB_AA)
records$p_XX_AB <- with(records, p_AA_AB + p_AB_AB + p_BB_AB)
records$p_XX_BB <- with(records, p_AA_BB + p_AB_BB + p_BB_BB)
records.l <-
    reshape(records, varying = v <- grep("p_", names(records), value = TRUE),
            v.names = "p", timevar = "GT", times = sub("p_", "", v), direction = "long")
```


```r
histogram(~ p_AA_XX + p_AB_XX + p_BB_XX, data = records, layout = c(3, 1), main = "Control sample genotype", xlab = "posterior probability")
```

<img src="figure/controlGT-1.png" title="plot of chunk controlGT" alt="plot of chunk controlGT" width="600px" />


```r
histogram(~ p_XX_AA + p_XX_AB + p_XX_BB, data = records, layout = c(3, 1), main = "Case sample genotype", xlab = "posterior probability")
```

<img src="figure/caseGT-1.png" title="plot of chunk caseGT" alt="plot of chunk caseGT" width="600px" />
