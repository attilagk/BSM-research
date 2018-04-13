---
layout: default
title: Dependence of Mutect2 LOD Score on Read Depth
tags: [ mutect2 ]
featimg: "lod-32-pass-curves-1.png"
---

Given some candidate somatic SNV the NLOD (normal log odds) score quantifies evidence for somatic mosaicism in the reference tissue---muscle or fibroblasts---(the "Normal" tissue) while the TLOD for that in NeuN+ nuclei ("Tumor" tissue).  Thus, the NLOD and TLOD score quantify two different subhypotheses: *mosaic reference tissue* and *mosaic NeuN+*.  How do these two LOD scores depend on the depth sequencing of the NeuN+ ("tumor") tissue when the depth for the reference tissue ("normal") is fixed?

## Preliminaries

To study the question the following analysis is carried out here:

1. take 1 or 32 Mb segments of reference tissue--NeuN+ BAM pairs produced by my me (from Chess lab)
1. subsample BAM for NeuN+ nuclei (but keep reference tissue BAM) at various sampling fractions
1. call variants with Mutect2 for each fraction
1. optionally filter candidate variants (VCF filter field)
1. inspect distribution of LOD scores

Notable differences between the current callset by me (Chess lab) and ([the one from the Park lab]({{ site.baseurl }}{% post_url /projects/bsm/2017-10-02-sept-freeze-snv %}))

|----------------------|-------------|----------|
|                      |  Chess lab  | Park lab |
|----------------------|-------------|----------|
|   caller             |  Mutect2    | Mutect   |
|----------------------|-------------|----------|
|   reference tissue   |  muscle     |fibroblast|
|----------------------|-------------|----------|
| decoy for GRCh37?    |   no decoy  |   yes    |
|----------------------|-------------|----------|
|analyzed segment      |  1 to 32 Mb |whole genome|
|----------------------|-------------|----------|
|   individual         |  MSSM179    | 5154 (Brain ID)   |
|----------------------|-------------|----------|

([As seen in a previous post]({{ site.baseurl }}{% post_url /projects/bsm/2017-05-03-strelka-mutect2-pilot %})) the 32 Mb segment analyzed by me (Chess lab) contains short burst of read-depth spikes suggesting repetitive region, which may complicate interpretation of variant calls.

Subsample BAM, run Mutect2 and parse NLOD and TLOD scores


```bash
indir="$HOME/projects/bsm/results/2017-10-05-readdepth-snvscore"
cd $indir
./doall.sh
```


```
## Loading required package: Matrix
```

```
## Loading required package: RColorBrewer
```

Read NLOD and TLOD scores for all levels of subsampling


```r
df <- lapply(list(seg1MB = "1MB", seg32MB = "32MB"), import.lod)
```

## Results

Recall ([see alignment stats]({{ site.baseurl }}{% post_url /projects/bsm/2017-05-24-alignment-stats %})) that NeuN+ nuclei were sequenced at $$\approx 140 \times$$ depth while muscle at $$\approx 35 \times$$ depth, which is $$\approx 0.25\times 140 \times$$ depth.  Thus subsampling the $$\approx 140 \times$$ BAM for NeuN+ at a sampling fraction of $$0.25$$ brings the depth for NeuN+ to the value of that for muscle.  Apart from 0.25 a number of smaller or larger fractions were also used to subsample the NeuN+ BAM.

### All candidates

The two plots below show how the LOD, separately for the *mosaic muscle* and *mosaic NeuN+* subhypothesis, depends on the fraction at which the NeuN+ BAM was subsampled.  The lower plot is just a finer $$y$$-scaling of the upper one.

When looking at the plots note that

* the total count of candidate variants increases with subsampling fraction
* both LOD distributions change with subsampling fraction but the way they change is the opposite:
  * while the distribution of NLOD (mosaic muscle) widens and shifts slightly right with decreasing fraction that of TLOD (mosaic NeuN+) gets narrower and left-shifted
* the distributions at fraction 1.00 show a somewhat similar pattern to the earlier result ([the one from the Park lab]({{ site.baseurl }}{% post_url /projects/bsm/2017-10-02-sept-freeze-snv %})) in that the *mosaic NeuN+* distribution has a stronger right tail (the extreme LOD scores are larger than those for *mosaic muscle* but this tendency is much stronger in the Park lab result than in the present case
  * that difference might be due to the different methodological details of generating BAM and VCF files (see table above)


```r
my.lab <- "fraction of reads in ~140x NeuN+ at constant ~35x muscle"
(tp <- useOuterStrips(histogram(~ LOD | frac * tissue, data = df$seg32MB, type = "count", between = list(y = 0.5), par.settings = list(par.main.text = list(font = 1, cex = 1.0)), scales = list(alternating = 1), main = my.lab)))
```

<img src="figure/lod-32-hist-1.png" title="plot of chunk lod-32-hist" alt="plot of chunk lod-32-hist" width="700px" />

```r
update(tp, ylim = c(0, 30))
```

<img src="figure/lod-32-hist-2.png" title="plot of chunk lod-32-hist" alt="plot of chunk lod-32-hist" width="700px" />

### Candidates `PASS`ing Mutect2 filters

A relatively small fraction (0.0308869) of calls have `FILTER = PASS` value.  The histograms of the corresponding LOD distributions (below) are therefore noisier than the ones above.  As much as the noise allows comparison, the LOD distribution for the two subhypotheses (mosaic muscle and mosaic NeuN+) are similar for `PASS`ing calls at all subsampling frequencies.


```r
(tp <- useOuterStrips(histogram(~ LOD | frac * tissue, data = df$seg32MB, type = "count", between = list(y = 0.5), par.settings = list(par.main.text = list(font = 1, cex = 1.0)), scales = list(alternating = 1), main = my.lab, subset = FILTER == "PASS")))
```

<img src="figure/lod-32-pass-hist-1.png" title="plot of chunk lod-32-pass-hist" alt="plot of chunk lod-32-pass-hist" width="700px" />

```r
update(tp, ylim = c(0, 10))
```

<img src="figure/lod-32-pass-hist-2.png" title="plot of chunk lod-32-pass-hist" alt="plot of chunk lod-32-pass-hist" width="700px" />

### A more detailed view

A more detailed view on the `PASS`ing calls in the 32Mb segment is presented below to look beyond the noisy histograms.  This shows that

1. higher depth for NeuN+ leads to more candidate SNVs
1. the candidate SNVs that are specific to higher NeuN+ depth tend to have larger TLOD score (LOD for *mosaic NeuN+*) than NLOD score (LOD for *mosaic muscle*) consistently with ([the result from the Park lab]({{ site.baseurl }}{% post_url /projects/bsm/2017-10-02-sept-freeze-snv %}))


```r
xyplot(LOD ~ as.numeric(as.character(frac)) | tissue, data = df$seg32MB, type = "b", groups = ID, subset = FILTER == "PASS", xlab = my.lab, ylim = c(0, 50))
```

<img src="figure/lod-32-pass-curves-1.png" title="plot of chunk lod-32-pass-curves" alt="plot of chunk lod-32-pass-curves" width="700px" />

The same representation for all calls in the 1MB segment shows a similar picture


```r
xyplot(LOD ~ as.numeric(as.character(frac)) | tissue, data = df$seg1MB, type = "b", groups = ID, xlab = my.lab, ylim = c(0, 35))
```

<img src="figure/lod-1MB-curves-1.png" title="plot of chunk lod-1MB-curves" alt="plot of chunk lod-1MB-curves" width="700px" />

## Conclusion

The result suggests that the right-shifted TLOD distribution relative to the NLOD distribution ([see the previous analysis]({{ site.baseurl }}{% post_url /projects/bsm/2017-10-02-sept-freeze-snv %})) is an artefact related to sequencing depth, rather than to a different frequency of true mosaic SNVs in the two tissues.  But the right shift is weaker for my results (from the Chess lab) than for those from the Park lab, which might be due to a number of differences in data analysis.
