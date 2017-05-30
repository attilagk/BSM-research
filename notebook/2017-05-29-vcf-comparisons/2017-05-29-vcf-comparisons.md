

## Goal

Assess concordance of call sets

1. between callers
   * compare call set by mutect2 and that by strelka for a given tissue-pair and variant type (snvs or indels)
1. between reference tissues
   * compare call set with muscle reference to that with NeuN_mn reference given a variant type
   * use call sets obtained by mutect2 $\cap$ strelka given a tissue pair and variant type

The first

## Results


```bash
./myisec.sh &&
    find results/ -name callset-sizes.tsv | xargs head
```

```
## ==> results/1_isec-callers/NeuN_mn-NeuN_pl/indels/callset-sizes.tsv <==
## 711	0000.vcf	private to	mutect2.bcf
## 81	0001.vcf	private to	strelka.bcf
## 14	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/NeuN_mn-NeuN_pl/snvs/callset-sizes.tsv <==
## 11252	0000.vcf	private to	mutect2.bcf
## 259	0001.vcf	private to	strelka.bcf
## 335	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/muscle-NeuN_pl/indels/callset-sizes.tsv <==
## 694	0000.vcf	private to	mutect2.bcf
## 91	0001.vcf	private to	strelka.bcf
## 10	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/1_isec-callers/muscle-NeuN_pl/snvs/callset-sizes.tsv <==
## 11275	0000.vcf	private to	mutect2.bcf
## 266	0001.vcf	private to	strelka.bcf
## 369	0002.vcf	shared by both	mutect2.bcf strelka.bcf
## 
## ==> results/2_cmp-reftissues/indels/callset-sizes.tsv <==
## 8	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 4	0001.vcf	private to	muscle-NeuN_pl.bcf
## 7	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
## 
## ==> results/2_cmp-reftissues/snvs/callset-sizes.tsv <==
## 110	0000.vcf	private to	NeuN_mn-NeuN_pl.bcf
## 144	0001.vcf	private to	muscle-NeuN_pl.bcf
## 226	0002.vcf	shared by both	NeuN_mn-NeuN_pl.bcf muscle-NeuN_pl.bcf
```


```r
#grid.draw(draw.pairwise.venn(100, 70, 30, c("First", "Second"), scaled = FALSE))
```
