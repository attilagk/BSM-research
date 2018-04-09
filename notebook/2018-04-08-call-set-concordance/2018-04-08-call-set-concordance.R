library(VariantAnnotation)
library(vcfR)

# For a given sample import VCFs across all segments, variant types
# (snvs/indels) and variant callers.
# The import function is readVcf or read.vcfR of the VariantAnnotation and
# vcfR packages, respectively.
import.across.segments <- function(sam = samples[1], sgs = segs, vart = vartypes, cal = callers) {
    across.vartypes <- function(seg) {
        across.callers <- function(var) {
            d <- paste0(sam, "/", seg, "/vcf/", var)
            fl <- paste0(d, "/", cal, ".vcf.gz")
            names(fl) <- sub("Tnseq", "Tnseq.Mutect2", cal)
            #lapply(fl, read.vcfR, verbose = FALSE)
            l <- lapply(fl, readVcf, "hg19")
            lapply(l, row.names)
        }
        lapply(vart, across.callers)
    }
    lapply(sgs, across.vartypes)
}

set.size.length <- function(vcf = cs.vcf, vtype = "snvs") {
    seg.len <- as.numeric(sub("wgs", "3235", sub("MB", "", names(vcf))))
    df <- cbind(data.frame(length.Mb = seg.len),
                d2 <- data.frame(t(sapply(vcf, function(x) sapply(x[[vtype]], length)))))
    reshape(df, varying = names(d2), v.names = "set.size", timevar = "caller",
            times = names(d2), direction = "long", idvar = "length.Mb")
}

partition.sizes <- function(vcfs = cs.vcf[["100MB"]]$snvs) {
    vp <- get.venn.partitions(vcfs)
    as.matrix(vp[ , seq.int(length(vcfs))])
}
