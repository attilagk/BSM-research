library(VariantAnnotation)

# calculates partition sizes based on the Venn diagram of a list of callsets
#
# vcfs
partition.sizes <- function(vcfs = cs.vcf[["100MB"]]$snvs) {
    vp <- get.venn.partitions(vcfs)
    callsets.containing.partition <- apply(as.matrix(vp[ , seq.int(length(vcfs))]), 1, sum)
    df <- data.frame(callsets.containing.partition = callsets.containing.partition, calls.in.partition = vp$..count..)
    df[with(df, order(callsets.containing.partition, calls.in.partition)), ]
}

