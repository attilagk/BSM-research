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



get.calls4indiv <-
    function(indiv = "MSSM_106", callers = c("lofreqSomatic", "somaticSniper", "strelka2Germline2s", "strelka2Somatic", "TNseq")) {
        import.vcfs <- function(vartype = "snvs") {
            fl <- paste0("../../results/calls/", indiv, "/", vartype, "/", callers, ".vcf.gz")
            names(fl) <- names(callers)
            l <- lapply(fl, readVcf, "hg19")
            lapply(l, row.names)
        }
        names(callers) <- sub("T(n|N)seq", "T\\1seq.Mutect2", callers)
        vartypes <- c("snvs", "indels")
        names(vartypes) <- vartypes
        vcf <- lapply(vartypes, import.vcfs)
        part <- rbind(cbind(data.frame(vartype = "snvs"), partition.sizes(vcf$snvs)),
                      cbind(data.frame(vartype = "indels"), partition.sizes(vcf$indels)))
        ssize <- data.frame(t(sapply(vartypes, function(vartype) sapply(vcf[[vartype]], length))))
        ssize <- reshape(ssize, varying = names(ssize), v.names = "set.size",
                         timevar = "caller", times = names(ssize), direction = "long", idvar = "var.type", ids = vartypes)
        ssize$var.type <- factor(ssize$var.type, ordered = TRUE, levels = vartypes)
        results <- list(vcf = vcf, part = part, ssize = ssize)
        return(results)
    }
