import.across.segments <- function(sam = samples[1], sgs = segs, vart = vartypes, cal = callers) {
    across.vartypes <- function(seg) {
        across.callers <- function(var) {
            d <- paste0(sam, "/", seg, "/vcf/", var)
            fl <- paste0(d, "/", cal, ".vcf.gz")
            names(fl) <- sub("Tnseq", "Tnseq.Mutect2", cal)
            lapply(fl, readVcf, "hg19")
        }
        lapply(vart, across.callers)
    }
    lapply(sgs, across.vartypes)
}
