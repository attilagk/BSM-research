
read1vcf <- function(indiv, subdir, names.only = TRUE) {
    postfix <-
        if (subdir == "S1")
            ""
        else
            paste0(".", subdir)
    fn <- paste0("greally-lab/", subdir, "/", indiv, "_S1", postfix, ".vcf.gz")
    rng <- GRanges(seqnames="chr22", IRanges(TRUE))
    #v <- readVcf(TabixFile(fn), "hg19", param = rng)
    v <- readVcf(fn, "hg19")
    if(names.only)
        return(row.names(v))
    else
        return(v)
}
