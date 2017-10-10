import.lod <- function(seg = "1MB") {
    frac <- rev(c("1.00", "0.50", "0.25", "0.12", "0.06", "0.03"))
    f <- paste0("~/projects/bsm/results/2017-10-05-readdepth-snvscore/MSSM179_NeuN_pl-", seg, "-", frac, "/mutect2/lod.csv")
    names(f) <- frac
    l <- lapply(f, read.csv)
    df <- do.call(rbind, lapply(names(l), function(x) cbind(data.frame(frac = x), l[[x]])))
    df <- cbind(data.frame(ID = factor(paste(df$CHROM, df$POS, sep = ":"))), df)
    df.long <- reshape(df, varying = c("NLOD", "TLOD"), v.names = "LOD", idvar = c("ID", "frac"), timevar = "tissue", times = paste("mosaic", c("muscle","NeuN+")), direction = "long")
    return(df.long)
}
