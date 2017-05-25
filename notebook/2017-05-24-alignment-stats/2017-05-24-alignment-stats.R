get.samstats <- function(statsfile, field = "COV") {
    read.delim(pipe(paste0("grep ^", field, " ", statsfile)), header = FALSE)
}

get.fai <- function(fastaix = "~/data/GRCh37/dna/Homo_sapiens.GRCh37.dna.fa.fai") {
    fai <- read.delim(fastaix , header = FALSE, stringsAsFactors = FALSE)[1:2]
    names(fai) <- c("contig", "length")
    return(fai)
}

get.read.depths <- function(fai, tissues = c("NeuN_pl", "NeuN_mn", "muscle"), suffices = "1000",
                            prefix = "../../data/MSSM_179/aln-stats/MSSM179_", extension = ".bam.depth.") {
    helper <- function(tis, suf) {
        fname <- paste0(prefix, tis, extension, suf)
        df <- read.delim(fname, header = FALSE, col.names = c("contig", "pos", "depth"))
        df$tissue <- factor(tis, levels = tissues, ordered = TRUE)
        df$suffix <- factor(suf, levels = suffices, ordered = TRUE)
        return(df)
    }
    rdepth <-
        lapply(suffices,
               function(suf)
                   do.call(rbind, lapply(tissues, helper, suf = suf)))
    rdepth <- do.call(rbind, rdepth)
    rdepth$contig <- factor(rdepth$contig, level = fai$contig, ordered = TRUE)
    return(rdepth)
}


depth.plot <- function(x, chromosomes, ...)
    xyplot(depth ~ pos / 1e6 | contig, groups = tissue, data = x,
           subset = x$contig %in% chromosomes & x$tissue != "muscle",
           layout = c(1, length(chromosomes)), strip = FALSE, strip.left = TRUE,
           type = "p", pch = ".", auto.key = list(column = 3),
           key = simpleKey(text = tissues[-3], lines = TRUE, points = FALSE, columns = 2),
           xlab = "position (Mb)", grid = TRUE, par.settings = list(superpose.symbol = list(alpha = 0.5)),
           ...)
