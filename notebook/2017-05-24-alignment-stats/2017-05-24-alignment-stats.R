library(latticeExtra)

get.samstats <- function(statsfile, field = "COV") {
    read.delim(pipe(paste0("grep ^", field, " ", statsfile)), header = FALSE)
}

get.fai <- function(fastaix = "~/data/GRCh37/karyotypic-order/Homo_sapiens.GRCh37.dna.fa.fai") {
    fai <- read.delim(fastaix , header = FALSE, stringsAsFactors = FALSE)[1:2]
    names(fai) <- c("contig", "length")
    return(fai)
}

get.read.depths <- function(fai, tissues = c("NeuN_pl", "NeuN_mn", "muscle"), suffices = "1000",
                            prefix = "~/results/bsm/2017-05-24-alignment-stats/MSSM179_", extension = ".bam.depth.") {
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


depth.plot <- function(x, sel.contigs, chromosomal = TRUE, ...) {
    tp <- xyplot(depth ~ pos / ifelse(chromosomal, 1e6, 1e3) | contig, groups = tissue, data = x,
                 subset = x$contig %in% sel.contigs & x$tissue != "muscle",
                 type = "p", pch = ".", auto.key = list(column = 3),
                 key = simpleKey(text = tissues[-3], lines = TRUE, points = FALSE, columns = 2),
                 xlab = "position (kb)", grid = TRUE,
                 ...)
    if(chromosomal)
        update(tp, par.settings = list(superpose.symbol = list(alpha = 0.5)),
               strip = FALSE, strip.left = TRUE,
               layout = c(1, length(sel.contigs)), xlab = "position (Mb)")
    else
        update(tp, par.strip.text = list(cex = 0.7))
}

horiz.depth.plot <- function(x, y = depth.hist, sel.contigs = c("22", "X", "Y"), fi = fai, histo = NULL, ...) {
    limits <- lapply(sel.contigs, function(x) c(0, subset(fi, subset = fi$contig %in% x, select = "length", drop = TRUE)))
    lengths <- subset(fi, subset = fi$contig %in% sel.contigs, select = "length", drop = TRUE)
    dp <- xyplot(depth ~ pos | contig, groups = tissue, data = x, subset = x$contig %in% sel.contigs,
                 scales = list(x = list(tick.number = 2, relation = "free", limits = limits), y = list()),
                 type = c("p"), lwd = 2, pch = ".", auto.key = list(column = 3),
                 key = simpleKey(text = tissues, lines = TRUE, points = FALSE, lwd = 2),
                 xlab = "base position", ylab = "read depth", between = list(x = 0.5, y = 0), ylim = c(0, 300),
                 layout = c(length(sel.contigs), 1),...)
    resizePanels(dp, w = lengths)
}

combined.depth.plot <- function(x, y = depth.hist, sel.contig.ix = c(22, 23, 24), fi = fai, histo = NULL, ...) {
    fi$length.Mb <- fi$length / 1e6
    limits <- lapply(sel.contig.ix, function(x) c(0, fi$length.Mb[x]))
    sel.contigs <- levels(x$contig)[sel.contig.ix]
    histo.w <- sum(fi$length.Mb[sel.contig.ix]) / 4
    dp <- xyplot(depth ~ pos / 1e6 | contig, groups = tissue, data = x, subset = x$contig %in% sel.contigs,
                 scales = list(x = list(relation = "free", limits = limits), y = list()),
                 type = c("p", "spline"), lwd = 2, pch = ".", auto.key = list(column = 3),
                 key = simpleKey(text = tissues, lines = TRUE, points = FALSE, lwd = 2),
                 xlab = "", between = list(x = 0.5, y = 0))
    histo <- xyplot(depth ~ frequency, groups = tissue, data = y, type = "l", lwd = 2, ylim = c(-10, 300))
    tp <- c(dp, histo, y.same = TRUE, layout = c(length(sel.contigs) + 1, 1))
    tp <- resizePanels(tp, w = c(fi$length.Mb[sel.contig.ix], histo.w))
    update(tp, ylim = c(0, 300), strip = strip.custom(factor.levels = c(sel.contigs, "histogram")),
           scales = list(x = list(tick.number = 2)), ...)
}
