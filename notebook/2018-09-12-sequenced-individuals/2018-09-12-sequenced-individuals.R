reshape.summ <- function(sm) {
    long <- expand.grid(row.names(sm), names(sm))
    long$done <- sapply(seq.int(dim(long)[1]),
                        function(y) sm[as.character(long[y, 1]), as.character(long[y, 2])])
    names(long)[1:2] <- c("workflow", "individuals")
    long$operation <- long$workflow
    long$operation <- sub(".*\\.fastq", "Sequencing", long$operation)
    long$operation <- sub(".*\\.bam", "Mapping", long$operation)
    long$operation <- sub(".*\\.vcf", "Calling", long$operation)
    long$operation <- factor(long$operation, ordered = TRUE, levels = c("Sequencing", "Mapping", "Calling"))
    long[long$operation == "Mapping" & long$done == 1 & ! is.na(long$done), "done"] <- 2
    long[long$operation == "Calling" & long$done == 1 & ! is.na(long$done), "done"] <- 3
    return(long)
}

my.levelplot <- function(df, my.col = c("gray", "darkgreen", "red", "blue"), indiv.as.columns = FALSE) {
    xs.tr <- rev # axis transformation function
    dfl <- lapply(df, reshape.summ)
    dfl <- do.call(rbind, lapply(names(dfl), function(y) data.frame(group = y, dfl[[y]])))
    ncols <- sapply(df, ncol)
    nrows <- sapply(df, nrow)
    nops <- length(levels(dfl$operation)) # number of operations
    rlim <- cumsum(ncols)
    llim <- rlim - ncols + 1
    limits <- lapply(seq_along(llim), function(y) c(llim[y], rlim[y]))
    at <- lapply(limits, function(l) seq(from = l[1], to = l[2], by = 1))
    labels <- lapply(at, function(y) levels(dfl$individuals)[y])
    #labels <- unlist(lapply(labels, function(y) c(list(y), rep(list(rep("", length(y))), nops - 1))), recursive = FALSE)
    lp.1 <- levelplot(done ~ workflow * individuals | group , data = dfl, col.regions = my.col, colorkey = FALSE,
              scales = list(x = list(rot = 90),
                            y = list(alternating = FALSE, relation = "free", limits = limits, labels = labels, at = at)),
              xlab = "", ylab = "", layout = c(1, length(df)), between = list(y = 0.5))
    lp.2 <- levelplot(done ~ individuals * xs.tr(workflow) | group, data = dfl, col.regions = my.col, colorkey = FALSE,
              scales = list(x = list(rot = 90, alternating = FALSE, relation = "free", limits = limits, labels = labels, at = at),
                            y = list(labels = levels(dfl$workflow), at = xs.tr(seq_len(nrows[1])))),
              xlab = "", ylab = "", layout = c(length(df), 1), between = list(x = 0.5))
    lp.1 <- resizePanels(lp.1, h = ncols)
    lp.2 <- resizePanels(lp.2, w = ncols)
    lp <- if(indiv.as.columns) lp.2 else lp.1
    lp <- update(lp, auto.key = list(text = levels(dfl$operation), rectangles = TRUE, points = FALSE),
                 par.settings = list(superpose.polygon = list(col = my.col[-1])))
    print(lp)
    return(dfl)
}
