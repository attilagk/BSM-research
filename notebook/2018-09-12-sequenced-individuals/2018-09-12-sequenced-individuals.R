reshape.summ <- function(sm) {
    long <- expand.grid(row.names(sm), names(sm))
    long$done <- sapply(seq.int(dim(long)[1]),
                        function(y) sm[as.character(long[y, 1]), as.character(long[y, 2])])
    names(long)[1:2] <- c("workflow", "individuals")
    return(long)
}
