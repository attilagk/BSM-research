# takes as input 'tsv', the output of mybcftools-isec-tsv
# output is a list of all sets whose partitinos are in 'tsv'
# While this function is fast, it depends on mybcftools-isec-tsv, which is
# limited to comparing 4 or 5 callsets (e.g. it cannot be applied to the case
# of only 2 callsets).
import.mybcftools.isec.tsv <- function(tsv, setnames = character(0)) {
    df <- read.delim(tsv, header = FALSE)
    cnames <- c("segment", "position", "ref", "alt")
    m <- as.matrix(df[ , seq_along(cnames)])
    v <- apply(m, 1, paste, collapse = ":")
    S <- df[ , seq(from = length(cnames) + 1, to = ncol(df))]
    names(S)<- if(length(setnames)) setnames else paste0("set", seq.int(ncol(S)))
    lapply(S, function(y) v[as.logical(y)])
}


# Similar to import.mybcftools.isec.tsv in its semantics and output but does
# not depend on mybcftools-isec-tsv and thus can be applied to comparisons of
# any number of callsets.
import.mybcftools.isec.tsv.2 <- function(tsv, setnames = character(0)) {
    df <- read.delim(tsv, header = FALSE, colClasses = c("integer", "integer", "factor", "factor", "character"))
    m <- as.matrix(df[ , 1:4])
    v <- apply(m, 1, paste, collapse = ":")
    S <- data.frame(t(sapply(data.frame(strsplit(df[[5]], split = ""), stringsAsFactors = FALSE), as.integer)))
    names(S)<- if(length(setnames)) setnames else paste0("set", seq.int(ncol(S)))
    lapply(S, function(y) v[as.logical(y)])
}
