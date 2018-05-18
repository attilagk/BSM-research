# takes as input 'tsv', the output of mybcftools-isec-tsv
# output is a list of all sets whose partitinos are in 'tsv'
import.mybcftools.isec.tsv <- function(tsv, setnames = character(0)) {
    df <- read.delim(tsv, header = FALSE)
    cnames <- c("segment", "position", "ref", "alt")
    m <- as.matrix(df[ , seq_along(cnames)])
    v <- apply(m, 1, paste, collapse = ":")
    S <- df[ , seq(from = length(cnames) + 1, to = ncol(df))]
    names(S)<- if(length(setnames)) setnames else paste0("set", seq.int(ncol(S)))
    lapply(S, function(y) v[as.logical(y)])
}
