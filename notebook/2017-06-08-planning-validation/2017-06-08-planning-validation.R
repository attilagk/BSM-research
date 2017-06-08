best.replicate <- function(d) {
    indiv <- unique(as.character(d$Individual.ID))
    indiv <- indiv[indiv != ""] # remove empty string
    best.within.indiv <- function(i) {
        di <- d[d$Individual.ID == i, ]
        bestrep <- which.max(di$Total.DNA..ng.)
        if(! length(bestrep)) NULL
        else di[bestrep, ]
    }
    df <- do.call(rbind, lapply(indiv, best.within.indiv))
    row.names(df) <- df$Individual.ID
    return(df)
}
