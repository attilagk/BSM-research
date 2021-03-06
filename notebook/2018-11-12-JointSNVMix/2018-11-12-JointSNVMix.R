# import records from a jsm.tsv file and optionally reshape (default)
get.records <- function(tsv, do.reshape = TRUE) {
    records <- read.delim(tsv)
    # aggregate probabilities
    # across genotypes for the case/tumor sample
    records$p_AA_UV <- with(records, p_AA_AA + p_AA_AB + p_AA_BB)
    records$p_AB_UV <- with(records, p_AB_AA + p_AB_AB + p_AB_BB)
    records$p_BB_UV <- with(records, p_BB_AA + p_BB_AB + p_BB_BB)
    # across genotypes for the control/normal sample
    records$p_UV_AA <- with(records, p_AA_AA + p_AB_AA + p_BB_AA)
    records$p_UV_AB <- with(records, p_AA_AB + p_AB_AB + p_BB_AB)
    records$p_UV_BB <- with(records, p_AA_BB + p_AB_BB + p_BB_BB)
    # all somatic mutations from AA
    records$p_AA_UB <- with(records, p_AA_AB + p_AA_BB)
    # all somatic mutations from AB
    records$p_AB_UU <- with(records, p_AA_AB + p_AA_BB)
    # all somatic mutations from BB
    records$p_BB_AU <- with(records, p_BB_AB + p_BB_AA)
    if (! do.reshape)
        return(records)
    records.l <-
        reshape(records, varying = v <- grep("p_", names(records), value = TRUE),
                v.names = "p", timevar = "GT", times = sub("p_", "", v), direction = "long")
    return(records.l)
}
