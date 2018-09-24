
get.vaf <- function(df) {
    get.allele.count <- function(allele = c("Ref", "Alt")[1], sample = c("Normal", "Tumor")[1])
        return(as.numeric(df[cbind(seq_len(nrow(df)), as.numeric(df[[allele]]) + ifelse(sample == "Normal", 6, 14))]))
    get.depth <- function(sample = c("Normal", "Tumor")[1])
        return(as.numeric(df[cbind(seq_len(nrow(df)), ifelse(sample == "Normal", 3, 11))]))
    Normal.DP <- get.depth("Normal")
    Tumor.DP <- get.depth("Tumor")
    cbind(data.frame(Ref.VAF.Normal = get.allele.count("Ref", "Normal") / Normal.DP,
                     Alt.VAF.Normal = get.allele.count("Alt", "Normal") / Normal.DP),
          data.frame(Ref.VAF.Tumor = get.allele.count("Ref", "Tumor") / Tumor.DP,
                     Alt.VAF.Tumor = get.allele.count("Alt", "Tumor") / Tumor.DP))
}

reshape.vaf <- function(v)
    reshape(v, varying = list(c("Ref.VAF.Normal", "Ref.VAF.Tumor"), c("Alt.VAF.Normal", "Alt.VAF.Tumor")),
            v.names = c("Ref.VAF", "Alt.VAF"), timevar = "Sample", times = c("Normal", "Tumor"), direction = "long")

vmc.precision <- function(svmprobs) {
    numer <- sort(svmprobs, decreasing = TRUE)
    denom <- seq_along(numer)
    return(numer / denom)
}
