
vmc.precision <- function(svmprobs) {
    numer <- cumsum(sort(svmprobs, decreasing = TRUE))
    denom <- seq_along(numer)
    return(numer / denom)
}

