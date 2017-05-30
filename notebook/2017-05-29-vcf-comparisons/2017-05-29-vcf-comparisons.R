my.grid.draw <- function(dir = "results/2_cmp-reftissues/snvs/", cls = clsets, ...) {
    ss <- with(cls, subset(cls, subset = directory == dir, select = set.size))$set.size
    grid.newpage()
    grid.draw(draw.pairwise.venn(ss[1] + ss[3], ss[2] + ss[3], ss[3],
                                 cex = 1.5, cat.cex = 1.5,
                                 ...))
}

#extractor <- function(pat = "indels|snvs")
#    factor(sub(paste0("^.*(", pat, ").*$"), "\\1", clsets$directory))
#clsets$var.type <- extractor("indels|snvs")
#clsets$ref.tissue <- extractor("NeuN_mn|muscle")
#clsets$comparison <- extractor("callers|reftissues")
