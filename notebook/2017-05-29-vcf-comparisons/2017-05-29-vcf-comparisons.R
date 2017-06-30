my.venn2 <- function(dir = "results/2_cmp-reftissues/snvs/", cls = clsets, ...) {
    ss <- with(cls, subset(cls, subset = directory == dir, select = set.size))$set.size
    grid.newpage()
    grid.draw(draw.pairwise.venn(ss[1] + ss[3], ss[2] + ss[3], ss[3],
                                 cex = 1.3, cat.cex = 1.3, cat.pos = 165 * c(-1, 1),
                                 ...))
}


my.venn3 <- function(dir = "1_isec-callers/NeuN_mn-NeuN_pl/snvs/", cls = clsets1, ...) {
    dir <- paste0("results/mutect2-", c("unfilt", "PASS"), "/", dir)
    ss <- with(cls, subset(cls, subset = directory %in% dir, select = set.size))$set.size
    area1 <- ss[1] + ss[3] # unfiltered mutect2
    area2 <- ss[2] + ss[3] # strelka
    area3 <- ss[4] + ss[6] # PASS + mutect2
    n12 <- ss[3] # both unfiltered mutect2 and strelka
    n23 <- ss[6] # both PASS + mutect2 and strelka
    n13 <- area3  # both unfiltered mutect2 and PASS + mutect2 = PASS + mutect2
    n123 <- ss[6] # all three = both PASS + mutect2 and strelka
    grid.newpage()
    grid.draw(draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                               category = c("mutect2", "strelka", "mutect2.PASS"),
                               cex = 1.3, cat.cex = 1.3, col = my.col <- c("cyan", "magenta", "darkcyan"),
                               fill = my.col,
                               ...))
}
