my.venn2 <- function(dir = "results/2_cmp-reftissues/snvs/", cls = clsets, ...) {
    ss <- with(cls, subset(cls, subset = directory == dir, select = set.size))$set.size
    grid.newpage()
    grid.draw(draw.pairwise.venn(ss[1] + ss[3], ss[2] + ss[3], ss[3],
                                 cex = 1.3, cat.cex = 1.3, cat.pos = 165 * c(-1, 1),
                                 ...))
}


my.venn3a <- function(dir = "1_isec-callers/NeuN_mn-NeuN_pl/snvs/", cls = clsets1, ...) {
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

my.venn3b <- function(dir = "results/mutect2-unfilt/2_cmp-reftissues/snvs/", ...) {
    df <- read.delim(paste0(dir, "callset-sizes.tsv"))
    l <-
        list(area1 = c(1, 4, 5, 7), area2 = c(2, 4, 6, 7), area3 = c(3, 5, 6, 7),
             n12 = 4, n23 = 6, n13 = 5, n123 = 7)
    arg <- lapply(l, function(i) sum(df$nrec[i]))
    arg$category <- as.character(unlist(df[1, 3:5]))
    arg$cex <- 1.3
    arg$cat.cex <- 1.3
    arg$col <- my.col <- c("skyblue", "green", "brown")
    arg$fill <- my.col
    grid.newpage()
    grid.draw(do.call(draw.triple.venn, c(arg, ...)))
    invisible(l)
}
