library(VennDiagram)
pos <- list(Tnseq.Mutect2 = read.delim("figures/benchmark-Tnseq-4.csv")$POS,
            strelka2Somatic = read.delim("figures/benchmark-strelka2Somatic-4.csv")$POS)
v <- venn.diagram(pos, NULL, col = "gray", fill = c("blue", "red"), cex = 2, cat.cex = 2, cat.pos = c(-30, 30))
png("figures/Tnseq-4-strelka2Somatic-4-venn.png")
grid.draw(v)
dev.off()
