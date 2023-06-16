library(circlize)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra)
library(vcfR)

sample_map <- data.frame(patient_ID = c("Sample_1"),
                         caller1 = c("example/DELLY/CW114.vcf"),
                         caller2 = c("example/Manta/CW114.vcf"))


VCF_List <- readVCFs(sample_map)

## Make the bed for translocation

set.seed(123)
bed1 = generateRandomBed(nr = 1000)

bed1_chr1 <- bed1 %>% filter(chr == "chr1")
bed1_chr1 <- bed1_chr1[ ,c(1,2,2)]
bed1_chr1 = bed1_chr1[sample(nrow(bed1_chr1), 20), ]



bed2 = generateRandomBed(nr = 1000)

bed2_chr2 <- bed2 %>% filter(chr == "chr2")
bed2_chr2 <- bed2_chr2[ ,c(1,2,2)]
bed2_chr2 = bed2_chr2[sample(nrow(bed2_chr2), 20), ]

translocation_example_table <- cbind(bed1_chr1[, c(1,2)], bed2_chr2[, c(1,2)])

fwrite(translocation_example_table, "./inst/translocation_example.csv")


pdf("./inst/extdata/Figures/translocation_example.pdf", width = 5, height = 5)
circos.initializeWithIdeogram(chromosome.index = c("chr1", "chr2"))
circos.genomicLink(bed1_chr1[, c(1,2)], bed2_chr2[, c(1,2)], lty = 1, lwd = 2)
circos.clear()
dev.off()

colnames(translocation_example_table) <- c("chr_1", "pos_1", "chr_2", "pos_2")

pdf("./inst/extdata/Figures/translocation_scatter_plot.pdf", width = 4, height = 4)

ggplot(translocation_example_table, aes(x = pos_1, y = pos_2)) +
    geom_point(size = 2) +
    theme_classic2() + theme(text = element_text(size = 10, face = "bold"),
                             axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust=0.5, colour = "black"),
                             axis.text.y = element_text(size = 10, colour = "black"),
                             strip.background = element_blank()) +
    xlab("Chromosome 1") +
    ylab("Chromosome 2") + xlim(c(0, 2.5e+08)) + ylim(c(0, 2.5e+08))

dev.off()

distance <- get_dist(translocation_example_table[, c(2,4)])
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

df <- translocation_example_table[, c(2,4)]

k2 <- kmeans(df, centers = 2, nstart = 25)
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df, ggtheme = theme_classic2()) + ggtitle("k = 2") + xlab("Chromosome 1") + ylab("Chromosome 2")+
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))+ xlim(c(-2, 2)) + ylim(c(-2, 2))

p2 <- fviz_cluster(k3, geom = "point",  data = df, ggtheme = theme_classic2()) + ggtitle("k = 3")+ xlab("Chromosome 1") + ylab("Chromosome 2")+
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))+ xlim(c(-2, 2)) + ylim(c(-2, 2))

p3 <- fviz_cluster(k4, geom = "point",  data = df, ggtheme = theme_classic2()) + ggtitle("k = 4")+ xlab("Chromosome 1") + ylab("Chromosome 2")+
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))+ xlim(c(-2, 2)) + ylim(c(-2, 2))

p4 <- fviz_cluster(k5, geom = "point",  data = df, ggtheme = theme_classic2()) + ggtitle("k = 5")+ xlab("Chromosome 1") + ylab("Chromosome 2")+
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))+ xlim(c(-2, 2)) + ylim(c(-2, 2))

pdf("./inst/extdata/Figures/translocation_scatter_plot.pdf", width = 5, height = 4)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
