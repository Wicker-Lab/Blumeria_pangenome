library(pheatmap)
library(RColorBrewer)

breaksList = seq(97, 100, by = 0.3)

orthoANI_chrom01_allisols_v2 <- read.delim("orthoANI_chrom01_allisols_v3", header =  T)
full_matrix_chr01 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom01_allisols_v2)
pheatmap_matrix_chr01 <- pheatmap(full_matrix_chr01, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom01_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr01
dev.off()

orthoANI_chrom02_allisols_v2 <- read.delim("orthoANI_chrom02_allisols_v3", header =  T)
full_matrix_chr02 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom02_allisols_v2)
pheatmap_matrix_chr02 <- pheatmap(full_matrix_chr02, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom02_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr02
dev.off()

orthoANI_chrom03_allisols_v2 <- read.delim("orthoANI_chrom03_allisols_v3", header =  T)
full_matrix_chr03 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom03_allisols_v2)
pheatmap_matrix_chr03 <- pheatmap(full_matrix_chr03, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom03_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr03
dev.off()

orthoANI_chrom04_allisols_v2 <- read.delim("orthoANI_chrom04_allisols_v3", header =  T)
full_matrix_chr04 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom04_allisols_v2)
pheatmap_matrix_chr04 <- pheatmap(full_matrix_chr04, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom04_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr04
dev.off()

orthoANI_chrom05_allisols_v2 <- read.delim("orthoANI_chrom05_allisols_v3", header =  T)
full_matrix_chr05 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom05_allisols_v2)
pheatmap_matrix_chr05 <- pheatmap(full_matrix_chr05, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom05_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr05
dev.off()

orthoANI_chrom06_allisols_v2 <- read.delim("orthoANI_chrom06_allisols_v3", header =  T)
full_matrix_chr06 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom06_allisols_v2)
pheatmap_matrix_chr06 <-pheatmap(full_matrix_chr06, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom06_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr06
dev.off()

orthoANI_chrom07_allisols_v2 <- read.delim("orthoANI_chrom07_allisols_v3", header =  T)
full_matrix_chr07 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom07_allisols_v2)
pheatmap_matrix_chr07 <- pheatmap(full_matrix_chr07, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom07_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr07
dev.off()

orthoANI_chrom08_allisols_v2 <- read.delim("orthoANI_chrom08_allisols_v3", header =  T)
full_matrix_chr08 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom08_allisols_v2)
pheatmap_matrix_chr08 <- pheatmap(full_matrix_chr08, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom08_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr08
dev.off()

orthoANI_chrom09_allisols_v2 <- read.delim("orthoANI_chrom09_allisols_v3", header =  T)
full_matrix_chr09 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom09_allisols_v2)
pheatmap_matrix_chr09 <- pheatmap(full_matrix_chr09, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom09_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr09
dev.off()

orthoANI_chrom10_allisols_v2 <- read.delim("orthoANI_chrom10_allisols_v3", header =  T)
full_matrix_chr10 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom10_allisols_v2)
pheatmap_matrix_chr10 <- pheatmap(full_matrix_chr10, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom10_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr10
dev.off()

orthoANI_chrom11_allisols_v2 <- read.delim("orthoANI_chrom11_allisols_v3", header =  T)
full_matrix_chr11 <- xtabs(orthoANI_value ~ ChromosomeA + ChromosomeB, data = orthoANI_chrom11_allisols_v2)
pheatmap_matrix_chr11 <- pheatmap(full_matrix_chr11, color = colorRampPalette(rev(brewer.pal(9, "Reds")))(10), breaks = breaksList, fontsize_row = 12, fontsize_col = 12)
pdf(file="orthoANI_chrom11_allisols_v4.pdf", width = 9, height = 9)
pheatmap_matrix_chr11
dev.off()
