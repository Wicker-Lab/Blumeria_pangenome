###Figure 1A
#Worldwide passport data for B. graminis isolates used
#Libraries
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggplot2)
library(ggrepel)

Input_Fig1A <- read.delim("Input_Fig1A")

world_map_data <- map_data("world")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
bg_worldwide <-ggplot(data = world)+
  geom_sf(fill = "white") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World Map")+
  geom_point(data = Input_Fig1A, aes(x=Longitude, y=Latitude, colour=Country, size=Size))+
  scale_size(range = c(2, 7)) + scale_color_manual(values=c("Argentina"="#88CCEE", "Australia"="#CC6677", "China"="#DDCC77", "Egypt"="#117733", "Switzerland"="#332288", "Poland"="#332288", "Iran"="#AA4499", "Israel/Palestine"="#44AA99", "Japan"="#999933",  "Russia"="#661100", "Turkey"="#6699CC", "USA"="#888888", "Kazakhstan"="#000000"))+
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=20), legend.title = element_text(size = 20)) + theme(plot.title = element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5))) + guides(shape = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(data=Input_Fig1A, aes(x=Longitude, y=Latitude, label=Genomes), size=5, box.padding = 0.5, max.overlaps = Inf, min.segment.length = unit(0, 'lines')) 
bg_worldwide

pdf(file="Output_Fig1A.pdf", width = 12, height = 9)
bg_worldwide
dev.off()











###Figure 1B
#PCA etc Blumeria graminis f.sp. tritici
#load tidyverse package and other ones
library(tidyverse)
library(ggplot2)
library(ggrepel)
#read in data
pca <- read_table("Input_Fig1B.eigenvec", col_names = FALSE)
eigenval <- scan("Input_Fig1B.eigenval")

#sort out the pca data
#remove nuisance column
pca <- pca[,-1]
#set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca <- pca[-1,]

#location
loc <- rep(NA, length(pca$ind))
loc[grep("ISR", pca$ind)] <- "Israel/Palestine"
loc[grep("EGY", pca$ind)] <- "Egypt"
loc[grep("CHN", pca$ind)] <- "China"
loc[grep("ARG", pca$ind)] <- "Argentina"
loc[grep("AUS", pca$ind)] <- "Australia"
loc[grep("CHE", pca$ind)] <- "Europe"
loc[grep("FRA", pca$ind)] <- "Europe"
loc[grep("GBR", pca$ind)] <- "Europe"
loc[grep("IRN", pca$ind)] <- "Iran"
loc[grep("RUS", pca$ind)] <- "Russia"
loc[grep("USA", pca$ind)] <- "USA"
loc[grep("TUR", pca$ind)] <- "Turkey"
loc[grep("KAZ", pca$ind)] <- "Kazakhstan"
loc[grep("POL", pca$ind)] <- "Europe"
loc[grep("JPN", pca$ind)] <- "Japan"

# remake data.frame
pca <- as_tibble(data.frame(pca, loc))

library(rcartocolor) #Safe palette, plus black colour
safe_colorblind_palette13 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#000000", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette13)

labels <- pca |> filter(ind %in% c("ARG_4_2", "CHE_96224", "CHN_17_40", "CHN_52_27", "IRN_GOR_2", "ISR_7", "JPN_CHIKA", "USA_2"))
labels <- as.tibble(data.frame(labels))
labels$PC1 <- as.numeric(as.character(labels$PC1))
labels$PC2 <- as.numeric(as.character(labels$PC2))

pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# plot pca
pca$PC1 <- as.numeric(as.character(pca$PC1))
pca$PC2 <- as.numeric(as.character(pca$PC2))
b <- ggplot(pca, aes(PC1, PC2, col = loc)) + geom_point(size = 2.1) 
b <- b + scale_colour_manual(values = safe_colorblind_palette13)
b <- b + theme_minimal()
b <- b + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(color='Country') + #theme(legend.text = element_text(face = "bold"))+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=25), legend.title = element_text(size = 25))+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_label_repel(data=labels, aes(x=PC1, y=PC2, label=ind, size=5), box.padding = 0.6, max.overlaps = Inf, min.segment.length = unit(0, 'lines'), 
                   nudge_y = .016, show.legend=F)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) 

#Blumeria graminis f.sp. tritici cohort with 388 isolate worldwide
png("Output_Fig1B.png")
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
dev.off()
pdf("Output_Fig1B.pdf", width = 12, height = 6)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
dev.off()










###Figure 1D
#For Nucleotide Diversity
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(multcomp)
library(patchwork)
library(palmerpenguins)
coh400_8perpop_allpops <- read.delim("Input_1D.windowed.pi", header=T)
ggplot(coh400_8perpop_allpops) +
  aes(x = Population, y = nuc_div_wins, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_8perpop_allpops$nuc_div_wins, coh400_8perpop_allpops$Population, p.adjust.method = "holm")

x <- which(names(coh400_8perpop_allpops) == "Population") # name of grouping variable
y <- which(names(coh400_8perpop_allpops) == "nuc_div_wins") # names of variables to test

anova_model <- aov(nuc_div_wins ~ Population, data = coh400_8perpop_allpops)
tukey_result <- TukeyHSD(anova_model)
tukey_result
tukey_df <- as.data.frame(tukey_result$Population)
group_names <- rownames(tukey_df)
group1 <- sapply(strsplit(group_names, "-"), `[`, 1)
group2 <- sapply(strsplit(group_names, "-"), `[`, 2)

# Add group columns
tukey_df$group1 <- sub(" -.*", "", rownames(tukey_df))
tukey_df$group2 <- sub(".*- ", "", rownames(tukey_df))

# Add y.position (you can adjust as needed)
tukey_df$y.position <- seq(max(coh400_8perpop_allpops$nuc_div_wins) * 1.05, 
                           length.out = nrow(tukey_df), by = 0.02)

# Add significance stars using p.adjusted values
tukey_df$p.adj.signif <- cut(tukey_df$`p adj`,
                             breaks = c(0, 0.001, 0.01, 0.05, 1),
                             labels = c("***", "**", "*", "ns"))
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("ARG", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISRdic"), c("ISRdic", "ISR"), c("ISR", "JPN"), c("JPN", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_8perpop_allpops,
                   x = colnames(coh400_8perpop_allpops[j]), y = colnames(coh400_8perpop_allpops[i]),
                   color = colnames(coh400_8perpop_allpops[j]), ylab = "Nucleotide Diversity in 10k windows", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_8perpop_allpops[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400 <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 0.00099) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 0.00088) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "ISR", label.y = 0.00083) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "ISRdic"="#000000"))
ggpar(p400, orientation = "vertical", ylim = c(min(coh400_8perpop_allpops$nuc_div_wins), 0.00099))

pdf("Output_Fig1D_windowed_pi.pdf", width=12, height=4)
ggpar(p400, orientation = "vertical", ylim = c(min(coh400_8perpop_allpops$nuc_div_wins), 0.00099))
dev.off()

png("Output_Fig1D_windowed_pi.png", width=9000, height=3000, res=600)
ggpar(p400, orientation = "vertical", ylim = c(min(coh400_8perpop_allpops$nuc_div_wins), 0.00099))
dev.off()









###S1 Fig
#Percentage of PCA etc Blumeria graminis f.sp. tritici variance explained from PCA in Figure 1B
eigenval <- scan("Input_Fig1B.eigenval")

# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# make plot
library(ggplot2)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a2 <- a + ylab("Percentage variance explained") + theme_light() + theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                                                                        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
                                                                        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
                                                                        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))
png("Output_S1Fig.png")
a2
dev.off()
pdf("Output_S1Fig.pdf", width = 12, height = 8)
a2
dev.off()










###S4 Fig. - Singletons (both SNPs and INDELs, or only one) for 8 random isolates or all isolates
# library
library(ggplot2)
library(viridis)
library(reshape2)

Input_S4A_Fig <- read.delim("Input_S4A_Fig")

###S4A Fig. - Singletons (both SNPs and INDELs) - 8 random isolates
#Singletons (both SNPs and INDELs) for a subset of 8 or less (when not available) random isolate
Input_S4A <- ggplot(Input_S4A_Fig, aes(fill=Population, y=Counts, x=interaction(ID))) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Singletons - Both SNPs and INDELs") +
  #facet_wrap(~Population) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", 
                             "ISR"="#44AA99", "JPN"="#999933", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "ISRdic"="#000000"))+
  xlab("Isolates")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(color='Species') + 
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 60, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size = 18))
Input_S4A

png("Output_S4A_Fig_SNPs_and_INDELs.png")
Input_S4A
dev.off()
pdf("Output_S4A_Fig_SNPs_and_INDELs.pdf", width = 12, height = 6)
Input_S4A
dev.off()

###S4B Fig. - Singletons (only SNPs) for all isolates
#Singletons (only SNPs) for all isolates
Input_S4BC_Fig <- read.delim("Input_S4BC_Fig")
Input_S4B <- ggplot(Input_S4BC_Fig, aes(fill=Population, y=Counts_SNPs, x=interaction(Isolate))) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Singletons - SNPs") +
  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", 
                             "ISR"="#44AA99", "JPN"="#999933", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "ISRdic"="#000000",
                             "KAZ"="#661100"))+
  xlab("Isolates")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(color='Species') + 
  theme(axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size = 25))
Input_S4B

png("Output_S4B_Fig_SNPs.png")
Input_S4B
dev.off()
pdf("Output_S4B_Fig_SNPs.pdf", width = 12, height = 4)
Input_S4B
dev.off()

###S4C Fig. - Singletons (only INDELs) for all isolates
#Singletons (only INDELs) for all isolates
Input_S4BC_Fig <- read.delim("Input_S4BC_Fig")
Input_S4C <- ggplot(Input_S4BC_Fig, aes(fill=Population, y=Counts_INDELs, x=interaction(Isolate))) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Singletons - INDELs") +
  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", 
                             "ISR"="#44AA99", "JPN"="#999933", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "ISRdic"="#000000",
                             "KAZ"="#661100"))+
  xlab("Isolates")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(color='Species') + 
  theme(axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size = 25))
Input_S4C

png("Output_S4C_Fig_INDELs.png")
Input_S4C
dev.off()
pdf("Output_S4C_Fig_INDELs.pdf", width = 12, height = 4)
Input_S4C
dev.off()




















###S5 Fig - Mantel Test
###For World - 366 Blumeria graminis isolates
library(vcfR)
vcf_coh366 <- read.vcfR("Input_S5A_Fig.vcf")   #read in all data 
aa.genlight_coh366 <- vcfR2genlight(vcf_coh366, n.cores=2)
library(adegenet)
locNames(aa.genlight_coh366) <- paste(vcf_coh366@fix[,1],vcf_coh366@fix[,2],sep="_")   # add real SNP.names 
pop(aa.genlight_coh366)<-substr(indNames(aa.genlight_coh366),1,3)               # add pop
library(StAMPP)
aa.D.ind_coh366 <- stamppNeisD(aa.genlight_coh366, pop = FALSE)  # Nei's 1972 distance between indivs 
aa.D.pop_coh366 <- stamppNeisD(aa.genlight_coh366, pop = TRUE)   # Nei's 1972 distance between pops 
colnames(aa.D.ind_coh366) <- rownames(aa.D.ind_coh366)
colnames(aa.D.ind_coh366) <- rownames(aa.D.ind_coh366)    
aa.D.ind.dist_coh366<-as.dist(aa.D.ind_coh366, diag=T, upper = T) 
attr(aa.D.ind.dist_coh366, "Labels")<-rownames(aa.D.ind_coh366)
coords_coh366 <- read.csv ("Input_S5A_Fig_coords", sep ="\t", header = T)     # tab-separated file for all pops 
xy.coords.only_coh366 <- subset(coords_coh366, select=c("lon","lat")) 
Dgeo_coh366 <- dist(xy.coords.only_coh366, diag = T, upper = T, method = "euclidean")
attr(Dgeo_coh366, "Labels")<-rownames(aa.D.ind_coh366)
IBD_coh366 <- mantel.randtest(Dgeo_coh366,aa.D.ind.dist_coh366,nrepet = 100000)
IBD_coh366
Dgeo_coh366 <- as.dist(Dgeo_coh366)
aa.D.ind.dist_coh366 <- as.dist(aa.D.ind.dist_coh366)
Dgeo_coh366 = as.matrix(Dgeo_coh366)
aa.D.ind.dist_coh366 = as.matrix(aa.D.ind.dist_coh366)
library(MASS)
dens_coh366 <- kde2d(Dgeo_coh366,aa.D.ind.dist_coh366, h=rep(1.5, 2), n=1000, lims=c(-1, 300, -1, 0.9)) 
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
pdf("Output_S5A_Fig.pdf", width = 10, height = 5)
plot(Dgeo_coh366, aa.D.ind.dist_coh366, pch=20,cex=.9, xlab="Geogrpahic distance (Euclidean)", ylab="Genetic Distance (Nei's 1972)", cex.lab=1.2, cex.axis=1.5)
image(dens_coh366, col=transp(myPal(300),.7), add=TRUE) 
abline(lm(as.numeric(as.character(aa.D.ind.dist_coh366))~as.numeric(as.character(Dgeo_coh366))))
title("Correlation of Genetic and Geographic distances (World - 366 Blumeria graminis isolates)")
dev.off()
png("Output_S5A_Fig.png", width = 5000, height = 3000, res = 600)
plot(Dgeo_coh366, aa.D.ind.dist_coh366, pch=20,cex=.9, xlab="Geogrpahic distance (Euclidean)", ylab="Genetic Distance (Nei's 1972)", cex.lab=1.2, cex.axis=1.5)
image(dens_coh366, col=transp(myPal(300),.7), add=TRUE) 
abline(lm(as.numeric(as.character(aa.D.ind.dist_coh366))~as.numeric(as.character(Dgeo_coh366))))
title("World - 366 Blumeria graminis isolates")
dev.off()

###For Middle Asia - 166 Blumeria graminis isolates
library(vcfR)
vcf_coh166 <- read.vcfR("Input_S5B_Fig.vcf")   #read in all data 
aa.genlight_coh166 <- vcfR2genlight(vcf_coh166, n.cores=2) 
library(adegenet)
locNames(aa.genlight_coh166) <- paste(vcf_coh166@fix[,1],vcf_coh166@fix[,2],sep="_")   # add real SNP.names 
pop(aa.genlight_coh166)<-substr(indNames(aa.genlight_coh166),1,3)               # add pop
library(StAMPP)
aa.D.ind_coh166 <- stamppNeisD(aa.genlight_coh166, pop = FALSE)  # Nei's 1972 distance between indivs 
aa.D.pop_coh166 <- stamppNeisD(aa.genlight_coh166, pop = TRUE)   # Nei's 1972 distance between pops 
colnames(aa.D.ind_coh166) <- rownames(aa.D.ind_coh166)    
aa.D.ind.dist_coh166<-as.dist(aa.D.ind_coh166, diag=T, upper = T) 
attr(aa.D.ind.dist_coh166, "Labels")<-rownames(aa.D.ind_coh166)
coords_coh166 <- read.csv ("Input_S5B_Fig_coords", sep ="\t", header = T)     # tab-separated file for all pops 
xy.coords.only_coh166 <- subset(coords_coh166, select=c("lon","lat")) 
Dgeo_coh166 <- dist(xy.coords.only_coh166, diag = T, upper = T, method = "euclidean")
attr(Dgeo_coh166, "Labels")<-rownames(aa.D.ind_coh166)
IBD_coh166 <- mantel.randtest(Dgeo_coh166,aa.D.ind.dist_coh166,nrepet = 100000)
IBD_coh166
Dgeo_coh166 <- as.dist(Dgeo_coh166)
aa.D.ind.dist_coh166 <- as.dist(aa.D.ind.dist_coh166)
Dgeo_coh166 = as.matrix(Dgeo_coh166)
aa.D.ind.dist_coh166 = as.matrix(aa.D.ind.dist_coh166)
library(MASS)
dens_coh166 <- kde2d(Dgeo_coh166,aa.D.ind.dist_coh166, h=rep(1.5, 2), n=1000, lims=c(-1, 300, -1, 0.9)) 
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
pdf("Output_S5B_Fig.pdf", width = 10, height = 5)
plot(Dgeo_coh166, aa.D.ind.dist_coh166, pch=20,cex=.9, xlab="Geogrpahic distance (Euclidean)", ylab="Genetic Distance (Nei's 1972)", cex.lab=1.2, cex.axis=1.5)
image(dens_coh166, col=transp(myPal(300),.7), add=TRUE) 
abline(lm(as.numeric(as.character(aa.D.ind.dist_coh166))~as.numeric(as.character(Dgeo_coh166))))
title("Correlation of Genetic and Geographic distances (Middle Asia - 166 Blumeria graminis isolates)")
dev.off()
png("Output_S5B_Fig.png", width = 5000, height = 3000, res = 600)
plot(Dgeo_coh166, aa.D.ind.dist_coh166, pch=20,cex=.9, xlab="Geogrpahic distance (Euclidean)", ylab="Genetic Distance (Nei's 1972)", cex.lab=1.2, cex.axis=1.5)
image(dens_coh166, col=transp(myPal(300),.7), add=TRUE) 
abline(lm(as.numeric(as.character(aa.D.ind.dist_coh166))~as.numeric(as.character(Dgeo_coh166))))
title("Middle Asia - 166 Blumeria graminis isolates")
dev.off()




