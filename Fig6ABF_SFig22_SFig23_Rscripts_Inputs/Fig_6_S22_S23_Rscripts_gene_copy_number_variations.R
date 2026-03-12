###Figure 6A - All isolates - Duplications
#All genes duplications with all isolates per population
library(ggplot2)
library(ggpubr)

normCovAllfinalRegionsDupSum_v2 <- read.delim("Input_6A", header = T, sep = " ")
coh400_allpops_dups <- normCovAllfinalRegionsDupSum_v2
ggplot(coh400_allpops_dups) +
  aes(x = Population, y = Dup_per_Isol, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_allpops_dups$Dup_per_Isol, coh400_allpops_dups$Population, p.adjust.method = "holm")
x <- which(names(coh400_allpops_dups) == "Population") # name of grouping variable
y <- which(names(coh400_allpops_dups) == "Dup_per_Isol") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Bgs", "Bgdic"), c("Bgdic", "ARG"), c("ARG", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "JPN"),c("JPN", "KAZ"), c("KAZ", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests
coh400_allpops_dups$Population <- factor(coh400_allpops_dups$Population , levels=c("Bgs", "Bgdic", "ARG", "CHN", "EGY", "EUR", "IRN", "ISR", "JPN", "KAZ", "RUS", "TUR", "USA"))

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_allpops_dups,
                   x = colnames(coh400_allpops_dups[j]), y = colnames(coh400_allpops_dups[i]),
                   color = colnames(coh400_allpops_dups[j]), ylab = "Gene Duplications", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_allpops_dups[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400 <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 300) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 275) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TUR", label.y = 270) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000", "KAZ"="#661100"))
ggpar(p400, orientation = "vertical")

pdf("Output_6A.pdf", width=12, height=4)
ggpar(p400, orientation = "vertical")
dev.off()
png("Output_6A.png", width=9000, height=3000, res=600)
ggpar(p400, orientation = "vertical")
dev.off()








###Figure 6B - PCA of effector gene presence/absence summary per isolate (including 5 B.g. secalis isolates)
#PCA of effector gene presence/absence summary per isolate (including 5 B.g. secalis isolates)

df <- read.delim("Input_Fig6B", header = T, row.names = 1)
cnv_matrix <- df

# Remove columns (genes) with any missing values
cnv_clean <- cnv_matrix[, colSums(is.na(cnv_matrix)) == 0]
# Remove columns (genes) with zero variance (all samples same CNV)
cnv_clean <- cnv_clean[, apply(cnv_clean, 2, var) != 0]

#Run PCA
pca_result <- prcomp(cnv_clean, center = TRUE, scale = TRUE)

library(factoextra)
eigenvalues <- get_eigenvalue(pca_result)
print(eigenvalues)

library(ggplot2)
pca <- as.data.frame(pca_result$x)
ggplot(pca, aes(x = PC1, y = PC2)) + geom_point() + theme_minimal()

library(tidyverse)
pca$ind <- rownames(pca)
# location #3 12colours
loc3 <- rep(NA, length(pca$ind))
loc3[grep("ISR", pca$ind)] <- "ISR"
loc3[grep("ISR_58", pca$ind)] <- "Bgdic"
loc3[grep("ISR_63", pca$ind)] <- "Bgdic"
loc3[grep("ISR_66", pca$ind)] <- "Bgdic"
loc3[grep("ISR_202_m", pca$ind)] <- "Bgdic"
loc3[grep("ISR_203", pca$ind)] <- "Bgdic"
loc3[grep("ISR_206_m", pca$ind)] <- "Bgdic"
loc3[grep("ISR_207", pca$ind)] <- "Bgdic"
loc3[grep("ISR_209", pca$ind)] <- "Bgdic"
loc3[grep("ISR_210_m", pca$ind)] <- "Bgdic"
loc3[grep("ISR_211_m", pca$ind)] <- "Bgdic"
loc3[grep("ISR_212", pca$ind)] <- "Bgdic"
loc3[grep("ISR_220", pca$ind)] <- "Bgdic"
loc3[grep("EGY", pca$ind)] <- "EGY"
loc3[grep("CHN", pca$ind)] <- "CHN"
loc3[grep("ARG", pca$ind)] <- "ARG"
loc3[grep("AUS", pca$ind)] <- "USA"
loc3[grep("CHE", pca$ind)] <- "EUR"
loc3[grep("FRA", pca$ind)] <- "EUR"
loc3[grep("GBR", pca$ind)] <- "EUR"
loc3[grep("IRN", pca$ind)] <- "IRN"
loc3[grep("RUS", pca$ind)] <- "RUS/KAZ"
loc3[grep("USA", pca$ind)] <- "USA"
loc3[grep("TUR", pca$ind)] <- "TUR"
loc3[grep("KAZ", pca$ind)] <- "RUS/KAZ"
loc3[grep("POL", pca$ind)] <- "EUR"
loc3[grep("JPN", pca$ind)] <- "JPN"
loc3[grep("Bgs_S_1201", pca$ind)] <- "Bgs"
loc3[grep("Bgs_S_1203", pca$ind)] <- "Bgs"
loc3[grep("Bgs_S_1391", pca$ind)] <- "Bgs"
loc3[grep("Bgs_S_1400", pca$ind)] <- "Bgs"
loc3[grep("Bgs_S_1459", pca$ind)] <- "Bgs"

# remake data.frame
pca <- as_tibble(data.frame(pca, loc3))

# PCA Plot
library(rcartocolor) #Safe palette, plus black colour
b <- ggplot(pca, aes(PC1, PC2, col = loc3)) + geom_point(size = 1.1) 
b <- b + theme_minimal()
b <- b + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(color='Country') +
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS/KAZ"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000"))+
  theme(legend.text=element_text(size=25), legend.title = element_text(size = 25))+
  guides(color = guide_legend(override.aes = list(size = 4))) 
b + xlab(paste0("Dim1 (", signif(eigenvalues$variance.percent[1], 3), "%)")) + ylab(paste0("Dim2 (", signif(eigenvalues$variance.percent[2], 3), "%)")) 

png("Output_Fig6B.png")
b + xlab(paste0("Dim1 (", signif(eigenvalues$variance.percent[1], 3), "%)")) + ylab(paste0("Dim2 (", signif(eigenvalues$variance.percent[2], 3), "%)")) 
dev.off()
pdf("Output_Fig6B.pdf", width = 12, height = 6)
b + xlab(paste0("Dim1 (", signif(eigenvalues$variance.percent[1], 3), "%)")) + ylab(paste0("Dim2 (", signif(eigenvalues$variance.percent[2], 3), "%)")) 
dev.off()









###Figure 6C - Expansion or shrinkage of effector families in the pangenome
#Bubble plot of the effector families expansion or shrinkage in the pangenome
library(ggfortify)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(reshape2)

CHN1740_inFam<-read.table("Untransformed_RUN_CHN1740.txt")
names(CHN1740_inFam)<-c("Orthogroup","Genename")
ARG42_inFam<-read.table("Untransformed_RUN_ARG42.txt")
names(ARG42_inFam)<-c("Orthogroup","Genename")
Bgs1459_inFam<-read.table("Untransformed_RUN_Bgs1459.txt")
names(Bgs1459_inFam)<-c("Orthogroup","Genename")
CHN5227_inFam<-read.table("Untransformed_RUN_CHN5227.txt")
names(CHN5227_inFam)<-c("Orthogroup","Genename")
IRNGOR2_inFam<-read.table("Untransformed_RUN_IRNGOR2.txt")
names(IRNGOR2_inFam)<-c("Orthogroup","Genename")
ISR211_inFam<-read.table("Untransformed_RUN_ISR211.txt")
names(ISR211_inFam)<-c("Orthogroup","Genename")
ISR7_inFam<-read.table("Untransformed_RUN_ISR7.txt")
names(ISR7_inFam)<-c("Orthogroup","Genename")
JPNCHIKA_inFam<-read.table("Untransformed_RUN_JPNCHIKARA.txt")
names(JPNCHIKA_inFam)<-c("Orthogroup","Genename")
THUN12_inFam<-read.table("Untransformed_RUN_THUN12.txt")
names(THUN12_inFam)<-c("Orthogroup","Genename")
USA2_inFam<-read.table("Untransformed_RUN_USA2.txt")
names(USA2_inFam)<-c("Orthogroup","Genename")


Eff_Fam_844<-read.table("Input_6C_Sorted_List_of_844_effectors_based_on_MuellerPraz_et_al_2019_withFamily")
head(Eff_Fam_844)
names(Eff_Fam_844)<-c("Genename","Fam")

FamIdentifier<-merge(CHN1740_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(CHN1740_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
head(Gene_per_Fam)
Gene_per_Fam_CHN1740=Gene_per_Fam[grep("^CHN1740",Gene_per_Fam$Genename),]

FamIdentifier<-merge(ARG42_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(ARG42_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_ARG42=Gene_per_Fam[grep("^BgtARG",Gene_per_Fam$Genename),]

FamIdentifier<-merge(CHN5227_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(CHN5227_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_CHN5227=Gene_per_Fam[grep("^BgtCHN52",Gene_per_Fam$Genename),]

FamIdentifier<-merge(ISR7_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(ISR7_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_ISR7=Gene_per_Fam[grep("^BgISR7",Gene_per_Fam$Genename),]

FamIdentifier<-merge(IRNGOR2_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(IRNGOR2_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_IRNGOR2=Gene_per_Fam[grep("^BgtIRNGOR2",Gene_per_Fam$Genename),]

FamIdentifier<-merge(JPNCHIKA_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(JPNCHIKA_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_JPNCHIKA=Gene_per_Fam[grep("^BgtJPNCHIKA",Gene_per_Fam$Genename),]

FamIdentifier<-merge(Bgs1459_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(Bgs1459_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_Bgs1459=Gene_per_Fam[grep("^Bgs",Gene_per_Fam$Genename),]

FamIdentifier<-merge(ISR211_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(ISR211_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_ISR211=Gene_per_Fam[grep("^BgtISR211",Gene_per_Fam$Genename),]

FamIdentifier<-merge(THUN12_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(THUN12_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_THUN12=Gene_per_Fam[grep("^BgTH12",Gene_per_Fam$Genename),]

FamIdentifier<-merge(USA2_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(USA2_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_USA2=Gene_per_Fam[grep("^BgUSA2",Gene_per_Fam$Genename),]

FamIdentifier<-merge(USA2_inFam,Eff_Fam_844,by=("Genename")) %>% select(Orthogroup,Fam)
Gene_per_Fam<-merge(USA2_inFam,FamIdentifier,by=c("Orthogroup")) %>% unique() 
Gene_per_Fam_96224=Gene_per_Fam[grep("^Bgt",Gene_per_Fam$Genename),]


fam96224=read.table("Input_6C_Genelist_Eff_Fam_Group_2018_39.txt",h=T) %>% mutate(Isolate="CHE_96224")# %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famCHN1740= Gene_per_Fam_CHN1740 %>% mutate(Isolate="CHN_17_40") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famCHN5227= Gene_per_Fam_CHN5227 %>% mutate(Isolate="CHN_52_27") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famARG42= Gene_per_Fam_ARG42 %>% mutate(Isolate="ARG_4_2") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famIRNGOR2= Gene_per_Fam_IRNGOR2 %>% mutate(Isolate="IRN_GOR_2") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famISR7= Gene_per_Fam_ISR7 %>% mutate(Isolate="ISR_7") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famJPNCHIKA= Gene_per_Fam_JPNCHIKA %>% mutate(Isolate="JPN_CHIKA") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famUSA2= Gene_per_Fam_USA2 %>% mutate(Isolate="USA_2") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famBgs1459= Gene_per_Fam_Bgs1459 %>% mutate(Isolate="Bgs_1459") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famISR211= Gene_per_Fam_ISR211 %>% mutate(Isolate="ISR_211") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
famTHUN12= Gene_per_Fam_THUN12 %>% mutate(Isolate="THUN12") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)
fam96224new= Gene_per_Fam_96224 %>% mutate(Isolate="CHE_96224") %>% group_by(Fam,Isolate) %>% summarise(NumMemb=length(Fam))  %>% select(Fam,NumMemb,Isolate)


famAllfsp=rbind(fam96224new,famCHN1740,famJPNCHIKA,famISR7,famIRNGOR2,famARG42,famCHN5227,famUSA2,famISR211,famBgs1459,famTHUN12)
famAllfsp
famAllfspdcast=dcast(famAllfsp,Fam~Isolate,value.var="NumMemb")
famAllfspdcast[,2:12]

scaledforcolofsp=as.data.frame(t(as.data.frame(apply(famAllfspdcast[,2:12],1,scale))))
scaledforcolofsp
head(famAllfspdcast[,2:12])

names(scaledforcolofsp)=c( "ARG_4_2","Bgs_1459","CHE_96224","CHN_17_40","CHN_52_27","IRN_GOR_2","JPN_CHIKA","ISR_211","ISR_7","THUN12","USA_2")
head(scaledforcolofsp)
scaledforcolofsp$Fam=famAllfspdcast$Fam
head(scaledforcolofsp)
scaledforcolofspm=melt(scaledforcolofsp)
head(scaledforcolofspm)
names(scaledforcolofspm)=c("Fam","Isolate","Scaled")

FinalMfsp=merge(famAllfsp,scaledforcolofspm) 
FinalMFilteredfsp=FinalMfsp %>% filter(Scaled != "NaN")
FinalMfsp %>% group_by(Isolate,Fam) %>% summarise(EffinFam=sum(NumMemb))

FinalMfsp %>% group_by(Fam) %>%select(Fam) %>% unique()
FinalMFilteredfsp %>% group_by(Fam) %>%select(Fam) %>% unique()
#62/171

FinalMFilteredfspPLOT <- ggplot(FinalMFilteredfsp) + geom_point(aes(x = Fam, y = Isolate,size=NumMemb, fill=Scaled),shape=21) + 
  theme_bw() + scale_fill_gradient2(high="red",mid="white",low="blue") + 
  xlab("Effector Family")+
  labs(size="# of Members", fill="Scaled")+
  theme(axis.text.x = element_text(color = "black", size = 19, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 19, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title=element_text(size=20))+
  theme(legend.text=element_text(size=rel(1.6)))+
  scale_size_area(breaks = c(1,5,10,30,60,90), max_size = 18)
FinalMFilteredfspPLOT

png("Output_Fig6C.png", width = 2000, height = 600)
FinalMFilteredfspPLOT
dev.off()
pdf("Output_Fig6C.pdf", width = 20, height = 6)
FinalMFilteredfspPLOT
dev.off()





































###S22A Fig. - Subset of isolates per population - Duplications
#All genes duplications, 8 random isolates per population

normCovAllfinalRegionsDupSum_v3_8perpop <- read.delim("Input_S22A_Fig", header = T, sep = " ")
coh400_8perpop_allpops_dups <- normCovAllfinalRegionsDupSum_v3_8perpop
ggplot(coh400_8perpop_allpops_dups) +
  aes(x = Population, y = Dup_per_Isol, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_8perpop_allpops_dups$Dup_per_Isol, coh400_8perpop_allpops_dups$Population, p.adjust.method = "holm")
x <- which(names(coh400_8perpop_allpops_dups) == "Population") # name of grouping variable
y <- which(names(coh400_8perpop_allpops_dups) == "Dup_per_Isol") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Bgs", "Bgdic"), c("Bgdic", "ARG"), c("ARG", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "JPN"),c("JPN", "KAZ"), c("KAZ", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_8perpop_allpops_dups,
                   x = colnames(coh400_8perpop_allpops_dups[j]), y = colnames(coh400_8perpop_allpops_dups[i]),
                   color = colnames(coh400_8perpop_allpops_dups[j]), ylab = "Gene Duplications", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_8perpop_allpops_dups[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400 <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 300) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 275) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TUR", label.y = 270) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000", "KAZ"="#661100"))
ggpar(p400, orientation = "vertical")

pdf("Output_S22A_Fig.pdf", width=12, height=4)
ggpar(p400, orientation = "vertical")
dev.off()
png("Output_S22A_Fig.png", width=9000, height=3000, res=600)
ggpar(p400, orientation = "vertical")
dev.off()



###S22B Fig. - All isolates - Duplications
#All genes duplications with all isolates per population
normCovAllfinalRegionsDupSum_v2 <- read.delim("Input_S22B_Fig", header = T, sep = " ")
coh400_allpops_dups <- normCovAllfinalRegionsDupSum_v2
ggplot(coh400_allpops_dups) +
  aes(x = Population, y = Dup_per_Isol, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_allpops_dups$Dup_per_Isol, coh400_allpops_dups$Population, p.adjust.method = "holm")
x <- which(names(coh400_allpops_dups) == "Population") # name of grouping variable
y <- which(names(coh400_allpops_dups) == "Dup_per_Isol") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Bgs", "Bgdic"), c("Bgdic", "ARG"), c("ARG", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "JPN"),c("JPN", "KAZ"), c("KAZ", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_allpops_dups,
                   x = colnames(coh400_allpops_dups[j]), y = colnames(coh400_allpops_dups[i]),
                   color = colnames(coh400_allpops_dups[j]), ylab = "Gene Duplications", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_allpops_dups[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400 <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 300) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 275) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TUR", label.y = 270) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000", "KAZ"="#661100"))
ggpar(p400, orientation = "vertical")

pdf("Output_S22B_Fig.pdf", width=12, height=4)
ggpar(p400, orientation = "vertical")
dev.off()
png("Output_S22B_Fig.png", width=9000, height=3000, res=600)
ggpar(p400, orientation = "vertical")
dev.off()



###S22C Fig. - Subset of isolates per population - Deletions
#All genes deletions, 8 random isolates per population
normCovAllfinalRegionsDelSum_v3_8perpop <- read.delim("Input_S22C_Fig", header = T, sep = " ")

coh400_8perpop_dels <- normCovAllfinalRegionsDelSum_v3_8perpop
ggplot(coh400_8perpop_dels) +
  aes(x = Population, y = Del_per_Isol, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_8perpop_dels$Del_per_Isol, coh400_8perpop_dels$Population, p.adjust.method = "holm")
x <- which(names(coh400_8perpop_dels) == "Population") # name of grouping variable
y <- which(names(coh400_8perpop_dels) == "Del_per_Isol") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("ARG", "Bgdic"), c("Bgdic", "Bgs"), c("Bgs", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "JPN"),c("JPN", "KAZ"), c("KAZ", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_8perpop_dels,
                   x = colnames(coh400_8perpop_dels[j]), y = colnames(coh400_8perpop_dels[i]),
                   color = colnames(coh400_8perpop_dels[j]), ylab = "Gene Deletions", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_8perpop_dels[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400_8perpop <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 300) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 275) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TUR", label.y = 270) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000", "KAZ"="#661100"))
#  scale_colour_manual(values = safe_colorblind_palette13)
ggpar(p400_8perpop, orientation = "vertical")

pdf("Output_S22C_Fig.pdf", width=12, height=4)
ggpar(p400_8perpop, orientation = "vertical")
dev.off()
png("Output_S22C_Fig.png", width=9000, height=3000, res=600)
ggpar(p400_8perpop, orientation = "vertical")
dev.off()



###S22D Fig. - All isolates - Deletions
#All genes deletions with all isolates per population
normCovAllfinalRegionsDelSum_v2 <- read.delim("Input_S22D_Fig", header = T, sep = " ")
coh400_allpops_dels <- normCovAllfinalRegionsDelSum_v2
ggplot(coh400_allpops_dels) +
  aes(x = Population, y = Del_per_Isol, color = Population) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh400_allpops_dels$Del_per_Isol, coh400_allpops_dels$Population, p.adjust.method = "holm")
x <- which(names(coh400_allpops_dels) == "Population") # name of grouping variable
y <- which(names(coh400_allpops_dels) == "Dup_per_Isol") # names of variables to test
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Bgs", "Bgdic"), c("Bgdic", "ARG"), c("ARG", "CHN"), c("CHN", "EGY"), c("EGY", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "JPN"),c("JPN", "KAZ"), c("KAZ", "RUS"), c("RUS", "TUR"), c("TUR", "USA")) # comparisons for post-hoc tests

for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh400_allpops_dels,
                   x = colnames(coh400_allpops_dels[j]), y = colnames(coh400_allpops_dels[i]),
                   color = colnames(coh400_allpops_dels[j]), ylab = "Gene Deletions", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh400_allpops_dels[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p400_allpops <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 300) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 275) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "TUR", label.y = 270) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "AUS"="#CC6677", "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "Bgs"="#882255", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgdic"="#000000", "KAZ"="#661100"))
ggpar(p400_allpops, orientation = "vertical")

pdf("Output_S22D_Fig.pdf", width=12, height=4)
ggpar(p400_allpops, orientation = "vertical")
dev.off()
png("Output_S22D_Fig.png", width=9000, height=3000, res=600)
ggpar(p400_allpops, orientation = "vertical")
dev.off()







###S23A Fig. - Heatmap of CNVs of effector genes in various isolates (Presence Abscence Variation)
#All Effector genes CNVs in various isolates
library(dplyr)
library(reshape2)
library(ggplot2)
library(pheatmap)

normCovAll=read.table("Input_S23_Fig_df_results_Bg_ALL_final_subset_for_pangenome_paper.csv",h=T,sep=",")
subset8isl=read.table("Input_S23_Fig_Isolates_8Random.txt",h=T)
islnamessubset8ils=as.vector(subset8isl$Isolate)
islnamessubset8ils=c("gene",islnamessubset8ils)
islnamessubset8ils

normCovAll8isl = normCovAll %>% select(any_of(islnamessubset8ils))

Names844=read.table("Input_S23_Fig_Names_of_844_effectors.txt")
vNames844=Names844$V1
Effectors8Isl=subset(normCovAll8isl, gene %in% vNames844)
row.names(Effectors8Isl)=Effectors8Isl$gene
Effectors8Isl$gene=NULL
length(normCovAll8isl)
col_names=read.table("Input_S23_Fig_Isolates_8Random.txt",h=T,row.names = "Isolate")
col_names=as.data.frame(col_names)
length(col_names$Population)
length(Effectors8Isl)
ann_colors = list(
  Population = c("ARG"="#88CCEE",  "CHN"="#DDCC77", "EGY"="#117733", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "RUS"="#661100", "TUR"="#6699CC", "USA"="#888888", "Bgd"="#000000"))

pheatmapEff=pheatmap(Effectors8Isl,show_rownames=FALSE,show_colnames=FALSE, annotation_col=col_names, annotation_colors =ann_colors, 
                     clustering_distance_cols="euclidean", cluster_rows =TRUE, annotation_names_col=FALSE, treeheight_col = FALSE, fontsize = 15)
ggsave("Output_S23A_Fig_Effectors.pdf",plot=pheatmapEff, width=10, height=5)
ggsave("Output_S23A_Fig_Effectors.png",plot=pheatmapEff, width=10, height=5)

###S23B Fig. - Heatmap of CNVs of non effector genes in various isolates (Presence Abscence Variation)
#All Non Effector genes CNVs in various isolates
NonEffectors8Isl=subset(normCovAll8isl, !(gene %in% vNames844))
vRandom844=as.vector(sample(NonEffectors8Isl$gene,844))
NonEffectors8Isl844=subset(NonEffectors8Isl, gene %in% vRandom844)
row.names(NonEffectors8Isl844)=NonEffectors8Isl844$gene
NonEffectors8Isl844$gene=NULL

pheatmapNonEff=pheatmap(NonEffectors8Isl844,show_rownames=FALSE,show_colnames=FALSE, annotation_col=col_names, annotation_colors =ann_colors, 
                        clustering_distance_cols="euclidean", cluster_rows =TRUE, annotation_names_col=FALSE, treeheight_col = FALSE, fontsize = 15)
ggsave("Output_S23B_Fig_NonEffectors.pdf",plot=pheatmapNonEff, width=10, height=5)
ggsave("Output_S23B_Fig_NonEffectors.png",plot=pheatmapNonEff, width=10, height=5)
