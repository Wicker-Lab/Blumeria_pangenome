###Figure 8A - TE insertions per isolate for populations
#Stacked barplot
#Create 'Relationships type of TEs' plot
library("ggplot2")
#load result file
file = read.csv(file = "Input_Fig8A.csv", sep = ',', header = TRUE)

g4 <- ggplot(data = file, aes(x = Type, fill = Isolate)) + geom_bar(position = "fill") + ylab("proportion") +
  stat_count(geom = "text", aes(label = stat(count)), position=position_fill(vjust=0.5), colour="white")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(title = "Relationships type of TEs") +
  scale_fill_manual(values=c("ARG_4_2"="#88CCEE", "CHN_17_40"="#DDCC77", "CHN_52_27"="#DDCC77", "CHE_96224"="#332288", "IRN_GOR_2"="#AA4499", 
                             "ISR_7"="#44AA99", "JPN_CHIKA"="#999933", "Bgs_1459"="#882255", "USA_2"="#888888", "Bgd_ISR_211"="#000000", "Bgtl_THUN_12"="#444444",
                             "KAZ"="#661100", "AUS"="#CC6677", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "EGY"="#117733"))+
  theme(axis.text.x=element_text(color = "black", size = 15, angle = 40, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(size=20))
g4

pdf("Output_Fig8A.pdf", width=12, height=6)
g4
dev.off()
png("Output_Fig8A.png", width=8000, height=4000, res=600)
g4
dev.off()





###Figure 8B - TE proximity to genes (Effectors and non-effectors)
#TE proximity to genes (Effectors and non-effectors)
library(ggpattern)
library(ggplot2)
tedens_pang <- read.delim("Input_Fig8B.csv", sep = ",")

tedens_pang_initial <- ggplot(tedens_pang, aes(fill=Isolate, y=Upstream, x=reorder(TE.family, -TIPs, sum))) + 
  geom_bar(position="stack", stat="identity")+
  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG_4_2"="#88CCEE", "CHN_17_40"="#DDCC77", "CHN_52_27"="#DDCC77", "CHE_96224"="#332288", "IRN_GOR_2"="#AA4499", 
                             "ISR_7"="#44AA99", "JPN_CHIKA"="#999933", "Bgs_1459"="#882255", "USA_2"="#888888", "Bgd_ISR_211"="#000000", "Bgtl_THUN_12"="#444444",
                             "KAZ"="#661100", "AUS"="#CC6677", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "EGY"="#117733"))+
  xlab("Isolate")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + #labs(color='Species') + 
  theme(axis.text.x=element_text(color = "black", size = 15, angle = 60, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size = 15))

tedens_pang$Window = factor(tedens_pang$Window, levels = c("0-500bp", "501-1000bp", "1001-2000bp", "2000+bp"), ordered = TRUE)

library(grid)
g.mid<-ggplot(tedens_pang,aes(x=1,y=Window))+geom_text(aes(label=Window))+
  geom_segment(aes(x=0.94,xend=0.96,yend=Window))+
  geom_segment(aes(x=1.04,xend=1.065,yend=Window))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(), panel.grid=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.background=element_blank(), axis.text.x=element_text(color=NA), axis.ticks.x=element_line(color=NA), plot.margin = unit(c(1,-1,1,-1), "mm"))+
  scale_fill_manual(values=c("ARG_4_2"="#88CCEE", "CHN_17_40"="#DDCC77", "CHN_52_27"="#DDCC77", "CHE_96224"="#332288", "IRN_GOR_2"="#AA4499", 
                             "ISR_7"="#44AA99", "JPN_CHIKA"="#999933", "Bgs_1459"="#882255", "USA_2"="#888888", "Bgd_ISR_211"="#000000", "Bgtl_THUN_12"="#444444",
                             "KAZ"="#661100", "AUS"="#CC6677", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "EGY"="#117733"))+
  theme(axis.text.y=element_text(size=18))

g1 <- ggplot(data = tedens_pang, aes(x = Window, y = Upstream, fill=Isolate, pattern=Type)) +
  ggtitle("Upstream") +
  geom_bar_pattern(position = position_dodge(preserve = "single"), stat="identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse() + coord_flip() +  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG_4_2"="#88CCEE", "CHN_17_40"="#DDCC77", "CHN_52_27"="#DDCC77", "CHE_96224"="#332288", "IRN_GOR_2"="#AA4499", 
                             "ISR_7"="#44AA99", "JPN_CHIKA"="#999933", "Bgs_1459"="#882255", "USA_2"="#888888", "Bgd_ISR_211"="#000000", "Bgtl_THUN_12"="#444444",
                             "KAZ"="#661100", "AUS"="#CC6677", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "EGY"="#117733"))+
  theme(legend.position="none")+
  scale_pattern_manual(values = c('Effector' = "stripe", 'Non-Effector' = "none"))+
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none")))+
  theme(axis.text.x=element_text(size=21))+
  theme(text=element_text(size=18))

g2 <- ggplot(data = tedens_pang, aes(x = Window, y = Downstream, fill=Isolate, pattern=Type)) +xlab(NULL)+
  ggtitle("Downstream") +
  geom_bar_pattern(position = position_dodge(preserve = "single"), stat="identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  coord_flip() +  #theme(legend.position="none") +
  scale_fill_manual(values=c("ARG_4_2"="#88CCEE", "CHN_17_40"="#DDCC77", "CHN_52_27"="#DDCC77", "CHE_96224"="#332288", "IRN_GOR_2"="#AA4499", 
                             "ISR_7"="#44AA99", "JPN_CHIKA"="#999933", "Bgs_1459"="#882255", "USA_2"="#888888", "Bgd_ISR_211"="#000000", "Bgtl_THUN_12"="#444444",
                             "KAZ"="#661100", "AUS"="#CC6677", "OTH"="#882255", "RUS"="#661100", "TUR"="#6699CC", "EGY"="#117733"))+
  scale_pattern_manual(values = c('Effector' = "stripe", 'Non-Effector' = "none"))+
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none")))+
  theme(axis.text.x=element_text(size=21))+
  theme(text=element_text(size=18))
g2
library(gridExtra)
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))
grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))

png("Output_Fig8B.png")
grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
dev.off()
pdf("Output_Fig8B.pdf", width = 12, height = 6)
grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
dev.off()
pdf("Output_Fig8B_legend.pdf", width = 12, height = 6)
g2
dev.off()










###Figure 8C - TE insertions per isolate for populations
#Boxplot of the number of TIPs per isolate for populations with eight or more isolates (isolates with less than 2000 TIPs are shown, extreme outliers were excluded, because they were interpreted as technical artifacts)
# library
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(multcomp)
library(patchwork)
library(palmerpenguins)
coh177_allpops_tips_counts <- read.delim("Input_Fig8C", sep = "\t")
ggplot(coh177_allpops_tips_counts) +
  aes(x = Group, y = TIPs, color = Group) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme_minimal()

pairwise.t.test(coh177_allpops_tips_counts$TIPs, coh177_allpops_tips_counts$Group, p.adjust.method = "holm")
x <- which(names(coh177_allpops_tips_counts) == "Group") # name of grouping variable
y <- which(names(coh177_allpops_tips_counts) == "TIPs") # names of variables to test

method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("CHN", "JPN"), c("JPN", "USA"), c("USA", "ARG"), c("ARG", "EUR"), c("EUR", "IRN"), c("IRN", "ISR"), c("ISR", "Bgdic")) # comparisons for post-hoc tests
coh177_allpops_tips_counts$Group <- factor(coh177_allpops_tips_counts$Group , levels=c("CHN", "JPN", "USA","ARG", "EUR", "IRN", "ISR", "Bgdic"))


for (i in y) {
  for (j in x) {
    p <- ggboxplot(coh177_allpops_tips_counts,
                   x = colnames(coh177_allpops_tips_counts[j]), y = colnames(coh177_allpops_tips_counts[i]),
                   color = colnames(coh177_allpops_tips_counts[j]), ylab = "TIPs", 
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(coh177_allpops_tips_counts[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

options(scipen = 999)
p177 <- p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)), method = method1, label.y = 1550) +
  stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format", label.y = 1300) +# remove if p-value of ANOVA or Kruskal-Wallis test >= alpha)
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Bgdic", label.x = 0, label.y = 1300) + 
  scale_color_manual(values=c("ARG"="#88CCEE", "CHN"="#DDCC77", "EUR"="#332288", "IRN"="#AA4499", "ISR"="#44AA99", "JPN"="#999933", "USA"="#888888", "Bgdic"="#000000"))
ggpar(p177, orientation = "horizontal", ylim = c(0, 1500))

pdf("Output_Fig8C.pdf", width=6, height=8)
ggpar(p177, ylim = c(0, 1500), orientation = "horizontal")
dev.off()
png("Output_Fig8C.png", width=9000, height=3000, res=600)
ggpar(p177, orientation = "horizontal")
dev.off()






###Figure 8D - TE insertions in populations
#Analyses of TE insertions (TIPs) in the different wheat powdery mildew populations normalized for sample size for each population
library("ggplot2")
library(dplyr)
library(data.table)

tips_all <- read.delim("Input_Fig8D")

tips_coh180_rnd8 <- ggplot(tips_all, aes(fill=Population, y=TIPs, x=reorder(TEfamily2, -TIPs, sum))) + 
  geom_bar(position="stack", stat="identity")+
  theme(legend.position="none") +
  scale_fill_manual(values=c("ARG"="#88CCEE", "CHN"="#DDCC77", "EUR"="#332288", "IRN"="#AA4499", 
                             "ISR"="#44AA99", "JPN"="#999933", "USA"="#888888", "Bgdic"="#000000"))+
  xlab("TE Family")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + #labs(color='Species') + 
  theme(axis.text.x=element_text(color = "black", size = 15, angle = 60, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=15), legend.title = element_text(size = 15))
tips_coh180_rnd8

png("Output_Fig8D.png")
tips_coh180_rnd8
dev.off()
pdf("Output_Fig8D.pdf", width = 6, height = 6)
tips_coh180_rnd8
dev.off()












###S10 Fig. - Genome statistics of the pangenome
#libraries
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(multcomp)
library(patchwork)
library(palmerpenguins)
library(tidyverse)

bg_chr_statistics_v2 <- read.delim("Input_S10_Fig", header = T)
bg_chr_statistics_v2
# make plot
pchr <- ggplot(bg_chr_statistics_v2, aes(x=factor(Chromosome), y=Bases, fill=Isolate)) + geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("Bgt_ARG_4_2_v1"="#88CCEE", "Bgt_CHN_17_40_v1"="#DDCC77", "Bgt_CHN_52_27_v1"="#DDCC77", 
                             "Bgt_CHE96224_v3_16"="#332288", "Bgt_IRN_GOR_2_v1"="#AA4499", "Bgt_ISR7_v1_4"="#44AA99", 
                             "Bgt_JPN_CHIKA_v1"="#999933", "Bgs_1459"="#882255", "Bgt_USA_2_v1"="#888888", "Bgd_ISR211_v1"="#000000", 
                             "Bgtl_THUN12_v1_2"="#444444"))
pchr
pchr2 <- pchr + ylab("Base pairs") + theme_light() + xlab("Chromosomes") + labs(fill='Isolates') + theme(legend.text=element_text(size=15),legend.title = element_text(size = 15))+
  theme(axis.text.x = element_text(color = "black", size = 17, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 17, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))
pchr2
options(scipen = 999)
ggpar(pchr2, orientation = "vertical")

pdf("Output_S10_Fig.pdf", width=12, height=4)
ggpar(pchr2, orientation = "vertical")
dev.off()
png("Output_S10_Fig.png", width=9000, height=3000, res=600)
ggpar(pchr2, orientation = "vertical")
dev.off()






###S25 Fig. - TE analyses of Blumeria graminis pangenome
#TE counts of TE superfamilies (and other TEs) for the pangenome
#Load ggplot2
library(ggplot2)

te_knowntes_bg_pangenome <- read.delim("Input_S25_Fig", header = T)

tes_pan <- ggplot(te_knowntes_bg_pangenome, aes(x = Isolate, y = Count, fill = Class)) + 
  geom_col() + coord_polar("y") +
  scale_fill_manual(values=c("XXX"="#000000", "DHH"="#924900", "DTT"="#920000", "DTX"="#ff6db6", "DTC"="#ffb6db", "DTM"="#b66dff", "RII"="#b6dbff", 
                            "RIJ"="#6db6ff", "RIX"="#490092", "RLC"="#004949", "RLG"="#009292", "RLX"="#24ff24", "RSX"="#ffff6d"))
ggpar(tes_pan, orientation = "vertical")

pdf("Output_S25_Fig.pdf", width=12, height=8)
ggpar(tes_pan, orientation = "vertical")
dev.off()
png("Output_S25_Fig.png", width=9000, height=6000, res=600)
ggpar(tes_pan, orientation = "vertical")
dev.off()

