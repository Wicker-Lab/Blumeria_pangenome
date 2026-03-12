library(ggplot2)

gene_stats_bg_all <- read.delim("Input_S9A_Fig_all_genes_stats_v3", header = T)
gene_stats_bg_all$Type <- factor(gene_stats_bg_all$Type, levels = c("Non-Effectors", "Grouped Eff", "Ungrouped Eff"))

b <- ggplot(gene_stats_bg_all, aes(fill=Type, x=Isolate_CDS, y=Gene_Number, label=Gene_Number)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  geom_text(size = 7, position = position_stack(vjust = 0.7))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + #theme(legend.text = element_text(face = "bold"))+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 60, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.text=element_text(size=25), legend.title = element_text(size = 25))+
  scale_fill_manual("legend", values = c("Non-Effectors"="#88CCEE", "Grouped Eff"="#CC6677", "All Effectors"="#DDCC77", "Ungrouped Eff"="#117733"))+
  xlab("Isolates")+ ylab("Number of Genes") + guides(fill=guide_legend(title="Type of genes"))
b  

png("Output_S9A_Fig.png", width = 1200, height = 800)
b
dev.off()
pdf("Output_S9A_Fig.pdf", width = 12, height = 8)
b
dev.off()















