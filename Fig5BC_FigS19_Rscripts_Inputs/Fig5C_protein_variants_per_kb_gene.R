library(ggplot2)
library(dplyr)
library(tidyr)


# set the path to the correct working directory
my_path <- "~/data/dir_Blumeria/dir_Blumeria_pan_genome_2024/dir_R_scripts/dir_for_github"

setwd(my_path)
list.files()


# plot for potein variants per kb gene, chr-11 long and short arms separate ------------------

infile <- "all_allele_profiles_summary_DE_chr11LR"
df <-read.table(infile, header=T) 
head(df)

title <-"Protein variants per kb gene"

colors <- c("#FF0000","#999999")


p <- ggplot(df,aes(x=norm_prot,y=group2,color=group))+
  geom_boxplot()+
  scale_color_manual(values=colors)+
  geom_jitter(width=0, height=0.2,size=0.1,alpha=0.3)+
  xlim(0,100)+
  ylab("Chromosome")+
  xlab("Protein variants per kb gene")+
  ggtitle(title)+
  theme_classic()+
  scale_y_discrete(limits = rev)+
  theme(axis.text.y = element_text(hjust = 0))

p


# write plot to output file
outpng <- "Fig5C_protein_variants_per_kb_CDS.png"
path <- getwd()
wide <- 6
high <- 8
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=300)


