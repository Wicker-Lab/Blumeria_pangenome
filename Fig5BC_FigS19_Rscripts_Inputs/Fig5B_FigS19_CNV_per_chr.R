library(ggplot2)
library(dplyr)
library(tidyr)


# set the path to the correct working directory
my_path <- "~/data/dir_Blumeria/dir_Blumeria_pan_genome_2024/dir_R_scripts/dir_for_github"

setwd(my_path)
list.files()


# gene presence-absence (GPA) plots --------
# CNV per gene per chromosome

infile <- "Bgt_399_isolates_CNV_per_gene"
infile
df <-read.table(infile, header=T) 
head(df)


df_long <- pivot_longer(df,names_to="CNV_type",values_to="counts",cols=8:9)
head(df_long)


PAP <-"CNV"


# make plot for single chromosome as example -----
chr <- 1
title <-paste(PAP," per gene ",chr,sep='')

colors = c("#FF0000","#999999")

p <- ggplot(df_long[df_long$chr==chr,],aes(x=be/1000000,y=counts,color=group))+
  geom_point(size=2,alpha=0.5,aes(shape=CNV_type))+
  ylab(PAP)+
  xlab(paste("Position chr",chr,"[Mb]",sep=''))+
  ylim(0,500)+
  ggtitle(title)+
  ylab(paste("Isolates with ",PAP,sep=''))+
  theme_classic()+
  scale_color_manual(values = colors)+
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 1))
  
p




# loop over all chromosomes for gene presence-absence (GPA) plots --------
for (i in 1:11) {

chr <- i

p <- ggplot(df_long[df_long$chr==chr,],aes(x=be/1000000,y=counts,color=group))+
  geom_point(size=2,alpha=0.5,aes(shape=CNV_type))+
  ylab(PAP)+
  xlab(paste("Position chr",chr,"[Mb]",sep=''))+
  ylim(0,400)+
  ylab(paste("Isolates with ",PAP,sep=''))+
  theme_classic()+
  scale_color_manual(values = colors)+
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 1))

p



# write plot to output file
outpng <- paste("CNV_per_gene_chr",chr,".png",sep='')
path <- getwd()
wide <- 14
high <- 1.8
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=300)


}


