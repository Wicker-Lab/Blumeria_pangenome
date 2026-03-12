library(gdsfmt)
library(ggplot2)
library(dplyr)
library(stringr)

bgd_isr211 <- read.delim("Input_Bgd_ISR_211_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgd_isr211_dist <- ggplot(bgd_isr211, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
  color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgd_ISR_211_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
#bgd_isr211_dist

bgt_usa2 <- read.delim("Input_Bgt_USA_2_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_usa2_dist <- ggplot(bgt_usa2, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_USA_2_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())
#bgt_usa2_dist

bgs_1459 <- read.delim("Input_Bgs_1459_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgs_1459_dist <- ggplot(bgs_1459, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgs_1459_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
#bgs_1459_dist

bgtl_THUN12 <- read.delim("Input_Bgt_THUN12_CDS_v1_2_aa_all_sizedistr_v3", header=TRUE)
bgtl_THUN12_dist <- ggplot(bgtl_THUN12, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgtl_THUN_12_CDS_v1_2")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
#bgtl_THUN12_dist

bgt_arg42 <- read.delim("Input_Bgt_ARG_4_2_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_arg42_dist <- ggplot(bgt_arg42, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_ARG_4_2_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())
#bgt_arg42_dist

bgt_che96224 <- read.delim("Input_Bgt_CDS_v4_23_aa_all_sizedistr_v3", header=TRUE)
bgt_che96224_dist <- ggplot(bgt_che96224, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_CHE_96224_v4_23_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
#bgt_che96224_dist

bgt_isr7 <- read.delim("Input_Bgt_ISR_7_CDS_v1_4_1_aa_all_sizedistr_v3", header=TRUE)
bgt_isr7_dist <- ggplot(bgt_isr7, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_ISR_7_CDS_v1_4")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.y=element_blank())
#bgt_isr7_dist

bgt_chn1740 <- read.delim("Input_Bgt_CHN_17_40_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_chn1740_dist <- ggplot(bgt_chn1740, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_CHN_17_40_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.y=element_blank())
#bgt_chn1740_dist

bgt_chn5227 <- read.delim("Input_Bgt_CHN_52_27_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_chn5227_dist <- ggplot(bgt_chn5227, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_CHN_52_27_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.y=element_blank())
#bgt_chn5227_dist

bgt_irngor2 <- read.delim("Input_Bgt_IRN_GOR_2_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_irngor2_dist <- ggplot(bgt_irngor2, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_IRN_GOR_2_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))
#bgt_irngor2_dist

bgt_jpnchika <- read.delim("Input_Bgt_JPN_CHIKA_CDS_v1_1_aa_all_sizedistr_v3", header=TRUE)
bgt_jpnchika_dist <- ggplot(bgt_jpnchika, aes(x=size)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  labs(y= "Number of proteins", x = "Size of proteins")+
  ggtitle("Bgt_JPN_CHIKA_CDS")+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 21, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 25, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
#bgt_jpnchika_dist


require(gridExtra)
grid.arrange(bgt_usa2_dist, bgd_isr211_dist, bgtl_THUN12_dist, bgs_1459_dist, bgt_arg42_dist, 
             bgt_jpnchika_dist, bgt_che96224_dist, bgt_isr7_dist, bgt_irngor2_dist, bgt_chn1740_dist, 
             bgt_chn5227_dist, ncol=4, nrow=3)

pdf("Output_S9B_Fig.pdf", width = 16, height = 12)
grid.arrange(bgt_usa2_dist, bgd_isr211_dist, bgtl_THUN12_dist, bgs_1459_dist, bgt_arg42_dist, 
             bgt_jpnchika_dist, bgt_che96224_dist, bgt_isr7_dist, bgt_irngor2_dist, bgt_chn1740_dist, bgt_chn5227_dist, ncol=4)
dev.off()

png("Output_S9B_Fig.png", width = 9900, height = 6000, res = 600)
grid.arrange(bgt_usa2_dist, bgd_isr211_dist, bgtl_THUN12_dist, bgs_1459_dist, bgt_arg42_dist, 
             bgt_jpnchika_dist, bgt_che96224_dist, bgt_isr7_dist, bgt_irngor2_dist, bgt_chn1740_dist, bgt_chn5227_dist, ncol=4)
dev.off()

