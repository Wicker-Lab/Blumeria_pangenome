library(gdsfmt)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)


setwd("~/data/dir_Blumeria/dir_Blumeria_pan_genome_2024/dir_R_scripts/dir_for_github/Fig6CDE_R_Perl_scripts_Inputs")
list.files()

colors <- c("#FF0000","#0088FF","#FF00FF")



# boxplot with regression curve non-effectors, exclude 1459 -----------------

infile <- "homologs_bi_directional_PAP_matrix_non_Eff_excl_THUN12_permutations"
df0 <-read.table(infile, sep='\t', header=T)
head(df0)
set <- "NE_ex1549"

# select only sample that include 1459
df <- df0[df0$Group=='excl_1459',]
head(df)


# caclulate  medians for each sample step
medians <- df %>%
  group_by(isolates) %>%
  summarise(core_genes_median = median(core_genes), .groups = "drop") %>%
  mutate(isolates_num = as.numeric(factor(isolates)))

medians_df <- as.data.frame(medians)
medians_df



# predict the convergence value of the curve, using the median values 
# for 8 medians
x <- 1:8
y <- medians$core_genes_median

# fit asymptotic regression 
fit <- nls(y ~ Asym + (R0 - Asym) * exp(-k * x), start = list(Asym = min(y), R0 = max(y), k = 0.1))

summary(fit)
# get convergence value 
convergence_val <-coef(fit)["Asym"]
convergence_val

title <- "Core genome size estimate isolates including Bgs1459\nNon-effector genes"
colors <- c("#FF0000","#0088FF","#FF00FF")


# draw box plot with regression curve and asymptotic value ++++++++++++++++++
p <- ggplot(df, aes(x = as.factor(isolates), y = core_genes)) +
  geom_boxplot(color="#FF0000") +
  geom_jitter(size=0.4,alpha=0.5,width=0.2,color="#FF0000")+
  ylim(6900,7800)+
  scale_x_discrete(limits = as.character(1:9))+
  #xlim(1,9)+
  # add median values
  geom_point(
    data = medians,
    aes(x = as.factor(isolates), y = core_genes_median),
    inherit.aes = FALSE,
    color = "black",
    size = 2
    ) +
  geom_text(
    data = medians,
    aes(x = as.factor(isolates), y = core_genes_median,label=round(core_genes_median),hjust=0.5,vjust=-5),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  ggtitle(title)+
  xlab("Number of Blumera isolates added to reference CHE96224")+
  ylab("Number of genes in core genome")+
  # add curve 
  geom_smooth(
    data = medians, linetype =2, linewidth=0.7,
    aes(x = isolates_num, y = core_genes_median),
    inherit.aes = FALSE,
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  # horizontal line at asymptotic value
  geom_hline(yintercept = convergence_val, linetype = "dashed", color = "#BB00FF")+

  # add asymptotic value 
  annotate("text",x = 8.5,y = convergence_val, label = round(convergence_val),color = "#BB00FF",hjust = 0,vjust=2)+
  theme_classic()

p


# make output image +++++++++++++++++++++++++++++++++++
outfile <- paste(infile,"_boxplot_curve",set,sep='')
outpng <- paste(outfile,".png",sep='')
high <- 6
wide  <- 5.5
path <- getwd()
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=400)








# boxplot with regression curve non-effectors, include 1459 ---------------------------

infile <- "homologs_bi_directional_PAP_matrix_non_Eff_excl_THUN12_permutations"
df0 <-read.table(infile, sep='\t', header=T)
head(df0)
set <- "NE_in1549"

# select only sample that include 1459
df <- df0[df0$Group=='incl_1459',]
head(df)


# caclulate  medians for each sample step
medians <- df %>%
  group_by(isolates) %>%
  summarise(core_genes_median = median(core_genes), .groups = "drop") %>%
  mutate(isolates_num = as.numeric(factor(isolates)))

medians_df <- as.data.frame(medians)
medians_df



# predict the convergence value of the curve, using the median values
# for medians
x <- 1:9
y <- medians$core_genes_median

# fit asymptotic regression 
fit <- nls(y ~ Asym + (R0 - Asym) * exp(-k * x), start = list(Asym = min(y), R0 = max(y), k = 0.1))

summary(fit)
# get convergence value 
convergence_val <-coef(fit)["Asym"]
convergence_val

title <- "Core genome size estimate isolates including Bgs1459\nNon-effector genes"
colors <- c("#FF0000","#0088FF","#FF00FF")


# draw box plot with regression curve and asymptotic value ++++++++++++++++++
p <- ggplot(df, aes(x = factor(isolates), y = core_genes)) +
  geom_boxplot(color="#0088FF") +
  geom_jitter(size=0.4,alpha=0.5,width=0.2,color="#0088FF")+
  ylim(6900,7800)+
  # add median values
  geom_point(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  geom_text(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median,label=round(core_genes_median),hjust=0.5,vjust=3.5),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  ggtitle(title)+
  xlab("Number of Blumera isolates added to reference CHE96224")+
  ylab("Number of genes in core genome")+
  # add curve 
  geom_smooth(
    data = medians, linetype =2, linewidth=0.7,
    aes(x = isolates_num, y = core_genes_median),
    inherit.aes = FALSE,
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  # horizontal line at asymptotic value
  geom_hline(yintercept = convergence_val, linetype = "dashed", color = "#BB00FF")+
  
  # add asymptotic value 
  annotate("text",x = 8.5,y = convergence_val, label = round(convergence_val),color = "#BB00FF",hjust = 0,vjust=2)+
  theme_classic()

p


# make output image +++++++++++++++++++++++++++++++++++
outfile <- paste(infile,"_boxplot_curve",set,sep='')
outpng <- paste(outfile,".png",sep='')
high <- 6
wide  <- 5.5
path <- getwd()
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=400)






# boxplot with regression curve effectors, exclude 1459 -----------------

infile <- "homologs_bi_directional_PAP_matrix_Eff_excl_THUN12_permutations"
df0 <-read.table(infile, sep='\t', header=T)
head(df0)
set <- "E_ex1549"

# select only sample that include 1459
df <- df0[df0$Group=='excl_1459',]
head(df)


# caclulate  medians for each sample step
medians <- df %>%
  group_by(isolates) %>%
  summarise(core_genes_median = median(core_genes), .groups = "drop") %>%
  mutate(isolates_num = as.numeric(factor(isolates)))

medians_df <- as.data.frame(medians)
medians_df



# predict the convergence value of the curve, using the median values 
# 8 medians
x <- 1:8
y <- medians$core_genes_median

# fit asymptotic regression 
fit <- nls(y ~ Asym + (R0 - Asym) * exp(-k * x), start = list(Asym = min(y), R0 = max(y), k = 0.1))

summary(fit)
# get convergence value 
convergence_val <-coef(fit)["Asym"]
convergence_val

title <- "Core genome size estimate isolates including Bgs1459\nEffector genes"
colors <- c("#FF0000","#0088FF","#FF00FF")


# draw box plot with regression curve and asymptotic value ++++++++++++++++++
p <- ggplot(df, aes(x = factor(isolates), y = core_genes)) +
  geom_boxplot(color="#FF0000") +
  geom_jitter(size=0.4,alpha=0.5,width=0.2,color="#FF0000")+
  ylim(600,850)+
  scale_x_discrete(limits = as.character(1:9))+
  # add median values
  geom_point(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  geom_text(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median,label=round(core_genes_median),hjust=0.5,vjust=-4),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  ggtitle(title)+
  xlab("Number of Blumera isolates added to reference CHE96224")+
  ylab("Number of genes in core genome")+
  # add curve 
  geom_smooth(
    data = medians, linetype =2, linewidth=0.7,
    aes(x = isolates_num, y = core_genes_median),
    inherit.aes = FALSE,
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  # horizontal line at asymptotic value
  geom_hline(yintercept = convergence_val, linetype = "dashed", color = "#BB00FF")+
  
  # add asymptotic value 
  annotate("text",x = 8.5,y = convergence_val, label = round(convergence_val),color = "#BB00FF",hjust = 0,vjust=2)+
  theme_classic()

p


# make output image 
outfile <- paste(infile,"_boxplot_curve",set,sep='')
outpng <- paste(outfile,".png",sep='')
high <- 6
wide  <- 5.5
path <- getwd()
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=400)








# boxplot with regression curve effectors, include 1459 ---------------------------

infile <- "homologs_bi_directional_PAP_matrix_Eff_excl_THUN12_permutations"
df0 <-read.table(infile, sep='\t', header=T)
head(df0)
set <- "E_in1549"

# select only sample that include 1459
df <- df0[df0$Group=='incl_1459',]
head(df)


# caclulate  medians for each sample step
medians <- df %>%
  group_by(isolates) %>%
  summarise(core_genes_median = median(core_genes), .groups = "drop") %>%
  mutate(isolates_num = as.numeric(factor(isolates)))

medians_df <- as.data.frame(medians)
medians_df



# predict the convergence value of the curve, using the median values
# x = 1:8 for 8 medians
x <- 1:9
#x <- 1:9
y <- medians$core_genes_median

# fit asymptotic regression 
fit <- nls(y ~ Asym + (R0 - Asym) * exp(-k * x), start = list(Asym = min(y), R0 = max(y), k = 0.1))

summary(fit)
# get convergence value 
convergence_val <-coef(fit)["Asym"]
convergence_val

title <- "Core genome size estimate isolates including Bgs1459\nEffector genes"
colors <- c("#FF0000","#0088FF","#FF00FF")


# draw box plot with regression curve and asymptotic value ++++++++++++++++++
p <- ggplot(df, aes(x = factor(isolates), y = core_genes)) +
  geom_boxplot(color="#0088FF") +
  geom_jitter(size=0.4,alpha=0.5,width=0.2,color="#0088FF")+
  ylim(600,850)+
  # add median values
  geom_point(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  geom_text(
    data = medians,
    aes(x = factor(isolates), y = core_genes_median,label=round(core_genes_median),hjust=0.5,vjust=-3.5),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  ggtitle(title)+
  xlab("Number of Blumera isolates added to reference CHE96224")+
  ylab("Number of genes in core genome")+
  # add curve 
  geom_smooth(
    data = medians, linetype =2, linewidth=0.7,
    aes(x = isolates_num, y = core_genes_median),
    inherit.aes = FALSE,
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  # horizontal line at asymptotic value
  geom_hline(yintercept = convergence_val, linetype = "dashed", color = "#BB00FF")+
  
  # add asymptotic value 
  annotate("text",x = 8.5,y = convergence_val, label = round(convergence_val),color = "#BB00FF",hjust = 0,vjust=2)+
  theme_classic()

p


# make output image +++++++++++++++++++++++++++++++++++
outfile <- paste(infile,"_boxplot_curve",set,sep='')
outpng <- paste(outfile,".png",sep='')
high <- 6
wide  <- 5.5
path <- getwd()
ggsave(p, filename=outpng, path = path, width=wide,height=high,dpi=400)







