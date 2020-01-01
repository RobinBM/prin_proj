###Project Principle statistics###
setwd("D:/Ugent/MaNaMa/principal_statistics/project")
library(reshape2)
library(ggplot2)
library(ggpubr)

armpit <- read.table("armpit.txt", sep = " ")
head(armpit)
summary(armpit)
#first 8 columns show the relative abundances, age, BMI and gender
armpit <- armpit[complete.cases(armpit),]

armpit$Gender <- factor(trimws(armpit[["Gender"]], which = c("both", "left", "right"), whitespace = "[ \t\r\n]"))
armpit$BMI <- factor(armpit$BMI)
#check if all relative abundances sum op to 100
rowSums(armpit[1:8]) >= 99.99

#Check visually with boxplots and numerically with the kolmogorov-smirnov test 
#to see wether the distributions from two species are the same.

armpit_bxp1 <- melt(armpit[1:8])
bxp1 <- ggplot(armpit_bxp1, aes(x=variable, y=value)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.title.x = element_blank())+
  theme_minimal() +
  ylab("Relative abundances (%)") + 
  geom_boxplot()

armpit_bxp2 <- melt(log(armpit[1:8]+1))
bxp2 <- ggplot(armpit_bxp2, aes(x=variable, y=value)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.title.x = element_blank())+
  theme_minimal() +
  ylab("Log-transformed relative abundances") + 
  geom_boxplot()

figure <- ggarrange(bxp1, bxp2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
annotate_figure(figure, top = "Relative abundance boxplots for the 8 different species")


ks_pvalues <- data.frame(matrix(nrow = dim(armpit)[2]-3, ncol = dim(armpit)[2] - 3))
rownames(ks_pvalues) <- colnames(armpit[1:8])
colnames(ks_pvalues) <- colnames(armpit[1:8])
for (x_index in 1:8){
  for (y_index in 1:8){
    x <- armpit[[x_index]]
    y <- armpit[[y_index]]
    ks_pvalues[x_index,y_index] <- ks.test(x,y)$p.value
  }
}

round(ks_pvalues,3)

#not all of the distributions are the same, therefore we permutate the order of the variables within each species

#Check correlation distribution
cor_permutation = c()
for (i in 1:10000){
  x = sample(armpit$Corynebacterium.1)
  y = sample(armpit$Staphylococcus.1)
  cor_permutation[i] = cor(x,y, method = "spearman")
}

qq <- ggplot(data.frame(c = cor_permutation), aes(sample = c)) +
  stat_qq() + 
  stat_qq_line() +
  theme_minimal() +
  ylab("Frequency") + 
  xlab("Spearman Correlation")

histo <- ggplot(data.frame(c = cor_permutation), aes(x = c)) +
  geom_histogram(color="black", fill="white") + 
  theme_minimal() +
  ylab("Spearman Correlation") + 
  xlab("Theoretical Quantiles Normal Distribution")

figure <- ggarrange(histo, qq,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure
annotate_figure(figure, top = "Histogram (A) and QQ-plot (B) for a permutation result \n between Corynebacterium 1 and Staphylococcus 1")

# Get upper triangle
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# Get lower triangle
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

library(psych)
#Create figure
a <- corr.test(armpit[1:8], method = "spearman", adjust = "BY")
#ccorrections possible: holm, hochberg, hommel, bonferroni, bh, by, fdr, none
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3045855/
cor_pvalues <- round(a$p,2)
cormat <- round(a$r,2)

#without adjusted p-values
lower_tri_p_values <- t(as.matrix(get_lower_tri(cor_pvalues)))
lower_tri <- t(as.matrix(get_lower_tri(cormat)))
# Melt the correlation matrix
melted_cormat <- melt(lower_tri, na.rm = TRUE)
melted_cormat_p_values <- melt(lower_tri_p_values, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()

not_adjusted <- ggheatmap + 
  geom_text( aes(melted_cormat_p_values$Var2, melted_cormat_p_values$Var1, label = melted_cormat_p_values$value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


#With adjusted p-values
upper_tri_p_values <- as.matrix(get_upper_tri(cor_pvalues))
upper_tri <- as.matrix(get_upper_tri(cormat))
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat_p_values <- melt(upper_tri_p_values, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()

adjusted <- ggheatmap + 
  geom_text( aes(melted_cormat_p_values$Var2, melted_cormat_p_values$Var1, label = melted_cormat_p_values$value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

figure <- ggarrange(not_adjusted, adjusted,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
annotate_figure(figure, top = "Correlation test results \n A: p-values not adjusted for multiple testing \n B: p-values adjusted with the Benjamini-Yekeutili method")

color_CI <- data.frame(t(a$ci.adj))
color_CI <- (colSums(color_CI > 0) == 2) | (colSums(color_CI < 0) == 2)
color_CI[color_CI == TRUE] <- "Red"
color_CI[color_CI == FALSE] <- "White"

rownames(a$ci.adj) <- rownames(a$ci)
par(mar=c(2,8,2,2))
boxplot(t(a$ci.adj), col = color_CI, horizontal = TRUE, main = "Adjusted confidence intervals", las = 1)
abline(v = 0)

