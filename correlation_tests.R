###Project Principle statistics###
setwd("D:/Ugent/MaNaMa/principal_statistics/project")
library(reshape2)
library(ggplot2)

armpit <- read.table("armpit.txt", sep = " ")
head(armpit)
summary(armpit)
#first 8 columns show the relative abundances, age, BMI and gender
armpit <- armpit[complete.cases(armpit),]

armpit$Gender <- factor(trimws(armpit[["Gender"]], which = c("both", "left", "right"), whitespace = "[ \t\r\n]"))
armpit$BMI <- factor(armpit$BMI)
#check if all relative abundances sum op to 100
rowSums(armpit[1:8]) >= 99.99
summary(armpit)

head(armpit)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cor_pvalues <- data.frame(matrix(nrow = dim(armpit)[2]-2, ncol = dim(armpit)[2] - 2))
rownames(cor_pvalues) <- colnames(armpit[1:9])
colnames(cor_pvalues) <- colnames(armpit[1:9])
for (x_index in 1:9){
  for (y_index in 1:9){
    x <- armpit[[x_index]]
    y <- armpit[[y_index]]
    cor_pvalues[x_index,y_index] <- cor.test(x,y, method = "spearman")$p.value
  }
}
cor_pvalues <- round(cor_pvalues,2)

cormat <- data.frame(matrix(nrow = dim(armpit)[2]-2, ncol = dim(armpit)[2] - 2))
rownames(cormat) <- colnames(armpit[1:9])
colnames(cormat) <- colnames(armpit[1:9])
for (x_index in 1:9){
  for (y_index in 1:9){
    x <- armpit[[x_index]]
    y <- armpit[[y_index]]
    cormat[x_index,y_index] <- cor.test(x,y, method = "spearman")$estimate
  }
}
cormat <- round(cormat,2)
cormat

#cormat <- reorder_cormat(cormat)
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
                                   size = 12, hjust = 1))+
  coord_fixed()

ggheatmap + 
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

Species_percentage <- data.frame(matrix(nrow = 39, ncol = 3))
rownames(Species_percentage) <- rownames(armpit)
colnames(Species_percentage) <- c("Corynebacterium","Staphylococus", "Age")
Species_percentage$Corynebacterium  <- rowSums(armpit[1:4])/100
Species_percentage$Staphylococus  <- rowSums(armpit[5:8])/100
Species_percentage$Age <- armpit$Age

Species_mean <- data.frame(matrix(nrow = 39, ncol = 3))
rownames(Species_mean) <- rownames(armpit)
colnames(Species_mean) <- c("Corynebacterium","Staphylococus", "Age")
Species_mean$Corynebacterium  <- rowMeans(armpit[1:4])/100
Species_mean$Staphylococus  <- rowMeans(armpit[5:8])/100
Species_mean$Age <- armpit$Age

#COrrelation between staph and coryne bact is obviously negative, if one increase the other decreases. Seems that there is a significant
#correlation between species and age. However this could be dominated by only one genus, try to code removal of geni and look at the effect
#on the correlation
#Means and sum are practically the same, just testing stuff out
cor.test(Species_percentage[["Corynebacterium"]], Species_percentage[["Staphylococus"]], method = "spearman" )
cor.test(Species_percentage[["Corynebacterium"]], Species_percentage[["Age"]], method = "spearman" )
cor.test(Species_percentage[["Staphylococus"]], Species_percentage[["Age"]], method = "spearman" )

cor.test(Species_mean[["Corynebacterium"]], Species_mean[["Staphylococus"]], method = "spearman" )
cor.test(Species_mean[["Corynebacterium"]], Species_mean[["Age"]], method = "spearman" )
cor.test(Species_mean[["Staphylococus"]], Species_mean[["Age"]], method = "spearman" )

