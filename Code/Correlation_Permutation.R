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

#To rely on a permutation null distribution, we need a strong formulation
#of H0 that states that the distributions are equal.

ggplot(melt(armpit[1:8]), aes(x=variable, y=value)) + 
  geom_boxplot()

ggplot(melt(log(armpit[1:8]+1)), aes(x=variable, y=value)) + 
  geom_boxplot()

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

#not all of the distributions are the same, therefore we cannot randomize the values, but permutate the order of the variables within each bacteria.
#test the assumptions
cor_permutation = c()
for (i in 1:10000){
  x = sample(armpit$Corynebacterium.1)
  y = sample(armpit$Staphylococcus.1)
  cor_permutation[i] = cor(x,y, method = "kendall")
}

cor_armpit = cor(armpit$Corynebacterium.1,armpit$Corynebacterium.2, method = "kendall")
cor_armpit
mean(abs(cor_permutation) > abs(cor_armpit))
mean(cor_permutation)
sd(cor_permutation)
hist(cor_permutation)
qqplot(y= cor_permutation, x =rnorm(10000))
#normality of the permutation correlation results can be assumed (needed for the confidence intervals)
#this can also be assumed even when the two sequences that we compare do not come form the same distribution

p_values <- data.frame(matrix(nrow = dim(armpit)[2]-3, ncol = dim(armpit)[2] - 3))
rownames(p_values) <- colnames(armpit[1:8])
colnames(p_values) <- colnames(armpit[1:8])

cor_values <- data.frame(matrix(nrow = dim(armpit)[2]-3, ncol = dim(armpit)[2] - 3))
rownames(cor_values) <- colnames(armpit[1:8])
colnames(cor_values) <- colnames(armpit[1:8])

#not correct i think?
confidence_intervals <- data.frame(matrix(nrow = dim(armpit)[2]-3, ncol = dim(armpit)[2] - 3))
rownames(confidence_intervals) <- colnames(armpit[1:8])
colnames(confidence_intervals) <- colnames(armpit[1:8])


for (x_index in 1:8){
  for (y_index in 1:8){
    x <- armpit[[x_index]]
    y <- armpit[[y_index]]
    cor_permutation = c()
    for (i in 1:2000){
      x = sample(x)
      y = sample(y)
      cor_permutation[i] = cor(x,y, method = "kendall")
    }
    print(colnames(armpit)[[x_index]])
    print(colnames(armpit)[[y_index]])
    print(cor(armpit[[x_index]], armpit[[y_index]]))
    print(mean( abs(cor_permutation) > abs(cor(armpit[[x_index]], armpit[[y_index]]))))
    p_values[x_index,y_index] <- mean( abs(cor_permutation) > abs(cor(armpit[[x_index]], armpit[[y_index]])))
    confidence_intervals[x_index,y_index] <- 2.3*sd(cor_permutation)/sqrt(36)
    cor_values[x_index,y_index]  <- cor(armpit[[x_index]], armpit[[y_index]])
  }
}


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

as.matrix(get_upper_tri(p_values))
as.matrix(get_upper_tri(cor_values))

upper_tri <- as.matrix(get_upper_tri(p_values))
dim(upper_tri) <- NULL
upper_tri = upper_tri[!is.na(upper_tri)]
upper_tri_sorted <- sort(upper_tri)
plot(y =upper_tri_sorted, x = 1:36 )

#Can't use benjamini hochberg as sequences are not independent and are negatively associated
#benjamini yekutieli
arbitrary_dependece <- 0
for (i in 1:36){
  arbitrary_dependece <- arbitrary_dependece + 1/i
}


k <- 0
for (value in upper_tri_sorted){
  if (value > (k/(36*arbitrary_dependece))*0.05){
    print(((k-1)/(36*arbitrary_dependece))*0.05)
    break
  }
  k <- k + 1
}


#holm bonferoni, less conservative
k <- 0
for (value in upper_tri_sorted){
  if (value > 0.05/(36-k+1)){
    print(0.05/(36-k))
    break
  }
  k <- k + 1
}


#confidence intervals
#the permutation distribution is normal so we can just use x_mean - t*s/sqrt(n)


