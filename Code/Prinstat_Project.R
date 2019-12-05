install.packages('corrgram')
install.packages('psych')
library(psych)
library(corrgram)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Data clean and tidy

#1.1 read in data
armpit <- read.csv(file.choose(), header = TRUE, sep = ' ')

#1.2 view data
str(armpit)
summary(armpit)

#1.3 Remove observation with missing value
armpit_cleaned <- na.omit(armpit)

#1.4 Tidy variables in suitable structure
armpit_tidied <- within(armpit_cleaned,{
                        BMI <- factor(BMI)
                        Gender[Gender == "M "] <- "M"
                        Gender[Gender == "F "] <- "F"
                        Gender <- droplevels(Gender)
                        })

#1.5 Create new variables
new_armpit <- transform(armpit_tidied, Agecat = Age)
new_armpit <- within(new_armpit,{
                    Agecat[Age <= 26] <- "Junior"
                    Agecat[Age >= 27] <- "Adult"
                    Agecat[Age >= 51] <- "Senior"
                    })
new_armpit[1:8] <- new_armpit[1:8] /  rowSums(armpit_tidied[, 1:8])
str(new_armpit)

#2. Descriptive analysis

#2.1 Univariate fAnalysis (age, BMI, Gender, organisms)

myvars <- c('Corynebacterium.1', 'Corynebacterium.2', 'Corynebacterium.3', 'Corynebacterium.4', 'Staphylococcus.1', 'Staphylococcus.2', 'Staphylococcus.3', 'Staphylococcus.4')
describe(new_armpit[myvars])

par(mfrow=c(2,4))
for(i in (1:8)){
  name <- names(new_armpit)[i]
  boxplot(new_armpit[i], main = name)
}

par(mfrow=c(1,3))
for(i in (10:12)){
  count <- table(new_armpit[,i])
  name <- names(new_armpit)[i]
  barplot(count, main = name)
}

#2.2 Bivariate Analysis

#2.2.1 Creat a long format dataframe for analysis
new_armpit <- transform(new_armpit, organisms = names(new_armpit[1]))
long_armpit <- data.frame(new_armpit[1], new_armpit[9:13])
names(long_armpit)[1] <- 'Relative Aboundance'
for (i in 2:8){
  tem_armpit <- transform(new_armpit, organisms = names(new_armpit[i]))
  sub_armpit <- data.frame(tem_armpit[i], tem_armpit[9:13])
  names(sub_armpit)[1] <- 'Relative Aboundance'
  long_armpit <- rbind(long_armpit, sub_armpit)
}
long_armpit$Agecat <- factor(long_armpit$Agecat, ordered = TRUE, levels = c('Junior', 'Adult', 'Senior'))
long_armpit$organisms <- factor(long_armpit$organisms)

#2.2.2 Relative aboundance relate to agecat
ggplot(long_armpit,
       aes(x = organisms, y = `Relative Aboundance`)) + 
       geom_bar(stat = 'identity') +
       labs(title = "Raletive Aboundanc of 8 Organisms")

#2.2.3 Relative aboundance relate to agecat

describeBy(new_armpit[myvars], list(Agecat=new_armpit$Agecat))

ggplot(long_armpit,
       aes(x = organisms, y = `Relative Aboundance`, fill = Agecat)) + 
       geom_bar(stat = 'identity', position = 'fill') +
       labs(title = "Raletive Aboundanc relate to Agecat")

#2.2.4 Relative aboundance relate to BMI

describeBy(new_armpit[myvars], list(BMI=new_armpit$BMI))

ggplot(long_armpit,
       aes(x = organisms, y = `Relative Aboundance`, fill = BMI)) + 
       geom_bar(stat = 'identity', position = 'fill') +
       labs(title = "Raletive Aboundanc relate to BMI")

#2.2.5 Relative aboundance relate to gender

describeBy(new_armpit[myvars], list(Gender=new_armpit$Gender))

ggplot(long_armpit,
       aes(x = organisms, y = `Relative Aboundance`, fill = Gender)) + 
       geom_bar(stat = 'identity', position = 'fill') +
       labs(title = "Raletive Aboundanc relate to Gender")

#2.3 Multivariate frequencies, supporting research question



#3.1 correlations in species
species <- new_armpit[1:8]
name <- c("S1", "S2", "S3", "S4", "C1", "C2", "C3", "C4")
cor_matrix <- matrix(cor(species), nrow = 8, ncol = 8,
                     dimnames = list(name, name))
cor_data <- data.frame(cor_matrix)
corrgram(cor_data, order = TRUE,
         lower.panel = panel.shade, upper.panel = panel.pie, text.panel = panel.txt,
         main = "Corrgram of species intercorrelations")


