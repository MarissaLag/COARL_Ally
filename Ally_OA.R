#ANOVA for Ally

##t-tests of SAR bacterial isolates from mb2021 - control and high salinity treatments
##sorce: https://www.statology.org/three-way-anova-in-r/

install.packages("datarium")
install.packages("rstatix")
install.packages("ggpubr")
install.packages("tidyverse")

library("datarium")
library(tidyverse)
library(ggpubr)
library(rstatix)

Data = X60_fert_Phenotype_Data_CSV

#Assumptions test ----
#source: https://www.datanovia.com/en/lessons/anova-in-r/#assumptions

ggboxplot(Data, x = "Family_Number", y = "Size")

#check for outliers ----
outliers <- Data %>% 
  group_by(Family_Number) %>%
  identify_outliers(Size)

View(outliers)

#outliers present in Family 10, 22, 49, 51, 77, 94, 95, 98 and more - most are just inputting errors

#To see size values
size_values <- Data$Size[Data$Family_Number == 10]
print(size_values)

#Remove outliers by manually checking (most were typos)


#check for distribution

# Build the linear model
model  <- lm(Size ~ Family_Number, data = Data)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

#check All the points fall approximately along the reference line, for each cell.
#looks like there are outliers resulting in non-normal dist - try square rooting dependent variable (sizes)

Data$Size <- sqrt(Data$Size)
View(Data)
model  <- lm(Size ~ Family_Number, data = Data)
ggqqplot(residuals(model))

#square-root did not help - try log
Data$Size <- log(Data$Size)
View(Data)
model  <- lm(Size ~ Family_Number, data = Data)
ggqqplot(residuals(model))

#Still did not help... May not be normally dist'd?

# Compute Shapiro-Wilk test of normality
#Note that, if your sample size is greater than 50, the normal QQ plot is preferred because at larger sample sizes the Shapiro-Wilk test becomes very sensitive even to a minor deviation from normality.
shapiro_test(residuals(model))

#if larger sample size:
ggqqplot(Data, "Size", facet.by = "Treatment")

#get means of groups ----
Data %>% sample_n_by(Treatment, size = 2)

Data %>% sample_n_by(Treatment, size = 3)

#trying z-score

Data$Size[Data$Size == ""] <- NA

# Remove NA values from the "Size" column
df_clean <- na.omit(Data$Size)

# Calculate the mean of the cleaned "Size" data
mean_size <- mean(df_clean)
#mean = 73.63
sd_size <- sd(df_clean)
#sd = 4.09

# Perform z-score transformation
Data$Size_zscore <- (Data$Size - mean_size) / sd_size
View(Data)
# Print the transformed dataframe
print(df)

#test normality with z-scores
model  <- lm(Size_zscore ~ Family_Number, data = Data)
ggqqplot(residuals(model))

#Data is not normal and transformations aren't helping
#Online suggests some deviations from line (especially at line ends) may not mean data is non-normal especially with a very large sample size
#since z-score is not helping, use untransformed data

#next assumption: Homogeneity of variance assumption
model  <- lm(Size_zscore ~ Family_Number, data = Data)
ggqqplot(residuals(model))

model  <- lm(Size ~ Family_Number, data = Data)
ggqqplot(residuals(model))
plot(model, 1)

#no evident relationships between residuals and fitted values (the mean of each groups), which is good. So, we can assume the homogeneity of variances.
#Can also use levene's test to test for homo variance

Data %>% levene_test(Size_zscore ~ Family_text)

#Very signficance, which indicates reject null hypothesis of homo variance - likely due to large sample size.
#May be able to ignore as residuals plot looks okay

#test for significance ----

#perform three-way ANOVA
model <- aov(Size ~ Treatment * Family_text, data=Data)
model

model <- lm(Size ~ Treatment * Family_text,
           data = Data
)

Anova(model)

#results

Anova(model)
Note: model has aliased coefficients
sums of squares computed by model comparison
Anova Table (Type II tests)

Response: Size
Sum Sq   Df  F value    Pr(>F)    
Treatment              18202    1 1871.774 < 2.2e-16 ***
  Family_text            23533  127   19.055 < 2.2e-16 ***
  Treatment:Family_text  13208  124   10.953 < 2.2e-16 ***
  Residuals              67868 6979                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Post-Hoc text - pairwise with p-adjustment

pairwise.t.test(Data$Size, Data$Treatment,
                p.adjust.method = "BH"
)

#pairwise wont test differences between families - too many families
pairwise <- pairwise.t.test(Data$Size, Data$Family_text,
                p.adjust.method = "BH"
)

View(pairwise)

mod <- aov(Size ~ Treatment,
           data = Data
)

TukeyHSD(mod,
         which = "Family_text"
)

par(mar = c(4.1, 13.5, 4.1, 2.1))

# create confidence interval for each comparison

ggplot(Data, aes(Family_text, Size, colour = Treatment)) + geom_boxplot()+ scale_color_manual(values = c("black","grey")) + theme_bw() + ylab("Number of mites")
