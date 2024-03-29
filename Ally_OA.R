#ANOVA for Ally

##t-tests of SAR bacterial isolates from mb2021 - control and high salinity treatments
##sorce: https://www.statology.org/three-way-anova-in-r/

install.packages("datarium")
install.packages("rstatix")
install.packages("ggpubr")
install.packages("tidyverse")
install.packages("ggResidpanel")
install.packages("DHARMa")
install.packages("lme4")
install.packages("fitdistrplus")

library("datarium")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggResidpanel)
library(DHARMa)
library(lme4)
library(fitdistrplus)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

Data = X60_fert_Phenotype_Data_CSV
#Data = Sizes_nomissing #tested with no missing values - made no difference in normality checks

#Assumptions test ----
#source: https://www.datanovia.com/en/lessons/anova-in-r/#assumptions

ggboxplot(Data, x = "Family_Number", y = "Size")

#check for outliers ----
outliers <- Data %>% 
  group_by(Family_text) %>%
  identify_outliers(Size)

View(outliers)

#outliers present in Family 10, 22, 49, 51, 77, 94, 95, 98 and more - most are just inputting errors

#To see size values
size_values <- Data$Size[Data$Family_Number == 10]
print(size_values)

#Remove outliers by manually checking (most were typos)


#check for distribution

#look at histogram of data dist
hist(Data$Size)

#As expected, size data looks normal but with right-skewdness

# Build the linear model
model  <- lm(Size ~ Family_text + Treatment, data = Data)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

#another way to look at residuals distribution

resid_panel(model)

#Can compare residuals plot when using a glm instead of a lm - let's go with a gamma dist'd (for right skewdness)
model_gamma <- glm(Size ~ Family_text + Treatment, data = Data, family = Gamma(link = "log"))
# Create a QQ plot of residuals
resid_panel(model_gamma)

#Try DHARMa package ----
#normal dist
model_2 <- glm(Size ~ Family_text + Treatment, data = Data, family = Gamma(link="log"))
resids <- simulateResiduals(model_2)
plot(resids)
#DHARMa package is giving odd results...

#gamma dist
model_2 <- glm(Size ~ Family_text + Treatment, data = Data, family = Gamma(link = "log"))
resids <- simulateResiduals(model_2)
plot(resids)

model_2 <- glm(Size ~ Family_text + Treatment, data = Data, family = poisson(link = "log"))
resids <- simulateResiduals(model_2)
plot(resids)


#Can also try data transformations
#try square rooting dependent variable (sizes)

Data$Size <- sqrt(Data$Size)
View(Data)
model_sqrt  <- lm(Size ~ Family_text + Treatment, data = Data)
resid_panel(model_sqrt)

#square-root did not help - try log
Data$Size <- log(Data$Size)
View(Data)
model  <- lm(Size ~ Family_Number, data = Data)
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
#Note that, if your sample size is greater than 50, the normal QQ plot is preferred because at larger sample sizes the Shapiro-Wilk test becomes very sensitive even to a minor deviation from normality.
shapiro_test(residuals(model))



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
#Can also use levene's test to test for homo variance when sample sizes are smaller (so not here)

Data %>% levene_test(Size_zscore ~ Family_text)

#Summary- assuming data passes all assumptions for ANOVA (see residual plot, ignoring DHARMa results)
#Run ANOVA ----

#test for significance ----

Data$Family_Number <- as.character(Data$Family_Number)

model <- lm(Size ~ Treatment * Family_Number,
           data = Data
)

Anova(model)

summary(model)

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



#Post-Hoc text - pairwise with p-adjustment ----

#pairwise wont test differences between families - too many families so giving NA 
pairwise <- pairwise.t.test(Data$Size, Data$Family_Number,
                p.adjust.method = "BH"
)

View(pairwise)

#summary: Need different method for post-hoc test for families -may need to group families into resilient or susceptible categories beforehand

#size plots ----

Data$Family_Number <- as.character(Data$Family_Number)

#NOTCH = If FALSE (default) make a standard box plot. If TRUE, make a notched box plot. Notches are used to compare groups; if the notches of two boxes do not overlap, this suggests that the medians are significantly different.

ggplot(Data, aes(x = Family_Number, y = Size, color = Treatment)) + 
  geom_boxplot(outlier.color = "grey", notch = TRUE) +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  ylab("Size (microns)") +
  xlab("Family") + 
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#remove families with no values

filtered_data <- Data %>%
  filter(!(Family_Number %in% c("133", "134", "135", "136", "137", "138")))


ggplot(filtered_data, aes(x = Family_Number, y = Size, color = Treatment)) + 
  geom_boxplot(outlier.color = "grey", notch = TRUE) +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8)) +
  ylab("Size (microns)") +
  xlab("Family") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#x-axis orders not correct

family_levels <- unique(Data$Family_Number)

#order as appears on data sheet (don't follow base R rules)
Data$Family_Number <- factor(Data$Family_Number, levels = family_levels)

filtered_data <- Data %>%
  filter(!(Family_Number %in% c("89", "90", "133", "134", "135", "136", "137", "138")))


ggplot(filtered_data, aes(x = Family_Number, y = Size, color = Treatment)) + 
  geom_boxplot(outlier.color = "grey", notch = TRUE) +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8)) +
  ylab("Size (microns)") +
  xlab("Family") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Look at average size of each family in each treatment

avg_data <- Data %>%
  group_by(Family_Number, Treatment) %>%
  summarise(
    Avg_Size = mean(Size),
    SD_Size = sd(Size),
    .groups = 'drop'
  )

View(avg_data)

#order families with largest-smallest size within DW treatment
ordered_avg_data_DW <- avg_data %>%
  filter(Treatment == "DW") %>%
  arrange(desc(Avg_Size))

View(ordered_avg_data_DW)

top_ten_largest <- ordered_avg_data_DW %>%
  slice(1:10)

top_ten_smallest <- ordered_avg_data_DW %>%
  slice(n() - 9:n())

#Top ten largest families in DW = 71, 1, 94, 45, 73, 140, 95, 199, 51, 122 (in order)
#Smallest 10 families in DW = 7, 6, 16, 4, 36, 102, 10, 31, 110, 9, 11 (in order from largest - smallest)

#order families with largest-smallest size within SW treatment
ordered_avg_data_SW <- avg_data %>%
  filter(Treatment == "SW") %>%
  arrange(desc(Avg_Size))
View(ordered_avg_data_SW)

top_ten_largest <- ordered_avg_data_SW %>%
  slice(1:10)

top_ten_smallest <- ordered_avg_data_SW %>%
  slice(n() - 9:n())

#Try plotting from largest to smallest in DW

# Ensure that the "Family_text" column is a factor with the levels set based on "Avg_Size"

avg_data <- Data %>%
  group_by(Family_Number, Treatment) %>%
  summarise(
    Avg_Size = mean(Size),
    SD_Size = sd(Size),
    .groups = 'drop'
  )

ordered_avg_data_DW <- avg_data %>%
  filter(Treatment == "DW") %>%
  arrange(desc(Avg_Size))

ordered_avg_data_DW$Family_Number <- factor(ordered_avg_data_DW$Family_Number, 
                                          levels = ordered_avg_data_DW$Family_Number)

#finalized largest-smallest SW and DW plots ----
# Plot the data with the families ordered by average size
#remove NA (larvae that did not develop in OA)

ordered_avg_data_DW <- na.omit(ordered_avg_data_DW)
View(ordered_avg_data_DW)

#no colour
ggplot(ordered_avg_data_DW, aes(x = Family_Number, y = Avg_Size)) + 
  geom_boxplot() +
  geom_errorbar(data = ordered_avg_data_DW, aes(x = Family_Number, ymin = Avg_Size - SD_Size, ymax = Avg_Size + SD_Size), width = 0.4, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Average Size") +
  xlab("Family")

#add colour for treatment (red = DW)
plot2 <- ggplot(ordered_avg_data_DW, aes(x = Family_Number, y = Avg_Size)) + 
  geom_boxplot(fill = "white", color = "red", size = 1) + 
  geom_errorbar(data = ordered_avg_data_DW, aes(x = Family_Number, ymin = Avg_Size - SD_Size, ymax = Avg_Size + SD_Size), width = 0.4, position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("Average Size") +
  xlab("Family") +
  ggtitle("Deep water")

#Do the same for SW treatment

ordered_avg_data_SW$Family_Number <- factor(ordered_avg_data_SW$Family_Number, 
                                            levels = ordered_avg_data_SW$Family_Number)

ordered_avg_data_SW <- na.omit(ordered_avg_data_SW)
View(ordered_avg_data_SW)

#no colour
ggplot(ordered_avg_data_SW, aes(x = Family_Number, y = Avg_Size)) + 
  geom_boxplot() +
  geom_errorbar(data = ordered_avg_data_DW, aes(x = Family_Number, ymin = Avg_Size - SD_Size, ymax = Avg_Size + SD_Size), width = 0.4, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Average Size") +
  xlab("Family")

#add colour for treatment (blue = SW) and remove gridlines
plot1 <- ggplot(ordered_avg_data_SW, aes(x = Family_Number, y = Avg_Size)) + 
  geom_boxplot(fill = "white", color = "blue", size = 1) + 
  geom_errorbar(data = ordered_avg_data_SW, aes(x = Family_Number, ymin = Avg_Size - SD_Size, ymax = Avg_Size + SD_Size), width = 0.4, position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("Average Size") +
  xlab("") +
  ggtitle("Surface water")

#arrange plots into a grid ----

library(gridExtra)

# Assuming you have two ggplot objects named plot1 and plot2
combined_plots <- grid.arrange(plot1, plot2, ncol = 1)

#PROBLEM TO FIX - can't figure out how to add a colour legend...


#overlapping density plots ----

p2 <- ggplot(data=Data, aes(x=Size, group=Treatment, fill=Treatment)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()
p2

p2 <- ggplot(data=Data, aes(x=Size, group=Treatment, fill=Treatment)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum() +
  scale_fill_manual(labels = c("Ambient", "OA"), values=c("blue", "red"))

p2



#Grouping families into R, MR, or S? ----

#Ideas:
#Based on data set - could make arbitrary size cut offs in OA - e.g., Resilient = 73-76um, moderate-resilient = 70-72um, and susceptible = 69 and below
#however not sure how to incorporate standard deviations....or even if we need to.
#Alternatively could use Z-score values to define R, MR, or S?
#Could try to correlate scoring stats to size - e.g., if sizes below 70 correlate with having more deformed larvae?
#Also should remove families that did not fertilize well in SW


#Below, Grouping families based on R, MR, or S to use for size visuals based on avg family size in OA seawater


#first average DW sizes for all families

avg_data_DW <- Data %>%
  filter(Treatment == "DW") %>%
  group_by(Family_text) %>%
  summarise(
    Avg_Size = mean(Size),
    SD_Size = sd(Size),
    .groups = 'drop'
  )


avg_data_classified_DW <- avg_data_DW %>%
  mutate(Size_Category = case_when(
    Avg_Size > 74 ~ "R",
    Avg_Size >= 70 & Avg_Size <= 73.99 ~ "MR",
    TRUE ~ "S"
  ))


ggplot(data = avg_data_classified_DW, aes(x = factor(Size_Category, levels = c("R", "MR", "S")), y = Avg_Size, fill = Size_Category)) + 
  geom_violin() +
  scale_fill_brewer(palette = "Spectral") +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  ylab("Average Size (microns)") +
  xlab("") +
  guides(fill = FALSE)

#make list of families classified as R, MR, or S

#R
result_R <- avg_data_classified_DW %>%
  filter(Size_Category == "R") %>%
  dplyr::select(Family_text) %>%
  distinct()
View(result_R)

#MR
result_MR <- avg_data_classified_DW %>%
  filter(Size_Category == "MR") %>%
  dplyr::select(Family_text) %>%
  distinct()
View(result_MR)

#S
result_S <- avg_data_classified_DW %>%
  filter(Size_Category == "S") %>%
  dplyr::select(Family_text) %>%
  distinct()
View(result_S)



#Combine this infor with Data object (i.e., denote families as R, MR, or S)

Data <- Data %>%
  mutate(Size_Category = case_when(
    Family_text %in% result_R$Family_text ~ "R",
    Family_text %in% result_MR$Family_text ~ "MR",
    Family_text %in% result_S$Family_text ~ "S",
    TRUE ~ "Unknown"
  ))

View(Data)



avg_data_DW <- Data %>%
  filter(Treatment == "DW") %>%
  group_by(Family_text) %>%
  summarise(
    Avg_Size = mean(Size),
    SD_Size = sd(Size),
    .groups = 'drop'
  )


#Size category density plots ----

plot1 <- ggplot(data = subset(Data, Treatment == "DW"), aes(x = Size, fill = Size_Category)) +
  geom_density(alpha = 0.5) +
  labs(title = "DW Treatment", x = "Size", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(55, 90)

#to add blue border to denote treatment

plot1 <- ggplot(data = subset(Data, Treatment == "DW"), aes(x = Size, fill = Size_Category)) +
  geom_density(alpha = 0.5) +
  labs(title = "", x = "Size", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(55, 90) +
  theme(panel.border = element_rect(color = "red", fill = NA, linewidth = 3)) +
  guides(fill=FALSE)


plot2 <- ggplot(data = subset(Data, Treatment == "SW"), aes(x = Size, fill = Size_Category)) +
  geom_density(alpha = 0.5) +
  labs(title = "SW Treatment", x = "Size", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(55, 90)

#Add border colour to denote treatment

plot2 <- ggplot(data = subset(Data, Treatment == "SW"), aes(x = Size, fill = Size_Category)) +
  geom_density(alpha = 0.5) +
  labs(title = "", x = "Size", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(55, 90) +
  theme(panel.border = element_rect(color = "blue", fill = NA, linewidth = 3))

combined_plots <- grid.arrange(plot2, plot1, ncol = 1)

#add dotted line indicating peak size - code not working below 

p <- ggplot(data = subset(Data, Treatment == "DW"), aes(x = Size, fill = Size_Category)) +
  geom_density(alpha = 0.5) +
  labs(title = "DW Treatment", x = "Size", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Calculate peak points for each Size_category
# Remove missing values and then calculate the peak density
peaks <- Data %>%
  filter(!is.na(Size), Treatment == "DW") %>%
  group_by(Size_Category) %>%
  summarise(peak_density = density(Size)$x[which.max(density(Size)$y)], .groups = "drop")

#peaks: MR=73.1, R=76.1, and S=69.9


# Add annotated peak lines to the plot
p + geom_vline(xintercept = c(73.1, 76.1, 69.9), linetype = "dotted")


#NOT part of Ally's project (ignore below)
#OA selection event----
#With Ally's size data, mimic an OA selection event
#Will need to do more analysis to determine surviving size in OA (using scoring and Sizes) but for now, lets say 
#Larvae smaller than 73um will not survive (filter out)
#And look at new density plot

filtered_data <- Data %>%
  
  filter(Size > 73)

# Create a density plot for the filtered data
plot_filtered <- ggplot(data = subset(filtered_data, Treatment =="DW"), aes(x = Size)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Density Plot of Sizes < 75", x = "Size", y = "Density") +
  theme_minimal() +
  xlim(50, 90)

plot_filtered

plot <- ggplot(data = subset(Data, Treatment =="DW"), aes(x = Size)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Sizes < 75", x = "Size", y = "Density") +
  theme_minimal() +
  xlim(50, 90)


# Create a single overlay plot

combined_data <- rbind(subset(Data, Treatment == "DW" & Size < 75), subset(Data, Treatment == "DW"))

overlay_plot <- ggplot(combined_data, aes(x = Size, fill = Treatment, color = as.factor(Size < 75))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Sizes < 75 Overlayed", x = "Size", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("red", "red")) +
  xlim(50, 90)

overlay_plot

plot <- ggplot(data = subset(Data, Treatment == "DW"), aes(x = Size, fill = cut(Size, c(-Inf, 73, Inf)))) +
  +     geom_density(alpha = 0.5) +
  +     labs(title = "Density Plot with Size < 75 Colored Red", x = "Size", y = "Density") +
  +     theme_minimal() +
  +     scale_fill_manual(values = c("blue", "red")) +
  +     xlim(50, 90)
> 
  > plot


#make a dotted line at 73um (geom_vline fxn)

plot <- ggplot(data = subset(Data, Treatment == "DW"), aes(x = Size)) +
  +     geom_density(fill = "blue", alpha = 0.5) +
  +     labs(title = "Density Plot with Dotted Line at x = 73", x = "Size", y = "Density") +
  +     theme_minimal() +
  +     geom_vline(xintercept = 73, linetype = "dotted") +
  +     xlim(50, 90)

 plot


 #Lets try to assess mortality rate of a given size in OA
 #Issue = scoring does not relate to size values for each individual... so can't use scoring as 
 #direct measurement for survival...will have to go back into samples and score and take measurements of same larva
 
 #From Johnson 2022:
 #Daily mortality rates for each family were calculated from the negative of the slope of a linear regression 
 #relating the natural log of the proportion of fish remaining to time. 
 #Average size of larvae within each family was measured by digitally photographing a random sample of 10–20 larvae under a microscope

#could look up score counts for each family in DW and then look at the family's average size in DW?
 
 str(Data)
 Data$Family_Number <- as.character(Data$Family_Number)
 
 scores_DW <- Data %>%
   filter(Treatment == "DW") %>%
   group_by(Family_Number, Score) %>%
   summarise(number_of_scores = n())
 
 View(Data)
 View(scores_DW)
 
 # Remove NA and "S" values under Score
 scores_DW <- scores_DW %>%
   filter(!is.na(Score), Score != "S")
 
 
 ggplot(scores_DW, aes(x = Family_Number, y = number_of_scores, fill = Score)) +
   geom_bar(stat = "identity") +
   coord_flip() +
   labs(title = "Proportion of Different Scores for Each Family",
        x = "Family Number",
        y = "Proportion",
        fill = "Score") +
   theme(axis.text.y = element_text(size = 8))  # Change the font size as needed
 
 #reorder families with highest number of "G" scores
 
 ordered_family <- scores_DW %>%
   filter(Score == "G") %>%
   group_by(Family_Number, Score) %>%
   arrange(desc(number_of_scores))
 #Make G count as percentage of total count (total count = 30 observations)
 
 ordered_family$number_of_scores <- (ordered_family$number_of_scores / 30) * 100
 
 str(ordered_family)
 View(ordered_family)
 
 ggplot(ordered_family, aes(x = Family_Number, y = number_of_scores)) +
   geom_bar(stat = "identity") +
   labs(title = "Proportion of Different Scores for Each Family",
        x = "Family Number",
        y = "Number of G's") +
   theme(axis.text.y = element_text(size = 8))  # Change the font size as needed
 
#to order in order specified by ordered_family (highest to lowest number of G's)
 
 library(forcats)
 
 # Reorder the factor levels of Family_Number
 ordered_family$Family_Number <- fct_inorder(ordered_family$Family_Number)

 ggplot(data = ordered_family, aes(x = Family_Number, y = number_of_scores)) +
   geom_bar(stat = "identity") +
   labs(title = "Scores for Each Family",
        x = "Family Number",
        y = "Number of G's") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   ylim(0, 100) +
   theme(axis.text.x = element_text(size = 8))
 
#get avergage size in SW and DW
 
 avg_data <- Data %>%
   group_by(Treatment) %>%
   summarise(
     Avg_Size = mean(Size, na.rm = TRUE),
     SD_Size = sd(Size, na.rm = TRUE),
     .groups = 'drop'
   )
 
 View(avg_data)
 
#get % difference
 DW_avg <- avg_data$Avg_Size[avg_data$Treatment == "DW"]
 SW_avg <- avg_data$Avg_Size[avg_data$Treatment == "SW"]
 
 change <- SW_avg - DW_avg
 averages <- (SW_avg + DW_avg)/2
 
 percent_difference <- change/averages *100

print(percent_difference) #4.155% difference between DW and SW average size

#Based off Bitter et al., 2019, size differences will disappear
#What size "cut-off" will make the DW mean the same as the SW mean?

#Plot histogram of DW and SW

p2 <- ggplot(data=final_filtered_data, aes(x=Size, group=Treatment, fill=Treatment)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum() +
  scale_fill_manual(labels = c("Ambient", "OA"), values=c("blue", "red"))

p2

#filter values smaller than 72.32um and check the mean of DW
filtered_data <- dplyr::filter(Data, Size >= 73.5)
View(filtered_data)

filtered_data_DW <- dplyr::filter(Data, Size >= 72.32, Treatment == "DW")
filtered_data_SW <- dplyr::filter(Data, Treatment == "SW")

final_filtered_data <- rbind(filtered_data_DW, filtered_data_SW)

View(filtered_data_SW)

avg_data <- final_filtered_data %>%
  group_by(Treatment) %>%
  summarise(
    Avg_Size = mean(Size, na.rm = TRUE),
    SD_Size = sd(Size, na.rm = TRUE),
    .groups = 'drop'
  )

View(avg_data)


sample_data <- rnorm(10000, mean = 63.8)
mean <- mean(sample_data)


# Plot a histogram
hist(sample_data, breaks = 10, col = "lightblue", main = "", xlab = "Values")



#In order for DW and SW averages to equal the same, all sizes below 72.32 must be filtered out from DW
#Paper by Bitter et al., 2019 proposed all slow-growing larvae will diminish over time until mean avg size of 
#ambient and OA is the same by the end of larval rearing
#eventually, may want to test this empirically (e.g., rear larvae in OA conditions and measure avg size overtime and compare to control)



