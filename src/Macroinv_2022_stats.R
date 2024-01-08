
####################################################### Thesis Paper 1 - Macroinvertebrates_Statistics #####################################################

library(dplyr)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(AICcmodavg)
library(Hmisc)
library(dendextend)
library(RColorBrewer)
library(usdm)
library(ggplot2)
library(visreg)
# library(car) # It changes vif(). Use as "car::" when needed. 
library(DFIT)
library(forcats) # to modify gray facet titles in facet_wrap(). 
# install.packages("r2glmm")
library(r2glmm)
library(MuMIn)
library(gridExtra)

##############  1. Preparing base data frames #################

#### 1.1. Working physicochemical data ####

physchem_2022 <- read.csv("data/Other_vars/Other_factors_CERESTRES_2022.csv", fileEncoding="latin1", na.strings=c("","NA"))
physchem_2022$Sampling_date <- as.Date(physchem_2022$Sampling_date)

# To have a single value for each physicochemical variable representing the conditions in which organisms were developed within plots, the average of the three previous
# physicochemical measurements to the biodiversity samplings is calculated:

## Assigning Sampling, Rep and Treat:
Sampling <- c("1", "2", "3", "4")
Treat <- c("AWD", "MSD", "CON", "MSD", "AWD", "CON", "MSD", "CON", "AWD", "AWD", "MSD", "CON", "MSD", "AWD", "CON")
Rep <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5")

## Subsetting physchem_2022, only considering three Sampling dates previous to each biodiversity samplings.
## Sampling column indicating to which biodiversity sampling each physicochemical sampling is assigned.
## Biodiversity sampling date 1: "2022-06-10" -> 3 previous dates: "2022-06-09", "2022-06-02", "2022-05-26"
## Biodiversity sampling date 2: "2022-07-15" -> 3 previous dates: "2022-07-14", "2022-07-07", "2022-07-05"
## Biodiversity sampling date 3: "2022-08-02" -> 3 previous dates: "2022-07-29", "2022-07-27", "2022-07-21"
## Biodiversity sampling date 4: "2022-08-31" -> 3 previous dates: "2022-07-29", "2022-07-27", "2022-07-21"

prev.dates <- c("2022-06-09", "2022-06-02", "2022-05-26", "2022-07-14", "2022-07-07", "2022-07-05", "2022-07-29", "2022-07-27", "2022-07-21", "2022-07-29", "2022-07-27", "2022-07-21")
sampling.prev <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4")
PrevDate_Samp <- data.frame(prev.dates, sampling.prev) # Creates data frame with previous dates to average
PrevDate_Samp$prev.dates <- as.Date(PrevDate_Samp$prev.dates)
physchem_avg_2022 <- merge(physchem_2022, PrevDate_Samp, by.x="Sampling_date", by.y="prev.dates") # Assigning sampling.prev values
colnames(physchem_avg_2022)[colnames(physchem_avg_2022) == "sampling.prev"] <- "Sampling" # Replacing column name "sampling.prev" for "Sampling"

## Creating sampdateID and averaging the three previous dates:
physchem_avg_2022$sampdateID <- paste0(physchem_avg_2022$Plot, "_" ,physchem_avg_2022$Sampling) # ID to group.by() and summarise(mean)

physchem_avg_2022 <- physchem_avg_2022 %>% 
                     group_by(sampdateID) %>% 
                     summarise(Conduct_microS_cm = mean(Conduct_microS_cm, na.rm = TRUE), Temp_10_cm = mean(Temp_10_cm, na.rm = TRUE), pH_soil = mean(pH_soil, na.rm = TRUE), 
                               Redox_pot = mean(Redox_pot, na.rm = TRUE), Water_temp = mean(Water_temp, na.rm = TRUE), O2_percent = mean(O2_percent, na.rm = TRUE), 
                               O2_mg_l = mean(O2_mg_l, na.rm = TRUE), Salinity = mean(Salinity, na.rm = TRUE), pH_water = mean(pH_water, na.rm = TRUE), 
                               across(Plot, ~., .names = "Plot"), across(Treat, ~., .names = "Treat"), 
                               across(Rep, ~., .names = "Rep"), across(Sampling, ~., .names = "Sampling")) %>% # Keeps previous Plot, Treat, Rep and Sampling data.
                               filter(row_number() == max(row_number())) # As the across() function triplicates each row (due to takingo Plot, Rep... from each of the three averaged values) this keeps only one row. 

Sampling_date <- c("2022-06-10", "2022-07-15", "2022-08-02", "2022-08-31") # Macroinvertebrate sampling dates
Sampling <- c("1", "2", "3", "4")
Sam.Date <- data.frame(Sampling_date, Sampling) # data frame to include "Sampling"
Sam.Date$Sampling_date <- as.Date(Sam.Date$Sampling_date)
physchem_avg_2022 <- merge(physchem_avg_2022, Sam.Date, by.x="Sampling", by.y="Sampling") # Assigning Sampling values
physchem_avg_2022$siteID <- paste0(physchem_avg_2022$Plot, "_", physchem_avg_2022$Sampling, "_", physchem_avg_2022$Treat) # Creates siteID to merge later with Hills_ColOdoHet data frame
physchem_avg_2022 <- physchem_avg_2022 %>% 
                     arrange(Sampling_date, Plot) %>% # Sorts by Sampling_date and then by Plot
                     select(Sampling_date, Sampling, Plot, Treat, Rep, Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water, sampdateID, siteID) %>%  # Re-orders.
                     mutate(across(c(Water_temp, O2_percent, O2_mg_l, Salinity, pH_water), ~ifelse(is.nan(.), NA, .))) # Replaces NaN for NA values.

#### 1.2. Creating base biodiv-physchem data frame ####

Hills_Physchem <- merge(Hills_ColOdoHet, physchem_avg_2022, by = "siteID", all = TRUE)
Hills_Physchem <- Hills_Physchem %>% 
                  mutate(across(c(q2.se), ~ifelse(is.nan(.), NA, .))) # Replaces NaN for NA values.

colnames(Hills_Physchem)[colnames(Hills_Physchem) == "Plot.x"] <- "Plot" # Renames ".x" columns
colnames(Hills_Physchem)[colnames(Hills_Physchem) == "Sampling.x"] <- "Sampling"
colnames(Hills_Physchem)[colnames(Hills_Physchem) == "Treat.x"] <- "Treat"

Hills_Physchem <- Hills_Physchem[, !(colnames(Hills_Physchem) %in% c("Sampling.y", "Plot.y", "Treat.y"))] # Removes duplicated ".y" columns

Hills_Physchem <- Hills_Physchem %>% 
                  arrange(Sampling_date, Plot) %>% # Sorts by Sampling_date and then by Plot
                  select(Sampling_date, Sampling, Plot, Treat, Rep, q0.obs, q0.est, q0.se, q1.obs, q1.est, q1.se, q2.obs, q2.est, q2.se,
                         Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water, sampdateID, siteID)  # Re-orders.

################ 2. STATS ###################

##### 2.1. Data validation ######
# According to "Analysing the impact of multiple stressors in aquatic biomonitoring data: A ‘cookbook’ with applications in R" - Feld et al., 2016. ##

#### 2.1.1. Outlier analysis ####

summary(Hills_Physchem$Conduct_microS_cm)
boxplot(Hills_Physchem$Conduct_microS_cm)
boxplot(Hills_Physchem$Conduct_microS_cm)$out
max(boxplot(Hills_Physchem$Conduct_microS_cm)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$Conduct_microS_cm)$out) # prints the minimum outlier value

summary(Hills_Physchem$Temp_10_cm)
boxplot(Hills_Physchem$Temp_10_cm)
boxplot(Hills_Physchem$Temp_10_cm)$out
max(boxplot(Hills_Physchem$Temp_10_cm)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$Temp_10_cm)$out) # prints the minimum outlier value

summary(Hills_Physchem$pH_soil)
boxplot(Hills_Physchem$pH_soil)
boxplot(Hills_Physchem$pH_soil)$out
max(boxplot(Hills_Physchem$pH_soil)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$pH_soil)$out) # prints the minimum outlier value

summary(Hills_Physchem$Redox_pot)
boxplot(Hills_Physchem$Redox_pot)
boxplot(Hills_Physchem$Redox_pot)$out
max(boxplot(Hills_Physchem$Redox_pot)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$Redox_pot)$out) # prints the minimum outlier value

summary(Hills_Physchem$Water_temp)
boxplot(Hills_Physchem$Water_temp)
boxplot(Hills_Physchem$Water_temp)$out
max(boxplot(Hills_Physchem$Water_temp)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$Water_temp)$out) # prints the minimum outlier value
## Outliers identified for var: Water_temp - 31.4 and 29.4, I decide to keep them as they are probable values and not results of sampling error.

summary(Hills_Physchem$O2_percent)
boxplot(Hills_Physchem$O2_percent)
boxplot(Hills_Physchem$O2_percent)$out
max(boxplot(Hills_Physchem$O2_percent)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$O2_percent)$out) # prints the minimum outlier value

summary(Hills_Physchem$O2_mg_l)
boxplot(Hills_Physchem$O2_mg_l)
boxplot(Hills_Physchem$O2_mg_l)$out
max(boxplot(Hills_Physchem$O2_mg_l)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$O2_mg_l)$out) # prints the minimum outlier value

summary(Hills_Physchem$Salinity)
boxplot(Hills_Physchem$Salinity)
boxplot(Hills_Physchem$Salinity)$out
max(boxplot(Hills_Physchem$Salinity)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$Salinity)$out) # prints the minimum outlier value
## Outliers identified for var: Salinity - 0.7466667 0.6600000 0.6600000, I decide to keep them as they are probable values and not results of sampling error.

summary(Hills_Physchem$pH_water)
boxplot(Hills_Physchem$pH_water)
boxplot(Hills_Physchem$pH_water)$out
max(boxplot(Hills_Physchem$pH_water)$out) # prints the maximum outlier value 
min(boxplot(Hills_Physchem$pH_water)$out) # prints the minimum outlier value

## From here onwards tests and models will be evaluated both considering outliers and removing them (suffix: "_nooutliers")

Hills_Physchem_nooutliers <-  Hills_Physchem %>%  # New data frame with a (1/0) Outliers column, then removing outlier rows
                              mutate(Outliers =case_when(Salinity > 0.74  | Water_temp > 29 ~ 1,TRUE ~ 0)) %>% 
                              filter(Outliers == 0)

#### 2.1.2. Testing correlations ####

#### i) Spearman rank correlation ####

## The Spearman rank correlation is a non-parametric measure of association between two variables. 
## It assesses the strength and direction of the monotonic relationship (whether it goes up or down) between two variables. 

### Considering outliers:

cor_matrix <- Hills_Physchem %>% 
              select(q0.obs, q1.obs, q2.obs,
                     Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water, Sampling) %>% 
              na.omit

corr <- round(cor(cor_matrix, method =  "spearman"), 1) # Calculates the Spearman rank correlation matrix

pdf("outputs/Plots/BIO/Corr_plot.pdf", width = 11)
corrplot::corrplot.mixed(corr, order = 'hclust', addrect = 2)
dev.off()

### Removing outliers:

cor_matrix_nooutliers <- Hills_Physchem_nooutliers %>% 
              select(q0.obs, q1.obs, q2.obs,
                     Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water, Sampling) %>% 
              na.omit

corr_nooutliers <- round(cor(cor_matrix_nooutliers, method =  "spearman"), 1) # Calculates the Spearman rank correlation matrix

pdf("outputs/Plots/BIO/Corr_plot_nooutliers.pdf", width = 11)
corrplot::corrplot.mixed(corr_nooutliers, order = 'hclust', addrect = 2)
dev.off()

## Comparing corrplots the effect of removing rows with outliers for Water_temp and Salinity can be seen in correlation between other variables, for having discarded as 
## well data regarding them that was contained in these rows.

# Data Histogram (q0) # Checking each variable distribution (all 9 physicochemical independent variables)
hist(Hills_Physchem$q0.obs)
hist(Hills_Physchem_nooutliers$q0.obs)
hist(Hills_Physchem$q1.obs)
hist(Hills_Physchem_nooutliers$q1.obs)
hist(Hills_Physchem$Conduct_microS_cm)
hist(Hills_Physchem_nooutliers$Conduct_microS_cm)
hist(Hills_Physchem$Water_temp) # Significant variation
hist(Hills_Physchem_nooutliers$Water_temp) # Significant variation
hist(Hills_Physchem$O2_percent) # Significant variation
hist(Hills_Physchem_nooutliers$O2_percent) # Significant variation
hist(Hills_Physchem$Redox_pot)
hist(Hills_Physchem_nooutliers$Redox_pot)
hist(Hills_Physchem$pH_water)
hist(Hills_Physchem_nooutliers$pH_water)
hist(Hills_Physchem$Temp_10_cm)
hist(Hills_Physchem_nooutliers$Temp_10_cm)
hist(Hills_Physchem$Salinity) # Significant variation
hist(Hills_Physchem_nooutliers$Salinity) # Significant variation
hist(Hills_Physchem$pH_soil)
hist(Hills_Physchem_nooutliers$pH_soil)
hist(Hills_Physchem$O2_mg_l)
hist(Hills_Physchem_nooutliers$O2_mg_l)

#### ii) Variance inflation factors (VIF) ####

### Method: Excluding variables with VIF > 5 stepwise, starting with the variable that has the highest VIF
### A VIF > 8 is applied in Feld et al., 2016. We lower the threshold as pH_water resulted in High collinearity in posterior tested GLMMs.

###  Considering outliers:

# Step 1a: all variables considered
selected_vars <- c('Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected <- (Hills_Physchem[, selected_vars])
vif(data_selected,) 

# Step 2a: Removing variable with highest VIF (and VIF>5) - O2_mg_l (VIF = 1309.181048)
selected_vars_2 <- c('Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_2 <- (Hills_Physchem[, selected_vars_2])
vif(data_selected_2) 

# Step 3a: Removing next variable with highest VIF (and VIF>5) - Temp_10_cm  (VIF = 15.981256)
selected_vars_3 <- c('Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_3 <- (Hills_Physchem[, selected_vars_3])
vif(data_selected_3) 

# Step 4a: Removing next variable with highest VIF (and VIF>5) - pH_water  (VIF = 5.385582)
selected_vars_4 <- c('Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity')
data_selected_4 <- (Hills_Physchem[, selected_vars_4])
vif(data_selected_4) 

## Removing O2_mg_l, Temp_10_cm and pH_water achieves already a set where all variables have VIF < 5.

### Removing outliers:

# Step 1b: all variables considered
selected_vars_5 <- c('Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_5 <- (Hills_Physchem_nooutliers[, selected_vars_5])
vif(data_selected_5,) 

# Step 2b: Removing variable with highest VIF (and VIF>5) - O2_mg_l (VIF = 1845.687620)
selected_vars_6 <- c('Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_6 <- (Hills_Physchem_nooutliers[, selected_vars_6])
vif(data_selected_6) 

# Step 3b: Removing next variable with highest VIF (and VIF>5) - Temp_10_cm  (VIF = 21.248961)
selected_vars_7 <- c('Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_7 <- (Hills_Physchem_nooutliers[, selected_vars_7])
vif(data_selected_7) 

# Step 4b Removing next variable with highest VIF (and VIF>5) - pH_water  (VIF = 5.754363)
selected_vars_8 <- c('Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity')
data_selected_8 <- (Hills_Physchem_nooutliers[, selected_vars_8])
vif(data_selected_8) 

# Extra step: Step 5b Removing next variable with highest VIF (and VIF>5) - pH_water  (VIF = 5.754363) and adding Sampling
selected_vars_9 <- c('Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'Sampling')
data_selected_9 <- (Hills_Physchem_nooutliers[, selected_vars_9])
vif(data_selected_9) 

## Removing O2_mg_l, Temp_10_cm and pH_water achieves already a set where all variables have VIF < 5.

#### iii) Data Dendrogram ####

###  Considering outliers:

# Dend 1: Considering all variables (without VIF removal)
# Variable clustering:
similarity="pearson"
vclust_1 <- varclus(x=as.matrix(data_selected),
                  similarity=similarity,
                  type="data.matrix", 
                  method="complete",
                  na.action=na.retain,trans="abs")
dend_1 <- as.dendrogram(vclust_1)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_1 <- color_labels(dend_1, h=1-0.7,col=col_vector)
dend_1 <- color_branches(dend_1, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/BIO/Cluster_variables_pearson_1.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_1,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Dend 2: Considering only remaining variables after VIF removal
# Variable clustering:
similarity="pearson"
vclust_2 <- varclus(x=as.matrix(data_selected_4),
                  similarity=similarity,
                  type="data.matrix", 
                  method="complete",
                  na.action=na.retain,trans="abs")
dend_2 <- as.dendrogram(vclust_2)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_2 <- color_labels(dend_2, h=1-0.7,col=col_vector)
dend_2 <- color_branches(dend_2, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/BIO/Cluster_variables_pearson_2_VIFvarsremoval.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_2,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Note: Dend 2 (removing variables after the VIF analysis) achieves already that all variables stay underneath a 0.75 threshold.

###  Removing outliers:

# Dend 3: Considering all variables (without VIF removal)
# Variable clustering:
similarity="pearson"
vclust_3 <- varclus(x=as.matrix(data_selected_5),
                    similarity=similarity,
                    type="data.matrix", 
                    method="complete",
                    na.action=na.retain,trans="abs")
dend_3 <- as.dendrogram(vclust_3)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_3 <- color_labels(dend_3, h=1-0.7,col=col_vector)
dend_3 <- color_branches(dend_3, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/BIO/Cluster_variables_pearson_3_nooutliers.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_3,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Dend 4: Considering only remaining variables after VIF removal
# Variable clustering:
similarity="pearson"
vclust_4 <- varclus(x=as.matrix(data_selected_8),
                    similarity=similarity,
                    type="data.matrix", 
                    method="complete",
                    na.action=na.retain,trans="abs")
dend_4 <- as.dendrogram(vclust_4)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_4 <- color_labels(dend_4, h=1-0.7,col=col_vector)
dend_4 <- color_branches(dend_4, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/BIO/Cluster_variables_pearson_4_VIFvarsremoval_nooutliers.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_4,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

## Dend 5  Considering Sampling (removing variables after the VIF analysis) achieves already that all variables stay underneath a 0.75 threshold.

# Dend 5: Considering only remaining variables after VIF removal
# Variable clustering:
similarity="pearson"
vclust_5 <- varclus(x=as.matrix(data_selected_9),
                    similarity=similarity,
                    type="data.matrix", 
                    method="complete",
                    na.action=na.retain,trans="abs")
dend_5 <- as.dendrogram(vclust_5)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_5 <- color_labels(dend_5, h=1-0.7,col=col_vector)
dend_5 <- color_branches(dend_5, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/BIO/Cluster_variables_pearson_4_VIFvarsremoval_nooutliers2.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_5,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()


#### 2.2. GLMMs ####

Hills_Physchem$Rep_Treat <- paste0(Hills_Physchem$Rep, "_", Hills_Physchem$Treat) # Creates the random effects variable, which is, in this case, controlled by the blocks design.
Hills_Physchem_nooutliers$Rep_Treat <- paste0(Hills_Physchem_nooutliers$Rep, "_", Hills_Physchem_nooutliers$Treat) # Creates the random effects variable, which is, in this case, controlled by the blocks design.

# Testing effect of editing vector types: # These were a trial to check if changing dataframe structure would impact the GLMMs result. They sopped fitting...
# Hills_Physchem$Sampling <- as.factor(Hills_Physchem$Sampling)   
# Hills_Physchem$Treat <- as.character(Hills_Physchem$Treat)

##### Checking q0 vs Samplings ####
# Testing for linear or quadratic relation

Hills_Physchem_summary_q0 <- Hills_Physchem %>%
                group_by(Sampling) %>%
                summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

q0_Sampling <- ggplot(Hills_Physchem, aes(Sampling, q0.obs)) +
                geom_point()  +
                geom_point(data = Hills_Physchem_summary_q0, aes(x = Sampling, y = mean_q0.obs), shape = 19, colour = "black", size = 6) 

print(q0_Sampling) # Relation between q0 and Sampling seems to be quadratic, "Treat*I(Sampling^2)" interaction should be tested.

##### Checking q1 vs Samplings ####
# Testing for linear or quadratic relation

Hills_Physchem_summary_q1 <- Hills_Physchem %>%
                group_by(Sampling) %>%
                summarise(mean_q1.obs = mean(q1.obs),se_q1.obs = sd(q1.obs) / sqrt(n()))

q1_Sampling <- ggplot(Hills_Physchem, aes(Sampling, q1.obs))  +
                geom_point() +
                geom_point(data = Hills_Physchem_summary_q1, aes(x = Sampling, y = mean_q1.obs), shape = 19, colour = "black", size = 6) 

print(q1_Sampling) # Relation between q1 and Sampling seems to be quadratic, "Treat*I(Sampling^2)" interaction should be tested.

##### Model 1: q0 - Poisson - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering all original variables (not considering all prev. correlation analysis)  ####

# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q0.pois1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois1, plot = T)
summary(glmm.q0.pois1)
car::Anova(glmm.q0.pois1)
performance::r2(glmm.q0.pois1)
performance::check_collinearity(glmm.q0.pois1)
performance::check_singularity(glmm.q0.pois1)

###  Removing outliers:

glmm.q0.pois1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois1.noout, plot = T)
summary(glmm.q0.pois1.noout)
car::Anova(glmm.q0.pois1.noout)
performance::r2(glmm.q0.pois1.noout)
performance::check_collinearity(glmm.q0.pois1.noout)
performance::check_singularity(glmm.q0.pois1.noout)

##### Model 2: q0 - Poisson - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering only remaining variables after correlation analysis  ####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity 
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q0.pois2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity  + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois2, plot = T)
summary(glmm.q0.pois2)
car::Anova(glmm.q0.pois2)
performance::r2(glmm.q0.pois2)
performance::check_collinearity(glmm.q0.pois2)
performance::check_singularity(glmm.q0.pois2)

###  Removing outliers:

glmm.q0.pois2.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity  + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois2.noout, plot = T)
summary(glmm.q0.pois2.noout)
car::Anova(glmm.q0.pois2.noout)
performance::r2(glmm.q0.pois2.noout)
performance::check_collinearity(glmm.q0.pois2.noout)
performance::check_singularity(glmm.q0.pois2.noout)

##### Model 3: q0 - Poisson - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis  ####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity 

###  Considering outliers:

glmm.q0.pois3 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                           (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois3, plot = T)
summary(glmm.q0.pois3)
car::Anova(glmm.q0.pois3)
performance::r2(glmm.q0.pois3)
performance::check_collinearity(glmm.q0.pois3)
performance::check_singularity(glmm.q0.pois3)
visreg(glmm.q0.pois3, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.pois3.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + 
                                 Salinity + (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois3.noout, plot = T)
summary(glmm.q0.pois3.noout)
car::Anova(glmm.q0.pois3.noout)
performance::r2(glmm.q0.pois3.noout)
performance::check_collinearity(glmm.q0.pois3.noout)
performance::check_singularity(glmm.q0.pois3.noout)
visreg(glmm.q0.pois3.noout, scale="response") # Plotting conditional residuals

##### Model 4: q0 - Poisson - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Removing High VIF vars from Model 3  ####
# Vars removed from Model 3 (High VIF): pH_soil + Redox_pot + Conduct_microS_cm
# Family: Poisson
# Interacting independent variables: Treat*Sampling and Treat*I(Sampling^2)
# Additional independent variables: Conduct_microS_cm + Water_temp + O2_percent + Salinity 

###  Considering outliers:

glmm.q0.pois4 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Water_temp + O2_percent + Salinity +
                           (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois4, plot = T)
summary(glmm.q0.pois4)
car::Anova(glmm.q0.pois4)
performance::r2(glmm.q0.pois4)
performance::check_collinearity(glmm.q0.pois4)
performance::check_singularity(glmm.q0.pois4)
visreg(glmm.q0.pois4, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.pois4.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Water_temp + O2_percent + Salinity +
                                 (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois4.noout, plot = T)
summary(glmm.q0.pois4.noout)
car::Anova(glmm.q0.pois4.noout)
performance::r2(glmm.q0.pois4.noout)
performance::check_collinearity(glmm.q0.pois4.noout)
performance::check_singularity(glmm.q0.pois4.noout)
visreg(glmm.q0.pois4.noout, scale="response") # Plotting conditional residuals

##### Model 4b: q0 - Poisson - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Removing High VIF vars from Model 3  ####
# Vars removed from Model 3 (High VIF): pH_soil + Redox_pot + Conduct_microS_cm
# Family: Poisson
# Interacting independent variables: Treat*Sampling and Treat*I(Sampling^2)
# Additional independent variables: Conduct_microS_cm + Water_temp + O2_percent + Salinity 

###  Considering outliers:

glmm.q0.pois5 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent + Salinity +
                           (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois5, plot = T)
summary(glmm.q0.pois5)
car::Anova(glmm.q0.pois5)
performance::r2(glmm.q0.pois5)
performance::check_collinearity(glmm.q0.pois5)
performance::check_singularity(glmm.q0.pois5)
visreg(glmm.q0.pois5, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.pois5.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent + Salinity +
                                 (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.pois5.noout, plot = T)
summary(glmm.q0.pois5.noout)
car::Anova(glmm.q0.pois5.noout)
performance::r2(glmm.q0.pois5.noout)
performance::check_collinearity(glmm.q0.pois5.noout)
performance::check_singularity(glmm.q0.pois5.noout)
visreg(glmm.q0.pois5.noout, scale="response") # Plotting conditional residuals

##### Model 5: q0 - nbinom2 - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: nbinom2
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q0.nbinom1 <-  glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                              (1|Rep) , family = "nbinom2")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.nbinom1, plot = T)
summary(glmm.q0.nbinom1)
car::Anova(glmm.q0.nbinom1)
performance::r2(glmm.q0.nbinom1)
performance::check_collinearity(glmm.q0.nbinom1)
performance::check_singularity(glmm.q0.nbinom1)
visreg(glmm.q0.nbinom1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.nbinom1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                                   (1|Rep) , family = "nbinom2")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.nbinom1.noout, plot = T)
summary(glmm.q0.nbinom1.noout)
car::Anova(glmm.q0.nbinom1.noout)
performance::r2(glmm.q0.nbinom1)
performance::check_collinearity(glmm.q0.nbinom1.noout)
performance::check_singularity(glmm.q0.nbinom1.noout)
visreg(glmm.q0.nbinom1.noout, scale="response") # Plotting conditional residuals

# # transf sqrt:
# Hills_Physchem$q0.sqrt <- sqrt(Hills_Physchem$q0.obs) # Test

##### Model 6: q0 - Gaussian - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering all original variables (not considering all prev. correlation analysis) ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q0.gaus1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus1, plot = T)
summary(glmm.q0.gaus1)
car::Anova(glmm.q0.gaus1)
performance::r2(glmm.q0.gaus1)
performance::check_collinearity(glmm.q0.gaus1)
performance::check_singularity(glmm.q0.gaus1)
visreg(glmm.q0.gaus1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus1.noout, plot = T)
summary(glmm.q0.gaus1.noout)
car::Anova(glmm.q0.gaus1.noout)
performance::r2(glmm.q0.gaus1.noout)
performance::check_collinearity(glmm.q0.gaus1.noout)
performance::check_singularity(glmm.q0.gaus1.noout)
visreg(glmm.q0.gaus1.noout, scale="response") # Plotting conditional residuals

##### Model 7: q0 - Gaussian - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q0.gaus2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus2, plot = T)
summary(glmm.q0.gaus2)
car::Anova(glmm.q0.gaus2)
performance::r2(glmm.q0.gaus2)
performance::check_collinearity(glmm.q0.gaus2)
performance::check_singularity(glmm.q0.gaus2)
visreg(glmm.q0.gaus2, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus2.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + Salinity + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus2.noout, plot = T)
summary(glmm.q0.gaus2.noout)
car::Anova(glmm.q0.gaus2.noout)
performance::r2(glmm.q0.gaus2.noout)
performance::check_collinearity(glmm.q0.gaus2.noout)
performance::check_singularity(glmm.q0.gaus2.noout)
visreg(glmm.q0.gaus2.noout, scale="response") # Plotting conditional residuals

##### Model 8: q0 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering all original variables (not considering all prev. correlation analysis) ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep

###  Considering outliers:

glmm.q0.gaus3 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus3, plot = T)
summary(glmm.q0.gaus3)
car::Anova(glmm.q0.gaus3)
performance::r2(glmm.q0.gaus3)
performance::check_collinearity(glmm.q0.gaus3)
performance::check_singularity(glmm.q0.gaus3)
visreg(glmm.q0.gaus3, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus3.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus3.noout, plot = T)
summary(glmm.q0.gaus3.noout)
car::Anova(glmm.q0.gaus3.noout)
performance::r2(glmm.q0.gaus3.noout)
performance::check_collinearity(glmm.q0.gaus3.noout)
performance::check_singularity(glmm.q0.gaus3.noout)
visreg(glmm.q0.gaus3.noout, scale="response") # Plotting conditional residuals

##### Model 9: q0 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep

###  Considering outliers:

glmm.q0.gaus4 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus4, plot = T)
summary(glmm.q0.gaus4)
car::Anova(glmm.q0.gaus4)
performance::r2(glmm.q0.gaus4)
performance::check_collinearity(glmm.q0.gaus4)
performance::check_singularity(glmm.q0.gaus4)
visreg(glmm.q0.gaus4, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus4.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus4.noout, plot = T)
summary(glmm.q0.gaus4.noout)
car::Anova(glmm.q0.gaus4.noout)
performance::r2(glmm.q0.gaus4.noout)
performance::check_collinearity(glmm.q0.gaus4.noout)
performance::check_singularity(glmm.q0.gaus4.noout)
visreg(glmm.q0.gaus4.noout, scale="response") # Plotting conditional residuals

##### Model 10: q0 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Only vars after correlation analysis and removing High VIF stepwise from Mod. 9 ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Water_temp + Salinity
# Independent variables removed from model 9: Redox_pot 
# Random effect: Rep

###  Considering outliers:

## Note: If I iterate this model and remove stepwise the highest VIF vars, finally I have to remove all variables, there's never a point without High correlations...

glmm.q0.gaus5 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent + pH_soil + Water_temp + Conduct_microS_cm +
                             Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus5, plot = T)
summary(glmm.q0.gaus5)
car::Anova(glmm.q0.gaus5)
performance::r2(glmm.q0.gaus5)
performance::check_collinearity(glmm.q0.gaus5)
performance::check_singularity(glmm.q0.gaus5)
visreg(glmm.q0.gaus5, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus5.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Water_temp +
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus5.noout, plot = T)
summary(glmm.q0.gaus5.noout)
car::Anova(glmm.q0.gaus5.noout)
performance::r2(glmm.q0.gaus5.noout)
performance::check_collinearity(glmm.q0.gaus5.noout)
performance::check_singularity(glmm.q0.gaus5.noout)
visreg(glmm.q0.gaus5.noout, scale="response") # Plotting conditional residuals

##### Model 11: q0 - Gaussian - All interactions bet. vars. - "Rep" Random factor - Only vars after correlation analysis and removing High VIF stepwise from Mod. 9 ####

# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Water_temp + Salinity
# Independent variables removed from model 9: Redox_pot 
# Random effect: Rep

###  Considering outliers:

## Note: If I iterate this model and remove stepwise the highest VIF vars, finally I have to remove all variables, there's never a point without High correlations...

glmm.q0.gaus6 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent*pH_soil + O2_percent*Water_temp + O2_percent*Conduct_microS_cm +
                           O2_percent*Salinity + pH_soil*Water_temp + pH_soil*Conduct_microS_cm + pH_soil*Salinity + Water_temp*Conduct_microS_cm  +  Water_temp*Salinity +
                           Conduct_microS_cm*Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus6, plot = T)
summary(glmm.q0.gaus6)
car::Anova(glmm.q0.gaus6)
performance::r2(glmm.q0.gaus6)
performance::check_collinearity(glmm.q0.gaus6)
performance::check_singularity(glmm.q0.gaus6)
visreg(glmm.q0.gaus6, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q0.gaus6.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent*pH_soil + O2_percent*Water_temp + O2_percent*Conduct_microS_cm +
                                O2_percent*Salinity + pH_soil*Water_temp + pH_soil*Conduct_microS_cm + pH_soil*Salinity + Water_temp*Conduct_microS_cm  +  Water_temp*Salinity +
                                Conduct_microS_cm*Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus6.noout, plot = T)
summary(glmm.q0.gaus6.noout)
car::Anova(glmm.q0.gaus6.noout)
performance::r2(glmm.q0.gaus6.noout)
performance::check_collinearity(glmm.q0.gaus6.noout)
performance::check_singularity(glmm.q0.gaus6.noout)
visreg(glmm.q0.gaus6.noout, scale="response") # Plotting conditional residuals

##### Model 11b: q0 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling + Treat*Sampling2
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + O2_percent + Salinity
# Random effect: Rep

###  Considering outliers:

Hills_Physchem$Sampling2 <- (Hills_Physchem$Sampling)^2 # Sampling^2 variable created for plotting purposes.

glmm.q0.gaus7 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*Sampling2 + Conduct_microS_cm + pH_soil + Redox_pot + 
                           O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Note: Model diagnostic results are the same when using Treat*I(Sampling^2) or Treat*Sampling2.

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus7, plot = T)
summary(glmm.q0.gaus7)
car::Anova(glmm.q0.gaus7)
r.squaredGLMM(glmm.q0.gaus7) # Calculates Pseudo-R-squared for Generalized Mixed-Effect models
performance::r2(glmm.q0.gaus7)
performance::check_collinearity(glmm.q0.gaus7)
performance::check_singularity(glmm.q0.gaus7)
visreg(glmm.q0.gaus7, xvar="Sampling2", by = "Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals
visreg(glmm.q0.gaus7, xvar="Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals

###  Removing outliers:

Hills_Physchem_nooutliers$Sampling2 <- (Hills_Physchem_nooutliers$Sampling)^2 # Sampling^2 variable created for plotting purposes.

glmm.q0.gaus7.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Redox_pot +  
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus7.noout, plot = T)
summary(glmm.q0.gaus7.noout)
car::Anova(glmm.q0.gaus7.noout) 
r.squaredGLMM(glmm.q0.gaus7.noout) # Calculates Pseudo-R-squared for Generalized Mixed-Effect models
performance::r2(glmm.q0.gaus7.noout)
performance::check_collinearity(glmm.q0.gaus7.noout)
performance::check_singularity(glmm.q0.gaus7.noout)
visreg(glmm.q0.gaus7.noout, scale="response") # Plotting conditional residuals

##### Model 12: q1 - Poisson - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering all original variables (not considering all prev. correlation analysis)  ####

# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q1.pois1 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois1, plot = T)
summary(glmm.q1.pois1)
car::Anova(glmm.q1.pois1)
performance::r2(glmm.q1.pois1)
performance::check_collinearity(glmm.q1.pois1)
performance::check_singularity(glmm.q1.pois1)
visreg(glmm.q1.pois1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.pois1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois1.noout, plot = T)
summary(glmm.q1.pois1.noout)
car::Anova(glmm.q1.pois1.noout)
performance::r2(glmm.q1.pois1.noout)
performance::check_collinearity(glmm.q1.pois1.noout)
performance::check_singularity(glmm.q1.pois1.noout)

##### Model 13: q1 - Poisson - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering only remaining variables after correlation analysis  ####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity 
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q1.pois2 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity  + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois2, plot = T)
summary(glmm.q1.pois2)
car::Anova(glmm.q1.pois2)
performance::r2(glmm.q1.pois2)
performance::check_collinearity(glmm.q1.pois2)
performance::check_singularity(glmm.q1.pois2)

###  Removing outliers:

glmm.q1.pois2.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity  + (1|Rep_Treat) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois2.noout, plot = T)
summary(glmm.q1.pois2.noout)
car::Anova(glmm.q1.pois2.noout)
performance::r2(glmm.q1.pois2.noout)
performance::check_collinearity(glmm.q1.pois2.noout)
performance::check_singularity(glmm.q1.pois2.noout)

##### Model 14: q1 - Poisson - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis  ####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity 

###  Considering outliers:

glmm.q1.pois3 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                           (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois3, plot = T)
summary(glmm.q1.pois3)
car::Anova(glmm.q1.pois3)
performance::r2(glmm.q1.pois3)
performance::check_collinearity(glmm.q1.pois3)
performance::check_singularity(glmm.q1.pois3)
visreg(glmm.q1.pois3, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.pois3.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q0.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + 
                                 Salinity + (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois3.noout, plot = T)
summary(glmm.q1.pois3.noout)
car::Anova(glmm.q1.pois3.noout)
performance::r2(glmm.q1.pois3.noout)
performance::check_collinearity(glmm.q1.pois3.noout)
performance::check_singularity(glmm.q1.pois3.noout)
visreg(glmm.q1.pois3.noout, scale="response") # Plotting conditional residuals

##### Model 15: q1- Poisson - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Removing High VIF vars from Model 14  ####
# Vars removed from Model 14 (High VIF): pH_soil + Redox_pot + Conduct_microS_cm
# Family: Poisson
# Interacting independent variables: Treat*Sampling and Treat*I(Sampling^2)
# Additional independent variables: Water_temp + O2_percent + Salinity 

###  Considering outliers:

glmm.q1.pois4 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Water_temp + O2_percent + Salinity +
                           (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois4, plot = T)
summary(glmm.q1.pois4)
car::Anova(glmm.q1.pois4)
performance::r2(glmm.q1.pois4)
performance::check_collinearity(glmm.q1.pois4)
performance::check_singularity(glmm.q1.pois4)
visreg(glmm.q1.pois4, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.pois4.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Water_temp + O2_percent + Salinity +
                                 (1|Rep) , family = "poisson")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.pois4.noout, plot = T)
summary(glmm.q1.pois4.noout)
car::Anova(glmm.q1.pois4.noout)
performance::r2(glmm.q1.pois4.noout)
performance::check_collinearity(glmm.q1.pois4.noout)
performance::check_singularity(glmm.q1.pois4.noout)
visreg(glmm.q1.pois4.noout, scale="response") # Plotting conditional residuals

##### Model 16: q1 - nbinom2 - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: nbinom2
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q1.nbinom1 <-  glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                              (1|Rep) , family = "nbinom2")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.nbinom1, plot = T)
summary(glmm.q1.nbinom1)
car::Anova(glmm.q1.nbinom1)
performance::r2(glmm.q1.nbinom1)
performance::check_collinearity(glmm.q1.nbinom1)
performance::check_singularity(glmm.q1.nbinom1)
visreg(glmm.q1.nbinom1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.nbinom1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity +
                                   (1|Rep) , family = "nbinom2")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.nbinom1.noout, plot = T)
summary(glmm.q1.nbinom1.noout)
car::Anova(glmm.q1.nbinom1.noout)
performance::r2(glmm.q1.nbinom1.noout)
performance::check_collinearity(glmm.q1.nbinom1.noout)
performance::check_singularity(glmm.q1.nbinom1.noout)
visreg(glmm.q1.nbinom1.noout, scale="response") # Plotting conditional residuals

# # transf sqrt:
# Hills_Physchem$q0.sqrt <- sqrt(Hills_Physchem$q0.obs) # Test

##### Model 17: q1 - Gaussian - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering all original variables (not considering all prev. correlation analysis) ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q1.gaus1 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus1, plot = T)
summary(glmm.q1.gaus1)
car::Anova(glmm.q1.gaus1)
performance::r2(glmm.q1.gaus1)
performance::check_collinearity(glmm.q1.gaus1)
performance::check_singularity(glmm.q1.gaus1)
visreg(glmm.q1.gaus1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus1.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus1.noout, plot = T)
summary(glmm.q1.gaus1.noout)
car::Anova(glmm.q1.gaus1.noout)
performance::r2(glmm.q1.gaus1.noout)
performance::check_collinearity(glmm.q1.gaus1.noout)
performance::check_singularity(glmm.q1.gaus1.noout)
visreg(glmm.q1.gaus1.noout, scale="response") # Plotting conditional residuals

##### Model 18: q1 - Gaussian - Treat*Sampling interaction - "Rep_Treat" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep_Treat

###  Considering outliers:

glmm.q1.gaus2 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus2, plot = T)
summary(glmm.q1.gaus2)
car::Anova(glmm.q1.gaus2)
performance::r2(glmm.q1.gaus2)
performance::check_collinearity(glmm.q1.gaus2)
performance::check_singularity(glmm.q1.gaus2)
visreg(glmm.q1.gaus2, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus2.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + Salinity + (1|Rep_Treat) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus2.noout, plot = T)
summary(glmm.q1.gaus2.noout)
car::Anova(glmm.q1.gaus2.noout)
performance::r2(glmm.q1.gaus2.noout)
performance::check_collinearity(glmm.q1.gaus2.noout)
performance::check_singularity(glmm.q1.gaus2.noout)
visreg(glmm.q1.gaus2.noout, scale="response") # Plotting conditional residuals

##### Model 19: q1 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering all original variables (not considering all prev. correlation analysis) ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep

###  Considering outliers:

glmm.q1gaus3 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1gaus3, plot = T)
summary(glmm.q1gaus3)
car::Anova(glmm.q1gaus3)
performance::r2(glmm.q1gaus3)
performance::check_collinearity(glmm.q1gaus3)
performance::check_singularity(glmm.q1gaus3)
visreg(glmm.q1gaus3, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus3.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus3.noout, plot = T)
summary(glmm.q1.gaus3.noout)
car::Anova(glmm.q1.gaus3.noout)
performance::r2(glmm.q1.gaus3.noout)
performance::check_collinearity(glmm.q1.gaus3.noout)
performance::check_singularity(glmm.q1.gaus3.noout)
visreg(glmm.q1.gaus3.noout, scale="response") # Plotting conditional residuals

##### Model 20: q1 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep

###  Considering outliers:

glmm.q1.gaus4 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus4, plot = T)
summary(glmm.q1.gaus4)
car::Anova(glmm.q1.gaus4)
performance::r2(glmm.q1.gaus4)
performance::check_collinearity(glmm.q1.gaus4)
performance::check_singularity(glmm.q1.gaus4)
visreg(glmm.q1.gaus4, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus4.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp +
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus4.noout, plot = T)
summary(glmm.q1.gaus4.noout)
car::Anova(glmm.q1.gaus4.noout)
performance::r2(glmm.q1.gaus4.noout)
performance::check_collinearity(glmm.q1.gaus4.noout)
performance::check_singularity(glmm.q1.gaus4.noout)
visreg(glmm.q1.gaus4.noout, scale="response") # Plotting conditional residuals

##### Model 20b: q1 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis - Removing Water_temp ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep

###  Considering outliers:

glmm.q1.gaus8 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + 
                           O2_percent + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus8, plot = T)
summary(glmm.q1.gaus8)
car::Anova(glmm.q1.gaus8)
performance::r2(glmm.q1.gaus8)
performance::check_collinearity(glmm.q1.gaus8)
performance::check_singularity(glmm.q1.gaus8)
visreg(glmm.q1.gaus8, scale="response") # Plotting conditional residuals

ggplot(data = Hills_Physchem, aes(Treat, Conduct_microS_cm)) +
  geom_point() + stat_summary(fun = "mean", geom = "point", size = 5)

# visreg(glmm.q1.gaus8, scale="response", "Sampling", by = "Treat") 
# visreg(glmm.q1.gaus8, scale="response", "Treat", by = "Sampling") 
# visreg(glmm.q1.gaus8, "Treat") 
# visreg(glmm.q1.gaus8, overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals
# visreg(glmm.q1.gaus8) # Plotting conditional residuals
# car::residualPlots(glmm.q1.gaus8, terms = "Treat")
# predict(glmm.q1.gaus8, type = "response")
# residuals(glmm.q1.gaus8) # Trying to plot partial residuals
# terms(glmm.q1.gaus8)

glmm.q1.gaus25 <- glm(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + 
                           O2_percent + Salinity + (1/Rep)+1 , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus25, plot = T)
summary(glmm.q1.gaus25)
car::Anova(glmm.q1.gaus25)
performance::r2(glmm.q1.gaus25)
performance::check_collinearity(glmm.q1.gaus25)
performance::check_singularity(glmm.q1.gaus25)
visreg(glmm.q1.gaus25, scale="response") # Plotting conditional residuals

predict(glmm.q1.gaus25, type = "terms")
residuals(glmm.q1.gaus8, type = "partials") # Trying to plot partial residuals
terms(glmm.q1.gaus8)


###  Removing outliers:

glmm.q1.gaus8.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Redox_pot +
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus8.noout, plot = T)
summary(glmm.q1.gaus8.noout)
car::Anova(glmm.q1.gaus8.noout)
performance::r2(glmm.q1.gaus8.noout)
performance::check_collinearity(glmm.q1.gaus8.noout)
performance::check_singularity(glmm.q1.gaus8.noout)
visreg(glmm.q1.gaus8.noout, scale="response") # Plotting conditional residuals

##### Model 21: q1 - Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Only vars after correlation analysis and removing High VIF stepwise from Mod. 20 ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Water_temp + Salinity
# Independent variables removed from model 9: Redox_pot 
# Random effect: Rep

###  Considering outliers:

## Note: If I iterate this model and remove stepwise the highest VIF vars, finally I have to remove all variables, there's never a point without High correlations...

glmm.q1.gaus5 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent + pH_soil + Water_temp + Conduct_microS_cm +
                             Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus5, plot = T)
summary(glmm.q1.gaus5)
car::Anova(glmm.q1.gaus5)
performance::r2(glmm.q1.gaus5)
performance::check_collinearity(glmm.q1.gaus5)
performance::check_singularity(glmm.q1.gaus5)
visreg(glmm.q1.gaus5, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus5.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling  + Treat*I(Sampling^2)  + Conduct_microS_cm + pH_soil + Water_temp +
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus5.noout, plot = T)
summary(glmm.q1.gaus5.noout)
car::Anova(glmm.q1.gaus5.noout)
performance::r2(glmm.q1.gaus5.noout)
performance::check_collinearity(glmm.q1.gaus5.noout)
performance::check_singularity(glmm.q1.gaus5.noout)
visreg(glmm.q1.gaus5.noout, scale="response") # Plotting conditional residuals

##### Model 22: q1 - Gaussian - All interactions bet. vars. - "Rep" Random factor - Only vars after correlation analysis and removing High VIF stepwise from Mod. 20 ####

# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + pH_soil + Water_temp + Salinity + O2_percent
# Independent variables removed from model 9: Redox_pot 
# Random effect: Rep_Treat

###  Considering outliers:

## Note: If I iterate this model and remove stepwise the highest VIF vars, finally I have to remove all variables, there's never a point without High correlations...

glmm.q1.gaus6 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent*pH_soil + O2_percent*Water_temp + O2_percent*Conduct_microS_cm +
                           O2_percent*Salinity + pH_soil*Water_temp + pH_soil*Conduct_microS_cm + pH_soil*Salinity + Water_temp*Conduct_microS_cm  +  Water_temp*Salinity +
                           Conduct_microS_cm*Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus6, plot = T)
summary(glmm.q1.gaus6)
car::Anova(glmm.q1.gaus6)
performance::r2(glmm.q1.gaus6)
performance::check_collinearity(glmm.q1.gaus6)
performance::check_singularity(glmm.q1.gaus6)
visreg(glmm.q1.gaus6, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.q1.gaus6.noout <- glmmTMB(data = Hills_Physchem_nooutliers, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + O2_percent*pH_soil + O2_percent*Water_temp + O2_percent*Conduct_microS_cm +
                                O2_percent*Salinity + pH_soil*Water_temp + pH_soil*Conduct_microS_cm + pH_soil*Salinity + Water_temp*Conduct_microS_cm  +  Water_temp*Salinity +
                                Conduct_microS_cm*Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q1.gaus6.noout, plot = T)
summary(glmm.q1.gaus6.noout)
car::Anova(glmm.q1.gaus6.noout)
performance::r2(glmm.q1.gaus6.noout)
performance::check_collinearity(glmm.q1.gaus6.noout)
performance::check_singularity(glmm.q1.gaus6.noout)
visreg(glmm.q1.gaus6.noout, scale="response") # Plotting conditional residuals
residuals()



##### Selected model for q0: Model 11b (not removing outliers) ####

## Selection criteria: 
# 1. Considering only remaining variables after initial correlation analyses (Rem. high correlated vars. if Dendrogram > 0.75 and VIF > 5).
# 2. Considers Treatments and Sampling^2 interaction.
# 3. Considers Rep (Blocks) as random effect.
# 4. Does not result in significant deviations or patterns for DHARMA residual diagnostics.
# 5. Removes Water_temp as ind. variable after plotting Dend 5 (incl. Sampling), corr > 0.75 with Sampling.
# As the model ANOVA indicates significant effect of the interaction  Treat*Sampling2 (a factorial and a continuous variable), post-hoc tests are not 
# the most adequate analyses. The effect is referred to from ANOVA results and "q0 vs Treat for each Sampling" plots (with empirical data).
# Selected model:
# glmm.q0.gaus7 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Treat*Sampling2 + Conduct_microS_cm + pH_soil + Redox_pot + 
#                  O2_percent + Salinity + (1|Rep) , family = "gaussian")

## i) q0 - Plot q0 vs Treat for each Sampling ####
# Using empirical (not predicted) data

ColOdoHet_summary_q0_Sampling <- Hills_ColOdoHet %>%
                                  group_by(Sampling, Treat) %>%
                                  summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

# 1st version: facet_wrap, 4 plots, each representing ine Sampling

q0.Treats_Sampling1 <- ggplot(Hills_Physchem, aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                                geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10)+
                                scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                                scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                                theme_bw() +
                                ylab("") +
                                ggtitle("") +
                                theme(plot.title = element_text(size=20, hjust=0.5)) +
                                theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                      axis.title.y =  element_blank(), axis.title.x = element_blank(), legend.position = "none", 
                                      axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                                geom_point(data = ColOdoHet_summary_q0_Sampling, aes(x = Treat, y = mean_q0.obs), shape = 19, colour = "black", size = 12) +
                                geom_point(data = ColOdoHet_summary_q0_Sampling, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 10) +
                                geom_errorbar(data = ColOdoHet_summary_q0_Sampling, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.7, size = 1) +
                                scale_y_continuous(limits = c(1, 12), breaks = seq(2, 12, by = 2)) +
                                facet_wrap(~Sampling, labeller = as_labeller(c("1" = "Sampling 1",
                                                                               "2" = "Sampling 2",
                                                                               "3" = "Sampling 3",
                                                                               "4" = "Sampling 4"))) # adds the "Sampling" word to each gray facet title.

print(q0.Treats_Sampling1)

ggsave("outputs/Plots/BIO/q0.Treats_Sampling.pdf", plot = q0.Treats_Sampling1 ,width = 10, height = 10)

# 2nd version: One plot q0 vs Sampling:

q0.Treats_Sampling2 <- ggplot(Hills_Physchem, aes(Sampling, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                                  geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.1,shape = 21,colour = "black",size = 10)+
                                  scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                                  scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                                  theme_bw() +
                                  ylab("") +
                                  ggtitle(expression("Species richness (q"[0]*")")) +
                                  theme(plot.title = element_text(size=20, hjust=0.5)) +
                                  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                        axis.title.y = element_text(size = 20, margin = margin(r = 12)), axis.title.x = element_blank(), legend.position = "none", 
                                        axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                                  # geom_point(data = ColOdoHet_summary_q0_Sampling, aes(x = Sampling, y = mean_q0.obs), position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), shape = 19, colour = "black", size = 12) +
                                  geom_point(data = ColOdoHet_summary_q0_Sampling, aes(x = Sampling, y = mean_q0.obs), position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), shape = 19, size = 10) +
                                  scale_y_continuous(limits = c(1, 12), breaks = seq(2, 12, by = 2))

print(q0.Treats_Sampling2)

# Note: with these plots, there seems to be Treat difference for Sampling 2 but not for 1, 3, 4. In ii, sub-models are fit per Sampling to test this.

# Arranging q0 individual plot and per-sampling plots (1st version):

ColOdoHet_plot_SpRich_indiv_arr <- readRDS("outputs/Plots/BIO/SpRich_indiv.rds") + # Calling individual plot created in Macroinv_2022 Script.
                                annotate(geom="text", x=0.55, y=1, label="(a)", color="black", size = 10, family = "serif", fontface = "bold")

q0.Treats_Sampling1_arr <- q0.Treats_Sampling1 + 
                            geom_text(data = subset(Hills_Physchem, Sampling == 3), 
                                      aes(label = "(b)", x = 0.7, y = 1.7), color="black", size = 8, family = "serif", fontface = "bold")

Arr.q0_ind.samp <- grid.arrange(ColOdoHet_plot_SpRich_indiv_arr, q0.Treats_Sampling1_arr, nrow = 1)

ggsave("outputs/Plots/BIO/Arr.q0_ind.samp.pdf", plot = Arr.q0_ind.samp ,width = 20, height = 10)

## ii) q0 - GLMM sub-models and ANOVA per Sampling ####

# Subsetting per Sampling:

Hills_Physchem_Samp1 <- subset(Hills_Physchem, Hills_Physchem$Sampling == "1")
Hills_Physchem_Samp2 <- subset(Hills_Physchem, Hills_Physchem$Sampling == "2")
Hills_Physchem_Samp3 <- subset(Hills_Physchem, Hills_Physchem$Sampling == "3")
Hills_Physchem_Samp4 <- subset(Hills_Physchem, Hills_Physchem$Sampling == "4")

# Sampling 1: 

glmm.q0.gaus7_Samp1 <- glmmTMB(data = Hills_Physchem_Samp1, q0.obs ~ Treat + Conduct_microS_cm + pH_soil + Redox_pot + 
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus7_Samp1, plot = T)
summary(glmm.q0.gaus7_Samp1)
car::Anova(glmm.q0.gaus7_Samp1)
r.squaredGLMM(glmm.q0.gaus7_Samp1) # Calculates Pseudo-R-squared for Generalized Mixed-Effect models
performance::r2(glmm.q0.gaus7_Samp1)
performance::check_collinearity(glmm.q0.gaus7_Samp1)
performance::check_singularity(glmm.q0.gaus7_Samp1)
visreg(glmm.q0.gaus7_Samp1, scale="response")
visreg(glmm.q0.gaus7_Samp1, xvar="Sampling2", by = "Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals
visreg(glmm.q0.gaus7_Samp1, xvar="Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals

# Note: No Treat effect detected.

# Sampling 2: 

# Note:For sub-model "glmm.q0.gaus7_Samp2" water physicochemical variables O2_percent and Salinity are removed due to lack of data for most of AWD and MSD
# plots as the Sampling was done after drainage for both treats, so the average of the three previous physicochemical samplings does not apply (no water 
# to assess)

glmm.q0.gaus7_Samp2 <- glmmTMB(data = Hills_Physchem_Samp2, q0.obs ~ Treat + Conduct_microS_cm + pH_soil + Redox_pot + 
                               (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus7_Samp2, plot = T)
summary(glmm.q0.gaus7_Samp2)
car::Anova(glmm.q0.gaus7_Samp2)
r.squaredGLMM(glmm.q0.gaus7_Samp2) # Calculates Pseudo-R-squared for Generalized Mixed-Effect models
performance::r2(glmm.q0.gaus7_Samp2)
performance::check_collinearity(glmm.q0.gaus7_Samp2)
performance::check_singularity(glmm.q0.gaus7_Samp2)
visreg(glmm.q0.gaus7_Samp2, scale="response")
visreg(glmm.q0.gaus7_Samp2, xvar="Sampling2", by = "Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals
visreg(glmm.q0.gaus7_Samp2, xvar="Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals

# Post-hoc test for sub-model glmm.q0.gaus7_Samp2:

# Sampling 3: 
# 
emmeans(glmm.q0.gaus7_Samp3, ~Treat , type = "response")
pairs(emmeans(glmm.q0.gaus7_Samp3, ~Treat , type = "response"))

glmm.q0.gaus7_Samp3 <- glmmTMB(data = Hills_Physchem_Samp3, q0.obs ~ Treat + Conduct_microS_cm + pH_soil + Redox_pot + 
                                 O2_percent + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.q0.gaus7_Samp2, plot = T)
summary(glmm.q0.gaus7_Samp2)
car::Anova(glmm.q0.gaus7_Samp2)
r.squaredGLMM(glmm.q0.gaus7_Samp2) # Calculates Pseudo-R-squared for Generalized Mixed-Effect models
performance::r2(glmm.q0.gaus7_Samp2)
performance::check_collinearity(glmm.q0.gaus7_Samp2)
performance::check_singularity(glmm.q0.gaus7_Samp2)
visreg(glmm.q0.gaus7_Samp2, scale="response")
visreg(glmm.q0.gaus7_Samp2, xvar="Sampling2", by = "Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals
visreg(glmm.q0.gaus7_Samp2, xvar="Treat", overlay = TRUE, type="conditional", scale = "response") # Plotting conditional residuals

# Post-hoc test for sub-model glmm.q0.gaus7_Samp2:

emmeans(glmm.q0.gaus7_Samp2, ~Treat , type = "response")
pairs(emmeans(glmm.q0.gaus7_Samp2, ~Treat , type = "response"))


# Note: Significant Treat effects on q0 identified for pairs: CON - MSD and CON - AWD

## iii) q0 - Comparing overall species richness averages per Treat (for Sampling 2) ####

q0.Sampling2.avg <- Hills_Physchem_Samp2 %>% 
                      group_by(Treat) %>% 
                      summarise(Avg_q0.obs = mean(q0.obs), se_q0.obs = sd(q0.obs) / sqrt(n()))

(q0.Sampling2.avg[1,2]-q0.Sampling2.avg[3,2])/q0.Sampling2.avg[1,2] # Emission % decrease CON vs AWD 
(q0.Sampling2.avg[1,2]-q0.Sampling2.avg[2,2])/q0.Sampling2.avg[1,2] # Emission % decrease CON vs MSD

##### Selected model for q1: Model 20b (not removing outliers) ####

# Same selection criteria as for q0.
# ANOVA results in significant Treat effect, we move forward with post-hoc tests.
# Selected model:
# glmm.q1.gaus8 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Conduct_microS_cm + pH_soil + Redox_pot + 
#                            O2_percent + Salinity + (1|Rep) , family = "gaussian")

emmeans(glmm.q1.gaus8, ~Treat , type = "response")
pairs(emmeans(glmm.q1.gaus8, ~Treat , type = "response"))

# Next steps, from example script:
# emmeans(m.rare, ~sampling_month, type = "response")
# pairs(emmeans(m.rare, ~sampling_month, type = "response"))
# 
# pairs(emmeans(m.rare, ~treatment|sampling_month, type = "response"))
# 
# df.rare <- as.data.frame(emmeans(m.rare, ~treatment|sampling_month, type = "response")) %>% 
#   mutate(divindex = rep("q0.obs",8),
#          emmean = rate) %>% 
#   select(-rate)


##### Abundance Models ####

# Including Rep to each Plot:
Plot <- c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P13", "P14", "P15")
Rep <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5")
Plot_Rep <- data.frame(Plot, Rep)
Abundance_2022 <- merge(Abundance_2022, Plot_Rep, by.x="Plot", by.y="Plot")

##### Model 1: Abund. - Gaussian - Treat*Order_SubOrder interaction - "Rep" Random factor ####
# Family: Gaussian
# Interacting independent variables: Treat*Order_SubOrder
# Random effect: Rep

glmm.abu.gaus1 <- glmmTMB(data = Abundance_2022, Abundance ~ Treat*Order_SubOrder + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.abu.gaus1, plot = T)
summary(glmm.abu.gaus1)
car::Anova(glmm.abu.gaus1)
performance::r2(glmm.abu.gaus1)
performance::check_collinearity(glmm.abu.gaus1)
performance::check_singularity(glmm.abu.gaus1)
visreg(glmm.abu.gaus1, scale="response") # Plotting conditional residuals

##### Model 2: Abund. - Gaussian - Treat*Order_SubOrder interaction - No Random factor ####
# Family: Gaussian
# Interacting independent variables: Treat*Order_SubOrder
# Random effect: -

glmm.abu.gaus2 <- glmmTMB(data = Abundance_2022, Abundance ~ Treat*Order_SubOrder , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.abu.gaus2, plot = T)
summary(glmm.abu.gaus2)
car::Anova(glmm.abu.gaus2)
performance::r2(glmm.abu.gaus2)
performance::check_collinearity(glmm.abu.gaus2)
performance::check_singularity(glmm.abu.gaus2)
visreg(glmm.abu.gaus2, scale="response") # Plotting conditional residuals

emmeans(glmm.abu.gaus2, ~Treat , type = "response")
pairs(emmeans(glmm.abu.gaus2, ~Treat , type = "response"))

emmeans(glmm.abu.gaus2, ~Treat|Order_SubOrder, type = "response")
pairs(emmeans(glmm.abu.gaus2, ~Treat|Order_SubOrder, type = "response"))

