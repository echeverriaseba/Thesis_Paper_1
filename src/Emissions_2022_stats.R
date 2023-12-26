
####################################################### Thesis Paper 1 - GHG_Emissions_Statistics #####################################################

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

################ 1. GHG_emissions-physchem data frame ###################
# "Master_GHG_2022_no_NA" dataframe already contains "CH4_flux_corrected" and physicochemical data and excludes all NA values (due to non GHG measurement dates, only physchem). 
# Dataframe called from the "Emissions_2022" Script.

# source("Emissions_2022.R")

################ 2. STATS ###################

##### 2.1. Data validation ######
# According to "Analyzing the impact of multiple stressors in aquatic biomonitoring data: A ‘cookbook’ with applications in R" - Feld et al., 2016. ##

#### 2.1.1. Outline analysis ####

summary(Master_GHG_2022_no_NA$Temp_soil)
boxplot(Master_GHG_2022_no_NA$Temp_soil)
boxplot(Master_GHG_2022_no_NA$Temp_soil)$out
max(boxplot(Master_GHG_2022_no_NA$Temp_soil)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Temp_soil)$out) # prints the minimum outlier value
## Outliers identified for var: Temp_soil - 0 and 38 

summary(Master_GHG_2022_no_NA$Rice_cover_prop)
boxplot(Master_GHG_2022_no_NA$Rice_cover_prop)
boxplot(Master_GHG_2022_no_NA$Rice_cover_prop)$out
max(boxplot(Master_GHG_2022_no_NA$Rice_cover_prop)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Rice_cover_prop)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Env_temp_initial)
boxplot(Master_GHG_2022_no_NA$Env_temp_initial)
boxplot(Master_GHG_2022_no_NA$Env_temp_initial)$out
max(boxplot(Master_GHG_2022_no_NA$Env_temp_initial)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Env_temp_initial)$out) # prints the minimum outlier value
## Outliers identified for var: Env_temp_initial - 232

summary(Master_GHG_2022_no_NA$Env_temp_final)
boxplot(Master_GHG_2022_no_NA$Env_temp_final)
boxplot(Master_GHG_2022_no_NA$Env_temp_final)$out
max(boxplot(Master_GHG_2022_no_NA$Env_temp_final)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Env_temp_final)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Water_level_corr)
boxplot(Master_GHG_2022_no_NA$Water_level_corr)
boxplot(Master_GHG_2022_no_NA$Water_level_corr)$out
max(boxplot(Master_GHG_2022_no_NA$Water_level_corr)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Water_level_corr)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Conduct_microS_cm)
boxplot(Master_GHG_2022_no_NA$Conduct_microS_cm)
boxplot(Master_GHG_2022_no_NA$Conduct_microS_cm)$out
max(boxplot(Master_GHG_2022_no_NA$Conduct_microS_cm)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Conduct_microS_cm)$out) # prints the minimum outlier value
## Outliers identified for var: Conduct_microS_cm - 0.46 and 2.24

summary(Master_GHG_2022_no_NA$Temp_10_cm)
boxplot(Master_GHG_2022_no_NA$Temp_10_cm)
boxplot(Master_GHG_2022_no_NA$Temp_10_cm)$out
max(boxplot(Master_GHG_2022_no_NA$Temp_10_cm)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Temp_10_cm)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$pH_soil)
boxplot(Master_GHG_2022_no_NA$pH_soil)
boxplot(Master_GHG_2022_no_NA$pH_soil)$out
max(boxplot(Master_GHG_2022_no_NA$pH_soil)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$pH_soil)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Redox_pot)
boxplot(Master_GHG_2022_no_NA$Redox_pot)
boxplot(Master_GHG_2022_no_NA$Redox_pot)$out
max(boxplot(Master_GHG_2022_no_NA$Redox_pot)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Redox_pot)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Water_temp)
boxplot(Master_GHG_2022_no_NA$Water_temp)
boxplot(Master_GHG_2022_no_NA$Water_temp)$out
max(boxplot(Master_GHG_2022_no_NA$Water_temp)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Water_temp)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$O2_percent)
boxplot(Master_GHG_2022_no_NA$O2_percent)
boxplot(Master_GHG_2022_no_NA$O2_percent)$out
max(boxplot(Master_GHG_2022_no_NA$O2_percent)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$O2_percent)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$O2_mg_l)
boxplot(Master_GHG_2022_no_NA$O2_mg_l)
boxplot(Master_GHG_2022_no_NA$O2_mg_l)$out
max(boxplot(Master_GHG_2022_no_NA$O2_mg_l)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$O2_mg_l)$out) # prints the minimum outlier value

summary(Master_GHG_2022_no_NA$Salinity)
boxplot(Master_GHG_2022_no_NA$Salinity)
boxplot(Master_GHG_2022_no_NA$Salinity)$out
max(boxplot(Master_GHG_2022_no_NA$Salinity)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$Salinity)$out) # prints the minimum outlier value
## Outliers identified for var: Salinity - 0.37 and 1.11

summary(Master_GHG_2022_no_NA$pH_water)
boxplot(Master_GHG_2022_no_NA$pH_water)
boxplot(Master_GHG_2022_no_NA$pH_water)$out
max(boxplot(Master_GHG_2022_no_NA$pH_water)$out) # prints the maximum outlier value 
min(boxplot(Master_GHG_2022_no_NA$pH_water)$out) # prints the minimum outlier value

## From here onwards tests and models will be evaluated both considering outliers and removing them (suffix: "_nooutliers")

Master_GHG_2022_no_NA_nooutliers <-  Master_GHG_2022_no_NA %>%  # New data frame with a (1/0) Outliers column, then removing outlier rows
  mutate(Outliers =case_when(Temp_soil > 37  | Temp_soil < 5 | Conduct_microS_cm < 0.5 |  Conduct_microS_cm > 2.2 | Env_temp_initial > 230 | 
                               Salinity < 0.4 | Salinity > 1.1 ~ 1,TRUE ~ 0)) %>% 
                               filter(Outliers == 0)

#### 2.1.2. Testing correlations ####

#### i) Spearman rank correlation ####

## The Spearman rank correlation is a non-parametric measure of association between two variables. 
## It assesses the strength and direction of the monotonic relationship (whether it goes up or down) between two variables. 

### Considering outliers:

cor_matrix_GHG <- Master_GHG_2022_no_NA %>% 
  select(CH4_flux_corrected, Water_level_corr, Temp_soil, Rice_cover_prop, Env_temp_initial, Env_temp_final,
         Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water) %>% 
  na.omit

corr_GHG <- round(cor(cor_matrix_GHG, method =  "spearman"), 1) # Calculates the Spearman rank correlation matrix

pdf("outputs/Plots/GHG/Corr_plot_GHG.pdf", width = 11)
corrplot::corrplot.mixed(corr_GHG, order = 'hclust', addrect = 2)
dev.off()

### Removing outliers:

cor_matrix_GHG_nooutliers <- Master_GHG_2022_no_NA_nooutliers %>% 
                              select(CH4_flux_corrected, Water_level_corr, Temp_soil, Rice_cover_prop, Env_temp_initial, Env_temp_final,
                                     Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water) %>% 
                              na.omit

corr_GHG_nooutliers <- round(cor(cor_matrix_GHG_nooutliers, method =  "spearman"), 1) # Calculates the Spearman rank correlation matrix

pdf("outputs/Plots/GHG/Corr_plot_GHG_nooutliers.pdf", width = 11)
corrplot::corrplot.mixed(corr_GHG_nooutliers, order = 'hclust', addrect = 2)
dev.off()

#### ii) Variance inflation factors (VIF) ####

### Method: Excluding variables with VIF > 5 stepwise, starting with the variable that has the highest VIF
### A VIF > 8 is applied in Feld et al., 2016. We lower the threshold as pH_water resulted in High collinearity in posterior tested GLMMs.

###  Considering outliers:

# Step 1a: all variables considered
selected_vars_GHG <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
                       'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_GHG  <- (Master_GHG_2022_no_NA[, selected_vars_GHG])
vif(data_selected_GHG,) 

# Step 2a: Removing variable with highest VIF (and VIF>5) - O2_mg_l (VIF = 22.642541)
selected_vars_GHG_2 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
                       'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_GHG_2  <- (Master_GHG_2022_no_NA[, selected_vars_GHG_2])
vif(data_selected_GHG_2,) 

# Step 3a: Removing next variable with highest VIF (and VIF>5) - Env_temp_initial  (VIF = 13.342621)
selected_vars_GHG_3 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_final',
                         'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'Salinity', 'pH_water')
data_selected_GHG_3  <- (Master_GHG_2022_no_NA[, selected_vars_GHG_3])
vif(data_selected_GHG_3,) 

#Note: Candidate variables for removal from posterior analyses:  O2_mg_l, Env_temp_initial.

### Removing outliers:

# Step 1b: all variables considered
selected_vars_GHG_4 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
                       'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_GHG_4  <- (Master_GHG_2022_no_NA_nooutliers[, selected_vars_GHG_4])
vif(data_selected_GHG_4,) 

# Step 2b: Removing variable with highest VIF (and VIF>5) - O2_percent (VIF = 44.548571)
selected_vars_GHG_5 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
                         'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_GHG_5  <- (Master_GHG_2022_no_NA_nooutliers[, selected_vars_GHG_5])
vif(data_selected_GHG_5,) 

# Step 3b: Removing next variable with highest VIF (and VIF>5) - Env_temp_initial   (VIF = 12.771963)
selected_vars_GHG_6 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_final',
                         'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_GHG_6  <- (Master_GHG_2022_no_NA_nooutliers[, selected_vars_GHG_6])
vif(data_selected_GHG_6,) 

# Step 4b Removing next variable with highest VIF (and VIF>5) - Temp_10_cm   (VIF = 5.357361)
selected_vars_GHG_7 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_final',
                         'Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity', 'pH_water')
data_selected_GHG_7  <- (Master_GHG_2022_no_NA_nooutliers[, selected_vars_GHG_7])
vif(data_selected_GHG_7,) 

#Note: Candidate variables for removal from posterior analyses:  O2_percent, Env_temp_initial, Temp_10_cm.

#### iii) Data Dendrogram ####

###  Considering outliers:

# Dend 1: Considering all variables (without VIF removal)
# Variable clustering:
similarity="pearson"
vclust_GHG_1 <- varclus(x=as.matrix(data_selected_GHG),
                    similarity=similarity,
                    type="data.matrix", 
                    method="complete",
                    na.action=na.retain,trans="abs")
dend_GHG_1 <- as.dendrogram(vclust_GHG_1)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_1 <- color_labels(dend_GHG_1, h=1-0.7,col=col_vector)
dend_GHG_1 <- color_branches(dend_GHG_1, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_1.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_1,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Dend 2: Considering only remaining variables after VIF removal
# Variable clustering:
similarity="pearson"
vclust_GHG_2 <- varclus(x=as.matrix(data_selected_GHG_3),
                    similarity=similarity,
                    type="data.matrix", 
                    method="complete",
                    na.action=na.retain,trans="abs")
dend_GHG_2 <- as.dendrogram(vclust_GHG_2)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_2 <- color_labels(dend_GHG_2, h=1-0.7,col=col_vector)
dend_GHG_2 <- color_branches(dend_GHG_2, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_2_VIFvarsremoval.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_2,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Note: Dend 2 (removing variables after the VIF analysis) does not achieve that all variables stay underneath a 0.75 threshold. 
# We proceed to remove Env_temp_final and pH_water and re-check with Dend 3

# Dend 3: Considering only remaining variables after VIF removal and removing correlated vars. acording to Dend 2

selected_vars_GHG_8 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop',
                         'Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity')
data_selected_GHG_8 <- (Master_GHG_2022_no_NA[, selected_vars_GHG_8])

# Variable clustering:
similarity="pearson"
vclust_GHG_3 <- varclus(x=as.matrix(data_selected_GHG_8),
                        similarity=similarity,
                        type="data.matrix", 
                        method="complete",
                        na.action=na.retain,trans="abs")
dend_GHG_3 <- as.dendrogram(vclust_GHG_3)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_3 <- color_labels(dend_GHG_3, h=1-0.7,col=col_vector)
dend_GHG_3 <- color_branches(dend_GHG_3, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_2_VIFvarsremoval_2.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_3,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Note: Remaining (non-correlated) variables after VIF and dendrogram analyses not removing outliers: 
# 'Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity'

###  Removing outliers:

# Dend 4: Considering all variables (without VIF removal)
# Variable clustering:
similarity="pearson"
vclust_GHG_4 <- varclus(x=as.matrix(data_selected_GHG_4),
                        similarity=similarity,
                        type="data.matrix", 
                        method="complete",
                        na.action=na.retain,trans="abs")
dend_GHG_4 <- as.dendrogram(vclust_GHG_4)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_4 <- color_labels(dend_GHG_4, h=1-0.7,col=col_vector)
dend_GHG_4 <- color_branches(dend_GHG_4, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_nooutliers.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_4,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Dend 5: Considering only remaining variables after VIF removal
# Variable clustering:
similarity="pearson"
vclust_GHG_5 <- varclus(x=as.matrix(data_selected_GHG_7),
                        similarity=similarity,
                        type="data.matrix", 
                        method="complete",
                        na.action=na.retain,trans="abs")
dend_GHG_5 <- as.dendrogram(vclust_GHG_5)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_5 <- color_labels(dend_GHG_5, h=1-0.7,col=col_vector)
dend_GHG_5 <- color_branches(dend_GHG_5, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_VIFvarsremoval_nooutliers.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_5,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Note: Dend 5 (removing variables after the VIF analysis) does not achieve that all variables stay underneath a 0.75 threshold. 
# We proceed to remove pH_water and re-check with Dend 6

# Dend 6: Considering only remaining variables after VIF removal and removing correlated vars. acording to Dend 5

selected_vars_GHG_9 <- c('Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_final',
                         'Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity')
data_selected_GHG_9  <- (Master_GHG_2022_no_NA_nooutliers[, selected_vars_GHG_9])

# Variable clustering:
similarity="pearson"
vclust_GHG_6 <- varclus(x=as.matrix(data_selected_GHG_9),
                        similarity=similarity,
                        type="data.matrix", 
                        method="complete",
                        na.action=na.retain,trans="abs")
dend_GHG_6 <- as.dendrogram(vclust_GHG_6)
# Random colors and plot dendrogram:
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- brewer.pal(n=8,"Paired")
dend_GHG_6 <- color_labels(dend_GHG_6, h=1-0.7,col=col_vector)
dend_GHG_6 <- color_branches(dend_GHG_6, h=1-0.7,col=col_vector)
cairo_pdf("outputs/Plots/GHG/Cluster_variables_pearson_GHG_VIFvarsremoval_nooutliers2.pdf",width=7,height=4)
par(mar=c(5,2,4,17)+0.1)
plot(dend_GHG_6,horiz = TRUE,xlab="",axes = FALSE)
axis(1,at=1-seq(0,1,0.2),labels=seq(0,1,0.2))
dev.off()

# Note: Remaining (non-correlated) variables after VIF and dendrogram analyses removing outliers: 
# 'Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_final', 'Conduct_microS_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_mg_l', 'Salinity'


#### 2.2. GLMMs ####

# Creating a "Sampling" column that assigns a number to each unique Sampling_date. This way we can have "Sampling" as a model variable. 
Master_GHG_2022_no_NA <- Master_GHG_2022_no_NA %>%
                          arrange(Sampling_date) %>%
                          mutate(Sampling = match(Sampling_date, unique(Sampling_date)))

Master_GHG_2022_no_NA_nooutliers <- Master_GHG_2022_no_NA_nooutliers %>%
                          arrange(Sampling_date) %>%
                          mutate(Sampling = match(Sampling_date, unique(Sampling_date)))

##### Model 1: Gaussian - Treat*Sampling interaction - "Rep" Random factor - Considering all original variables (not considering all prev. correlation analysis) ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: 'Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
# 'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water'
# Random effect: Rep

###  Considering outliers:

glmm.CH4.gaus1 <- glmmTMB(data = Master_GHG_2022_no_NA, CH4_flux_corrected ~ Treat*Sampling + Water_level_corr + Temp_soil + Rice_cover_prop + 
                          Env_temp_initial + Env_temp_final + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + 
                          O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus1, plot = T)
summary(glmm.CH4.gaus1)
car::Anova(glmm.CH4.gaus1)
performance::r2(glmm.CH4.gaus1)
performance::check_collinearity(glmm.CH4.gaus1)
performance::check_singularity(glmm.CH4.gaus1)
visreg(glmm.CH4.gaus1, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.CH4.gaus1.noout <- glmmTMB(data = Master_GHG_2022_no_NA_nooutliers, CH4_flux_corrected ~ Treat*Sampling + Water_level_corr + Temp_soil + Rice_cover_prop + 
                            Env_temp_initial + Env_temp_final + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + 
                            O2_mg_l + Salinity + pH_water + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus1.noout, plot = T)
summary(glmm.CH4.gaus1.noout)
car::Anova(glmm.CH4.gaus1.noout)
performance::r2(glmm.CH4.gaus1.noout)
performance::check_collinearity(glmm.CH4.gaus1.noout)
performance::check_singularity(glmm.CH4.gaus1.noout)
visreg(glmm.CH4.gaus1.noout, scale="response") # Plotting conditional residuals

##### Model 2: Gaussian - Treat*Sampling interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: 'Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
# 'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water'
# Random effect: Rep

###  Considering outliers:

glmm.CH4.gaus2 <- glmmTMB(data = Master_GHG_2022_no_NA, CH4_flux_corrected ~ Treat*Sampling + Water_level_corr + Temp_soil + Rice_cover_prop + 
                            Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_mg_l + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus2, plot = T)
summary(glmm.CH4.gaus2)
car::Anova(glmm.CH4.gaus2)
performance::r2(glmm.CH4.gaus2)
performance::check_collinearity(glmm.CH4.gaus2)
performance::check_singularity(glmm.CH4.gaus2)
visreg(glmm.CH4.gaus2, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.CH4.gaus2.noout <- glmmTMB(data = Master_GHG_2022_no_NA_nooutliers, CH4_flux_corrected ~ Treat*Sampling + Water_level_corr + Temp_soil + 
                                  Env_temp_final + Rice_cover_prop + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_mg_l + Salinity + 
                                  (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus2.noout, plot = T)
summary(glmm.CH4.gaus2.noout)
car::Anova(glmm.CH4.gaus2.noout)
performance::r2(glmm.CH4.gaus2.noout)
performance::check_collinearity(glmm.CH4.gaus2.noout)
performance::check_singularity(glmm.CH4.gaus2.noout)
visreg(glmm.CH4.gaus2.noout, scale="response") # Plotting conditional residuals

##### Model 3: Gaussian - Treat*Sampling and Treat*I(Sampling^2) interaction - "Rep" Random factor - Considering only remaining variables after correlation analysis ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: 'Water_level_corr', 'Temp_soil', 'Rice_cover_prop', 'Env_temp_initial', 'Env_temp_final',
# 'Conduct_microS_cm', 'Temp_10_cm', 'pH_soil', 'Redox_pot', 'Water_temp', 'O2_percent', 'O2_mg_l', 'Salinity', 'pH_water'
# Random effect: Rep

###  Considering outliers:

glmm.CH4.gaus3 <- glmmTMB(data = Master_GHG_2022_no_NA, CH4_flux_corrected ~ Treat*Sampling + Treat*I(Sampling^2) + Water_level_corr + Temp_soil + Rice_cover_prop + 
                            Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_mg_l + Salinity + (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus3, plot = T)
summary(glmm.CH4.gaus3)
car::Anova(glmm.CH4.gaus3)
performance::r2(glmm.CH4.gaus3)
performance::check_collinearity(glmm.CH4.gaus3)
performance::check_singularity(glmm.CH4.gaus3)
visreg(glmm.CH4.gaus3, scale="response") # Plotting conditional residuals

###  Removing outliers:

glmm.CH4.gaus3.noout <- glmmTMB(data = Master_GHG_2022_no_NA_nooutliers, CH4_flux_corrected ~ Treat*Sampling + Treat*I(Sampling^2) + Water_level_corr + Temp_soil + 
                                  Env_temp_final + Rice_cover_prop + Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_mg_l + Salinity + 
                                  (1|Rep) , family = "gaussian")

# Model diagnostics:
DHARMa::simulateResiduals(glmm.CH4.gaus3.noout, plot = T)
summary(glmm.CH4.gaus3.noout)
car::Anova(glmm.CH4.gaus3.noout)
performance::r2(glmm.CH4.gaus3.noout)
performance::check_collinearity(glmm.CH4.gaus3.noout)
performance::check_singularity(glmm.CH4.gaus3.noout)
visreg(glmm.CH4.gaus3.noout, scale="response") # Plotting conditional residuals


##### Selected model: Model 2 (removing outliers) ####

emmeans(glmm.CH4.gaus2.noout, ~Treat , type = "response")
pairs(emmeans(glmm.CH4.gaus2.noout, ~Treat , type = "response"))

# Check if the following makes sense:

emmeans(glmm.CH4.gaus2.noout, ~Sampling, type = "response")
pairs(emmeans(glmm.CH4.gaus2.noout, ~Sampling, type = "response"))

pairs(emmeans(glmm.CH4.gaus2.noout, ~Treat|Sampling, type = "response"))