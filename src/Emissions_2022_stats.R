
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
