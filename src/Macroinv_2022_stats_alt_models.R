
####################################################### Thesis Paper 1 - Macroinvertebrates_Statistics #####################################################

library(dplyr)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(AICcmodavg)

##############  1. Preparing base data frames #################

#### 1.1. Working physicochemical data ####

physchem_2022 <- read.csv("data/Other_factors_CERESTRES_2022.csv", fileEncoding="latin1", na.strings=c("","NA"))
physchem_2022$Sampling_date <- as.Date(physchem_2022$Sampling_date)

# To have a single value for each physicochemical variable representing the conditions in which organisms were developed within plots, the average of the three previous
# physicochemical measurements to the biodiversity samplings is calculated:

## Assigning Sampling, Rep and Treat:
Sampling <- c("1", "2", "3", "4")
Treat <- c("AWD", "MSD", "CON", "MSD", "AWD", "CON", "MSD", "CON", "AWD", "AWD", "MSD", "CON", "MSD", "AWD", "CON")
Rep <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5")

## Subseting physchem_2022, only considering three Sampling dates previous to each biodiversity samplings.
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

##############  2. Applying GLMMs #################

#### 2.1. Testing correlations ####

cor_matrix <- Hills_Physchem %>% 
  select(q0.obs, q1.obs, q2.obs,
         Conduct_microS_cm, Temp_10_cm, pH_soil, Redox_pot, Water_temp, O2_percent, O2_mg_l, Salinity, pH_water, Sampling) %>% 
  na.omit

corr <- round(cor(cor_matrix, method =  "spearman"), 1) # Calculates the Spearman rank correlation matrix
## The Spearman rank correlation is a non-parametric measure of association between two variables. 
## It assesses the strength and direction of the monotonic relationship (whether it goes up or down) between two variables. 

pdf("outputs/Plots/BIO/Corr_plot.pdf", width = 11)
corrplot::corrplot.mixed(corr, order = 'hclust', addrect = 2)
dev.off()

# Model for rare species (q0)
hist(Hills_Physchem$q0.obs)

#### 2.2. GLMM - All variables (fixed effects) included ####

Hills_Physchem$Rep_Treat <- paste0(Hills_Physchem$Rep, "_", Hills_Physchem$Treat) # Creates the random effects variable, which is, in this case, controlled by the blocks design.

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
