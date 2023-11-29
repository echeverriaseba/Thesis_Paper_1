
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

# Testing effect of editing vector types:
# Hills_Physchem$Sampling <- as.factor(Hills_Physchem$Sampling)
# Hills_Physchem$Treat <- as.character(Hills_Physchem$Treat)

##### Model 1: Poisson w/interaction ####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.pois1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

DHARMa::simulateResiduals(glmm.q0.pois1, plot = T)
summary(glmm.q0.pois1)
car::Anova(glmm.q0.pois1)
performance::r2(glmm.q0.pois1)
performance::check_collinearity(glmm.q0.pois1)

##### Model 1.2: Poisson w/interaction (VIF>5 vars. removed)####
# Family: Poisson
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.pois1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

DHARMa::simulateResiduals(glmm.q0.pois1, plot = T)
summary(glmm.q0.pois1)
car::Anova(glmm.q0.pois1)
performance::r2(glmm.q0.pois1)
performance::check_collinearity(glmm.q0.pois1)


##### Model 2: Poisson wo/interaction ####
# Family: Poisson
# Interacting independent variables: -
# Additional independent variables: Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.pois2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "poisson")

DHARMa::simulateResiduals(glmm.q0.pois2, plot = T)
summary(glmm.q0.pois2)
car::Anova(glmm.q0.pois2)
performance::r2(glmm.q0.pois2)
performance::check_collinearity(glmm.q0.pois2)

##### Model 3: Negative binomial 1 w/interaction ####
# Family: Negative binomial1
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.nb1.1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "nbinom1")

DHARMa::simulateResiduals(glmm.q0.nb1.1, plot = T)
summary(glmm.q0.nb1.1)
car::Anova(glmm.q0.nb1.1)
performance::r2(glmm.q0.nb1.1)
performance::check_collinearity(glmm.q0.nb1.1)

##### Model 4: Negative binomial 1 wo/interaction ####
# Family: Negative binomial1
# Interacting independent variables: -
# Additional independent variables: Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.nb1.2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "nbinom1")

DHARMa::simulateResiduals(glmm.q0.nb1.2, plot = T)
summary(glmm.q0.nb1.2)
car::Anova(glmm.q0.nb1.2)
performance::r2(glmm.q0.nb1.2)
performance::check_collinearity(glmm.q0.nb1.2)

##### Model 5: Negative binomial 2 w/interaction ####
# Family: Negative binomial1
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.nb2.1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "nbinom2")

DHARMa::simulateResiduals(glmm.q0.nb2.1, plot = T)
summary(glmm.q0.nb2.1)
car::Anova(glmm.q0.nb2.1)
performance::r2(glmm.q0.nb2.1)
performance::check_collinearity(glmm.q0.nb2.1)

##### Model 6: Negative binomial 2 wo/interaction ####
# Family: Negative binomial1
# Interacting independent variables: -
# Additional independent variables: Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.nb2.2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "nbinom2")

DHARMa::simulateResiduals(glmm.q0.nb2.2, plot = T)
summary(glmm.q0.nb2.2)
car::Anova(glmm.q0.nb2.2)
performance::r2(glmm.q0.nb2.2)
performance::check_collinearity(glmm.q0.nb2.2)

##### Model 7: Gaussian w/interaction ####
# Family: Gaussian
# Interacting independent variables: Treat*Sampling
# Additional independent variables: Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.gaus1 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat*Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

DHARMa::simulateResiduals(glmm.q0.gaus1, plot = T)
summary(glmm.q0.gaus1)
car::Anova(glmm.q0.gaus1)
performance::r2(glmm.q0.gaus1)
performance::check_collinearity(glmm.q0.gaus1)

##### Model 8: Gaussian wo/interaction ####
# Family: Gaussian
# Interacting independent variables: -
# Additional independent variables: Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp + O2_percent + O2_mg_l + Salinity + pH_water
# Random effect: Rep_Treat

glmm.q0.gaus2 <- glmmTMB(data = Hills_Physchem, q0.obs ~ Treat + Sampling + Conduct_microS_cm + Temp_10_cm + pH_soil + Redox_pot + Water_temp +
                           O2_percent + O2_mg_l + Salinity + pH_water + (1|Rep_Treat) , family = "gaussian")

DHARMa::simulateResiduals(glmm.q0.gaus2, plot = T)
summary(glmm.q0.gaus2)
car::Anova(glmm.q0.gaus2)
performance::r2(glmm.q0.gaus2)
performance::check_collinearity(glmm.q0.gaus2)

##### Model comparison using AIC ####
models.glmm.q0 <- list(glmm.q0.pois1, glmm.q0.pois2, glmm.q0.nb1.1, glmm.q0.nb1.2, glmm.q0.nb2.1, glmm.q0.nb2.2, glmm.q0.gaus1, glmm.q0.gaus2)
mod.names <- c('glmm.q0.pois1', 'glmm.q0.pois2', 'glmm.q0.nb1.1', 'glmm.q0.nb1.2', 'glmm.q0.nb2.1', 'glmm.q0.nb2.2', 'glmm.q0.gaus1', 'glmm.q0.gaus2')
aictab(cand.set = models.glmm.q0, modnames = mod.names)
