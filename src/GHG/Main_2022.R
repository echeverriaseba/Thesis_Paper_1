##############################################################################################################################################
################################################### Chromatography Results - 2022 ############################################################
##############################################################################################################################################

library(tidyverse)
library(dplyr)
library(ggrepel)
library(ggpmisc)
library(ggplot2)
library(ggpubr)
library(writexl)

##### 1. Working raw data #######

# First, a csv file merging Joan's raw data and sample references has to be created as initial input (see "Field_sheet_chrom_2022").

Chrom_results_2022 <- read.csv("data/Field_sheet_chrom_2022.csv", fileEncoding="latin1", na.strings=c("","NA"))
summary(Chrom_results_2022)
head(Chrom_results_2022)

##### 2. Replicating "GHG_matrix" (book from previous "Post-harvest field experiment data base" file) calculations ####

Chrom_results_2022$Chamber_temp_k <- Chrom_results_2022$Chamber_temp + 273 ## Adding column "Chamber_temp_k".
Chrom_results_2022$Surface_Area <- 0.129 ## Adding column "Surface_Area" : Be careful to edit this in case different chamber types are used.
Chrom_results_2022$Chamber_Height <- 0.72 ## Adding column "Chamber_Height" : Be careful to edit this in case different chamber types are used.
Chrom_results_2022$Volume <- Chrom_results_2022$Surface_Area * Chrom_results_2022$Chamber_Height ## Adding column "Volume".
Chrom_results_2022$C_density_g_m3 <- (12 / (82.0575 * Chrom_results_2022$Chamber_temp_k)) * 1000000
Chrom_results_2022$N_density_g_m3 <- (14 / (82.0575 * Chrom_results_2022$Chamber_temp_k)) * 1000000
Chrom_results_2022$CH4_byMass_mgm3 <- (Chrom_results_2022$C_density_g_m3 * Chrom_results_2022$CCH4_ppm) / 1000
Chrom_results_2022$CH4_byMass_mgm2 <- (Chrom_results_2022$CH4_byMass_mgm3 * Chrom_results_2022$Volume) / Chrom_results_2022$Surface_Area
Chrom_results_2022$N2O_byMass_mgm3 <- (Chrom_results_2022$N_density_g_m3 * Chrom_results_2022$NN2O_ppm) / 1000
Chrom_results_2022$N2O_byMass_mgm2 <- (Chrom_results_2022$N2O_byMass_mgm3 * Chrom_results_2022$Volume) / Chrom_results_2022$Surface_Area
Chrom_results_2022$CO2_byMass_mgm3 <- (Chrom_results_2022$C_density_g_m3 * Chrom_results_2022$CCO2_ppm) / 1000
Chrom_results_2022$CO2_byMass_mgm2 <- (Chrom_results_2022$CO2_byMass_mgm3 * Chrom_results_2022$Volume) / Chrom_results_2022$Surface_Area

#### 3. Creating new data frame to calculate rates and R2 ####

Emission_rates_2022 <- Chrom_results_2022 %>% distinct(Chrom_results_2022$Code, Chrom_results_2022$Sampling_date, Chrom_results_2022$Plot, Chrom_results_2022$Tr1, Chrom_results_2022$Rep) ## Creates data frame with unique values for certain columns
Emission_rates_2022 <- cbind(empty_column1=NA,Emission_rates_2022,empty_column2=NA, empty_column3=NA, empty_column4=NA, empty_column5=NA, empty_column6=NA, empty_column7=NA, empty_column8=NA, empty_column9=NA, empty_column10=NA) 
                             
#                              empty_column11=NA, empty_column12=NA, empty_column13=NA, empty_column14=NA, empty_column15=NA, empty_column16=NA, empty_column17=NA, empty_column18=NA, empty_column19=NA, empty_column20=NA, empty_column21=NA, empty_column22=NA, empty_column23=NA, empty_column24=NA, empty_column25=NA, empty_column26=NA, empty_column27=NA, empty_column28=NA, empty_column29=NA, empty_column30=NA, empty_column31=NA, empty_column23=NA, empty_column33=NA, empty_column34=NA, empty_column35=NA) ## Adds empty columns

colnames(Emission_rates_2022) = c("Code_Nr", "Code", "Sampling_date", "Plot", "Treat", "Rep", "CH4_flux_mgm2h", "R2_CH4", "p_CH4", "N2O_flux_mgm2h", "R2_N2O", "p_N2O", "CO2_flux_mgm2h", "R2_CO2", "p_CO2")
  
  # "Code_Nr", "Code", "Sampling_date", "Tr1", "Tr2", "Rep", "CH4_flux_mgm2h", "R2_CH4", "p_CH4", "Case_CH4", "CH4_flux_Alt1", "R2_CH4_Alt1", "p_CH4_Alt1", "CH4_flux_Alt2", "R2_CH4_Alt2", "p_CH4_Alt2", "CH4_flux_Alt3", "R2_CH4_Alt3", "p_CH4_Alt3", "CH4_flux_Alt4", "R2_CH4_Alt4", "p_CH4_Alt4", "CH4_flux_Log", "R2_CH4_Log", "p_CH4_Log", "CH4_model", "CH4_flux_corrected", "R2_CH4_corrected", "N2O_flux_mgm2h", "R2_N2O", "p_N2O", "Case_N2O", "N2O_flux_corrected", "R2_N2O_corrected", "CO2_flux_mgm2h", "R2_CO2", "p_CO2", "Case_CO2", "CO2_flux_corrected", "R2_CO2_corrected") ## Renames empty columns
Emission_rates_2022$Code_Nr <- 1:nrow(Emission_rates_2022)

#### 4. Loops to fill up new data frame (from calculations in Script "PH_GHG_test") ####

#### 4.1 Loops for CH4 ####

pdf('outputs/2022/CH4_Plots_original.pdf')

for (i in 1:length(Emission_rates_2022$Code)) {
  Code_i <- Emission_rates_2022$Code[i]
  Filt_i <- filter(Chrom_results_2022,Chrom_results_2022$Code == Code_i) # if returned as Time-Series, re-run library(dplyr)
  
  ## Loop section 1: Rate calculation.
  lm_i <- lm(CH4_byMass_mgm2~Sample_time_min, data=Filt_i)
  Emission_rates_2022$CH4_flux_mgm2h[i] <- coef(lm_i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  Emission_rates_2022$R2_CH4[i] <- summary(lm_i)$r.squared # Returns R2_CH4 for each Code.
  lmp_i <- function (modelobject) {  # Function created to call later the model's p-value
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  Emission_rates_2022$p_CH4[i] <- lmp_i(lm_i) # Returns p-value for each Code.
  
  # ## Loop section 2: Rate correction.
  # #### Loop section 2.1: Identifying Graph Cases:
  # Emission_rates_2022$Case_CH4[i] <- if_else((Filt_i$CH4_byMass_mgm2[1] > Filt_i$CH4_byMass_mgm2[2]) & (Filt_i$CH4_byMass_mgm2[2] < Filt_i$CH4_byMass_mgm2[3]) & (Filt_i$CH4_byMass_mgm2[3] >                                  Filt_i$CH4_byMass_mgm2[4]), 1, ## Case 1
  #                                       if_else((Filt_i$CH4_byMass_mgm2[1] < Filt_i$CH4_byMass_mgm2[2]) & (Filt_i$CH4_byMass_mgm2[2] > Filt_i$CH4_byMass_mgm2[3]) & (Filt_i$CH4_byMass_mgm2[3] >                                  Filt_i$CH4_byMass_mgm2[4]), 2,  ## Case 2
  #                                               0))
  # 
  # ## Loop section 2.2: Fitting 4 alternative "3-values" models (each one removing one time-step) and a Log model:
  # #### Alt_1: excluding concentration from time step 0 mins
  # Filt_Alt1i <- if(is.na(Filt_i$CH4_byMass_mgm2[1]) == TRUE) {Filt_i} else {Filt_i[-1,]}# Filters excluding one concentration. In case there are NA it doesn't exclude values.
  # lm_Alt1i <- lm(CH4_byMass_mgm2~Sample_time_min, data=Filt_Alt1i) # Linear model of these 3 values
  # Emission_rates_2022$CH4_flux_Alt1[i] <- coef(lm_Alt1i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  # Emission_rates_2022$R2_CH4_Alt1[i] <- summary(lm_Alt1i)$r.squared # Returns R2_CH4 for each Code.
  # Emission_rates_2022$p_CH4_Alt1[i] <- lmp_i(lm_Alt1i) # Returns p-value for each Code.
  # 
  # #### Alt_2: excluding concentration from time step 10 mins
  # Filt_Alt2i <- if(is.na(Filt_i$CH4_byMass_mgm2[2]) == TRUE) {Filt_i} else {Filt_i[-2,]}# Filters excluding one concentration. In case the value to be excluded is NA it doesn't exclude values
  # lm_Alt2i <- lm(CH4_byMass_mgm2~Sample_time_min, data=Filt_Alt2i) # Linear model of these 3 values
  # Emission_rates_2022$CH4_flux_Alt2[i] <- coef(lm_Alt2i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  # Emission_rates_2022$R2_CH4_Alt2[i] <- summary(lm_Alt2i)$r.squared # Returns R2_CH4 for each Code.
  # Emission_rates_2022$p_CH4_Alt2[i] <- lmp_i(lm_Alt2i) # Returns p-value for each Code.
  # 
  # #### Alt_3: excluding concentration from time step 20 mins
  # Filt_Alt3i <- if(is.na(Filt_i$CH4_byMass_mgm2[3]) == TRUE) {Filt_i} else {Filt_i[-3,]}# Filters excluding one concentration. In case the value to be excluded is NA it doesn't exclude values
  # lm_Alt3i <- lm(CH4_byMass_mgm2~Sample_time_min, data=Filt_Alt3i) # Linear model of these 3 values
  # Emission_rates_2022$CH4_flux_Alt3[i] <- coef(lm_Alt3i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  # Emission_rates_2022$R2_CH4_Alt3[i] <- summary(lm_Alt3i)$r.squared # Returns R2_CH4 for each Code.
  # Emission_rates_2022$p_CH4_Alt3[i] <- lmp_i(lm_Alt3i) # Returns p-value for each Code.
  # 
  # #### Alt_4: excluding concentration from time step 30 mins
  # Filt_Alt4i <- if(is.na(Filt_i$CH4_byMass_mgm2[4]) == TRUE) {Filt_i} else {Filt_i[-4,]}# Filters excluding one concentration. In case the value to be excluded is NA it doesn't exclude values
  # lm_Alt4i <- lm(CH4_byMass_mgm2~Sample_time_min, data=Filt_Alt4i) # Linear model of these 3 values
  # Emission_rates_2022$CH4_flux_Alt4[i] <- coef(lm_Alt4i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  # Emission_rates_2022$R2_CH4_Alt4[i] <- summary(lm_Alt4i)$r.squared # Returns R2_CH4 for each Code.
  # Emission_rates_2022$p_CH4_Alt4[i] <- lmp_i(lm_Alt4i) # Returns p-value for each Code.
  # 
  # #### Log Model:
  # Log_mod_i<- lm(CH4_byMass_mgm2~log(Sample_time_min+1), data=Filt_i)
  # Emission_rates_2022$CH4_flux_Log[i] <- coef(Log_mod_i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
  # Emission_rates_2022$R2_CH4_Log[i] <- summary(Log_mod_i)$r.squared # Returns R2_CH4 for each Code.
  # Emission_rates_2022$p_CH4_Log[i] <- lmp_i(Log_mod_i) # Returns p-value for each Code.
  # 
  # 
  # ## Loop section 2.2: Including restrictions in the loop with nested if_else:
  # Emission_rates_2022$CH4_flux_corrected[i] <- if((Emission_rates_2022$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {0}
  # else if(Emission_rates_2022$R2_CH4[i] > 0.7) {Emission_rates_2022$CH4_flux_mgm2h[i]}
  # else if((Emission_rates_2022$R2_CH4[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt1[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt2[i] < 0.7) &                                                             (Emission_rates_2022$R2_CH4_Alt3[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt4[i] < 0.7)) {0}
  # else if((Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &
  #         (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt4[i]))                                                               {Emission_rates_2022$CH4_flux_mgm2h[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt2[i]) & (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {Emission_rates_2022$CH4_flux_Alt1[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {Emission_rates_2022$CH4_flux_Alt2[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {Emission_rates_2022$CH4_flux_Alt3[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {Emission_rates_2022$CH4_flux_Alt4[i]}
  # else {0}
  # 
  # ## Column with chosen model:
  # Emission_rates_2022$CH4_model[i] <- if((Emission_rates_2022$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {"No flux"}
  # else if(Emission_rates_2022$R2_CH4[i] > 0.7) {"Original"}
  # else if((Emission_rates_2022$R2_CH4[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt1[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt2[i] < 0.7) &                                                             (Emission_rates_2022$R2_CH4_Alt3[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt4[i] < 0.7)) {"No flux"}
  # else if((Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &
  #         (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt4[i]))                                                               {"Original"}
  # else if ((Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt2[i]) & (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {"Alt. 1"}
  # else if ((Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {"Alt. 2"}
  # else if ((Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {"Alt. 3"}
  # else if ((Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {"Alt. 4"}
  # else {"No flux"}
  # 
  # ## Loop section 2.3: Calculating R2 according to the applied correction (if any):
  # Emission_rates_2022$R2_CH4_corrected[i] <- if((Emission_rates_2022$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {0}
  # else if(Emission_rates_2022$R2_CH4[i] > 0.7) {Emission_rates_2022$R2_CH4[i]}
  # else if((Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &
  #         (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt4[i]))                                                               {Emission_rates_2022$R2_CH4[i]}
  # else if((Emission_rates_2022$R2_CH4[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt1[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt2[i] < 0.7) &                                                             (Emission_rates_2022$R2_CH4_Alt3[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt4[i] < 0.7)) {0}
  # else if ((Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt2[i]) & (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {Emission_rates_2022$R2_CH4_Alt1[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {Emission_rates_2022$R2_CH4_Alt2[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {Emission_rates_2022$R2_CH4_Alt3[i]}
  # else if ((Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {Emission_rates_2022$R2_CH4_Alt4[i]}
  # else {0}
  # 
  # ## Add column qith method and rate selection logic:
  # Emission_rates_2022$Logic_CH4[i] <- if((Emission_rates_2022$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {"Original model has negative rate and R2 < 0.7"}
  # else if(Emission_rates_2022$R2_CH4[i] > 0.7) {"Original model has R2 > 0.7"}
  # else if((Emission_rates_2022$R2_CH4[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt1[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt2[i] < 0.7) &                                                             (Emission_rates_2022$R2_CH4_Alt3[i] < 0.7) & (Emission_rates_2022$R2_CH4_Alt4[i] < 0.7)) {"No model achieves  R2 > 0.7"}
  # else if((Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &
  #         (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (Emission_rates_2022$R2_CH4[i] > Emission_rates_2022$R2_CH4_Alt4[i]))                                                               {"Original model achieves the highest R2"}
  # else if ((Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt2[i]) & (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt1[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {"Original model has R2 < 0.7 and Alt. 1 achieves the highest R2 (> 0.7) w/positive rate"}
  # else if ((Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt3[i]) &                                           (Emission_rates_2022$R2_CH4_Alt2[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {"Original model has R2 < 0.7 and Alt. 2 achieves the highest R2 (> 0.7) w/positive rate"}
  # else if ((Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt3[i] > Emission_rates_2022$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {"Original model has R2 < 0.7 and Alt. 3 achieves the highest R2 (> 0.7) w/positive rate"}
  # else if ((Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt1[i]) & (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt2[i]) &                                           (Emission_rates_2022$R2_CH4_Alt4[i] > Emission_rates_2022$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {"Original model has R2 < 0.7 and Alt. 4 achieves the highest R2 (> 0.7) w/positive rate"}
  # else {"Alternative achieves higher R2 but negative rate"}
  # 
  ## Plot - Original values (before corrections):
  Plot_i <- ggplot(data = Filt_i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
    geom_point() +
    xlab("Sample time (min)") +
    ylab("CH4 by mass (mgm2)") +
    ggtitle(paste("Code = ", i, "; Case = ", Emission_rates_2022$Case_CH4[i],"; Original")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks=c(0, 10, 20, 30)) +
    stat_poly_line() +
    stat_poly_eq() +
    annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(Emission_rates_2022$CH4_flux_mgm2h[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)
  
print(Plot_i)
  
  # 
  # ## Plot Alt_1:
  # Plot_Alt_1 <- ggplot(data = Filt_Alt1i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CH4 by mass (mgm2)") +
  #   ggtitle(paste("Alt. Model 1")) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   scale_x_continuous(breaks=c(0, 10, 20, 30)) +
  #   stat_poly_line() +
  #   stat_poly_eq()  +
  #   annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(Emission_rates_2022$CH4_flux_Alt1[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)
  # 
  # ## Plot Alt_2:
  # Plot_Alt_2 <- ggplot(data = Filt_Alt2i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CH4 by mass (mgm2)") +
  #   ggtitle(paste("Alt. Model 2")) +
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   scale_x_continuous(breaks=c(0, 10, 20, 30)) +
  #   stat_poly_line() +
  #   stat_poly_eq()  +
  #   annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(Emission_rates_2022$CH4_flux_Alt2[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)
  # 
  # ## Plot Alt_3:
  # Plot_Alt_3 <- ggplot(data = Filt_Alt3i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CH4 by mass (mgm2)") +
  #   ggtitle(paste("Alt. Model 3")) +
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   scale_x_continuous(breaks = c(0, 10, 20, 30)) +
  #   stat_poly_line() +
  #   stat_poly_eq()  +
  #   annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(Emission_rates_2022$CH4_flux_Alt3[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)
  # 
  # ## Plot Alt_4:
  # Plot_Alt_4 <- ggplot(data = Filt_Alt4i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CH4 by mass (mgm2)") +
  #   ggtitle(paste("Alt. Model 4")) +
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   scale_x_continuous(breaks=c(0, 10, 20, 30)) +
  #   stat_poly_line() +
  #   stat_poly_eq()  +
  #   annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(Emission_rates_2022$CH4_flux_Alt4[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)
  # 
  # ## Log Model:
  # lm_eqn <- function(Filt_i){                                                    ## Function defined to include R2 in log fitted plot.
  #   log_coef <- lm(CH4_byMass_mgm2~log(Sample_time_min+1), data=Filt_i);
  #   eq <- substitute(italic(y) == ~~italic(R)^2~"="~r2, 
  #                    list(r2 = format(summary(log_coef)$r.squared, digits = 3)))
  #   as.character(as.expression(eq));
  # }
  # 
  # Plot_Log <- ggplot(data = Filt_i, aes(x=Sample_time_min, y=CH4_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CH4 by mass (mgm2)") +
  #   ggtitle(paste("Log Model: y~log(x+1)")) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   stat_poly_line(data = Filt_i, method="lm",formula=y~log(x+1),fill="red") +
  #   stat_poly_eq(data = Filt_i, method="lm",formula=y~log(x+1))
  # 
  # ## Arrange plots:
  # CH4_arrange <- ggarrange(Plot_i, Plot_Alt_1, Plot_Alt_2, Plot_Alt_3, Plot_Alt_4, Plot_Log,
  #                          ncol = 2, nrow = 3)   
  # print(CH4_arrange)
}  

dev.off()

# ## Data frame comparing original with chosen alternative model (if any):
# Results_CH4 <- Emission_rates_2022[c("Code_Nr", "Code", "CH4_flux_mgm2h", "R2_CH4", "CH4_model","CH4_flux_corrected", "R2_CH4_corrected", "Logic_CH4")] ## Creates data frame with certain columns form Emission_rates_2022
# 
# ## Plot original vs modeled rates:
# Plot_compare_CH4 <- Results_CH4 %>% filter(Code_Nr != 43) %>% filter(Code_Nr != 45) %>%  # Compares original and modeled rates, excluding codes 43 and 45 with outlier rates
#   ggplot(aes(x=CH4_flux_mgm2h, y=CH4_flux_corrected)) +
#   geom_point() +
#   xlab("Original rate (mgm2h)") +
#   ylab("Chosen model rate (mgm2h)") +
#   ggtitle(paste("Chosen Model vs Original Rates")) +
#   theme(plot.title = element_text(hjust = 0.5))+
#   stat_poly_line() +
#   stat_poly_eq() +
#   stat_regline_equation(aes(label = paste(after_stat(eq.label))), label.y.npc = 'top')
# print(Plot_compare_CH4)
# 
# Results_CH4 %>% count(CH4_model) ## Codes per chosen model

# Create data frame with selected method per Code.
# Plot original versus selected rates.

#### 4.2 Loops for N2O ####

# pdf('outputs/2022/N2O_Plots_original.pdf')

for (j in 1:length(Emission_rates_2022$Code)) {
  Code_j <- Emission_rates_2022$Code[j]
  Filt_j <- filter(Chrom_results_2022,Chrom_results_2022$Code == Code_j) # if returned as Time-Series, ru-run library(dplyr)
  
  ## Loop section 1: Rate calculation.
  lm_j <- lm(N2O_byMass_mgm2~Sample_time_min, data=Filt_j)
  Emission_rates_2022$N2O_flux_mgm2h[j] <- coef(lm_j)[2]*60 # Returns N2O_flux_mgm2h for each Code.
  Emission_rates_2022$R2_N2O[j] <- summary(lm_j)$r.squared # Returns R2_N2O for each Code.
  lmp_j <- function (modelobject) {  # Function created to call later the model's p-value
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  Emission_rates_2022$p_N2O[j] <- lmp_j(lm_j) # Returns p-value for each Code.
  
  # ## Loop section 2: Rate correction.
  # #### Loop section 2.1: Parameters for Outliers restriction:
  # 
  # Filt_out_j <- subset(Chrom_results_2022, Code == unique(Code)[j], select=c(N2O_byMass_mgm2, Sample_time_min)) 
  # Filt_out_j <- Filt_out_j %>% drop_na()
  # summary(Filt_out_j$N2O_byMass_mgm2)
  # IQR(Filt_out_j$N2O_byMass_mgm2)
  # Tmin_j = (summary(Filt_out_j$N2O_byMass_mgm2)[2]) - (1.5*IQR(Filt_out_j$N2O_byMass_mgm2)) ## T = Threshold values
  # Tmax_j = (summary(Filt_out_j$N2O_byMass_mgm2)[5]) + (1.5*IQR(Filt_out_j$N2O_byMass_mgm2)) 
  # Out_j <- Filt_out_j$N2O_byMass_mgm2[which(Filt_out_j$N2O_byMass_mgm2 < Tmin_j | Filt_out_j$N2O_byMass_mgm2 > Tmax_j)] # Identifying outliers
  # RemOut_j <- Filt_out_j$N2O_byMass_mgm2[which(Filt_out_j$N2O_byMass_mgm2 > Tmin_j & Filt_out_j$N2O_byMass_mgm2 < Tmax_j)] # Returns NOT outlier values
  # GHG_rates_RemOut_j <- subset(Filt_out_j, N2O_byMass_mgm2 %in% RemOut_j)
  # lm_RemOut_j <- lm(N2O_byMass_mgm2~Sample_time_min, data=GHG_rates_RemOut_j)
  # N2O_flux_OutRem_j <- coef(lm_RemOut_j)[2]*60 # Returns N2O_flux_corrected for each Code.
  # N2O_R2_OutRem_j <- summary(lm_RemOut_j)$r.squared # Returns R2_N2O_corrected for each Code.
  # 
  # ## Loop section 2.3: Including restrictions in the loop with nested if_else:
  # Emission_rates_2022$N2O_flux_corrected[j] <- if_else(Emission_rates_2022$R2_N2O[j]<0.3, 0, ## Restriction 1: R2 < 0.3 are removed.
  #                                                 if_else(Emission_rates_2022$R2_N2O[j]<0.7 & (length(RemOut_j) != length(Filt_out_j$N2O_byMass_mgm2)), N2O_flux_OutRem_j, ## Restriction 2: Outliers removal
  #                                                         Emission_rates_2022$N2O_flux_mgm2h[j]))
  # 
  # ## Loop section 2.4: Calculating j2 according to the applied correction (if any):
  # Emission_rates_2022$R2_N2O_corrected[j] <- if_else(Emission_rates_2022$R2_N2O[j]<0.3, 0,
  #                                               if_else(Emission_rates_2022$R2_N2O[j]<0.7 & (length(RemOut_j) != length(Filt_out_j$N2O_byMass_mgm2)), N2O_R2_OutRem_j,
  #                                                       Emission_rates_2022$R2_N2O[j]))
  # ## Plots before corrections:
  # 
  # Plot_j <- ggplot(data = Filt_j, aes(x=Sample_time_min, y=N2O_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("N2O by mass (mgm2)") +
  #   ggtitle(paste("Tr1 = ", Filt_j$Tr1[1], "; Tr2 = ", Filt_j$Tr2[1], "; Rep = ", Filt_j$Rep[1], ";", "Samp. date = ", Filt_j$Sampling_date[1])) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # 
  # print(Plot_j)
  
} 

# dev.off()

##### 4.3 Loops for CO2 ####

# pdf('outputs/2022/CO2_Plots_original.pdf')

for (k in 1:length(Emission_rates_2022$Code)) {
  Code_k <- Emission_rates_2022$Code[k]
  Filt_k <- filter(Chrom_results_2022,Chrom_results_2022$Code == Code_k) # if returned as Time-Series, re-run library(dplyr)
  
  ## Loop section 1: Rate calculation.
  lm_k <- lm(CO2_byMass_mgm2~Sample_time_min, data=Filt_k)
  Emission_rates_2022$CO2_flux_mgm2h[k] <- coef(lm_k)[2]*60 # Returns CO2_flux_mgm2h for each Code.
  Emission_rates_2022$R2_CO2[k] <- summary(lm_k)$r.squared # Returns R2_CO2 for each Code.
  lmp_k <- function (modelobject) {  # Function created to call later the model's p-value
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  Emission_rates_2022$p_CO2[k] <- lmp_k(lm_k) # Returns p-value for each Code.
  
  # ## Loop section 2: Rate correction.
  # #### Loop section 2.1: Parameters for Outliers restriction:
  # 
  # Filt_out_k <- subset(Chrom_results_2022, Code == unique(Code)[k], select=c(CO2_byMass_mgm2, Sample_time_min)) 
  # Filt_out_k <- Filt_out_k %>% drop_na()
  # summary(Filt_out_k$CO2_byMass_mgm2)
  # IQR(Filt_out_k$CO2_byMass_mgm2)
  # Tmin_k = (summary(Filt_out_k$CO2_byMass_mgm2)[2]) - (1.5*IQR(Filt_out_k$CO2_byMass_mgm2)) ## T = Threshold values
  # Tmax_k = (summary(Filt_out_k$CO2_byMass_mgm2)[5]) + (1.5*IQR(Filt_out_k$CO2_byMass_mgm2)) 
  # Out_k <- Filt_out_k$CO2_byMass_mgm2[which(Filt_out_k$CO2_byMass_mgm2 < Tmin_k | Filt_out_k$CO2_byMass_mgm2 > Tmax_k)] # Identifying outliers
  # RemOut_k <- Filt_out_k$CO2_byMass_mgm2[which(Filt_out_k$CO2_byMass_mgm2 > Tmin_k & Filt_out_k$CO2_byMass_mgm2 < Tmax_k)] # Returns NOT outlier values
  # GHG_rates_RemOut_k <- subset(Filt_out_k, CO2_byMass_mgm2 %in% RemOut_k)
  # lm_RemOut_k <- lm(CO2_byMass_mgm2~Sample_time_min, data=GHG_rates_RemOut_k)
  # CO2_flux_OutRem_k <- coef(lm_RemOut_k)[2]*60 # Returns CO2_flux_corrected for each Code.
  # CO2_R2_OutRem_k <- summary(lm_RemOut_k)$r.squared # Returns R2_CO2_corrected for each Code.
  # 
  # ## Loop section 2.3: Including restrictions in the loop with nested if_else:
  # Emission_rates_2022$CO2_flux_corrected[k] <- if_else(Emission_rates_2022$R2_CO2[k]<0.3, 0, ## Restriction 1: R2 < 0.3 are removed.
  #                                                 if_else(Emission_rates_2022$R2_CO2[k]<0.7 & (length(RemOut_k) != length(Filt_out_k$CO2_byMass_mgm2)), CO2_flux_OutRem_k, ## Restriction 2: Outliers removal
  #                                                         Emission_rates_2022$CO2_flux_mgm2h[k]))
  # 
  # ## Loop section 2.4: Calculating j2 according to the applied correction (if any):
  # Emission_rates_2022$R2_CO2_corrected[k] <- if_else(Emission_rates_2022$R2_CO2[k]<0.3, 0,
  #                                               if_else(Emission_rates_2022$R2_CO2[k]<0.7 & (length(RemOut_k) != length(Filt_out_k$CO2_byMass_mgm2)), CO2_R2_OutRem_k,
  #                                                       Emission_rates_2022$R2_CO2[k]))
  # 
  # ## Plots before corrections:
  #     
  # Plot_k <- ggplot(data = Filt_k, aes(x=Sample_time_min, y=CO2_byMass_mgm2)) +
  #   geom_point() +
  #   xlab("Sample time (min)") +
  #   ylab("CO2 by mass (mgm2)") +
  #   ggtitle(paste("Tr1 = ", Filt_k$Tr1[1], "; Tr2 = ", Filt_k$Tr2[1], "; Rep = ", Filt_k$Rep[1], ";", "Samp. date = ", Filt_k$Sampling_date[1])) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # 
  # print(Plot_k)
  # 
} 

# dev.off()

write_xlsx(Emission_rates_2022, "outputs/2022/Emission_rates_original_2022.xlsx")
