
####################################################### Thesis Paper 1 - Emissions #####################################################

library(dplyr)
library(ggplot2)
library("purrr")
library(writexl)
library(reshape2)
library(plotly)
library(hrbrthemes)
library(RColorBrewer)
library(ggraph)
library(gridExtra)

############### 1. Create Master_GHG_2022 Dataframe #########################

load("C:/Users/SECHEVERRIA/R_git/Thesis_Paper_1_git/outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData") # Load emissions 2022 with corrections:  

## Master_GHG_2022 taking corrected CH4. Original data for N2O and CO2: 
Emissions_2022 <- select(Emission_rates_w_corrections_2022, Sampling_date, Plot, Treat, Rep, CH4_flux_corrected, N2O_flux_corrected)

Water_level_2022 <- read.csv("data/Other_vars/Piezo_2022.csv", fileEncoding="latin1", na.strings=c("","NA")) %>% # Import water level (piezometer) data:
                    rename(Water_level_piezo = Water_level_cm)  

Master_GHG_2022  <- merge(Emissions_2022, Water_level_2022, by.x = c("Sampling_date", "Treat", "Plot", "Rep"), by.y = c("Date", "Treat", "Plot", "Rep"), all.x= TRUE, all.y = TRUE) %>%  # with all.x= TRUE, all.y = TRUE, the resulting dataframe contains as well dates where either only water level or emissions were recorded.
                        select(Sampling_date, Treat, Plot, Rep, CH4_flux_corrected, N2O_flux_corrected, Water_level_piezo)

Field_sheet_chrom_2022 <- read.csv("C:/Users/SECHEVERRIA/R_git/Thesis_Paper_1_git/data/GHG/Field_sheet_chrom_2022.csv", fileEncoding="latin1", na.strings=c("","NA")) %>% # Load field sheet 2022:
                          rename(Water_level_piezo = Water_level_cm)

Field_sheet_other_factors_2022 <- Field_sheet_chrom_2022 %>% 
                          group_by(Sampling_date, Plot) %>% 
                          summarise(Water_level_ruler = mean(Water_level_piezo), Temp_soil = mean(Temp_soil), Rice_cover_prop = mean(Rice_cover_prop), Env_temp_initial = mean(Env_temp_initial), Env_temp_final = mean(Env_temp_final)) # creates dataframe with all factors measured during chromatography campaign dates, besides GHG emissions.

#######  Water level correction:

Master_GHG_2022  <- Master_GHG_2022  %>% 
                        left_join(Field_sheet_other_factors_2022, by = c("Sampling_date", "Plot"))
Master_GHG_2022 $Sampling_date <- as.Date(Master_GHG_2022 $Sampling_date)
Master_GHG_2022  <- Master_GHG_2022 [!(Master_GHG_2022 $Plot %in% c("P10","P11", "P12", "P13", "P14", "P15")), ] #Removes Rep 4 and Rep 5
Master_GHG_2022 $Row_Nr <- 1:nrow(Master_GHG_2022 ) # Adding a row number column
Master_GHG_2022 $Water_level_corr <- NA # add column "Water_level_corr"
Master_GHG_2022  <- Master_GHG_2022 [, c(13, 1, 2, 3, 4, 5, 6, 7, 8, 14, 9, 10, 11, 12)] # Reorder columns
Master_GHG_2022  <- Master_GHG_2022 [!(is.na(Master_GHG_2022$CH4_flux_corrected) & is.na(Master_GHG_2022$N2O_flux_corrected) & is.na(Master_GHG_2022$Water_level_piezo)), ] # Removes rows without emissions and piezometer records.

# The created water level column must fulfill the following conditions:
  # a) If there is a valid value in Water_level_ruler and Water_level_piezo, take the Water_level_ruler value, unless the Water_level_ruler value is 0, then take the Water_level_piezo value; 
  # b) If there is only a Water_level_piezo value, take it.
  # c) If there is only a Water_level_ruler value > 0 then take it. 
  # d) If there is no Water_level_piezo and Water_level_ruler is either absent or = 0:
      #  d.i) AWD & MSD: look for the closest negative Water_level_piezo to the Sampling_date of this row that also has the same Plot value of this row.
      #  d.ii) CON: look for the closest negative Water_level_piezo to this Sampling_date of this row, corresponding to the MSD plot within this repetition (block):

## Create Water_level_corr (corrected) column considering piezometer and field ruler measurements (taken on sampling dates):

## This "for in loop" works already to apply conditions to the Master_GHG_2022, but it takes time to run:

for (i in 1:length(Master_GHG_2022 $Row_Nr)) {
  current_plot <- Master_GHG_2022 $Plot[i]
  current_rep <- Master_GHG_2022 $Rep[i]
  neg_piezo_indices <- which(Master_GHG_2022 $Water_level_piezo <= 0 & Master_GHG_2022 $Plot == current_plot)
  neg_piezo_indices_CON <- which(Master_GHG_2022 $Water_level_piezo <= 0 & Master_GHG_2022 $Rep == current_rep & Master_GHG_2022 $Treat == "MSD")
  closest_neg_index <- ifelse(Master_GHG_2022 $Treat[i] %in% c("AWD", "MSD"), neg_piezo_indices[which.min(abs(as.numeric(Master_GHG_2022 $Sampling_date[i] - Master_GHG_2022 $Sampling_date[neg_piezo_indices])))], 1) # the ifelse (... , ... , 1) solves the "replacement of length zero" for cases CON / piezo: NA / ruler: 0
  closest_neg_index_CON <- neg_piezo_indices_CON[which.min(abs(as.numeric(Master_GHG_2022 $Sampling_date[i] - Master_GHG_2022 $Sampling_date[neg_piezo_indices_CON])))]

  # Applying conditions for all cases (Opening each ifelse() in different LHS ~ RHS functions):

  Master_GHG_2022 $Water_level_corr[i] <- case_when(
    !is.na(Master_GHG_2022 $Water_level_ruler[i]) & !is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Water_level_ruler[i] == 0 ~ Master_GHG_2022 $Water_level_piezo[i], # a) part I
    !is.na(Master_GHG_2022 $Water_level_ruler[i]) & !is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Water_level_ruler[i] != 0 ~ Master_GHG_2022 $Water_level_ruler[i], # a) part II
    is.na(Master_GHG_2022 $Water_level_ruler[i]) & !is.na(Master_GHG_2022 $Water_level_piezo[i]) ~ Master_GHG_2022 $Water_level_piezo[i], # b)
    !is.na(Master_GHG_2022 $Water_level_ruler[i]) & is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Water_level_ruler[i] != 0 ~ Master_GHG_2022 $Water_level_ruler[i], # c)
    (is.na(Master_GHG_2022 $Water_level_ruler[i]) | Master_GHG_2022 $Water_level_ruler[i] == 0) & is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Treat[i] %in% c("AWD", "MSD") & is.na(closest_neg_index) ~ Master_GHG_2022 $Water_level_ruler[i], # d.i) part I
    (is.na(Master_GHG_2022 $Water_level_ruler[i]) | Master_GHG_2022 $Water_level_ruler[i] == 0) & is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Treat[i] %in% c("AWD", "MSD") & !is.na(closest_neg_index) ~ Master_GHG_2022 $Water_level_piezo[closest_neg_index], # d.i) part II
    (is.na(Master_GHG_2022 $Water_level_ruler[i]) | Master_GHG_2022 $Water_level_ruler[i] == 0) & is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Treat[i] %in% "CON" & is.na(closest_neg_index_CON) ~ Master_GHG_2022 $Water_level_ruler[i], # d.ii) part I
    (is.na(Master_GHG_2022 $Water_level_ruler[i]) | Master_GHG_2022 $Water_level_ruler[i] == 0) & is.na(Master_GHG_2022 $Water_level_piezo[i]) & Master_GHG_2022 $Treat[i] %in% "CON" & !is.na(closest_neg_index_CON) ~ Master_GHG_2022 $Water_level_piezo[closest_neg_index_CON], # d.ii) part II
    TRUE ~ NA_real_
  )
}

# write_xlsx(Master_GHG_2022 , "outputs/Master_GHG_2022.xlsx")

save(Master_GHG_2022, file = "outputs/Master_GHG_2022.RData")

# Including physicochemical parameters (recorded for water and soils during sampling dates):
physchem_2022 <- read.csv("data/Other_vars/Other_factors_CERESTRES_2022.csv", fileEncoding="latin1", na.strings=c("","NA"))
physchem_2022$Sampling_date <- as.Date(physchem_2022$Sampling_date)
Master_GHG_2022  <- Master_GHG_2022  %>% 
  left_join(physchem_2022, by = c("Sampling_date", "Plot", "Rep", "Treat"))

############### 2. Plot Analysis #########################

# Averaged GHG emission rates and other variables:
Avged_Master_GHG_2022 <- Master_GHG_2022 %>% 
                          group_by(Sampling_date, Treat) %>% 
                          summarize(across(
                            .cols = where(is.numeric),  # Exclude grouping columns
                            .fns = mean,  # Summary function (mean in this case)
                            .names = "avg_{.col}"  # New column names for the summarized values
                          )) %>% 
                          select(-c(avg_Rep, avg_Row_Nr)) # Excludes these non logical outcome columns.

### Plotting emission fluxes and water level against time: ###

#### I. CH4 Plots: ####

# Creating a dataframe in plotting format:
CH4_melted_Avged_Master_GHG_2022 <- melt(Avged_Master_GHG_2022, id.vars = c("Sampling_date", "Treat"),
                                     measure.vars = c("avg_CH4_flux_corrected", "avg_Water_level_corr"))
CH4_melted_Avged_Master_GHG_2022  <- CH4_melted_Avged_Master_GHG_2022 [!(is.na(CH4_melted_Avged_Master_GHG_2022$value)), ] # Removes rows without emissions and piezometer records.

# I.i. All in one plot:

CH4_flux_water_A <- ggplot(CH4_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_CH4_flux_corrected", color = Treat),
                              data = subset(CH4_melted_Avged_Master_GHG_2022, variable == "avg_CH4_flux_corrected")) +
                    geom_line(aes(y = value/10, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(CH4_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr")) +
                    scale_y_continuous(name = expression('CH4 flux (mg m'^"-2"*' h'^"-1"*')'), limits = c(-5, 10), breaks = seq(-6, 10, by = 2), 
                                       sec.axis = sec_axis(~ .*10, name = "Water level (cm)", breaks = seq(-60, 100, by = 20), labels = seq(-60, 100, by = 20))) +
                    scale_colour_brewer(palette = "Set2", direction = - 1) +
                    theme_bw() +
                    labs(x = "Sampling Date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_CH4_flux_corrected" = "solid", "avg_Water_level_corr" = "dashed"), labels = c("CH4 flux (mg/m2*h)", "Water level (cm)")) +
                    guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"))

print(CH4_flux_water_A)

ggsave("outputs/Plots/GHG/CH4_flux_water_A.pdf", width = 10, height = 5)

# I.ii. Faceted plots:

# Already working code, averaging plots per treatment:

CH4_flux_water_B <- ggplot(CH4_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                     geom_line(aes(y = value, linetype = "avg_CH4_flux_corrected", color = Treat),
                               data = subset(CH4_melted_Avged_Master_GHG_2022, variable == "avg_CH4_flux_corrected")) +
                     geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                               data = subset(CH4_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr")) +
                     # scale_y_continuous(name = expression('CH4 flux (mg m'^"-2"*' h'^"-1"*')'), limits = c(-5, 10), breaks = seq(-6, 10, by = 2),
                     #                    sec.axis = sec_axis(~ .*10, name = "Water level (cm)", breaks = seq(-60, 100, by = 20), labels = seq(-60, 100, by = 20))) +
                     scale_colour_brewer(palette = "Set2", direction = - 1) +
                     theme_bw() +
                     labs(x = "Sampling Date") +
                     geom_hline(yintercept=0, color = "grey") +
                     scale_linetype_manual(name = "Variable",
                                           values = c("avg_CH4_flux_corrected" = "solid", "avg_Water_level_corr" = "dashed"), labels = c("CH4 flux (mg/m2*h)", "Water level (cm)")) +
                     guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                     theme(axis.title.y = element_text(color = "black"),
                           axis.text.y = element_text(color = "black"),
                           axis.title.y.right = element_text(color = "black"),
                           axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                           strip.placement = "outside")+
                     facet_wrap(~variable,  ncol=1, scales="free_y", strip.position = "left", labeller = as_labeller(c(avg_CH4_flux_corrected = "CH4 flux (mg/m2*h)", avg_Water_level_corr = "Water level (cm)") ))+
                     ylab(NULL)

print(CH4_flux_water_B)

ggsave("outputs/Plots/GHG/CH4_flux_water_B.pdf", width = 10, height = 5)

# I.iii. Including all repetitions, not just their average:

Master_GHG_2022_no_NA <- subset(Master_GHG_2022, is.na(CH4_flux_corrected) == FALSE) 
Master_GHG_2022_NA <- subset(Master_GHG_2022, is.na(CH4_flux_corrected) == TRUE)
 
CH4_flux_C <- ggplot(Master_GHG_2022_no_NA, aes(x = Sampling_date, y = CH4_flux_corrected, color = Treat, group = Plot)) +
                      geom_line(alpha = 0.5, linetype = "dotted") +  # Adjust transparency by setting alpha
                      scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                      theme_bw() +
                      labs(y = expression(paste(CH[4], " flux (mg ", m^-2, " ", h^-1, ")"))) +
                      geom_hline(yintercept = 0, color = "grey") +
                      guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                      theme(
                        axis.title.y = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        axis.title.y.right = element_text(color = "black"),
                        axis.text.y.right = element_text(color = "black"),
                        strip.background = element_blank(),
                        strip.placement = "outside",
                        axis.text.x = element_blank(),
                        legend.position="top",
                        plot.margin = unit(c(1, 1, 0, 1), "lines")) +
                      xlab(NULL) +
                      scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01"))) 

Avg_Master_GHG_2022_no_NA <- Master_GHG_2022_no_NA %>% 
                              group_by(Sampling_date, Treat) %>% 
                              summarise(CH4_flux_corrected = mean(CH4_flux_corrected), 
                                        N2O_flux_corrected = mean(N2O_flux_corrected))

CH4_flux_C <- CH4_flux_C +
                      geom_line(data = Avg_Master_GHG_2022_no_NA, aes(x = Sampling_date, y = CH4_flux_corrected, group = Treat))

print(CH4_flux_C) # CH4 emissions with plt rates in dotted lines

CH4_water_C <- ggplot(CH4_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(CH4_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr"), show.legend = FALSE) +
                    scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                    theme_bw() +
                    labs(y = "Water level (cm)") +
                    labs(x = "Sampling date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                          values = c("avg_Water_level_corr" = "dashed"), labels = "Water level (cm)") +
                    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                          strip.placement = "outside",
                          plot.margin = unit(c(0, 1, 1, 1), "lines")) +
                    scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01")))

print(CH4_water_C) # Water level all treats

# Arrange plots:

CH4_flux_water_C <- grid.arrange(arrangeGrob(CH4_flux_C, CH4_water_C, nrow = 2, ncol = 1))

ggsave("outputs/Plots/GHG/CH4_flux_water_C.pdf", width = 7, height = 5, plot = CH4_flux_water_C) # Arrange CH4 emissions + Water level all treats

# I.iv. CH4 Plots - One Treat's water level per arrange - BES 2023 Ppt 

Water_CON <- CH4_melted_Avged_Master_GHG_2022 %>% 
              filter(Treat == "CON")
Water_MSD <- CH4_melted_Avged_Master_GHG_2022 %>% 
              filter(Treat == "MSD")
Water_AWD <- CH4_melted_Avged_Master_GHG_2022 %>% 
              filter(Treat == "AWD")

CH4_water_CON <- ggplot(Water_CON, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(Water_CON, variable == "avg_Water_level_corr"), show.legend = FALSE) +
                    scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                    theme_bw() +
                    labs(y = "Water level (cm)") +
                    labs(x = "Sampling date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_Water_level_corr" = "dashed"), labels = "Water level (cm)") +
                    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                          strip.placement = "outside",
                          plot.margin = unit(c(0, 1, 1, 1), "lines")) +
                    scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01")))

print(CH4_water_CON) # Water level for just CON treat

CH4_water_MSD <- ggplot(Water_MSD, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(Water_MSD, variable == "avg_Water_level_corr"), show.legend = FALSE) +
                    scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                    theme_bw() +
                    labs(y = "Water level (cm)") +
                    labs(x = "Sampling date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_Water_level_corr" = "dashed"), labels = "Water level (cm)") +
                    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                          strip.placement = "outside",
                          plot.margin = unit(c(0, 1, 1, 1), "lines")) +
                    scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01")))

print(CH4_water_MSD) # Water level for just MSD treat


CH4_water_AWD <- ggplot(Water_AWD, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(Water_AWD, variable == "avg_Water_level_corr"), show.legend = FALSE) +
                    scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                    theme_bw() +
                    labs(y = "Water level (cm)") +
                    labs(x = "Sampling date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_Water_level_corr" = "dashed"), labels = "Water level (cm)") +
                    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                          strip.placement = "outside",
                          plot.margin = unit(c(0, 1, 1, 1), "lines")) +
                    scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01")))

print(CH4_water_AWD) # Water level for just AWD treat

# Arrange plots:

# CH4 emissions + CON water level:
CH4_flux_water_CON <- grid.arrange(arrangeGrob(CH4_flux_C, CH4_water_CON, nrow = 2, ncol = 1))
ggsave("outputs/Plots/GHG/CH4_flux_water_CON.pdf", width = 7, height = 5, plot = CH4_flux_water_CON) # Arrange CH4 emissions + CON water level 

# CH4 emissions + MSD water level:
CH4_flux_water_MSD <- grid.arrange(arrangeGrob(CH4_flux_C, CH4_water_MSD, nrow = 2, ncol = 1))
ggsave("outputs/Plots/GHG/CH4_flux_water_MSD.pdf", width = 7, height = 5, plot = CH4_flux_water_MSD) # Arrange CH4 emissions + CON water level 

# CH4 emissions + AWD water level:
CH4_flux_water_AWD <- grid.arrange(arrangeGrob(CH4_flux_C, CH4_water_AWD, nrow = 2, ncol = 1))
ggsave("outputs/Plots/GHG/CH4_flux_water_AWD.pdf", width = 7, height = 5, plot = CH4_flux_water_AWD) # Arrange CH4 emissions + CON water level 

#### II. N2O Plots: ####

# Creating a dataframe in plotting format:
N2O_melted_Avged_Master_GHG_2022 <- melt(Avged_Master_GHG_2022, id.vars = c("Sampling_date", "Treat"),
                                         measure.vars = c("avg_N2O_flux_corrected", "avg_Water_level_corr"))
N2O_melted_Avged_Master_GHG_2022  <- N2O_melted_Avged_Master_GHG_2022 [!(is.na(N2O_melted_Avged_Master_GHG_2022$value)), ] # Removes rows without emissions and piezometer records.

# II.i. All in one plot:

N2O_flux_water_A <- ggplot(N2O_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_N2O_flux_corrected", color = Treat),
                              data = subset(N2O_melted_Avged_Master_GHG_2022, variable == "avg_N2O_flux_corrected")) +
                    geom_line(aes(y = value/10, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(N2O_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr")) +
                    scale_y_continuous(name = expression('N2O flux (mg m'^"-2"*' h'^"-1"*')'), limits = c(-5, 10), breaks = seq(-6, 10, by = 2), 
                                       sec.axis = sec_axis(~ .*10, name = "Water level (cm)", breaks = seq(-60, 100, by = 20), labels = seq(-60, 100, by = 20))) +
                    scale_colour_brewer(palette = "Set2", direction = - 1) +
                    theme_bw() +
                    labs(x = "Sampling Date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_N2O_flux_corrected" = "solid", "avg_Water_level_corr" = "dashed"), labels = c("N2O flux (mg/m2*h)", "Water level (cm)")) +
                    guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"))

print(N2O_flux_water_A)

ggsave("outputs/Plots/GHG/N2O_flux_water_A.pdf", width = 10, height = 5)

# I.ii. Facetted plots:

N2O_flux_water_B <- ggplot(N2O_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                    geom_line(aes(y = value, linetype = "avg_N2O_flux_corrected", color = Treat),
                              data = subset(N2O_melted_Avged_Master_GHG_2022, variable == "avg_N2O_flux_corrected")) +
                    geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                              data = subset(N2O_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr")) +
                    # scale_y_continuous(name = expression('N2O flux (mg m'^"-2"*' h'^"-1"*')'), limits = c(-5, 10), breaks = seq(-6, 10, by = 2), 
                    #                    sec.axis = sec_axis(~ .*10, name = "Water level (cm)", breaks = seq(-60, 100, by = 20), labels = seq(-60, 100, by = 20))) +
                    scale_colour_brewer(palette = "Set2", direction = - 1) +
                    theme_bw() +
                    labs(x = "Sampling Date") +
                    geom_hline(yintercept=0, color = "grey") +
                    scale_linetype_manual(name = "Variable",
                                          values = c("avg_N2O_flux_corrected" = "solid", "avg_Water_level_corr" = "dashed"), labels = c("N2O flux (mg/m2*h)", "Water level (cm)")) +
                    guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                    theme(axis.title.y = element_text(color = "black"),
                          axis.text.y = element_text(color = "black"),
                          axis.title.y.right = element_text(color = "black"),
                          axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                          strip.placement = "outside")+
                    facet_wrap(~variable,  ncol=1, scales="free_y", strip.position = "left", labeller = as_labeller(c(avg_N2O_flux_corrected = "N2O flux (mg/m2*h)", avg_Water_level_corr = "Water level (cm)") ))+
                    ylab(NULL) 

print(N2O_flux_water_B)

ggsave("outputs/Plots/GHG/N2O_flux_water_B.pdf", width = 10, height = 5)

# I.iii. Including all repetitions, not just their average:

Master_GHG_2022_no_NA_N2O <- subset(Master_GHG_2022, is.na(N2O_flux_corrected) == FALSE) 
Master_GHG_2022_NA_N2O <- subset(Master_GHG_2022, is.na(N2O_flux_corrected) == TRUE)

N2O_flux_C <- ggplot(Master_GHG_2022_no_NA_N2O, aes(x = Sampling_date, y = N2O_flux_corrected, color = Treat, group = Plot)) +
                      geom_line(alpha = 0.5, linetype = "dotted") +  # Adjust transparency by setting alpha
                      scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                      theme_bw() +
                      labs(y = expression(paste(N[2], "O flux (mg ", m^-2, " ", h^-1, ")"))) +
                      geom_hline(yintercept = 0, color = "grey") +
                      guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                      theme(
                        axis.title.y = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        axis.title.y.right = element_text(color = "black"),
                        axis.text.y.right = element_text(color = "black"),
                        strip.background = element_blank(),
                        strip.placement = "outside",
                        axis.text.x = element_blank(),
                        legend.position="top",
                        plot.margin = unit(c(1, 1, 0, 1), "lines")) +
                      xlab(NULL) +
                      scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01"))) 

Avg_Master_GHG_2022_no_NA_N2O <- Master_GHG_2022_no_NA_N2O %>% 
                             group_by(Sampling_date, Treat) %>% 
                             summarise(N2O_flux_corrected = mean(N2O_flux_corrected), 
                                      N2O_flux_corrected = mean(N2O_flux_corrected))

N2O_flux_C <- N2O_flux_C +
              geom_line(data = Avg_Master_GHG_2022_no_NA_N2O, aes(x = Sampling_date, y = N2O_flux_corrected, group = Treat))

print(N2O_flux_C)

N2O_water_C <- ggplot(N2O_melted_Avged_Master_GHG_2022, aes(x = Sampling_date, color = Treat, linetype = variable)) +
                      geom_line(aes(y = value, linetype = "avg_Water_level_corr", color = Treat),
                                data = subset(N2O_melted_Avged_Master_GHG_2022, variable == "avg_Water_level_corr"), show.legend = FALSE) +
                      scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                      theme_bw() +
                      labs(y = "Water level (cm)") +
                      labs(x = "Sampling date") +
                      geom_hline(yintercept=0, color = "grey") +
                      scale_linetype_manual(name = "Variable",
                                            values = c("avg_Water_level_corr" = "dashed"), labels = "Water level (cm)") +
                      guides(linetype = guide_legend(override.aes = list(color = "black"))) +
                      theme(axis.title.y = element_text(color = "black", margin = margin(r = 6)),
                            axis.text.y = element_text(color = "black"),
                            axis.title.y.right = element_text(color = "black"),
                            axis.text.y.right = element_text(color = "black"), strip.background = element_blank(),
                            strip.placement = "outside",
                            plot.margin = unit(c(0, 1, 1, 1), "lines")) +
                      scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01")))

print(N2O_water_C)

# Arrange plots:

N2O_flux_water_C <- grid.arrange(arrangeGrob(N2O_flux_C, N2O_water_C, nrow = 2, ncol = 1))

ggsave("outputs/Plots/GHG/N2O_flux_water_C.pdf", width = 7, height = 5, plot = N2O_flux_water_C)

# CH4, N2O and water arranged: 

CH4_flux_D <- ggplot(Master_GHG_2022_no_NA, aes(x = Sampling_date, y = CH4_flux_corrected, color = Treat, group = Plot)) +
                      geom_line(alpha = 0.5, linetype = "dotted") +  # Adjust transparency by setting alpha
                      scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                      theme_bw() +
                      labs(y = expression(paste(CH[4], " flux (mg ", m^-2, " ", h^-1, ")"))) +
                      geom_hline(yintercept = 0, color = "grey") +
                      guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                      theme(
                        axis.title.y = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        axis.title.y.right = element_text(color = "black"),
                        axis.text.y.right = element_text(color = "black"),
                        strip.background = element_blank(),
                        strip.placement = "outside",
                        axis.text.x = element_blank(),
                        legend.position = c(0.92, 0.7),
                        legend.background = element_rect(fill="white", size = 0.4, linetype="solid", colour = "black"),
                        plot.margin = unit(c(1, 1, 0, 1), "lines")) +
                      xlab(NULL) +
                      scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01"))) +
                      geom_line(data = Avg_Master_GHG_2022_no_NA, aes(x = Sampling_date, y = CH4_flux_corrected, group = Treat))

print(CH4_flux_D)

N2O_flux_D <- ggplot(Master_GHG_2022_no_NA_N2O, aes(x = Sampling_date, y = N2O_flux_corrected, color = Treat, group = Plot)) +
                      geom_line(alpha = 0.5, linetype = "dotted") +  # Adjust transparency by setting alpha
                      scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                      theme_bw() +
                      labs(y = expression(paste(N[2], "O flux (mg ", m^-2, " ", h^-1, ")"))) +
                      geom_hline(yintercept = 0, color = "grey") +
                      guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                      theme(
                        axis.title.y = element_text(color = "black"),
                        axis.text.y = element_text(color = "black", margin = margin(r = 0)),
                        axis.title.y.right = element_text(color = "black"),
                        axis.text.y.right = element_text(color = "black"),
                        strip.background = element_blank(),
                        strip.placement = "outside",
                        axis.text.x = element_blank(),
                        legend.position = "none",
                        plot.margin = unit(c(0, 1, 0, 1), "lines")) +
                      xlab(NULL) +
                      scale_x_date(limits = c(as.Date("2022-05-15"), as.Date("2022-10-01"))) +
                      geom_line(data = Avg_Master_GHG_2022_no_NA_N2O, aes(x = Sampling_date, y = N2O_flux_corrected, group = Treat))

print(N2O_flux_D)

# Arrange plots:

CH4_N2O_flux_water <- gridExtra::grid.arrange(gridExtra::arrangeGrob(CH4_flux_D, N2O_flux_D, CH4_water_C, nrow = 3, ncol = 1))

ggsave("outputs/Plots/GHG/CH4_N2O_flux_water.pdf", width = 7, height = 7, plot = CH4_N2O_flux_water)
