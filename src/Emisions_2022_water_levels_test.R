
# i = 158 # i value to test cases

for (i in 1:length(Emissions_water_2022$Row_Nr)) {
  current_plot <- Emissions_water_2022$Plot[i]
  current_rep <- Emissions_water_2022$Rep[i]
  neg_piezo_indices <- which(Emissions_water_2022$Water_level_piezo <= 0 & Emissions_water_2022$Plot == current_plot)  
  neg_piezo_indices_CON <- which(Emissions_water_2022$Water_level_piezo <= 0 & Emissions_water_2022$Rep == current_rep & Emissions_water_2022$Treat == "MSD")
  closest_neg_index <- ifelse(Emissions_water_2022$Treat[i] %in% c("AWD", "MSD"), neg_piezo_indices[which.min(abs(as.numeric(Emissions_water_2022$Sampling_date[i] - Emissions_water_2022$Sampling_date[neg_piezo_indices])))], 1) # the ifelse (... , ... , 1) solves the "replacement of length zero" for cases CON / piezo: NA / ruler: 0
  closest_neg_index_CON <- neg_piezo_indices_CON[which.min(abs(as.numeric(Emissions_water_2022$Sampling_date[i] - Emissions_water_2022$Sampling_date[neg_piezo_indices_CON])))]

# Applying conditions for all cases (Opening each ifelse() in different LHS ~ RHS functions):   
   
   Emissions_water_2022$Water_level_corr[i] <- case_when(
     !is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Water_level_ruler[i] == 0 ~ Emissions_water_2022$Water_level_piezo[i], # a) part I
     !is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Water_level_ruler[i] != 0 ~ Emissions_water_2022$Water_level_ruler[i], # a) part II
     is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i]) ~ Emissions_water_2022$Water_level_piezo[i], # b)
     !is.na(Emissions_water_2022$Water_level_ruler[i]) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Water_level_ruler[i] != 0 ~ Emissions_water_2022$Water_level_ruler[i], # c)
     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% c("AWD", "MSD") & is.na(closest_neg_index) ~ Emissions_water_2022$Water_level_ruler[i], # d.i) part I
     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% c("AWD", "MSD") & !is.na(closest_neg_index) ~ Emissions_water_2022$Water_level_piezo[closest_neg_index], # d.i) part II
     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% "CON" & is.na(closest_neg_index_CON) ~ Emissions_water_2022$Water_level_ruler[i], # d.ii) part I
     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% "CON" & !is.na(closest_neg_index_CON) ~ Emissions_water_2022$Water_level_piezo[closest_neg_index_CON], # d.ii) part II    
     TRUE ~ NA_real_
     )
}   
   
# Original:     
   
 #    Emissions_water_2022$Water_level_corr[i] <- case_when(
 #     !is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i]) ~ ifelse(Emissions_water_2022$Water_level_ruler[i] == 0, Emissions_water_2022$Water_level_piezo[i], Emissions_water_2022$Water_level_ruler[i]), # a)
 #     is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i]) ~ Emissions_water_2022$Water_level_piezo[i], # b)
 #     !is.na(Emissions_water_2022$Water_level_ruler[i]) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Water_level_ruler[i] != 0 ~ Emissions_water_2022$Water_level_ruler[i], # c)
 #     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% c("AWD", "MSD") ~ { 
 #     ifelse(is.na(closest_neg_index), Emissions_water_2022$Water_level_ruler[i], Emissions_water_2022$Water_level_piezo[closest_neg_index])}, #d.i)
 #     (is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% "CON" ~ {
 #     ifelse(is.na(closest_neg_index_CON), Emissions_water_2022$Water_level_ruler[i], Emissions_water_2022$Water_level_piezo[closest_neg_index_CON])}, #d.ii)
 #     TRUE ~ NA_real_
 #   )
 #   
    
# Written as if(){} else {}:
   
  # Emissions_water_2022$Water_level_corr[i] <- case_when(
  #   if(!is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i])) {
  #     if(Emissions_water_2022$Water_level_ruler[i] == 0) {
  #       Emissions_water_2022$Water_level_piezo[i]
  #     } else {Emissions_water_2022$Water_level_ruler[i]}
  #     } else {NA}, # a)
  #   if(is.na(Emissions_water_2022$Water_level_ruler[i]) & !is.na(Emissions_water_2022$Water_level_piezo[i])) {
  #     Emissions_water_2022$Water_level_piezo[i]
  #   } else {NA}, # b)
  #   if(!is.na(Emissions_water_2022$Water_level_ruler[i]) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Water_level_ruler[i] != 0) {
  #     Emissions_water_2022$Water_level_ruler[i]
  #   } else {NA}, # c)
  #   if((is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% c("AWD", "MSD")) {
  #     if(is.na(closest_neg_index)) {
  #       Emissions_water_2022$Water_level_ruler[i]
  #     } else {Emissions_water_2022$Water_level_piezo[closest_neg_index]}
  #   } else {NA}, # d.i)
  #   if((is.na(Emissions_water_2022$Water_level_ruler[i]) | Emissions_water_2022$Water_level_ruler[i] == 0) & is.na(Emissions_water_2022$Water_level_piezo[i]) & Emissions_water_2022$Treat[i] %in% "CON") {
  #     if(is.na(closest_neg_index_CON)) {
  #       Emissions_water_2022$Water_level_ruler[i]
  #     } else {Emissions_water_2022$Water_level_piezo[closest_neg_index_CON]}
  #   } else {NA}, # d.ii)
  #   TRUE ~ NA_real_
  # )
  # 
  
  
  Emissions_water_2022$Water_level_ruler[i]
  Emissions_water_2022$Water_level_piezo[i]
  Emissions_water_2022$Treat[i]