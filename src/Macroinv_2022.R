
####################################################### Thesis Paper 1 - Macroinvertebrates #####################################################

library(dplyr)
library(rlang)
library(ggplot2)
library(readr)
library(iNEXT)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library('unikn') 
library(rsvg)
library(png)

Macroinv_2022 <- read.csv("data/BIO/Macroinv_2022.csv", fileEncoding="latin1", na.strings=c("","NA"))

##############  1. Working initial sampling data #################

## Assigning repetitions to plots
Plot <- c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P13", "P14", "P15")
Rep <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5")
Plot_Rep <- data.frame(Plot, Rep)
Macroinv_2022 <- merge(Macroinv_2022, Plot_Rep, by.x="Plot", by.y="Plot")
Macroinv_2022  <- Macroinv_2022 [, c(2, 3, 1, 4, 19, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)] # Reorder columns

#### 1.1. Modifications regarding to UdG stay with Dani Boix. ####

Macroinv_2022 <- Macroinv_2022 %>% 
  mutate(Species = case_when(Genus == "Microvelia" ~ "pygmaea", TRUE ~ Species),  # i. All Microvelia found are Microvelia pygmaea
         Species = case_when(Genus == "Anisops" ~ "sardeus", TRUE ~ Species), # ii. All Anisops found are Anisops sardeus
         Genus = case_when(Family_SubFamily == "Curculionidae" ~ "Lissorhoptrus", TRUE ~ Genus),
         Species = case_when(Family_SubFamily == "Curculionidae" ~ "oryzophilus", TRUE ~ Species),  # iii. All Curculionidae found are Lissorhoptrus oryzophilus
         Genus = case_when(Family_SubFamily == "Planorbidae" ~ "Planorbis", TRUE ~ Genus),
         Species = case_when(Family_SubFamily == "Planorbidae" ~ "carinatus", TRUE ~ Species), # iv. All Planorbidae found are Planorbis carinatus
         Genus = case_when(Genus == "Dytiscus" ~ "Rhantus", TRUE ~ Genus),
         Species = case_when(Genus == "Rhantus" ~ "suturalis", TRUE ~ Species),# v. All Dytiscus found are Rhantus suturalis
         Genus = case_when(is.na(Genus) & Sampling == 1 & Plot == "P10" & Family_SubFamily == "Dytiscidae" ~ "Rhantus", TRUE ~ Genus),
         Species = case_when(is.na(Species) & Sampling == 1 & Plot == "P10" & Family_SubFamily == "Dytiscidae" ~ "suturalis", TRUE ~ Species),  # vi. P10 / Samp. 1: Dytiscidae is Rhantus suturalis
         Genus = case_when(is.na(Genus) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "Rhantus", TRUE ~ Genus),
         Species = case_when(is.na(Species) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "suturalis", TRUE ~ Species),
         Phylum_SubPhylum = case_when(is.na(Phylum_SubPhylum) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "Arthropoda", TRUE ~ Phylum_SubPhylum),
         Class_SubClass = case_when(is.na(Class_SubClass) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "Insecta", TRUE ~ Class_SubClass),
         Order_SubOrder = case_when(is.na(Order_SubOrder) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "Coleoptera", TRUE ~ Order_SubOrder),
         Family_SubFamily = case_when(is.na(Family_SubFamily) & Sampling == 1 & Plot == "P13" & Unclassified  == "Dytiscus sp" ~ "Dytiscidae", TRUE ~ Family_SubFamily),
         Unclassified  = case_when(Unclassified == "Dytiscus sp" & Sampling == 1 & Plot == "P13" ~ NA, TRUE ~ Unclassified), # vii. P13 / Samp. 1: Unclassified Dytiscus sp is Rhantus suturalis
         Stadium  = case_when(is.na(Stadium) & Genus == "Rhantus" ~ "Larvae", TRUE ~ Stadium), # viii. All remaining Rhantus suturalis are larvae
         Genus = case_when(is.na(Genus) & Family_SubFamily == "Culicinae" ~ "Culex", TRUE ~ Genus),
         Species = case_when(is.na(Species) & Family_SubFamily == "Culicinae" ~ "sp", TRUE ~ Species), # ix. All remaining Culicinae are Culex sp.
         Order_SubOrder = case_when(Family_SubFamily == "Daphniidae/ Moinidae" ~ "Cladocera", TRUE ~ Order_SubOrder),
         Class_SubClass = case_when(Family_SubFamily == "Daphniidae/ Moinidae" ~ "Branchiopoda", TRUE ~ Class_SubClass),
         Phylum_SubPhylum = case_when(Family_SubFamily == "Daphniidae/ Moinidae" ~ "Crustacea", TRUE ~ Phylum_SubPhylum), # x. Filling incomplete Daphniidae/ Moinidae row 1271
         Genus = case_when(Genus == "Bidessus/ Hydroglyphus" ~ "Hydroglyphus", TRUE ~ Genus), # xi. All Bidessus/ Hydroglyphus were Hydroglyphus sp.
         Genus = case_when(Genus == "Berosus " ~ "Enochrus", TRUE ~ Genus), # xi. Individual identified as Berosus was Enochrus sp.
         Genus = case_when(Genus == "Orthetrum" ~ "Sympetrum", TRUE ~ Genus),
         Species = case_when(Genus == "Orthetrum" ~ "fonscolombii", TRUE ~ Species), # xii. Individual identified as Orthetrum was Sympetrum fonscolombii (can still be Crocotemis, waiting for Dani)
         Genus = case_when(Genus == "Diplacodes" ~ "Sympetrum", TRUE ~ Genus),
         Species = case_when(Genus == "Species" ~ "fonscolombii", TRUE ~ Species), # xiii. Individual identified as Diplacodes lefrebvrii was Sympetrum fonscolombii (can still be Crocotemis, waiting for Dani)
         Species = case_when(Genus == "Radix" ~ "peregra", TRUE ~ Species), # xiv. All Radix found are Radix peregra
         Genus = case_when(Family_SubFamily == "Corixinae" ~ "Sigara", TRUE ~ Genus),
         Species = case_when(Family_SubFamily == "Corixinae" ~ "nigrolineata", TRUE ~ Species) # xv. All Corixinae found are Sigara nigrolineata (If Dani identifies some as a different species, discuss with him how to distribute, as all adults found were left with him)
  )

Macroinv_2022 <- subset(Macroinv_2022, Macroinv_2022$Unclassified != "Brichius elevatus" & 
                        Macroinv_2022$Unclassified != "Nysius sp" | is.na(Unclassified)) # xvi. Individuals named as "Brichius elevatus" and "Nysius sp" removed, they are terrestrial.
  
#### 1.2. Other modifications to original database ####

Macroinv_2022 <- Macroinv_2022 %>% 
                 mutate(Stadium = case_when(Order_SubOrder == "Coleoptera" & is.na(Stadium) ~ "Larvae", TRUE ~ Stadium), # i. Coleoptera without recorded stadium assigned as larvae.
                         Species = case_when(Species == "sp" ~ NA, TRUE ~ Species) # ii. All Species labelled as sp left as NA so assignations (1.3) apply to these rows as well. 
                  ) 

#### 1.3. Taxa Assignation ####

## - To observations with empty Family/Order/Species according to those in its same Sampling_date/Treat, weighted and keeping track of Stadium.
## - Developed in RScript Assignation_test.

Macroinv_2022$distID <- paste0(Macroinv_2022$Sampling, Macroinv_2022$Treat, Macroinv_2022$Order_SubOrder, Macroinv_2022$Stadium) # create distID

######### Loop for Family_SubFamily #####

Macroinv_2022_assign_Fam <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                       Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                       Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                       Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                       Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)
                                       )# Initialize an empty output DataFrame

for (i in 1:nrow(Macroinv_2022)) { # Iterate through rows in the original DataFrame
  row <- Macroinv_2022[i, ]
  
  distribution_values <- Macroinv_2022[Macroinv_2022$distID == row$distID & !is.na(Macroinv_2022$Family_SubFamily), ] %>%    # Get the distribution values
                          group_by(Family_SubFamily) %>%
                          summarise(Weighted_Abundance = sum(Weighted_Abundance))
                        
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    Macroinv_2022_assign_Fam <- rbind(Macroinv_2022_assign_Fam, row)
  } 
  
  else if (!is.na(row$Family_SubFamily)) {   # If Species is not NA, add the row to the output DataFrame
      Macroinv_2022_assign_Fam <- rbind(Macroinv_2022_assign_Fam, row)
      
  } else {

    total_weight <- sum(distribution_values$Weighted_Abundance)    # Distribute Weighted_Abundance based on weights
    weights <- distribution_values$Weighted_Abundance / total_weight
    new_Weighted_Abundance_values <- round(row$Weighted_Abundance * weights)
    
    new_rows <- data.frame(    # Create new rows in the output DataFrame
                          Date = rep(row$Date, nrow(distribution_values)), 
                          Sampling = rep(row$Sampling, nrow(distribution_values)), 
                          Plot = rep(row$Plot, nrow(distribution_values)),
                          Treat = rep(row$Treat, nrow(distribution_values)), 
                          Rep = rep(row$Rep, nrow(distribution_values)),
                          Sub_sample = rep(row$Sub_sample, nrow(distribution_values)), 
                          Weight = rep(row$Weight, nrow(distribution_values)), 
                          Phylum_SubPhylum = rep(row$Phylum_SubPhylum, nrow(distribution_values)),
                          Class_SubClass = rep(row$Class_SubClass, nrow(distribution_values)), 
                          Order_SubOrder = rep(row$Order_SubOrder, nrow(distribution_values)), 
                          Genus = rep(row$Genus, nrow(distribution_values)),
                          Species = rep(row$Species, nrow(distribution_values)), 
                          Stadium = rep(row$Stadium, nrow(distribution_values)), 
                          Unclassified = rep(row$Unclassified, nrow(distribution_values)),
                          Trophic_level = rep(row$Trophic_level, nrow(distribution_values)), 
                          Abundance = rep(row$Abundance, nrow(distribution_values)), 
                          distID = rep(row$distID, nrow(distribution_values)),
                          Obs = rep(row$Obs, nrow(distribution_values)),
                          Family_SubFamily = distribution_values$Family_SubFamily, # Applies distribution criteria
                          Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
                          )
    
    Macroinv_2022_assign_Fam <- rbind(Macroinv_2022_assign_Fam, new_rows)
  }
}

rownames(Macroinv_2022_assign_Fam) <- NULL # Reset row names and arrange output DataFrame
Macroinv_2022_assign_Fam <- Macroinv_2022_assign_Fam[order(Macroinv_2022_assign_Fam$distID), ]

######### Loop for Genus #####

Macroinv_2022_assign_Genus <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                       Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                       Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                       Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                       Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)) # Initialize an empty output DataFrame

for (i in 1:nrow(Macroinv_2022_assign_Fam)) { # Iterate through rows in the original DataFrame
  row <- Macroinv_2022_assign_Fam[i, ]
  
  distribution_values <- Macroinv_2022_assign_Fam[Macroinv_2022_assign_Fam$Family_SubFamily == row$Family_SubFamily & Macroinv_2022_assign_Fam$distID == row$distID & !is.na(Macroinv_2022_assign_Fam$Genus), ] %>%    # Get the distribution values # IMPORTANT: AGGREGATES FOR FAMILY EQUAL TO THAT OF THE ITERATING ROW
                                       group_by(Genus) %>%
                                       summarise(Weighted_Abundance = sum(Weighted_Abundance))
                        
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    Macroinv_2022_assign_Genus <- rbind(Macroinv_2022_assign_Genus, row)
  } 
  
  else if (!is.na(row$Genus)) {   # If Species is not NA, add the row to the output DataFrame
    Macroinv_2022_assign_Genus <- rbind(Macroinv_2022_assign_Genus, row)
  
  } else {
    
    total_weight <- sum(distribution_values$Weighted_Abundance)    # Distribute Weighted_Abundance based on weights
    weights <- distribution_values$Weighted_Abundance / total_weight
    new_Weighted_Abundance_values <- round(row$Weighted_Abundance * weights)
    
    new_rows <- data.frame(    # Create new rows in the output DataFrame
                          Date = rep(row$Date, nrow(distribution_values)), 
                          Sampling = rep(row$Sampling, nrow(distribution_values)), 
                          Plot = rep(row$Plot, nrow(distribution_values)),
                          Treat = rep(row$Treat, nrow(distribution_values)),
                          Rep = rep(row$Rep, nrow(distribution_values)),
                          Sub_sample = rep(row$Sub_sample, nrow(distribution_values)), 
                          Weight = rep(row$Weight, nrow(distribution_values)), 
                          Phylum_SubPhylum = rep(row$Phylum_SubPhylum, nrow(distribution_values)),
                          Class_SubClass = rep(row$Class_SubClass, nrow(distribution_values)), 
                          Order_SubOrder = rep(row$Order_SubOrder, nrow(distribution_values)), 
                          Family_SubFamily = rep(row$Family_SubFamily, nrow(distribution_values)),
                          Species = rep(row$Species, nrow(distribution_values)), 
                          Stadium = rep(row$Stadium, nrow(distribution_values)), 
                          Unclassified = rep(row$Unclassified, nrow(distribution_values)),
                          Trophic_level = rep(row$Trophic_level, nrow(distribution_values)), 
                          Abundance = rep(row$Abundance, nrow(distribution_values)), 
                          distID = rep(row$distID, nrow(distribution_values)),
                          Obs = rep(row$Obs, nrow(distribution_values)),
                          Genus = distribution_values$Genus,  # Applies distribution criteria
                          Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
                        )
    
    Macroinv_2022_assign_Genus <- rbind(Macroinv_2022_assign_Genus, new_rows)
  }
}

rownames(Macroinv_2022_assign_Genus) <- NULL # Reset row names and arrange output DataFrame
Macroinv_2022_assign_Genus <- Macroinv_2022_assign_Genus[order(Macroinv_2022_assign_Genus$distID), ]

######### Loop for Species #####

Macroinv_2022_assign_Species <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                         Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                         Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                         Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                         Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)
                                         
)# Initialize an empty output DataFrame

for (i in 1:nrow(Macroinv_2022_assign_Genus)) { # Iterate through rows in the original DataFrame
  row <- Macroinv_2022_assign_Genus[i, ]
  
  distribution_values <- Macroinv_2022_assign_Genus[Macroinv_2022_assign_Genus$Family_SubFamily == row$Family_SubFamily & Macroinv_2022_assign_Genus$Genus == row$Genus & Macroinv_2022_assign_Genus$distID == row$distID & !is.na(Macroinv_2022_assign_Genus$Species), ] %>%    # Get the distribution values # IMPORTANT: AGGREGATES FOR FAMILY EQUAL TO THAT OF THE ITERATING ROW
                          group_by(Species) %>%
                          summarise(Weighted_Abundance = sum(Weighted_Abundance))
  
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    Macroinv_2022_assign_Species <- rbind(Macroinv_2022_assign_Species, row)
  } 
  
  else if (!is.na(row$Species)) {   # If Species is not NA, add the row to the output DataFrame
    Macroinv_2022_assign_Species <- rbind(Macroinv_2022_assign_Species, row)
    
  } else {
    
    total_weight <- sum(distribution_values$Weighted_Abundance)    # Distribute Weighted_Abundance based on weights
    weights <- distribution_values$Weighted_Abundance / total_weight
    new_Weighted_Abundance_values <- round(row$Weighted_Abundance * weights)
    
    new_rows <- data.frame(    # Create new rows in the output DataFrame
                          Date = rep(row$Date, nrow(distribution_values)), 
                          Sampling = rep(row$Sampling, nrow(distribution_values)), 
                          Plot = rep(row$Plot, nrow(distribution_values)),
                          Treat = rep(row$Treat, nrow(distribution_values)),
                          Rep = rep(row$Rep, nrow(distribution_values)),
                          Sub_sample = rep(row$Sub_sample, nrow(distribution_values)), 
                          Weight = rep(row$Weight, nrow(distribution_values)), 
                          Phylum_SubPhylum = rep(row$Phylum_SubPhylum, nrow(distribution_values)),
                          Class_SubClass = rep(row$Class_SubClass, nrow(distribution_values)), 
                          Order_SubOrder = rep(row$Order_SubOrder, nrow(distribution_values)), 
                          Family_SubFamily = rep(row$Family_SubFamily, nrow(distribution_values)),
                          Genus = rep(row$Genus, nrow(distribution_values)), 
                          Stadium = rep(row$Stadium, nrow(distribution_values)), 
                          Unclassified = rep(row$Unclassified, nrow(distribution_values)),
                          Trophic_level = rep(row$Trophic_level, nrow(distribution_values)), 
                          Abundance = rep(row$Abundance, nrow(distribution_values)), 
                          distID = rep(row$distID, nrow(distribution_values)),
                          Obs = rep(row$Obs, nrow(distribution_values)),
                          Species = distribution_values$Species,  # Applies distribution criteria
                          Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
                        )
    
    Macroinv_2022_assign_Species <- rbind(Macroinv_2022_assign_Species, new_rows)
  }
}

rownames(Macroinv_2022_assign_Species) <- NULL # Reset row names and arrange output DataFrame
Macroinv_2022_assign_Species <- Macroinv_2022_assign_Species[order(Macroinv_2022_assign_Species$distID), ]

#### 1.4. Merge abundances of duplicated rows to avoid false Hill Number results. ####
## In case there are any, these must be merged into one row, adding up Weigthed_Abundance.

Macroinv_2022_assign_ALL <- Macroinv_2022_assign_Species %>%
                group_by(Date, Plot, Rep, Treat, Family_SubFamily, Genus, Species, Unclassified, Stadium) %>%
                summarise(Weighted_Abundance = sum(Weighted_Abundance)) %>% ungroup() %>%
                left_join(Macroinv_2022_assign_Species %>% 
                     select(Date, Sampling, Plot, Rep, Treat, Phylum_SubPhylum, Sub_sample, Weight, Class_SubClass, Order_SubOrder, 
                            Family_SubFamily, Genus, Species, Stadium, Unclassified, Trophic_level, Obs), 
                            by = c("Date","Plot", "Rep", "Treat", "Family_SubFamily", "Genus", "Species", "Unclassified", "Stadium"), multiple = "all") %>% 
                            distinct(Date, Plot, Rep, Treat, Family_SubFamily, Genus, Species, Unclassified, Stadium, .keep_all = TRUE)

Macroinv_2022_assign_ALL  <- Macroinv_2022_assign_ALL [, c(1, 11, 2, 4, 3, 13, 14, 12, 15, 16, 5, 6, 7, 8, 9, 10, 17, 18)] # Reorder columns

## Review: 

# - Check if distributions are correct.
# - Check cases in which there is only Family (All Genus and Species are NA), e.g. Curculionidae
# - Find review tables to compare assignations of each Dist_example_2 ... .csv in C:\Users\SECHEVERRIA\Documents\R\Thesis_Paper_1\data\Assignation_test\Dist_example.xlsx

sum(Macroinv_2022$Weighted_Abundance) # Check 1       
sum(Macroinv_2022_assign_Fam$Weighted_Abundance) # Check 2
sum(Macroinv_2022_assign_Genus$Weighted_Abundance) # Check 3
sum(Macroinv_2022_assign_Species$Weighted_Abundance) # Check 4 
sum(Macroinv_2022_assign_ALL$Weighted_Abundance) # Check 5 

# Output: .csv with all the previous corrections and assignations 

#### 1.5. Create Taxres_max before applying Hill Numbers ####

Macroinv_2022_assign_ALL <- Macroinv_2022_assign_ALL %>% 
  mutate(Taxres_max = 
           case_when(
             !is.na(Species) ~ Species,
             is.na(Species) & !is.na(Genus) ~ Genus,
             is.na(Species) & is.na(Genus) & !is.na(Family_SubFamily) ~ Family_SubFamily,
             is.na(Species) & is.na(Genus) & is.na(Family_SubFamily) & !is.na(Order_SubOrder) ~ Order_SubOrder,
             is.na(Species) & is.na(Genus) & is.na(Family_SubFamily) & is.na(Order_SubOrder) & !is.na(Class_SubClass) ~ Class_SubClass,
             is.na(Species) & is.na(Genus) & is.na(Family_SubFamily) & is.na(Order_SubOrder) & is.na(Class_SubClass) & !is.na(Phylum_SubPhylum) ~ Phylum_SubPhylum, 
             TRUE ~ NA
         ))

write_csv2(Macroinv_2022_assign_ALL, "outputs/csv/BIO/Macroinv_2022_assign_ALL.csv")

##############  2. Diversity Analysis - Hill Numbers #################

#### 2.1. Define extrar_divmetrics function ####

extrar_divmetrics <- function (datafile){
  
  list.sites <- datafile %>%
    group_by(siteID) %>%
    summarise(list(Abundance))
  
  list.sites2 <- list.sites[[2]]
  
  names(list.sites2) <- list.sites[[1]]
  
  x <- iNEXT::iNEXT(list.sites2, q=0, datatype="abundance",size=NULL)
  #ggiNEXT(x, type=1, se=TRUE, facet.var="none", color.var="none", grey=FALSE)
  
  inextparams <- as.data.frame(x[3]) %>%
    dplyr::select(siteID = AsyEst.Assemblage,
                  divmetric = AsyEst.Diversity,
                  obsmetric = AsyEst.Observed,
                  estmetric = AsyEst.Estimator,
                  se.metric = AsyEst.s.e.)
  
  rich.df <- inextparams %>%
    filter(divmetric=="Species richness") %>%
    dplyr::select(siteID,
                  q0.obs = obsmetric,
                  q0.est = estmetric,
                  q0.se = se.metric)
  
  shannon.df <- inextparams %>%
    filter(divmetric=="Shannon diversity") %>%
    dplyr::select(siteID,
                  q1.obs = obsmetric,
                  q1.est = estmetric,
                  q1.se = se.metric)
  
  simpson.df <- inextparams %>%
    filter(divmetric=="Simpson diversity") %>%
    dplyr::select(siteID,
                  q2.obs = obsmetric,
                  q2.est = estmetric,
                  q2.se = se.metric)
  
  div.params <- rich.df %>%
    left_join(.,shannon.df, by="siteID") %>%
    left_join(.,simpson.df,by="siteID")
  
  dataname <-  deparse(substitute(datafile))
  write_csv(div.params,paste0("data/",dataname,"_inext_params.csv"))

}

#### 2.2. Select and prepare dataset for diversity metrics ####
## Checking Order_SubOrder to work at diversity level

Tax_lev <- Macroinv_2022_assign_ALL %>% 
             group_by(Order_SubOrder) %>% 
             summarise(Family_SubFamily = length(unique(Family_SubFamily)), Genus = length(unique(Genus)), Species = length(unique(Species)))

ColOdoHet <- subset(Macroinv_2022_assign_ALL, Order_SubOrder %in% c("Coleoptera", "Heteroptera", "Odonata") & (is.na(Stadium) | Stadium != "Adult")) %>%  # Mobile adults out of analysis
             group_by(Sampling, Treat, Plot, Taxres_max) %>% 
             summarise(Abundance = mean(Weighted_Abundance))

ColOdoHet <- ColOdoHet %>% 
             mutate(siteID = paste0(Plot, "_", Sampling, "_", Treat)) # Create siteID column to use as variable within extrar_divmetrics

write_csv2(ColOdoHet, "outputs/csv/BIO/ColOdoHet.csv")

#### 2.3 Calculate Hill numbers through extrar_divmetrics function ####

ColOdoHet_inext_params <- extrar_divmetrics(ColOdoHet)

write_csv(ColOdoHet_inext_params, "outputs/csv/BIO/ColOdoHet_inext_params.csv")

Hills_ColOdoHet <- ColOdoHet %>% 
                    group_by(Plot, Sampling, Treat, siteID) %>% 
                    summarise(siteID = first(siteID),
                              Plot = first(Plot),
                              Sampling = first(Sampling),
                              Treat = first(Treat)) %>% 
                    left_join(ColOdoHet_inext_params, by = "siteID")

Hills_ColOdoHet$Treat <- factor(Hills_ColOdoHet$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable

write_csv(Hills_ColOdoHet, "outputs/csv/BIO/Hills_ColOdoHet.csv")

#  To plot Annexes with plots for Col, Het and Odo separate (see: ColOdoHet_q0q1_2022_all), see script:  ColOdoHet_separate

#### 2.4 Plot Hill Numbers ####

ColOdoHet_summary_q0 <- Hills_ColOdoHet %>%
                      group_by(Treat) %>%
                      summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

ColOdoHet_summary_q1 <- Hills_ColOdoHet %>%
                      group_by(Treat) %>%
                      summarise(mean_q1.obs = mean(q1.obs),se_q1.obs = sd(q1.obs) / sqrt(n()))

## Individual plots:

ColOdoHet_plot_SpRich_indiv <- ggplot(Hills_ColOdoHet, 
                                aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                            geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10)+
                            scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                            scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                            theme_bw() +
                            ylab("") +
                            ggtitle(expression("Species richness (q"[0]*")")) +
                            theme(plot.title = element_text(size=20, hjust=0.5)) +
                            theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                  axis.title.y = element_text(size = 20, margin = margin(r = 12)), axis.title.x = element_blank(), legend.position = "none", 
                                  axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                            geom_point(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs), shape = 19, colour = "black", size = 12) +
                            geom_point(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 10) +
                            geom_errorbar(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.3, size = 1) +
                            scale_y_continuous(limits = c(1, 12), breaks = seq(2, 12, by = 2)) 

print(ColOdoHet_plot_SpRich_indiv)

ggsave("outputs/Plots/BIO/SpRich_indiv.pdf", plot = ColOdoHet_plot_SpRich_indiv ,width = 10, height = 10)

ColOdoHet_plot_Shannon_indiv <- ggplot(Hills_ColOdoHet, 
                                 aes(Treat, q1.obs, group = Treat, colour = Treat, fill = Treat)) +
                            geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10)+
                            scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                            scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                            # scale_y_sqrt() +
                            theme_bw() +
                            # xlab("Treatment") +
                            ylab("") +
                            ggtitle(expression("Shannon diversity (q"[1]*")")) +
                            theme(plot.title = element_text(size=20, hjust=0.5)) +
                            theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                  axis.title.y = element_text(size = 20, margin = margin(r = 12)), axis.title.x = element_blank(), legend.position = "none", 
                                  axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                            geom_point(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs), shape = 19, colour = "black", size = 12) +
                            geom_point(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs), shape = 19, size = 10) +
                            geom_errorbar(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs, ymin = mean_q1.obs - se_q1.obs, ymax = mean_q1.obs + se_q1.obs), width = 0.3, size = 1) 

print(ColOdoHet_plot_Shannon_indiv)

ggsave("outputs/Plots/BIO/Shannon_indiv.pdf", plot = ColOdoHet_plot_Shannon_indiv ,width = 10, height = 10)

## Arranging individual plots:

ColOdoHet_inddivplots <- grid.arrange(arrangeGrob(ColOdoHet_plot_SpRich_indiv, ColOdoHet_plot_Shannon_indiv, nrow = 1, ncol = 2))

ggsave("outputs/Plots/BIO/ColOdoHet_indiv_arrange.pdf", plot = ColOdoHet_inddivplots ,width = 20, height = 10)

## Plots for arrange:

ColOdoHet_plot_SpRich <- ggplot(Hills_ColOdoHet, 
                                aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                            geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10)+
                            scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                            scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                            # scale_y_sqrt() +
                            theme_bw() +
                            # xlab("Treatment") +
                            ylab(expression("Species richness (q"[0]*")")) +
                            ggtitle("") +
                            theme(plot.title = element_text(size=20, hjust=0.5)) +
                            theme(axis.title.x = element_blank(), axis.text = element_text(size = 14), strip.text = element_text(size = 14), 
                                  axis.title.y = element_text(size = 20, margin = margin(r = 0)), legend.position = "none", axis.text.x = element_blank(),
                                  axis.text.y = element_text(size = 20, margin = margin(r = 0)), panel.border = element_rect(size = 1)) +
                            geom_point(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs), shape = 19, colour = "black", size = 12) +
                            geom_point(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 10) +
                            geom_errorbar(data = ColOdoHet_summary_q0, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.3, size = 1) +
                            scale_y_continuous(limits = c(1, 12), breaks = seq(2, 12, by = 2)) 

print(ColOdoHet_plot_SpRich)

ColOdoHet_plot_Shannon <- ggplot(Hills_ColOdoHet, 
                           aes(Treat, q1.obs, group = Treat, colour = Treat, fill = Treat)) +
                            geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10)+
                            scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#FF5D5D")) +
                            scale_fill_manual(values = c("#002B5B", "#03C988", "#FF5D5D"), guide = "none") +
                            # scale_y_sqrt() +
                            theme_bw() +
                            # xlab("Treatment") +
                            ylab(expression("Shannon index (q"[1]*")")) +
                            theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                  axis.title.y = element_text(size = 20, margin = margin(r = 12)), axis.title.x = element_blank(), legend.position = "none", 
                                  axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                            geom_point(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs), shape = 19, colour = "black", size = 12) +
                            geom_point(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs), shape = 19, size = 10) +
                            geom_errorbar(data = ColOdoHet_summary_q1, aes(x = Treat, y = mean_q1.obs, ymin = mean_q1.obs - se_q1.obs, ymax = mean_q1.obs + se_q1.obs), width = 0.3, size = 1) 

print(ColOdoHet_plot_Shannon)

# Arrange plots:

ColOdoHet_divplots <- grid.arrange(arrangeGrob(ColOdoHet_plot_SpRich, ColOdoHet_plot_Shannon, nrow = 2, ncol = 1))

ggsave("outputs/Plots/BIO/ColOdoHet_arrange.pdf", plot = ColOdoHet_divplots ,width = 10, height = 20)

#### 2.5 Prepare data for abundance plots ####

Sam.Date <- unique(Macroinv_2022[, c("Sampling", "Date")])
Order.taxres_max <- unique(Macroinv_2022_assign_ALL[, c("Taxres_max", "Order_SubOrder")])
 
ColOdoHet_merged <- subset(Macroinv_2022_assign_ALL, Order_SubOrder %in% c("Coleoptera", "Heteroptera", "Odonata") & (is.na(Stadium) | Stadium != "Adult"))  %>%  # Mobile adults out of analysis
                    select(Sampling, Treat, Plot, Taxres_max, Weighted_Abundance) 

colnames(ColOdoHet_merged)[colnames(ColOdoHet_merged) == "Weighted_Abundance"] <- "Abundance" 

ColOdoHet_merged <- merge(ColOdoHet_merged, Sam.Date[, c("Sampling", "Date")], by = "Sampling", all.x = TRUE)
ColOdoHet_merged <- merge(ColOdoHet_merged, Order.taxres_max[, c("Taxres_max", "Order_SubOrder")], by = "Taxres_max", all.x = TRUE) %>% 
                    select(Date, Plot, Treat, Order_SubOrder, Abundance)

Abundance_2022 <- read.csv("data/BIO/Macrofauna_2022.csv", fileEncoding="latin1", na.strings=c("","NA")) %>% # Import Macrofauna_2022 data
                  filter(!(c(Organism == "Pelophylax perezi" & Stadium == "Adult"))) %>% # Removes frog adults (leaving only tadpoles) and rows with empty traps
                  filter(!(c(Type == "Coleoptera" & Stadium == "Adult"))) %>%  # Removes Coleoptera adults
                  filter(!(Type == "Cricket")) %>%  # Removes Crickets 
                  filter(!(Type == "Diptera")) # Removes Diptera

colnames(Abundance_2022)[colnames(Abundance_2022) == "Type"] <- "Order_SubOrder" 

Abundance_2022 <- Abundance_2022 %>% 
                  select(Date, Plot, Treat, Order_SubOrder, Abundance)  
Abundance_2022 <- rbind(Abundance_2022, ColOdoHet_merged) %>% 
                  group_by(Date , Plot, Treat, Order_SubOrder) %>% 
                  summarise(Abundance = sum(Abundance)) %>% 
                  ungroup()

Abundance_2022$Treat <- factor(Abundance_2022$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable

write_csv2(Abundance_2022, "outputs/csv/BIO/Abundance_2022.csv")

acc_Abundance_2022 <- Abundance_2022 %>% # Creates dataframe with accumulated abundances per Treat and Order_SubOrder
                      group_by(Treat, Order_SubOrder) %>% 
                      summarise(Abundance = sum(Abundance))

#### 2.6 Plot abundance ####

Abundance_2022_plot <- ggplot(acc_Abundance_2022, aes(Treat, Abundance, group = Order_SubOrder, colour = Order_SubOrder, fill = Order_SubOrder, shape = Order_SubOrder)) +
                              geom_point(colour = "black", size = 13) +
                              geom_point(aes(colour = Order_SubOrder), size = 11) +
                              scale_colour_manual(values = c("#071952", "#E9B824", "#5BC0F8", "#D83F31", "#35A29F", "#F6FA70")) +
                              scale_fill_manual(values = c("#071952", "#E9B824", "#5BC0F8", "#D83F31", "#35A29F", "#F6FA70"), guide = "none") +
                              geom_line(aes(colour = Order_SubOrder), linetype = "dashed", linewidth = 1.5) +
                              scale_shape_manual(values = c(20, 23, 15, 16, 17, 18)) +
                              theme_bw() +
                              labs(title = "Accumulated abundance", y = NULL) +
                              theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                    axis.title.y = element_text(size = 20, margin = margin(r = 8)), axis.title.x = element_blank(), 
                                    axis.text.y = element_text(size = 20, margin = margin(r = 0), angle = 90), legend.position = c(0.9, 0.87), legend.title = element_blank(),
                                    legend.background = element_rect(fill="white", size = 0.7, linetype="solid", colour = "black"), 
                                    legend.text = element_text(colour="black", size = 15),  axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1),
                                    plot.title = element_text(hjust = 0.5, size = 20)) +  
                              # scale_y_continuous(breaks = seq(0, 30, by = 5)) +
                              guides(color = guide_legend(override.aes = list(size = 6, vjust = 10), byrow = TRUE)) +
                              scale_y_sqrt()  
                              # annotation_custom(rasterGrob(Hydrophilidae_png), xmin = 0.05, xmax = 0.1, ymin = 0.4, ymax = 0.5)

print(Abundance_2022_plot)  

ggsave("outputs/Plots/BIO/Abundance_2022.pdf", plot = Abundance_2022_plot, width = 10, height = 10)

## Arrange oh diversity plots and abundance plot:

Abu.div_2022_plots <- grid.arrange(arrangeGrob(ColOdoHet_divplots, Abundance_2022_plot, ncol = 2, nrow = 1))

# ggsave("outputs/Plots/BIO/Abu.div_2022_plots.pdf", plot = Abu.div_2022_plots, width = 20, height = 10)

