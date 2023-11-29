
####################################################### Thesis Paper 1 - Assignation_test #####################################################

library(dplyr)

### - Assigns Taxa to empty Family/Order/Species according to those non NA with same Sampling_date/Treat/Order/Stadium (assigns Taxa to all observations with empty fields from Family onward, only if there are non NA with same Sampling_date/Treat/Order/Stadium). 
### - Abundance is weighted according to these non NA rows. 

######### I. Assignation through loops #####
## Chain of loops: (1) Family, (2) Genus and (3) Species. Each loop creates an empty dataframe and uses the previous loop output to apply conditions through a "row" element.

## Trying with a simple example: THIS CODE WORKS ALREADY (for "data/Assignation_test/Dist_example_2.csv", "...Dist_example_2.2.csv", "...2.3" ), examples were built under C:\Users\SECHEVERRIA\Documents\R\Thesis_Paper_1\data\Assignation_test\Dist_example.xlsx

dist_eg_2 <-  read.csv("data/Assignation_test/Dist_example_2.4.csv", fileEncoding="latin1", na.strings=c("","NA"))

######### Loop for Family_SubFamily #####

dist_eg_2_assigned <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                 Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                 siteID = character(0), Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                 Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                 Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)
                                 )# Initialize an empty output DataFrame

for (i in 1:nrow(dist_eg_2)) { # Iterate through rows in the original DataFrame
  row <- dist_eg_2[i, ]
  
  distribution_values <- dist_eg_2[dist_eg_2$distID == row$distID & !is.na(dist_eg_2$Family_SubFamily), ] %>%    # Get the distribution values
                          group_by(Family_SubFamily) %>%
                          summarise(Weighted_Abundance = sum(Weighted_Abundance))
  
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    dist_eg_2_assigned <- rbind(dist_eg_2_assigned, row)
  } else if (!is.na(row$Family_SubFamily)) {   # If Species is not NA, add the row to the output DataFrame
    dist_eg_2_assigned <- rbind(dist_eg_2_assigned, row)
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
                          siteID = rep(row$siteID, nrow(distribution_values)),
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
                          Family_SubFamily = distribution_values$Family_SubFamily, # Applies distribution criteria
                          Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
    )
    
    dist_eg_2_assigned <- rbind(dist_eg_2_assigned, new_rows)
  }
}

rownames(dist_eg_2_assigned) <- NULL # Reset row names and arrange output DataFrame
dist_eg_2_assigned <- dist_eg_2_assigned[order(dist_eg_2_assigned$distID), ]

######### Loop for Genus #####

dist_eg_2_assigned_Genus <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                       Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                       siteID = character(0), Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                       Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                       Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)
                                       )# Initialize an empty output DataFrame

for (i in 1:nrow(dist_eg_2_assigned)) { # Iterate through rows in the original DataFrame
  row <- dist_eg_2_assigned[i, ]
  
  distribution_values <- dist_eg_2_assigned[dist_eg_2_assigned$Family_SubFamily == row$Family_SubFamily & dist_eg_2_assigned$distID == row$distID & !is.na(dist_eg_2_assigned$Genus), ] %>%    # Get the distribution values # IMPORTANT: AGGREGATES FOR FAMILY EQUAL TO THAT OF THE ITERATING ROW
                        group_by(Genus) %>%
                        summarise(Weighted_Abundance = sum(Weighted_Abundance))
  
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    dist_eg_2_assigned_Genus <- rbind(dist_eg_2_assigned_Genus, row)
  } else if (!is.na(row$Genus)) {   # If Species is not NA, add the row to the output DataFrame
    dist_eg_2_assigned_Genus <- rbind(dist_eg_2_assigned_Genus, row)
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
                        siteID = rep(row$siteID, nrow(distribution_values)),
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
                        Genus = distribution_values$Genus,  # Applies distribution criteria
                        Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
    )
    
    dist_eg_2_assigned_Genus <- rbind(dist_eg_2_assigned_Genus, new_rows)
  }
}

rownames(dist_eg_2_assigned_Genus) <- NULL # Reset row names and arrange output DataFrame
dist_eg_2_assigned_Genus <- dist_eg_2_assigned_Genus[order(dist_eg_2_assigned_Genus$distID), ]

######### Loop for Species #####

dist_eg_2_assigned_Species <- data.frame(distID = character(0), Family_SubFamily = character(0), Weighted_Abundance = numeric(0), # Original
                                         Date = numeric(0), Sampling = numeric(0), Plot = character(0), Treat = character(0), Rep = numeric(0),
                                         siteID = character(0), Sub_sample = character(0), Phylum_SubPhylum = character(0), Class_SubClass = character(0),
                                         Order_SubOrder = character(0), Genus = character(0), Species = character(0), Stadium = character(0), 
                                         Unclassified = character(0),Trophic_level = character(0), Abundance = numeric(0)
                                         )# Initialize an empty output DataFrame

for (i in 1:nrow(dist_eg_2_assigned_Genus)) { # Iterate through rows in the original DataFrame
  row <- dist_eg_2_assigned_Genus[i, ]

  distribution_values <- dist_eg_2_assigned_Genus[dist_eg_2_assigned_Genus$Family_SubFamily == row$Family_SubFamily & dist_eg_2_assigned_Genus$Genus == row$Genus & dist_eg_2_assigned_Genus$distID == row$distID & !is.na(dist_eg_2_assigned_Genus$Species), ] %>%    # Get the distribution values # IMPORTANT: AGGREGATES FOR FAMILY and GENUS EQUAL TO THAT OF THE ITERATING ROW
                        group_by(Species) %>%
                        summarise(Weighted_Abundance = sum(Weighted_Abundance))
                      
  if (nrow(distribution_values)==0) { # If no assigned Species for this Genus, add the row to the output DataFrame
    dist_eg_2_assigned_Species <- rbind(dist_eg_2_assigned_Species, row)
  } else if (!is.na(row$Species)) {   # If Species is not NA, add the row to the output DataFrame
    dist_eg_2_assigned_Species <- rbind(dist_eg_2_assigned_Species, row)
  } else {
    
    total_weight <- sum(distribution_values$Weighted_Abundance)    # Distribute Weighted_Abundance based on weights
    weights <- distribution_values$Weighted_Abundance / total_weight
    new_Weighted_Abundance_values <- round(row$Weighted_Abundance * weights)
    
    new_rows <- data.frame(    # Creates new rows in the output DataFrame
                          Date = rep(row$Date, nrow(distribution_values)),
                          Sampling = rep(row$Sampling, nrow(distribution_values)),
                          Plot = rep(row$Plot, nrow(distribution_values)),
                          Treat = rep(row$Treat, nrow(distribution_values)),
                          Rep = rep(row$Rep, nrow(distribution_values)),
                          siteID = rep(row$siteID, nrow(distribution_values)),
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
                          Species = distribution_values$Species,  # Applies distribution criteria
                          Weighted_Abundance = new_Weighted_Abundance_values # Applies weights
              )
    
    dist_eg_2_assigned_Species <- rbind(dist_eg_2_assigned_Species, new_rows)
  }
 }  

rownames(dist_eg_2_assigned_Species) <- NULL # Reset row names and arrange output DataFrame
dist_eg_2_assigned_Species <- dist_eg_2_assigned_Species[order(dist_eg_2_assigned_Species$distID), ]

######### II. Merge abundances of duplicated rows to avoid false Hill Number and assignation results.: #####
## In case there are, these must be merged into one row, so these observations are aggregated into a distinct siteID and Hill Numbers are clean of duplicates.

dist_eg_2_assigned_Final <- dist_eg_2_assigned_Species %>%
                            group_by(Date, Plot, Rep, siteID, Treat, Family_SubFamily, Genus, Species, Unclassified, Stadium) %>%
                            summarize(Weighted_Abundance = sum(Weighted_Abundance)) %>% ungroup() %>%
                            left_join(dist_eg_2_assigned_Species %>% select(Date, Sampling, Plot, Rep, siteID, Treat, Phylum_SubPhylum, Sub_sample, Weight, 
                            Class_SubClass, Order_SubOrder, Family_SubFamily, Genus, Species, Stadium, Unclassified, Trophic_level), 
                            by = c("Date","Plot", "Rep", "siteID", "Treat", "Family_SubFamily", "Genus", "Species", "Unclassified", "Stadium"), multiple = "all") %>% 
                            distinct(Date, Plot, Rep, siteID, Treat, Family_SubFamily, Genus, Species, Unclassified, Stadium, .keep_all = TRUE)

## Review: 

# - Check if distributions are correct.
# - Check cases in which there is only Family (All Genus and Species are NA), e.g. Curculionidae
# - Find review tables to compare assignations of each Dist_example_2 ... .csv in C:\Users\SECHEVERRIA\Documents\R\Thesis_Paper_1\data\Assignation_test\Dist_example.xlsx

sum(dist_eg_2$Weighted_Abundance) # Check 1       
sum(dist_eg_2_assigned$Weighted_Abundance) # Check 2
sum(dist_eg_2_assigned_Genus$Weighted_Abundance) # Check 3
sum(dist_eg_2_assigned_Species$Weighted_Abundance) # Check 4 
sum(dist_eg_2_assigned_Final$Weighted_Abundance) # Check 5 

subset(dist_eg_2_assigned_Final, dist_eg_2_assigned_Final$Sampling == 2, dist_eg_2_assigned_Final$Treat == "AWD", dist_eg_2_assigned_Final$Family_SubFamily == "Coenagrionidae")

check <- subset(dist_eg_2_assigned_Species, Sampling == 2 & Treat == "CON" & Order_SubOrder == "Odonata")

check <- check %>% 
          group_by(Family_SubFamily, Genus, Species) %>%
          summarize(Weighted_Abundance = sum(Weighted_Abundance)) %>% ungroup()
