rm(list = ls())
library(tidyverse)

#Loading data
data.macroinv <- read_csv2("data-raw/macroinvertebrates_def.csv") 



#Visualize data structure
str(data.macroinv)

#Change cell names of the trophic level and treatment and create a new column
#that give us the maximal taxonomical resolution of each taxa
data.predators.ont <- data.macroinv %>% 
  mutate(ontogeny=case_when(ontogeny=="Adulto" ~ "Adult",
                                 ontogeny == "Juvenil" ~ "Larvae",
                                 TRUE ~ ontogeny)) %>% 
  filter(trophic_level2 == "Aquatic predators",
         order %in% c("Diptera", "Hemiptera", "Coleoptera", "Odonata"),
         ! taxres_max %in% c("Dytiscidae", "Hydrophilidae")) %>% 
  mutate(taxres_max = case_when(taxres_max == "Mesoveliidae" ~ "Mesovelia",
                                taxres_max == "Corixidae" ~ "Sigara",
                                taxres_max == "Micronecta" ~ "Sigara",
                                taxres_max == "Libellulidae" ~ "Orthetrum",
                                taxres_max == "Aeshnidae" ~ "Anax",
                                taxres_max == "Coenagrionidae" ~ "Ischnura",
                                TRUE ~ taxres_max))
  
write_csv2(data.predators.ont, "data/predabund.ont.csv")


# Summarizing data by taxres_max without taking account the ontogeny

data.predators <- data.predators.ont %>% 
  group_by(sampling_month, treatment, plot, taxres_max) %>% 
  summarise(abundance = mean(abundance, na.rm = T))

write_csv2(data.predators, "data/predabund.csv")

#######DIVERSITY DATA #############
source("R_functions/functions.R")


data.preddiv <- data.predators %>% 
  mutate(siteID = paste0(plot,sampling_month,treatment))

extrar_divmetrics(data.preddiv)  
hills <- read_csv("data-raw/data.preddiv_inext_params.csv")
merge.hills <- data.preddiv %>% 
  group_by(plot,sampling_month,treatment) %>% 
  summarise(siteID = first(siteID),
            plot = first(plot),
            sampling_month = first(sampling_month),
            treatment = first(treatment)) %>% 
  left_join(hills, by = "siteID")

write_csv2(merge.hills, "data/predhill.csv")






