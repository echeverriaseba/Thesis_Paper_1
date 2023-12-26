
####################################################### Thesis Paper 1 - Production Statistics #####################################################

library(dplyr)
library(ggplot2)

#### 1. Prod data ####

Prod_2022 <- read.csv("data/PROD/Prod_2022.csv", fileEncoding="latin1", na.strings=c("","NA")) %>%  # Import Prod_2022 data
  subset(Experiment == "Cerestres") 

Prod_2022_means <- Prod_2022 %>% 
  group_by(Treat) %>% 
  summarise(Yield_mean = mean(Yield_kg_ha), SE = sd(Yield_kg_ha) / sqrt(n())) %>% 
  ungroup()

Prod_2022$Treat <- factor(Prod_2022$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable

#### 2. Prod stats ####

lm.prod1 <- lm(Yield_kg_ha~Treat, data = Prod_2022)

# Model diagnostics:
DHARMa::simulateResiduals(lm.prod1, plot = T)
summary(lm.prod1)
car::Anova(lm.prod1)
performance::r2(lm.prod1)
performance::check_collinearity(lm.prod1)
performance::check_singularity(lm.prod1)
visreg(lm.prod1, scale="response") # Plotting conditional residuals

emmeans(lm.prod1, ~Treat , type = "response")
pairs(emmeans(lm.prod1, ~Treat , type = "response"))
