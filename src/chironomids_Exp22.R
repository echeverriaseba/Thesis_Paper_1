rm(list = ls())

library(tidyverse)
library(glmmTMB)

#Loading data
chiro <- read_csv2("data/chironomids_clean_22.csv") 

#GLMM with the abundance of chironomids as response, treatment as fixed factor
#and site identity as random factor
m.chiro <- glmmTMB(data = chiro, Abundance ~ Treatment +
                          (1|siteID) , 
                        family = "poisson")
summary(m.chiro)
car::Anova(m.chiro)
#checking of model residuals
DHARMa::simulateResiduals(m.chiro, plot =T)

#estimating the R2 of the model
performance::r2(m.chiro)

#Obtaining the estimated means per treatment and the SE
emm <- emmeans::emmeans(
  m.chiro, ~ Treatment, type = "response"
  )
emm

#Obtaining the p-values
pairs(emmeans::emmeans(m.chiro, ~ Treatment, type = "response"))

#Creating a new dataframe with the output of the model (i.e., the mean and SE per treatment)
df.chiro <- as.data.frame(emmeans::emmeans(
  m.chiro, ~ Treatment, type = "response")
  )


#Plotting the observed values as a dotplot 
chiroplot1 <- ggplot(chiro, 
                     aes(
                       Treatment, 
                       Abundance, 
                       group = Treatment, 
                       colour = Treatment, 
                       fill = Treatment)
                     ) +
  #geom_violin(aes(fill = treatment), alpha = 0.2) +
  geom_point(
    position = position_jitterdodge (0.45, 
                                     jitter.width = 0.07, 
                                     jitter.height = 0), 
    alpha = 0.2,
    shape = 21,
    colour = "black",
    size = 5
    ) +
  #stat_summary(fun.data = mean_se, position = position_dodge(0.65),geom = "errorbar", width = 0.1) +
  #stat_summary(fun = mean, position = position_dodge(0.65), size = 1.8) +
  scale_colour_manual(
    name = "Management", 
    values = c("#F2A104", "#00743F")
    ) +
  scale_fill_manual(
    values = c("#F2A104", "#00743F"), 
    guide = "none"
    ) +
  scale_y_sqrt() +
  theme_bw() +
  xlab("Management") +
  ylab ("Chironomids density (individuals/core)") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14)
        )


#Overlaying the means and SE estimated by the model
chiroplot2 <- chiroplot1 + 
  geom_errorbar(
    mapping = aes(
      ymax =  rate + SE,
      ymin = rate - SE,
      x = Treatment,
      group = Treatment
    ),
    position = position_dodge(0.45),
    colour = "black",
    width = 0.05,
    data = df.chiro,
    inherit.aes = FALSE
  ) +
  geom_point(
    mapping = aes(
      y =  rate,
      x = Treatment,
      fill = Treatment),
    position = position_dodge(0.45),
    shape = 21,
    colour = "black",
    alpha = 1,
    size = 7,
    data = df.chiro,
    inherit.aes = FALSE
  ) 

chiroplot2

pdf("figures/chironomids_exp22.pdf", width = 4.13, height = 6.79)
chiroplot2
dev.off()

