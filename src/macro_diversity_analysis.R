rm(list = ls())

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
#########################################################################################
##################################ANALYSIS PREDATOR DIVERSITY############################
#########################################################################################

parameters <- read_csv2("data-raw/parameters_data.csv") %>% 
  mutate(siteID = paste0(plot,sampling_month,treatment))

hills <- read_csv2("data/predhill.csv") %>% 
  mutate(sampling_order = case_when(sampling_month == "May" ~ 1,
                                    sampling_month == "June" ~ 2,
                                    sampling_month == "July" ~ 3,
                                    sampling_month == "August" ~ 4,
                                    sampling_month == "September" ~ 5),
         siteID = paste0(plot,sampling_month,treatment)) %>% 
  filter(sampling_month %in% c("June", "July", "August", "September"))
  

hills$sampling_month <- factor(hills$sampling_month, levels = c("May", "June", "July", "August", "September"))

hills2 <- hills %>% 
  left_join(., parameters, by = "siteID", suffix = c("", "duplicated")) %>%
  select(-ends_with(".duplicated")) %>% 
  mutate(Plot_trat = paste0(plot, treatment))

cormatrix <- hills2 %>% 
  select(q0.obs, q1.obs, q2.obs, ph,
         temperature,prop_oxi, conc_oxi, conductivity, salinity, size, sampling_order) %>% 
  na.omit

corr <- round(cor(cormatrix, method =  "spearman"), 1)

pdf("figures/Corr_plot.pdf", width = 11)
corrplot::corrplot.mixed(corr, order = 'hclust', addrect = 2)
dev.off()

# 
# pdf("figures/richness_alltrophicgroups.pdf", width = 11)
# ggplot(hills, aes(x = sampling_month, y = q1.obs, colour = treatment)) + 
#   geom_jitter(position = position_jitterdodge(0.2), size = 4, alpha = 0.5) +
#   #facet_wrap ( ~ trophic_level2, scales = "free") +
#   stat_summary(fun.data = mean_se, position = position_dodge(0.75),geom = "errorbar", width = 0.2) +
#   stat_summary(fun = mean, position = position_dodge(0.75), size = 1.2) +
#   scale_colour_manual(name = "Management", values = c("#F2A104", "#00743F")) +
#   theme_bw() +
#   xlab("Month") +
#   ylab ("Macroinvertebrate richness") +
#   theme(axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14),
#         strip.text = element_text(size = 14))
# dev.off()


#MODEL FOR RARE SPECIES (q0)
hist(hills2$q0.obs)

#########FULL MODEL##########
m.rare.full <- glmmTMB(data = hills2, q0.obs ~ treatment*sampling_month + ph + 
                    salinity + temperature + conc_oxi + (1|Plot_trat) , family = "poisson")

DHARMa::simulateResiduals(m.rare.full, plot = T)
summary(m.rare.full)
car::Anova(m.rare.full)
performance::r2(m.rare.full)
performance::check_collinearity(m.rare.full)


###########FINAL MODEL######## I remove ph and temperature as VIF > 5
m.rare <- glmmTMB(data = hills2, q0.obs ~ treatment*sampling_month + 
                    salinity  + conc_oxi + (1|Plot_trat) , family = "poisson")

DHARMa::simulateResiduals(m.rare, plot = T)
summary(m.rare)
car::Anova(m.rare)
performance::r2(m.rare)
performance::check_collinearity(m.rare)


emmeans(m.rare, ~treatment , type = "response")
pairs(emmeans(m.rare, ~treatment, type = "response"))

emmeans(m.rare, ~sampling_month, type = "response")
pairs(emmeans(m.rare, ~sampling_month, type = "response"))

pairs(emmeans(m.rare, ~treatment|sampling_month, type = "response"))

df.rare <- as.data.frame(emmeans(m.rare, ~treatment|sampling_month, type = "response")) %>% 
  mutate(divindex = rep("q0.obs",8),
         emmean = rate) %>% 
  select(-rate)


# #MODEL FOR COMMON SPECIES (q1)
# hist(hills2$q1.obs)
# car::leveneTest(hills2$q1.obs, hills2$treatment)
# bartlett.test(hills2$q1.obs ~ hills2$treatment)
# summary(hills2$q1.obs)
# m.common <- glmmTMB(data = hills2, q1.obs ~ treatment*sampling_month + ph + 
#                       salinity , family = "gaussian")
# DHARMa::simulateResiduals(m.common, plot = T)
# car::Anova(m.common)
# performance::r2(m.common)
# performance::check_collinearity(m.common)
# emmeans(m.common, ~treatment, type = "response")
# pairs(emmeans(m.common, ~treatment, type = "response"))
# 
# df.common <- as.data.frame(emmeans(m.common, ~treatment, type = "response")) %>% 
#   mutate(divindex = rep("q1.obs", 2))

#FULL MODEL FOR HILL EVENNESS (q2)

m.dom.full <- glmmTMB(data = hills2, q2.obs ~ treatment*sampling_month + ph +  
                   salinity + temperature + conc_oxi + (1|Plot_trat), family = "gaussian")
simulateResiduals(m.dom.full, plot = T)
hist(resid(m.dom.full))
#car::leveneTest(hills2$q2.obs, hills2$treatment)
#bartlett.test(hills2$q2.obs ~ hills2$treatment)
car::Anova(m.dom.full)
performance::r2(m.dom.full)
performance::check_collinearity(m.dom.full)



####REDUCED MODEL FOR HILL EVENNESS (q2)### I removed ph as VIF >5
m.dom <- glmmTMB(data = hills2, q2.obs ~ treatment*sampling_month  +  
                   salinity +  conc_oxi + (1|Plot_trat), family = "gaussian")
simulateResiduals(m.dom, plot = T)
hist(resid(m.dom))
#car::leveneTest(hills2$q2.obs, hills2$treatment)
#bartlett.test(hills2$q2.obs ~ hills2$treatment)
car::Anova(m.dom)
performance::r2(m.dom)
performance::check_collinearity(m.dom)
emmeans(m.dom, ~treatment, type = "response")
pairs(emmeans(m.dom, ~treatment, type = "response"))

emmeans(m.dom, ~sampling_month, type = "response")
pairs(emmeans(m.dom, ~sampling_month, type = "response"))


df.dom <- as.data.frame(emmeans(m.dom, ~treatment|sampling_month, type = "response")) %>% 
  mutate(divindex = rep("q2.obs",8))



hills.long <- hills2 %>% 
  pivot_longer(cols = c(q0.obs, q2.obs), names_to = "divindex", values_to = "value") %>% 
  select(plot, treatment, sampling_month, temperature, sampling_order,divindex,value) %>% 
  filter(divindex %in% c("q0.obs", "q2.obs")) %>% 
  mutate(label = if_else(divindex == "q0.obs", "Taxa richness", "Hill evenness"))

hills.long$label <- factor(hills.long$label, levels = c("Taxa richness", "Hill evenness"))

estindex <- bind_rows(df.rare, df.dom) %>% 
  filter(divindex %in% c("q0.obs", "q2.obs")) %>% 
  mutate(label = if_else(divindex == "q0.obs", "Taxa richness", "Hill evenness"))
estindex$label <- factor(estindex$label, levels = c("Taxa richness", "Hill evenness"))

pdivindex <- ggplot(hills.long, aes(sampling_month, value, group = treatment, colour = treatment, fill = treatment)) +
  geom_point(position = position_jitterdodge (0.45, jitter.width = 0.1, jitter.height = 0.1), alpha = 0.2,  
             shape = 21,colour = "black", size = 5) +
  #stat_summary(fun.data = mean_se, position = position_dodge(0.65),geom = "errorbar", width = 0.1) +
  #stat_summary(fun = mean, position = position_dodge(0.65), size = 1.8) +
  scale_colour_manual(name = "Management", values = c("#F2A104", "#00743F")) +
  scale_fill_manual(values = c("#F2A104", "#00743F"), guide = "none") +
  #scale_x_discrete(labels = c("Taxa richness (q0)", "Hill evenness (q2)")) +
  theme_bw() +
  xlab("") +
  ylab ("Index value") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  facet_wrap(~ label, scales = "free_y")

pdivindex2 <- pdivindex + 
  geom_errorbar(
    mapping = aes(
      ymax =  emmean + SE,
      ymin = emmean - SE,
      x = sampling_month,
      group = treatment
    ),
    position = position_dodge(0.45),
    colour = "black",
    width = 0.1,
    data = estindex,
    inherit.aes = FALSE
  ) +
  geom_point(
    mapping = aes(
      y =  emmean,
      x = sampling_month,
      fill = treatment),
    position = position_dodge(0.45),
    shape = 21,
    colour = "black",
    alpha = 1, 
    size = 8,
    data = estindex,
    inherit.aes = FALSE
  ) 

pdf("figures/Fig2_Hill_profiles.pdf", width = 8, height = 9)
pdivindex2
dev.off()
