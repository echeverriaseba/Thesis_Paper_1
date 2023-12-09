##### Analysis from Stats meeting with Nestor and Carles (29.11.23) ####

# Example model: q1 - gaussian w/interaction (Sampling^2) - Removing high correlation variables from model 6  #
# Dependent variable: q1
# Variables removed from model 6: Treat*Sampling + Redox_pot + O2_percent 
# "Sampling" variable: I(Sampling^2)
# Family: Gaussian
# Interacting independent variables: Treat*Sampling + Treat*I(Sampling^2)
# Additional independent variables: Conduct_microS_cm + pH_soil + Redox_pot + Water_temp + O2_percent + Salinity
# Random effect: Rep

# glmm.q1.gaus7 <- glmmTMB(data = Hills_Physchem, q1.obs ~ Treat*Sampling + Treat*I(Sampling^2) + Treat*Conduct_microS_cm + Treat*Salinity + 
#                            (1|Rep) , family = "gaussian")
# 
# # Model diagnostics:
# DHARMa::simulateResiduals(glmm.q1.gaus7, plot = T)
# summary(glmm.q1.gaus7)
# car::Anova(glmm.q1.gaus7)
# performance::r2(glmm.q1.gaus7)
# performance::check_collinearity(glmm.q1.gaus7)
# performance::check_singularity(glmm.q1.gaus7)
# visreg(glmm.q1.gaus7, scale="response") # Plotting conditional residuals
# 
## Notes from meeting with Carles Alcaraz:
# plot(Hills_Physchem$Treat, Hills_Physchem$q1.obs)
# plot(Hills_Physchem$Salinity, Hills_Physchem$Conduct_microS_cm)
# plot(Hills_Physchem$Water_temp, Hills_Physchem$Temp_10_cm)
# 
# glm.Sal_Cond <- glm(data = Hills_Physchem, q1.obs ~ Conduct_microS_cm + Salinity , family = "gaussian")
# dffits(glm.Sal_Cond)
# plot(row.names(na.omit(Hills_Physchem)), dffits(glm.Sal_Cond))
# identify(row.names(na.omit(Hills_Physchem)), dffits(glm.Sal_Cond))
# 
# glm.Sal_Cond_nooutl <- glm(data = Hills_Physchem_nooutliers, q1.obs ~ Conduct_microS_cm + Salinity , family = "gaussian")
# dffits(glm.Sal_Cond_nooutl)
# plot(row.names(na.omit(Hills_Physchem_nooutliers)), dffits(glm.Sal_Cond_nooutl))
# identify(row.names(na.omit(Hills_Physchem_nooutliers)), dffits(glm.Sal_Cond_nooutl))