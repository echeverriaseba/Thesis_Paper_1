library(dplyr)

# Load emissions 2022 with corrections (calculated previously in script "GHG_rates_2022_w_corrections") 
load("outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData")

# 1. Data frame removing all R2 < 0 ####

GHG_2022_noRule <- Emission_rates_w_corrections_2022 %>% 
                    mutate(CH4_noRule = if_else(R2_CH4 < 0.7, NA_real_, CH4_flux_mgm2h)) %>% 
                    filter(!is.na(CH4_noRule))
                    

GHG_2022_noRule$Sampling_date <- as.Date(GHG_2022_noRule$Sampling_date)

# 2. CH4 flux Plot ####

CH4_noRule_plot <- ggplot(GHG_2022_noRule, aes(x = Sampling_date, y = CH4_noRule, color = Treat, group = Plot)) +
                      geom_line(alpha = 0.5, linetype = "dotted") +  # Adjust transparency by setting alpha
                      scale_colour_manual(name = "Irrigation strategies", values = c("#002B5B", "#03C988", "#FF5D5D"), breaks=c('CON', 'MSD', 'AWD')) +
                      theme_bw() +
                      labs(y = expression(paste(C-CH[4], " flux (mg ", m^-2, " ", h^-1, ")"))) +
                      geom_hline(yintercept = 0, color = "grey") +
                      guides(linetype = guide_legend(override.aes = list(color = c("black", "black")))) +
                      theme(
                        axis.title.y = element_text(color = "black"), legend.margin=margin(0,0,0,0),
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

Avg_GHG_2022_noRule <- GHG_2022_noRule %>% 
                      group_by(Sampling_date, Treat) %>% 
                      summarise(CH4_noRule = mean(CH4_noRule))

CH4_noRule_plot <- CH4_noRule_plot +
                      geom_line(data = Avg_GHG_2022_noRule, aes(x = Sampling_date, y = CH4_noRule, group = Treat))

print(CH4_noRule_plot) 

# 3. Calculating global cumulative C-CH4 emissions (kg C-CH4 ha-1) ####

GHG_2022_noRule$Row_Nr <- 1:nrow(GHG_2022_noRule) # Adding a row number column

Accum_CH4_noRule <- GHG_2022_noRule %>% 
                      select("Row_Nr", "Sampling_date", "Treat", "Plot", "Rep", "CH4_noRule") %>% 
                      arrange(Plot, Sampling_date) %>% 
                      group_by(Plot) %>% 
                      mutate(Days_passed = as.integer(Sampling_date - lag(Sampling_date))) %>% 
                      mutate(Hours_passed = Days_passed * 24) %>% 
                      mutate(CH4_mgm2 = CH4_noRule * Hours_passed) %>% 
                      mutate(CH4_kgha = CH4_mgm2 / 100) # *10.000 from m-2 to ha-1 ; /1.000.000 from mg to kg

Accum_CH4_noRule$CH4_kgha <- ifelse(is.na(Accum_CH4_noRule$CH4_kgha), 0, Accum_CH4_noRule$CH4_kgha)

Accum_CH4_noRule_tot <- Accum_CH4_noRule %>% 
                  group_by(Treat, Plot) %>%
                  summarise(CH4_kgha_tot = sum(CH4_kgha)) 

Accum_CH4_noRule_tot$Treat <- factor(Accum_CH4_noRule_tot$Treat, levels = c('CON', 'MSD', 'AWD')) # Treat to factor to reorder below's ggplot x axis

Accum_CH4_noRule_avg <- Accum_CH4_noRule_tot %>% 
                      group_by(Treat) %>%
                      summarise(CH4_kgha_avg = mean(CH4_kgha_tot)) %>% 
                      mutate(dif_vs_CON = (CH4_kgha_avg[1] - CH4_kgha_avg)/CH4_kgha_avg[1])