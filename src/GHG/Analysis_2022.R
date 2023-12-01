# Load source script:
source("src/Main_w_corrections_2022.R")

#### 1. Create a dataframe with average rates per treatment and Sampling_date.

# Averaged rates from original and corrected models:

Avg_rates_w_corr_2022 <- Emission_rates_w_corrections_2022 %>% 
                group_by(Sampling_date, Treat) %>% 
                summarize(avg_CH4_flux_mgm2h = mean(CH4_flux_mgm2h), avg_CH4_flux_corrected_mgm2h = mean(CH4_flux_corrected))

save(Avg_rates_w_corr_2022, file = "outputs/2022/Rates_corrected/Avg_rates_w_corr_2022.RData")

# Plotting rates across dates:

pdf('outputs/2022/Rates_corrected/CH4_2022_rates_time_original_vs_corr.pdf', width = 12)

Rates_vs_time_CH4 <- ggplot(data = Avg_rates_w_corr_2022, aes(color = Treat, x = Sampling_date, y = CH4_flux_mgm2h, group = Treat)) +
                geom_line(aes(y = avg_CH4_flux_mgm2h, linetype = "Average CH4 Flux")) +
                geom_line(aes(y = avg_CH4_flux_corrected_mgm2h, linetype = "Average Corrected CH4 Flux")) +
                scale_linetype_manual(values = c("dashed", "solid")) +  # Manually set line types
                theme_minimal() +
                xlab("Time") +
                ylab("CH4 flux (mgm2h)") +
                ggtitle("CH4 Emission rates") +
                theme(plot.title = element_text(hjust = 0.5)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(Rates_vs_time_CH4)

dev.off()
