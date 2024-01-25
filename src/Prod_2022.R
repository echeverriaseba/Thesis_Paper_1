
####################################################### Thesis Paper 1 - Production #####################################################

library(dplyr)
library(ggplot2)

Prod_2022 <- read.csv("data/PROD/Prod_2022.csv", fileEncoding="latin1", na.strings=c("","NA")) %>%  # Import Prod_2022 data
              subset(Experiment == "Cerestres") 

Prod_2022_means <- Prod_2022 %>% 
                    group_by(Treat) %>% 
                    summarise(Yield_mean = mean(Yield_kg_ha), SE = sd(Yield_kg_ha) / sqrt(n())) %>% 
                    ungroup()

Prod_2022$Treat <- factor(Prod_2022$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable

Prod_2022_plot <- ggplot(Prod_2022, aes(Treat, Yield_kg_ha, group = Treat, colour = Treat, fill = Treat)) +
                          geom_point(position = position_jitterdodge (0.80, jitter.width = 0.4, jitter.height = 0), alpha = 0.8, shape = 21, colour = "black",size = 10) +
                          geom_bar(data = Prod_2022_means, aes(x = Treat, y = Yield_mean, fill = Treat), alpha = 0.1, stat = "identity", width = 0.5) +
                          geom_errorbar(data = Prod_2022_means, aes(y = Yield_mean , ymin = Yield_mean - SE, ymax = Yield_mean + SE, color = Treat), width = 0.4, position = position_dodge(width = 0.9), size = 1) +
                          scale_colour_manual(name = "Treatment", values = c(CON = "#002B5B", MSD = "#03C988", AWD = "#FF5D5D")) +
                          scale_fill_manual(values = c(CON = "#002B5B", MSD = "#03C988", AWD = "#FF5D5D"), guide = "none") +
                          theme_bw() +
                          # scale_y_sqrt() +
                          # scale_y_log() +
                          # scale_y_continuous(trans = "log10") +
                          # ylim(6000, 10000) +
                          ylab(expression("Yield (kg ha"^-1*")")) +
                          ggtitle("") +
                          theme(plot.title = element_text(size=20, hjust=0.5)) +
                          theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                axis.title.y = element_text(size = 20, margin = margin(r = 12)), axis.title.x = element_blank(), legend.position = "none", 
                                axis.text.y = element_text(size = 20, margin = margin(r = 10), angle = 90), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
                          scale_y_continuous(expand = c(0, 0), limits = c(0, 9100), breaks = seq(1000, 9000, by = 1000))
                          # scale_y_continuous(expand = c(0, 0), limits = c(6000, 9100), breaks = seq(7000, 9000, by = 1000))

print(Prod_2022_plot)

ggsave("outputs/Plots/PROD/Prod_2022_plot.pdf", plot = Prod_2022_plot, width = 8, height = 10)
