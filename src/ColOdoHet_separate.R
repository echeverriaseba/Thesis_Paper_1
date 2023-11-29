
# To plot Annexes with plots for Col, Het and Odo separate (see: ColOdoHet_q0q1_2022_all), then siteID has to be changed to (Sampling&Rep&Treat&Order), and the Hills_ColOdoHet from script Macroinv_2022 merge has to be:

# Hills_ColOdoHet <- ColOdoHet %>% 
#   group_by(Plot, Sampling, Treat, siteID, ) %>% 
#   summarise(siteID = first(siteID),
#             Plot = first(Plot),
#             Sampling = first(Sampling),
#             Treat = first(Treat)) %>% 
#   left_join(ColOdoHet_inext_params, by = "siteID")

## Plotting iNEXT parameters:

Col_inext_params <- subset(Hills_ColOdoHet, Order_SubOrder == "Coleoptera")
Odo_inext_params <- subset(Hills_ColOdoHet, Order_SubOrder == "Odonata")
Het_inext_params <- subset(Hills_ColOdoHet, Order_SubOrder == "Heteroptera")

Col_inext_params$Treat <- factor(Col_inext_params$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable
Odo_inext_params$Treat <- factor(Odo_inext_params$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable
Het_inext_params$Treat <- factor(Het_inext_params$Treat, levels = c("CON", "MSD", "AWD")) # Reorder the Treat variable

# 1. Coleoptera:

summary_q0_Col <- Col_inext_params %>%
                  group_by(Treat) %>%
                  summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

summary_q1_Col <- Col_inext_params %>%
                  group_by(Treat) %>%
                  summarise(mean_q1.obs = mean(q1.obs),se_q1.obs = sd(q1.obs) / sqrt(n()))

Plot_SpRich_Col <- ggplot(Col_inext_params, 
                          aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                          geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
                          scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
                          scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
                          # scale_y_sqrt() +
                          theme_bw() +
                          # xlab("Treatment") +
                          ylab(expression("Species richness (q"[0]*")")) +
                          ggtitle("Coleoptera") +
                          theme(plot.title = element_text(size=20, hjust=0.5)) +
                          theme(axis.title.x = element_blank(), axis.text = element_text(size = 14), strip.text = element_text(size = 14), 
                                axis.title.y = element_text(size = 20), axis.text.x = element_blank(), legend.position = "none")+
                          geom_point(data = summary_q0_Col, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 5, color = "black") +
                          geom_errorbar(data = summary_q0_Col, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.2) +
                          scale_y_continuous(limits = c(0.8, 7.2), breaks = seq(1, 7, by = 1))

print(Plot_SpRich_Col)

Plot_Shannon_Col <- ggplot(Col_inext_params, 
                           aes(Treat, q1.obs, group = Treat, colour = Treat, fill = Treat)) +
                            geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
                            scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
                            scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
                            scale_y_sqrt() +
                            theme_bw() +
                            # xlab("Treatment") +
                            ylab(expression("Shannon index (q"[1]*")")) +
                            theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
                                  axis.title.x = element_blank(), legend.position = "none")+
                            geom_point(data = summary_q1_Col, aes(x = Treat, y = mean_q1.obs), shape = 19, size = 5, color = "black") +
                            geom_errorbar(data = summary_q1_Col, aes(x = Treat, y = mean_q1.obs, ymin = mean_q1.obs - se_q1.obs, ymax = mean_q1.obs + se_q1.obs), width = 0.2) +
                            scale_y_continuous(limits = c(0.8, 6.2), breaks = seq(1, 6, by = 1))


print(Plot_Shannon_Col)

# 2. Odonata:

summary_q0_Odo <- Odo_inext_params %>%
  group_by(Treat) %>%
  summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

summary_q1_Odo <- Odo_inext_params %>%
  group_by(Treat) %>%
  summarise(mean_q1.obs = mean(q1.obs),se_q1.obs = sd(q1.obs) / sqrt(n()))

Plot_SpRich_Odo <- ggplot(Odo_inext_params, 
                          aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
                          geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
                          scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
                          scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
                          #scale_y_sqrt() +
                          theme_bw() +
                          xlab("Treatment") +
                          # ylab ("Odonata (q0)") +
                          ggtitle("Odonata") +
                          theme(plot.title = element_text(size=20, hjust=0.5)) +
                          theme(axis.title.x = element_blank(), axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
                                strip.text = element_text(size = 14), axis.text.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), 
                                axis.text.y = element_blank())+
                          geom_point(data = summary_q0_Odo, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 5, color = "black") +
                          geom_errorbar(data = summary_q0_Odo, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.2) +
                          scale_y_continuous(limits = c(0.8, 7.2), breaks = seq(1, 7, by = 1))

print(Plot_SpRich_Odo)

Plot_Shannon_Odo <- ggplot(Odo_inext_params, 
                           aes(Treat, q1.obs, group = Treat, colour = Treat, fill = Treat)) +
  geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
  scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
  scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
  scale_y_sqrt() +
  theme_bw() +
  # xlab("Treatment") +
  # ylab ("Odonata (q1)") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
        axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())+
  geom_point(data = summary_q1_Odo, aes(x = Treat, y = mean_q1.obs), shape = 19, size = 5, color = "black") +
  geom_errorbar(data = summary_q1_Odo, aes(x = Treat, y = mean_q1.obs, ymin = mean_q1.obs - se_q1.obs, ymax = mean_q1.obs + se_q1.obs), width = 0.2) +
  scale_y_continuous(limits = c(0.8, 6.2), breaks = seq(1, 6, by = 1))

print(Plot_Shannon_Odo)

# Heteroptera:

summary_q0_Het <- Het_inext_params %>%
  group_by(Treat) %>%
  summarise(mean_q0.obs = mean(q0.obs),se_q0.obs = sd(q0.obs) / sqrt(n()))

summary_q1_Het <- Het_inext_params %>%
  group_by(Treat) %>%
  summarise(mean_q1.obs = mean(q1.obs),se_q1.obs = sd(q1.obs) / sqrt(n()))

Plot_SpRich_Het <- ggplot(Het_inext_params, 
                          aes(Treat, q0.obs, group = Treat, colour = Treat, fill = Treat)) +
  geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
  scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
  scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
  #scale_y_sqrt() +
  theme_bw() +
  xlab("Treatment") +
  # ylab ("Heteroptera (q0)") +
  ggtitle("Heteroptera") +
  theme(plot.title = element_text(size=20, hjust=0.5)) +
  theme(axis.title.x = element_blank(), axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        strip.text = element_text(size = 14),axis.text.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_blank())+
  geom_point(data = summary_q0_Het, aes(x = Treat, y = mean_q0.obs), shape = 19, size = 5, color = "black") +
  geom_errorbar(data = summary_q0_Het, aes(x = Treat, y = mean_q0.obs, ymin = mean_q0.obs - se_q0.obs, ymax = mean_q0.obs + se_q0.obs), width = 0.2) +
  scale_y_continuous(limits = c(0.8, 7.2), breaks = seq(1, 7, by = 1))

print(Plot_SpRich_Het)

Plot_Shannon_Het <- ggplot(Het_inext_params, 
                           aes(Treat, q1.obs, group = Treat, colour = Treat, fill = Treat)) +
  geom_point(position = position_jitterdodge (0.45, jitter.width = 0.07, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 5)+
  scale_colour_manual(name = "Treatment", values = c("#5FA7FF", "#90FF37", "#FF5D5D")) +
  scale_fill_manual(values = c("#5FA7FF", "#90FF37", "#FF5D5D"), guide = "none") +
  scale_y_sqrt() +
  theme_bw() +
  # xlab("Treatment") +
  # ylab ("Heteroptera (q1)") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
        axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())+
  geom_point(data = summary_q1_Het, aes(x = Treat, y = mean_q1.obs), shape = 19, size = 5, color = "black") +
  geom_errorbar(data = summary_q1_Het, aes(x = Treat, y = mean_q1.obs, ymin = mean_q1.obs - se_q1.obs, ymax = mean_q1.obs + se_q1.obs), width = 0.2) +
  scale_y_continuous(limits = c(0.8, 6.2), breaks = seq(1, 6, by = 1))

print(Plot_Shannon_Het)

# Arrange plots:

ColOdoHet_arrange_separate <- grid.arrange(arrangeGrob(Plot_SpRich_Col, Plot_SpRich_Odo, Plot_SpRich_Het, Plot_Shannon_Col, Plot_Shannon_Odo, Plot_Shannon_Het,nrow = 2, ncol = 3))

ggsave("outputs/Plots/BIO/ColOdoHet_q0q1_2022_alarvae.pdf", plot = ColOdoHet_arrange ,width = 20, height = 10)
