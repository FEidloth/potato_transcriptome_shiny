rm(list = ls())
#
#### LOAD LIBRARY ####
library(tidyverse)
library(rstatix)
library(shiny)
library(shinythemes)
library(viridis)
library(viridisLite)
library(ggpubr)
#### LOAD DATA ####
Meta_Data <- read.csv("../META_DATA/Hoopes_et_al.,2022/Meta_Data.csv")

leaf <- read.delim("../DATA/PLD3-6-e425-s001_leaf_LD.txt")
tuber <- read.delim("../DATA/PLD3-6-e425-s006_tuber_LD.txt")

#### PLOT THEME ####
my_theme <- 
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "#595959", fill = NA, linewidth = .6),
        panel.grid.major.y = element_line(colour = "#a5a5a5", linewidth = .1),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#a5a5a5", linewidth = .1),
        panel.grid.minor.x = element_blank(),
        axis.ticks = element_line(linewidth = .3, color = "#595959"), 
        axis.ticks.length = unit(.2, "cm"),
        axis.title.x = element_text(color = "#595959", vjust = 0, size = 10, face = 1),
        axis.title.y = element_text(color = "#595959", vjust = 2, size = 10, face = 1),
        axis.text = element_text(color="#595959", size=10, face=1),
        plot.title = element_text(margin = margin(10, 0, 10, 0),
                                  size = 14, hjust = 0.5,
                                  colour = "#595959"),
        strip.text.x = element_text(color = "#595959", size=12), 
        strip.background = element_rect(colour="transparent", fill="transparent"),
        legend.position = "none"
  )
#
#### DATA PREPROCESSING ####
log_leaf <- leaf %>% 
  select(-c(Module, BH_P_value,meta2d_AMP, meta2d_rAMP, Amplitude.Change.Coefficient,
            Oscillation.Type, Ave_Period, Ave_Phase)) %>% 
  rename_with(~ paste0(., "_leaf"), starts_with("T")) %>% 
  pivot_longer(cols = T0_1_leaf:T24_3_leaf,
               names_to = "Sample",
               values_to = "log") %>% 
  left_join(., Meta_Data, by = "Sample")

##
log_tuber <- tuber %>% 
  select(-c(Module, BH_P_value,meta2d_AMP, meta2d_rAMP, Amplitude.Change.Coefficient,
            Oscillation.Type, Ave_Period, Ave_Phase)) %>%
  rename_with(~ paste0(., "_tuber"), starts_with("T")) %>% 
  pivot_longer(cols = T0_1_tuber:T24_3_tuber,
               names_to = "Sample",
               values_to = "log") %>% 
  left_join(., Meta_Data, by = "Sample")

#### FILTER ####
gene_leaf <- log_leaf %>% 
  filter(gene_id == "Soltu.DM.01G000490")
gene_tuber <- log_tuber %>% 
  filter(gene_id == "Soltu.DM.01G000490")

gene_combined <- rbind(gene_leaf, gene_tuber)
gene_combined <- gene_combined %>% 
  mutate(Time = fct_relevel(Time, "ZT0", "ZT4", "ZT8",
                            "ZT12", "ZT16", "ZT20", 
                            "ZT24"))
#### STATS ####

stat_test <- compare_means(log ~ Time, gene_combined, group.by = "Tissue", 
                           ref.group = ".all.", method = "wilcox.test")

Max_log <- gene_combined %>% 
  group_by(Tissue, Time) %>% 
  summarise(max_log = max(log)) %>% 
  mutate(y.position = max_log + 0.5)

y_position <- as.vector(Max_log$y.position)

#### PLOT ####

gene_combined %>%
  ggplot(aes(x = Time, y = log, color = Tissue, fill = Tissue)) +
  geom_point() +
  scale_y_continuous(expand = expansion(mult = c(0, .2)), limits = c(0, NA)) +
  stat_summary(
    fun = mean, geom = "point", 
    shape = 95, size = 10, alpha = 0.8
  ) +
  stat_compare_means(aes(group = Time),ref.group = ".all.", label = "p.signif", 
                     #label.y = y_position,
                     size = 6, colour = "#595959",
                     hide.ns = T) +
  #stat_pvalue_manual(
   # stat_test, label = "p.signif", hide.ns = TRUE
  #) +
  facet_grid(~Tissue) +
  scale_color_viridis_d(option = "viridis", begin = .7, end = .3) +
  scale_fill_viridis_d(option = "viridis", begin = .7, end = .3) +
  my_theme










