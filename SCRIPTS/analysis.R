rm(list = ls())
#
#### LOAD LIBRARY ####
library(tidyverse)
library(rstatix)
library(shiny)
library(shinythemes)
#### LOAD DATA ####
Meta_Data <- read.csv("../META_DATA/Hoopes_et_al.,2022/Meta_Data.csv")

leaf <- read.delim("../DATA/PLD3-6-e425-s001_leaf_LD.txt")
tuber <- read.delim("../DATA/PLD3-6-e425-s006_tuber_LD.txt")

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
#### STATS####

stat_test <- gene_combined %>%
  group_by(Time, Tissue) %>%
  t_test(log ~ Time)
stat_test <- stat_test %>% add_y_position()


#### PLOT ####

gene_combined %>% 
  ggplot(aes(x = Time, y = log, color = Tissue, fill = Tissue)) +
  facet_grid(~Tissue)+
  geom_point()
















