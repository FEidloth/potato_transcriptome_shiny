#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
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

gene_combined <- rbind(log_leaf, log_tuber)
gene_combined <- gene_combined %>% 
  mutate(Time = fct_relevel(Time, "ZT0", "ZT4", "ZT8",
                            "ZT12", "ZT16", "ZT20", 
                            "ZT24"))



# Define UI for application that draws a histogram
# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                titlePanel("Potato Expression Atlas"), # title of App, optional
                
                navbarPage(
                  "Potato transcriptome", # Navbar side header
          
                  tabPanel("Expression Data",
                           sidebarPanel(
                             tags$h3("Gene ID:"),
                             textInput("genename", "Soltu.DM.:", "")), # sidebarPanel
                           
                           mainPanel(
                             # Header
                             h1("Expression Data"),
                             h3(textOutput("output_genename")),
                             
                             # Plot
                             plotOutput("gene_plot")
                             
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Background", "Background information will be displayed here")
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$output_genename <- renderText({
    paste("Expression of", input$genename, 
             sep = " ")
  })
  
  output$gene_plot <- renderPlot({
    gene_combined %>% 
      filter(gene_id == input$genename) %>% 
      ggplot(aes(x = Time, y = log, color = Tissue, fill = Tissue)) +
      facet_grid(~Tissue)+
      geom_point() +
      labs(title = paste("Log count of", input$genename, sep = " "))
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
