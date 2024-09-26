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
#
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

gene_combined <- rbind(log_leaf, log_tuber)
gene_combined <- gene_combined %>% 
  mutate(Time = fct_relevel(Time, "ZT0", "ZT4", "ZT8",
                            "ZT12", "ZT16", "ZT20", 
                            "ZT24"))


ZT_comparisons <- list(c(c("ZT0", "ZT4"), c("ZT0", "ZT8"), c("ZT0", "ZT12"),
                         c("ZT0", "ZT16"), c("ZT0", "ZT20"), c("ZT0", "ZT24"),
                         c("ZT4", "ZT8"), c("ZT4", "ZT12"), c("ZT4", "ZT16"),
                         c("ZT4", "ZT20"), c("ZT4", "ZT24"), c("ZT8", "ZT12"),
                         c("ZT8", "ZT16"), c("ZT8", "ZT20"), c("ZT8", "ZT24"),
                         c("ZT12", "ZT16"), c("ZT12", "ZT20"), c("ZT12", "ZT24"),
                         c("ZT16", "ZT20"), c("ZT16", "ZT24"), c("ZT20", "ZT24")))

#### CREATE APP ####
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
                             plotOutput("gene_plot"),
                             
                             # Download Button
                             downloadButton("download_plot", "Download Plot")
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
  
  gene_plot <- reactive({
    gene_combined %>%        
      filter(gene_id == input$genename) %>%        
      ggplot(aes(x = Time, y = log, color = Tissue, fill = Tissue)) +       
      facet_grid(~Tissue) +       
      geom_point() +
      scale_y_continuous(expand = expansion(mult = c(0, .2)), limits = c(0, NA)) +
      stat_summary(
        fun = mean, geom = "point", 
        shape = 95, size = 10, alpha = 0.8
      ) +
      stat_compare_means(
        method = "t.test", 
        aes(group = Time), 
        comparisons = ZT_comparisons,
        label = "p.signif", 
        hide.ns = TRUE,
        bracket.size = 0.2,
        tip.length = 0.01
      ) +
      labs(title = paste("Log count of", input$genename, sep = " ")) +
      scale_color_viridis_d(option = "viridis", begin = .7, end = .3) +
      scale_fill_viridis_d(option = "viridis", begin = .7, end = .3) +
      my_theme  
  })
  
  output$gene_plot <- renderPlot({
    gene_plot()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("gene_expression_plot_", input$genename, ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = gene_plot(), device = "pdf", width = 8, height = 4)
    }
  )
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
