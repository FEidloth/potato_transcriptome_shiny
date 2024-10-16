#
#### LOAD LIBRARY ####
library(tidyverse)
library(rstatix)
library(shiny)
library(shinythemes)
library(ggpubr)
library(pheatmap)
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


#### CREATE APP ####
# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                titlePanel("Potato Expression Atlas"), # title of App, optional
                
                navbarPage(
                  "Potato transcriptome", # Navbar side header
                  
                  tags$head(
                    tags$style(HTML("
                    #genename {
                    color: #a5a5a5;
                    background-color: #f0f0f0;
                    border: 1px solid #595959;
                    }"))
                  ), #change colour of textInput field genename
          
                  tabPanel("Expression Data",
                           sidebarPanel(
                             tags$h3("Gene ID:"),
                             textInput("genename", "Soltu.DM.:", "Soltu.DM.01G000010")), # sidebarPanel
                           
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
                  
                  tabPanel("Heatmap", 
                           sidebarPanel(
                             tags$h3("Enter Gene IDs (comma seperated):"),
                             textInput("gene_ids", "Gene_IDs:", ""),
                             actionButton("heatmap_button", "Generate Heatmap")
                           ),
                           mainPanel(
                             h1("Heatmap"),
                             plotOutput("heatmap_plot")
                           )),
                  
                  tabPanel("Background", "Background information will be displayed here")
                  # about paper
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$output_genename <- renderText({
    paste("Expression of", input$genename, 
             sep = " ")
  })
  
  stats <- reactive({
    gene_combined %>%
    filter(gene_id == input$genename) %>%
    group_by(Tissue) %>%
    t_test(log ~ Time)})
  
  stat_test <- reactive({stats %>% add_y_position()})
 
  
  gene_plot <- reactive({
    gene_combined %>%        
      filter(gene_id == input$genename) %>%        
      ggplot(aes(x = Time, y = log, color = Tissue, fill = Tissue)) +       
      facet_grid(~Tissue) +       
      geom_point() +
      scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, NA)) +
      stat_summary(
        fun = mean, geom = "point",
        shape = 95, size = 10, alpha = 0.8
      ) +
      stat_compare_means(aes(group = Time),ref.group = ".all.", label = "p.signif", 
                         size = 6, colour = "#595959", hide.ns = T) +
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
      paste0("gene_expression_plot_", input$genename, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = gene_plot(), device = "pdf", width = 8, height = 4)
    })
  
  output$heatmap_plot <- renderPlot({
      req(input$heatmap_button) #wait for click on button
    
    # Put input genes in list
    gene_list <- strsplit(input$gene_ids, ",")[[1]] # split input at ","
    gene_list <- trimws(gene_list) # remove spaces
    
    # filter input genes from gene_combined
    gene_subset <- gene_combined %>% 
      filter(gene_id %in% gene_list) %>% 
      pivot_wider(names_from = Time,
                  values_from = log) %>% 
      column_to_rownames("gene_id")
    
    # Generate Heatmap
    pheatmap::pheatmap(
      as.matrix(gene_subset),
      color = viridis(100),
      cluster_rows = T,
      cluster_cols = T,
      show_rownames = T,
      show_colnames = T
    )
    })
  
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
