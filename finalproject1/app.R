#import libraries
library(shiny)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DT)
library(bslib)
library(RColorBrewer)
library(gplots)
library(shinythemes)
library(colourpicker)
options(shiny.maxRequestSize=30*1024^2)

ui <- navbarPage(
  title = 'BF591 Final Project - Vrinda Jethalia', 
  theme = shinytheme("flatly"), 
  
  #------------------------------------First Tab Layout--------------------------------------#
  
  tabPanel("About the Application",
           h3(tags$b("This is a Shiny application that can be used for some quick Bioinformatics Analysis"), style = "text-align: center;"),
           h4(tags$b("Want to analyze the Huntington's Disease data? Try exploring the next 4 tabs using the various files in the data folder!"), style = "text-align: center;"),
           p("Developed by Vrinda Jethalia as part of the final project for BF591 at Boston University.", style = "text-align: center;"),
           img(src = "dna.jpg", width = 500, height = 200, style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);")
           ),
  
  #-------------------------------------Second Tab Layout-------------------------------------#
  
  tabPanel("Sample Exploration",
           sidebarLayout(
             sidebarPanel(
               h4(tags$b('Please load the Sample Information Matrix')),
               fileInput('file_samples', 'Choose file to upload',placeholder = "metadata.csv")),
             mainPanel (
               tabsetPanel(
                 tabPanel("Summary", tableOutput(outputId = "sample_summary")),
                 tabPanel("Data", dataTableOutput(outputId = "sample_data")),
                 tabPanel("Plot", plotOutput(outputId = "sample_den_aod"), plotOutput(outputId = "sample_den_rin"), 
                          plotOutput(outputId = "sample_den_pmi"), plotOutput(outputId = "sample_den_aoo"))
               )
             ))),
  
  #-------------------------------------Third Tab Layout-------------------------------------#
  
  tabPanel("Counts Exploration",
           sidebarLayout(
             sidebarPanel(
               h4(tags$b('Please load the Normalized Counts Information Matrix')),
               fileInput('file_counts', 'Choose file to upload', placeholder = "normalized_counts.csv"),
               
               sliderInput('variance_slider', 'Select the minimum percentile of variance genes', 
                           value = 50,
                           min = 0,
                           max = 100, 
                           step = 1),
               sliderInput('nonzero_slider', 'Select the minimum number of non-zero samples',
                           value = 22,
                           min = 0, 
                           max = 69),
               div(style = "text-align: center;",
               submitButton("Submit changes", icon = icon("thumbs-up", style = "color:white;"))
               )
               ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Summary", tableOutput(outputId = "counts_summary")),
                 tabPanel("Diagnostic Scatter Plots", plotOutput(outputId = "scatter_median_var"),
                          plotOutput(outputId = "scatter_median_nonzero")),
                 tabPanel("Heatmap", plotOutput(outputId = "counts_heatmap")),
                 tabPanel("PCA", selectInput(inputId = "component_1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                          selectInput(inputId = "component_2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                          plotOutput(outputId = "counts_pca"))
               )
             ))),
  
  #-------------------------------------Fourth Tab Layout-------------------------------------#
  
  tabPanel("Differential Expression",
           sidebarLayout(
             sidebarPanel(
               h4(tags$b("Please load the Differential Expression Results")),
               fileInput('deg_results', 'Choose file to upload', placeholder = 'deseq2_results.csv'),
               fluidRow(
                 column(
                   width = 12,
                   p("A volcano plot can be generated with ",
                     tags$b("log2 fold-change"), " on the x-axis and ",
                     tags$b("p-adjusted"), " on the y-axis.")
                 )),
               radioButtons("x_axis", "Choose the column for the x-axis",
                            choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                            selected = "log2FoldChange"),
               
               radioButtons("y_axis", "Choose the column for the y-axis",
                            choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                            selected = "padj"),
               
               colourInput("base", 'Base point color', "#22577A"),
               
               colourInput("highlight", 'Highlight point color', "#FFCF56"),
               
               sliderInput(inputId = "slider", min = -50, max = 0,
                           label = "Select the magnitude of p adjusted coloring:", value = -6, step = 1),
               
               div(style = "text-align: center;",
                   submitButton("Plot", icon("check", style = "color:white;"))),
                 ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Results Table", dataTableOutput(outputId = "deg_results")),
                 tabPanel("Plot", plotOutput(outputId = "deg_volcano")),
                 tabPanel("Plot Table", dataTableOutput(outputId = "deg_plot_data"))
                  )
               )
             )),
  
  #-------------------------------------Fifth Tab Layout-------------------------------------#
  
  tabPanel("Gene Set Enrichment Analysis",
           sidebarLayout(
             sidebarPanel(
               h4(tags$b("Please load the fgsea results")),
               fileInput("fgsea_results", "Choose file to upload", placeholder = "fgsea_results.csv"),
               h4(tags$b("1: GSEA NES top pathways filter")),
               sliderInput(inputId = 'gsea_pathways_slider', "Filter by number of pathways",
                           value = 10,
                           min = 1,
                           max = 40
                           ),
               h4(tags$b("2: Data table filter")),
               sliderInput('gsea_data_slider', 'Filter by p-adjusted value',
                           value = -2,
                           min = -20,
                           max = 0),
               radioButtons(inputId ='fgsea_direction', 'NES Direction:',
                           choices = c('all', 'positive', 'negative'),
                           selected = 'positive'),
               div(style = "text-align: center;",
               downloadButton('fgsea_download', 'Download Filtered Results', icon("file-download"), style = "background-color: #2c3e50; border-color: #2c3e50;"),
               ),
               dataTableOutput('gsea_output_download'),
               h4(tags$b("3: Scatter plot filter")),
               sliderInput("gsea_scatter_slider", 'Filter by p-adjusted value',
                           value = -1,
                           min = -20,
                           max = 0,
                           step = 1),
               div(style = "text-align: center;",
               submitButton("Submit changes", icon("check", style = "color:white;"))),
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Top Pathways", plotOutput(outputId = 'gsea_barplot')),
                 tabPanel("Data", dataTableOutput(outputId = 'gsea_data')),
                 tabPanel("NES vs padj Plot", plotOutput(outputId = 'gsea_scatter'))
               )
             )
             )
           ))

#Now it is the server code 
server <- function(input, output) {
  
#--------------------------Tab 2: Functions for Sample Information Exploration---------------------------#
  
  #loading the meta data 
  load_sample_data <- reactive({
    infile <- input$file_samples
    if (!is.null(input$file_samples$datapath)){
      meta_data <- read_csv(input$file_samples$datapath)
      return(meta_data)}
    else{
      return(NULL)
      }
  })
  
  #function to create the summary table
  summary_table <- function(dataf){
    if(!is.null(input$file_samples)){
      summ_df <- data.frame("Column Name" = colnames(dataf), "Type" = sapply(dataf, typeof),
                         "Mean" = sapply(dataf, mean, na.rm = TRUE), "SD" = sapply(dataf, sd, na.rm = TRUE))
      return(summ_df)}
    else{
      return(NULL)
      }
  }
  
  #function to create density plots for aod, pmi, rin and age of onset (aoo)
  aod_den_plot <- function(dataf){
    if(!is.null(input$file_samples))
    {
     density_plot <- ggplot(dataf, aes(Age_of_death))+
       geom_density(color = "black", fill = "coral2")+
       labs(title = 'Density Plot of Age of Death')+
       xlab('Age of Death')+
       ylab('Density')+
       theme_bw()
     return(density_plot)
    }
    else {
      return(NULL)
    }
  }
  
  pmi_den_plot <- function(dataf){
    if(!is.null(input$file_samples))
    {
      density_plot <- ggplot(dataf, aes(PMI))+
        geom_density(color = "black", fill = "indianred")+
        labs(title = 'Density Plot of PMI')+
        xlab('PMI')+
        ylab('Density')+
        theme_bw()
      return(density_plot)
    }
    else {
      return(NULL)
    }
  }
  
  rin_den_plot <- function(dataf){
    if(!is.null(input$file_samples))
    {
      density_plot <- ggplot(dataf, aes(RIN))+
        geom_density(color = "black", fill = "indianred4")+
        labs(title = 'Density Plot of RIN')+
        xlab('RIN')+
        ylab('Density')+
        theme_bw()
      return(density_plot)
    }
    else {
      return(NULL)
    }
  }
  
  aoo_den_plot <- function(dataf){
    if(!is.null(input$file_samples))
    {
      density_plot <- ggplot(dataf, aes(Age_of_onset))+
        geom_density(color = "black", fill = "violetred4")+
        labs(title = 'Density Plot of Age of onset of HD')+
        xlab('Age of Onset')+
        ylab('Density')+
        theme_bw()
      return(density_plot)
    }
    else {
      return(NULL)
    }
  }

  
  #----------------Tab 3: Functions for Counts Matrix Exploration---------------------#
  
  #Load the counts file
  load_counts <- reactive({
    if (!is.null(input$file_counts$datapath)){
      counts_data <- read_csv(input$file_counts$datapath)
      return(counts_data)}
    else{
      return(NULL)
    }
  })
  
  #Summary of the counts file 
  summary_counts <- function(dataf, slider_var, slider_nonzero){
    if(!is.null(input$file_counts))
    {
      no_of_samples <- ncol(dataf) - 1
      tot_no_of_genes <- nrow(dataf)
      dataf <- dataf %>% mutate(variance = apply(dataf[-1], MARGIN = 1, FUN = var))
      percentile_var <- quantile(dataf$variance, probs = slider_var/100)
      dataf <- na_if(dataf, 0)    
      dataf$non_zero <- no_of_samples-rowSums(is.na(dataf))  
      dataf <- filter(dataf, non_zero >= slider_nonzero)  #filter by non-zero samples
      filt_genes <- nrow(dataf)    
      perc_pass_genes1 <- filt_genes/tot_no_of_genes*100
      fail_genes <- tot_no_of_genes-filt_genes
      perc_fail_genes1 <- fail_genes/tot_no_of_genes*100
      #produce the summary tibble
      summary_tibb <- tibble('Parameter' = c('Number of Samples', 'Number of Genes', 'No of genes passing filters', "% of gene passing filters", 'No of genes not passing filters', '% of genes not passing filters'),
                         'Value' = c(no_of_samples, tot_no_of_genes, filt_genes, perc_pass_genes1, fail_genes, perc_fail_genes1))
      return(summary_tibb)
    }
    else {
      return(NULL)
    }
    
  }
  
  #Diagnostic scatter plots
  variance_vs_median <- function(dataf, slider_var){
    plotf <- dataf%>%
      mutate(Median = apply(dataf[-1], MARGIN = 1, FUN = median), 
             Variance = apply(dataf[-1], MARGIN = 1, FUN = var))
    perc_val <- quantile(plotf$Variance, probs = slider_var/100)   #calculate percentile
    plotf <- plotf %>% mutate(threshold = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE")) #sort
    #plot the scatter plot
    cols <- c("FALSE" = "lightpink", "TRUE" = "darkblue")
    scatter <- ggplot(plotf, aes(Median, Variance))+
      geom_point(aes(color=threshold), alpha=0.75)+
      scale_color_manual(values = cols)+
      labs(title = 'Plot of Median Count vs Variance', subtitle = "The filtered out genes are in light pink. X and Y axes are log-scaled")+
      scale_y_log10()+
      scale_x_log10()+
      theme_bw()+
      theme(legend.position = 'bottom')
    return(scatter)
  }
  
  median_vs_nonzeros <- function(dataf, slider_nonzero){
    no_tot_samples <- ncol(dataf)-1  #store original number of samples
    #make a plot tibble
    plotf <- dataf %>%   
      mutate(Median = apply(dataf[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  #calculate median, convert 0 to NA
    plotf$no_zeros <- rowSums(is.na(plotf))  
    plotf <- plotf %>% mutate(threshold = case_when(no_zeros <= slider_nonzero ~ "TRUE", TRUE ~ "FALSE")) #sort genes based on number of zeroes
    #plot scatter plot
    cols <- c("FALSE" = "lightpink", "TRUE" = "darkblue")
    scatter <- ggplot(plotf, aes(Median, no_zeros))+
      geom_point(aes(color=threshold), alpha=0.75)+
      scale_color_manual(values = cols)+
      scale_x_log10()+
      labs(title = 'Plot of Median Count vs Number of Non-Zero genes', subtitle = "Filtered out genes are in light pink. X-axis is log scaled.")+
      theme_bw()+
      ylab('Number of samples with zero count')+
      theme(legend.position = 'bottom')
    return(scatter)
  }
  
  
  #Function for generating a clustered heatmap of counts after filtering
  plot_heatmap <- function(dataf, slider_var, slider_nonzero){
      dataf <- na_if(dataf, 0)
      dataf$no_zeros <- rowSums(is.na(dataf))  #make new col, with counts.
      dataf <- filter(dataf, no_zeros <= slider_nonzero)
      dataf <- log10(dataf[,!colnames(dataf) %in% c("gene", "no_zeros")]) #exclude the gene names column and log scale the values  
      #produce plot
      plotf <- dataf %>% 
        mutate(variance = apply(dataf, MARGIN = 1, FUN = var)) #compute variance to filter the data
      perc_val <- quantile(plotf$variance, probs = slider_var/100, na.rm = TRUE)   #calculate percentile
      plotf <- filter(plotf, variance >= perc_val) #filter the tibble
      hmap <- heatmap.2(as.matrix(plotf[-ncol(plotf)]), scale = "row", col = brewer.pal(9, "Purples"), main="Clustered heatmap of the counts that remain after filtering")
      return(hmap)
  }
  
  #Function for PCA plot
  plot_pca <- function(dataf, slider_var, comp1, comp2){
      #make plot tib-
      filt_tib <- dataf %>% 
        mutate(variance = apply(dataf[-1], MARGIN = 1, FUN = var), .after = gene) 
      perc_val <- quantile(filt_tib$variance, probs = slider_var/100, na.rm = TRUE)   #calculate percentile
      filt_tib <- filter(filt_tib, variance >= perc_val) 
      pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) 
      #extract variance
      variance <- summary(pca_res)$importance[2,]
      x <- round(variance[comp1]*100, 2)
      y <- round(variance[comp2]*100, 2)
      #produce PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Principal Component Analysis")+
        xlab(str_c(comp1, x, "% variance", sep=" "))+
        ylab(str_c(comp2, y, "% variance", sep=" "))+
        theme_bw()
      return(pca)
      }
  
  
  #---------------------Tab 4: Functions for Differential Expression-------------------------#
  
  #Load data
  load_deg_results <- reactive({
    if (!is.null(input$deg_results$datapath)){
      deg_data <- read_csv(input$deg_results$datapath)
      return(deg_data)}
    else{
      return(NULL)
    }
  })
  
  #Function for Volcano Plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
      y_label = sprintf("-log10(%s)", y_name)
      col_label = sprintf("padj < (1*(10^%s))", slider)
      volc_plot <-  ggplot(dataf, aes(x = get(x_name), y = -log10(get(y_name)), color = padj < (1*(10**slider)))) +
        geom_point() +
        scale_colour_manual(name = col_label, values = setNames(c(color2, color1, "grey"),c(T, F, NA))) +
        theme(legend.position = "bottom", 
              legend.box = "horizontal") +
        xlab(x_name) +
        ylab(y_label) +
        theme_linedraw()
      return(volc_plot)
  }
  
  #Draw and Filter Table
  draw_table <- function(dataf, slider) {
    dataf <- dplyr::filter(dataf, dataf$padj < (1*(10**slider)))
    dataf$pvalue <- lapply(dataf$pvalue, formatC, digits = 9)
    dataf$padj <- lapply(dataf$padj, formatC, digits = 9)
    return(dataf)
  }
  
  
  #---------------------Tab 5: Functions for Gene Set Enrichment Analysis-------------------------#
  
  
  #reactive function -- load data
  load_gsea<- reactive({
      deg_data <- read_csv(input$fgsea_results$datapath)
      return(deg_data)
  })
  
  #barplot of top pathways
  top_pathways <- function(dataf, slider){
    top_pos <- dataf %>% 
      top_n(slider, dataf$NES) %>% 
      arrange(desc(NES)) %>% 
      pull(pathway)
    
    top_neg <- dataf %>% 
      top_n(slider, -dataf$NES) %>% 
      arrange(NES) %>% 
      pull(pathway)
    
    subset <- dataf %>% 
      filter(pathway %in% c(top_pos, top_neg)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(plot_name = str_replace_all(pathway, '_', ' '))
    
    plot <- subset %>% 
      mutate(plot_name = forcats::fct_reorder(factor(plot_name), subset$NES)) %>%
      ggplot() +
      geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'cadetblue4', 'FALSE' = 'darkseagreen1')) + 
      theme_minimal(base_size = 8) +
      ggtitle('fgsea results (Hallmark MSigDB Gene Sets)') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
      coord_flip()
    return(plot)
  }
  
  
  #filtered fgsea table
  draw_gsea_table<- function(dataf, slider, direction){
    data <- dataf %>% mutate(status = ifelse(NES > 0, 'positive', 'negative')) %>% filter(padj<=10**slider)
    if (direction != 'all') {
      filtered <- data %>% filter(data$status == direction)
      return(filtered)
    }
    else {
      return(data)
    }
  }
    

  #fgsea NES scatter plot  
  gsea_scatterplot<- function(dataf, slider){
    filtered <- dataf %>% mutate(filter_status=ifelse(padj<10**slider, "passed filter", "filtered out"))
    my_colors <- c("passed filter" = "turquoise", "filtered out" = "gray53")
    plot <- filtered %>% ggplot(aes(x=NES, y=-log10(padj), color=filter_status)) +
      geom_point()+
      scale_color_manual(values = my_colors)
    return(plot)
  }
  

  
#-----------------------------------------------Output--------------------------------------------#

  #-------------------Outputs from the second tab---------------------#
  
  output$sample_summary <- renderTable({
    summary_table(load_sample_data())
  })
  
  output$sample_data <- renderDataTable({
    load_sample_data()
  }, options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100))) #options to view 25 rows at a time by default
  
  
  output$sample_den_aod <- renderPlot({
    aod_den_plot(load_sample_data())
  })
  
  output$sample_den_pmi <- renderPlot({
    pmi_den_plot(load_sample_data())
  })
  
  output$sample_den_rin <- renderPlot({
    rin_den_plot(load_sample_data())
  })
  
  output$sample_den_aoo <- renderPlot({
    aoo_den_plot(load_sample_data())
  })
  
  #------------------Outputs from the third tab----------------------#
  
  output$counts_summary <- renderTable({
    summary_counts(load_counts(), input$variance_slider, input$nonzero_slider)
  })
  
  output$scatter_median_var <- renderPlot({
    variance_vs_median(load_counts(), input$variance_slider)
  })
  
  output$scatter_median_nonzero <- renderPlot({
    median_vs_nonzeros(load_counts(), input$nonzero_slider)
  })
  
  output$counts_heatmap <- renderPlot({
    plot_heatmap(load_counts(), input$variance_slider, input$nonzero_slider)
  })
  
  output$counts_pca <- renderPlot({
    plot_pca(load_counts(),input$variance_slider, input$component_1, input$component_2)
  })
  
  #------------------------------Outputs from the fourth tab--------------------------------#
  
  output$deg_results <- renderDataTable({
    load_deg_results()
  })
  
  output$deg_volcano <- renderPlot({
    volcano_plot(load_deg_results(), input$x_axis, input$y_axis, input$slider, input$base, input$highlight)
  }, height = 700, width = 900)
  
  
  output$deg_plot_data <- renderDataTable({
    draw_table(load_deg_results(), input$slider)
  })
  
  #------------------------------Outputs from the fifth tab--------------------------------#
  
  # bar chart, tab 1 
  output$gsea_barplot <- renderPlot({
    top_pathways(load_gsea(), input$gsea_pathways_slider)
  })
  
  
  # data table with filtered results
  output$gsea_data<-renderDataTable({
    draw_gsea_table(load_gsea(),input$gsea_data_slider, input$fgsea_direction)
  })

  
  #fgsea  filtered results table download
  output$fgsea_download<-downloadHandler(
    filename= function() {'NES_results2.csv'},
    content = function(file){
      write_csv(draw_gsea_table(load_gsea(),input$gsea_data_slider, input$fgsea_direction), file)
    })
  
  #fgsea normalized expression score -- scatter plot 
  output$gsea_scatter <- renderPlot({
    gsea_scatterplot(load_gsea(), input$gsea_scatter_slider)
    
  })
  
}

shinyApp(ui = ui, server = server)


