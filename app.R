library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(png)
library(grid)
library(pals)
library(paletteer)
library(zip)


ui = fluidPage(
  useShinyjs(),  # Initialize shinyjs
  # Custom banner/logo at the top
  fluidRow(
    column(2, tags$div(
      style = "text-align: left; margin-bottom: 20px;",
      tags$img(src = "vib_rf_cancer_biology_rgb_pos.png", height = "100px", style = "max-width: 100%; height: auto;")
    )),
    column(8, tags$div(
      style = "display: flex; align-items: center; justify-content: center; height: 100px; margin-bottom: 20px;",
      h1("Metabolite Data Analysis")
    )),
    column(2, tags$div(
      style = "text-align: right; margin-bottom: 20px; display: flex; align-items: center; justify-content: flex-end;",
      tags$img(src = "BIOINFORMATICS.png", height = "100px", style = "max-width: 100%; height: auto; margin-left: 10px;"),
      tags$img(src = "SpMet_pos_Vandamme.png", height = "100px", style = "max-width: 100%; height: auto; margin-left: 10px;")
    ))
  ),
  tabsetPanel(
    tabPanel("Input Data",
             sidebarLayout(
               sidebarPanel(
                 fileInput("raw_data", "Choose raw data (CSV file)", accept = ".csv"),
                 fileInput("metabolites_data", "Choose metabolites data (CSV file)", accept = ".csv"),
                 fileInput("pathway_annotation", "Choose pathway annotation (CSV file)", accept = ".csv"),
                 fileInput("annotation_files", "Choose annotation files (CSV files)", accept = ".csv", multiple = TRUE),
                 fileInput("msi_images", "Choose MSI images (PNG files)", accept = ".png", multiple = TRUE),
                 fluidRow(
                   column(6,
                          fluidRow(
                            column(12, actionButton("load_data", "Process all data"))
                          ),
                          fluidRow(
                            column(12, downloadButton("downloadData", "Download all processed data", disabled = TRUE))
                          )
                   ),
                   column(6,
                          fluidRow(
                            column(12, actionButton("process_msi", "Process MSI images"))
                          ),
                          fluidRow(
                            column(12, downloadButton("downloadAllMSI", "Download MSI images", disabled = TRUE))
                          )
                   )
                 )
               ),
               mainPanel(
                 column(2, tags$div(
                   height = "1000px", style = "max-width: 100%; height: auto; margin-left: 10px;",
                   tags$img(src = "lung_fendt_cropped.png")
                 )),
                 div(id = "loading_spinner", withSpinner(textOutput("loading_text")), style = "display:none;")
               )
             )
    ),
    tabPanel("Heatmaps",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("pathwayUI")
               ),
               mainPanel(
                 fluidRow(
                   column(4, selectInput("heatmapFileType", "File Type", choices = c("png", "pdf"), selected = "png")),
                   column(4, numericInput("heatmapWidth", "PDF / PNG width:", min = 4, max = 16, value = 8, step = 0.5)),
                   column(4, numericInput("heatmapHeight", "PDF / PNG height:", min = 4, max = 12, value = 6, step = 0.5))
                 ),
                 fluidRow(
                   column(12, downloadButton("downloadHeatmap", "Download heatmap"))
                 ),
                 plotOutput("heatmapPlot", width = "60vw", height = "60vh")
               )
             )
    ),
    tabPanel("Boxplots and MSI images",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("pathwayUI_boxplot"),
                 uiOutput("metaboliteUI")
               ),
               
               mainPanel(
                 fluidRow(
                   column(4, selectInput("boxplotFileType", "File Type", choices = c("png", "pdf"), selected = "png")),
                   column(4, numericInput("boxplotWidth", "PDF / PNG width:", min = 4, max = 16, value = 8, step = 0.5)),
                   column(4, numericInput("boxplotHeight", "PDF / PNG height:", min = 4, max = 12, value = 6, step = 0.5))
                 ),
                 fluidRow(
                   column(12, downloadButton("downloadBoxplot", "Download boxplot"))
                 ),
                 plotOutput("boxplotPlot", width = "60vw", height = "60vh"),
                 fluidRow(
                   column(12, downloadButton("downloadMSI", "Download MSI image"))
                 ),
                 plotOutput("msiImagePlot", width = "60vw", height = "60vh")
               )
             )
    )
  )
)

    
    


server = function(input, output, session) {
  # Increase file size limit
  options(shiny.maxRequestSize = 200*1024^2)

  data = reactiveValues(
    raw = NULL,
    metabolites = NULL,
    pathways = NULL,
    annotation = NULL,
    combined_dt = NULL,
    combined_dt_median = NULL,
    msi_images_dt = NULL
  )
  
  # Sanitize naming
  sanitize_dir = function(name) {
    temp = gsub("[^A-Za-z0-9]", "_", name)
    return(gsub("_+", "_", temp))
  }
  
  observeEvent(input$load_data, {
    req(input$raw_data, input$metabolites_data, input$pathway_annotation, input$annotation_files, input$msi_images)
    shinyjs::show("loading_spinner")
    output$loading_text = renderText("Loading data, please wait...")

    data$raw = fread(input$raw_data$datapath)
    data$metabolites = fread(input$metabolites_data$datapath, fill = TRUE, skip = "m/z", header = TRUE, sep = ";")
    data$pathways = fread(input$pathway_annotation$datapath)
    
    #basedir = "~/Downloads/Result_melanie"
    #data$raw = fread(paste0(basedir, "/raw_data/raw_data.csv"))
    #data$metabolites = fread(paste0(basedir, "/raw_data/metabolites-final.csv"))
    #data$pathways = fread(paste0(basedir, "/metabolites_pathway_annotation.csv"))
    #annots = paste0(basedir, "/annotation/", c("1_metastasis-1_control.csv", "2_Healthy-1_control.csv"))
    #msi_paths = list.files(basedir, pattern = ".png")
    data$annotation = lapply(1:nrow(input$annotation_files), function(i) {
    #data$annotation = lapply(annots, function(file) {
      dt = fread(input$annotation_files$datapath[i])
      dt[, annot_1 := gsub("-.*", "", input$annotation_files$name[i])]
      dt[, annot_2 := gsub(".*-|.csv", "", input$annotation_files$name[i])]
      return(dt)
    }) |> rbindlist()
    
    dt = data$raw[data$metabolites[, .(`m/z`, Name)], on = "m/z", roll = "nearest"]
    dt_m = data.table::melt(dt, id.vars = c("m/z", "Name"), variable.name = "Spot index")
    dt_m[, `Spot index` := as.integer(gsub("Spot ", "", `Spot index`))]
    
    data$combined_dt = merge(dt_m, data$annotation, by = "Spot index")
    data$combined_dt = merge(data$combined_dt, data$pathways, by.x = "Name", by.y = "Scils - Name", allow.cartesian = TRUE)
    
    setorder(data$combined_dt, annot_1, annot_2)
    data$combined_dt[, annot_1 := factor(annot_1, levels = unique(annot_1))]
    data$combined_dt[, annot_2 := factor(annot_2, levels = unique(annot_2))]
    
    # Process heatmap data
    data$combined_dt[, z_score := scale(value), by = Name]
    
    # Pathway specific
    data$combined_dt_median = data$combined_dt[, .(median_zscore = median(z_score)),
                                                   by = .(Name, Pathway, annot_1, annot_2)]
    
    setorder(data$combined_dt_median, -Name, annot_1, annot_2)
    data$combined_dt_median[, Name := factor(Name, levels = unique(Name))]
    data$combined_dt_median[, annot_1 := factor(annot_1, levels = unique(annot_1))]
    data$combined_dt_median[, annot_2 := factor(annot_2, levels = unique(annot_2))]
    
    # All metabolites
    data$combined_dt_nodups = unique(data$combined_dt, by = c("Name", "Spot index", "annot_1", "annot_2", "z_score"))
    data$combined_dt_median_all = data$combined_dt_nodups[, .(median_zscore = median(z_score)),
                                                          by = .(Name, annot_1, annot_2)]
    
    setorder(data$combined_dt_median_all, -Name, annot_1, annot_2)
    data$combined_dt_median_all[, Name := factor(Name, levels = unique(Name))]
    data$combined_dt_median_all[, annot_1 := factor(annot_1, levels = unique(annot_1))]
    data$combined_dt_median_all[, annot_2 := factor(annot_2, levels = unique(annot_2))]
    
    # Load MSI images
    data$msi_images_dt = lapply(1:nrow(input$msi_images), function(i) {
      mz_value = as.numeric(sub(" .*", "", input$msi_images$name[i]))
      data.table(name = input$msi_images$name[i], "m/z" = mz_value, imagepath = input$msi_images$datapath[i])
    }) |> rbindlist()

    # Match MSI image with name
    data$msi_images_dt = data$msi_images_dt[data$metabolites[, .(`m/z`, Name)], on = "m/z", roll = "nearest"]

    #print(merge(data$msi_images_dt, data$metabolites[, .(`m/z`, Name), by = "m/z"]))
    
    updateSelectInput(session, "pathway", choices = unique(data$pathways$Pathway))
    updateSelectInput(session, "pathway_boxplot", choices = unique(data$pathways$Pathway))
    
    # Create a temporary directory to store the files
    temp_dir = tempdir()
    full_download_dir = file.path(temp_dir, "full_download")
    
    # Loop through pathways and save heatmap for each pathway
    for (pathway in unique(data$combined_dt$Pathway)) {
      sanitized_pathway = sanitize_dir(pathway)
      pathway_dir = file.path(full_download_dir, sanitized_pathway)
      dir.create(paste0(pathway_dir, "/boxplots"), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(pathway_dir, "/msi-images"), showWarnings = FALSE, recursive = TRUE)
      
      
      heatmap_file = file.path(pathway_dir, paste("heatmap-", sanitized_pathway, ".png", sep = ""))
      heatmap_plot = ggplot(data$combined_dt_median[Pathway == pathway], aes(x = annot_1, y = Name, fill = median_zscore)) +
        facet_wrap(~annot_2) + 
        geom_tile() +
        paletteer::scale_fill_paletteer_c("pals::coolwarm") +
        labs(fill = "Median Z-score", x = "", title = input$pathway) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      ggsave(heatmap_file, plot = heatmap_plot, device = "png", width = 10, height = 8)
      
      # Loop through metabolites in the current pathway and save boxplot for each metabolite
      for (metabolite in unique(data$combined_dt[Pathway == pathway, Name])) {
        boxplot_file = file.path(pathway_dir, paste("boxplots/", sanitize_dir(metabolite), ".png", sep = ""))
        
        # Make sure it's unique
        boxplot_dt = unique(data$combined_dt[Name == metabolite], by = "Spot index")
        
        boxplot_plot = ggplot(boxplot_dt, aes(x = annot_1, y = value)) +
          facet_wrap(~annot_2) +
          geom_violin() +
          geom_boxplot(width = 0.1) +
          labs(y = "Raw Intensity", x = "", title = metabolite) +
          scale_color_brewer(palette = "Set1") + 
          theme_cowplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        ggsave(boxplot_file, plot = boxplot_plot, device = "png", width = 10, height = 8)
        
        msi_image = data$msi_images_dt[Name == metabolite, imagepath]
        if (length(msi_image) > 0) {
          dir.create(paste0(pathway_dir, "/msi_images"), showWarnings = FALSE)
          file.copy(msi_image, paste0(pathway_dir, "/msi-images/", sanitize_dir(metabolite), ".png", sep = ""))
        }
      }
    }
    
    # Create a ZIP file containing all the processed files
    zip_file = file.path(temp_dir, "processed_data.zip")
    zip::zip(
      zipfile = zip_file,
      files = list.files(full_download_dir, full.names = FALSE, recursive = TRUE),
      root = full_download_dir)
    
    # Provide the ZIP file for download
    output$downloadData = downloadHandler(
      filename = function() {
        "processed_data.zip"
      },
      content = function(file) {
        file.copy(zip_file, file)
      },
      contentType = "application/zip"
    )
    
    shinyjs::enable("downloadData")
    shinyjs::hide("loading_spinner")
    showNotification("All data processed successfully!", type = "message")
    
  })
  
  output$pathwayUI = renderUI({
    selectInput("pathway", "Select Pathway", choices =  c("All", unique(data$combined_dt$Pathway)))
  })
  
  output$pathwayUI_boxplot = renderUI({
    selectInput("pathway_boxplot", "Select Pathway", choices = c("Any", unique(data$combined_dt$Pathway)))
  })
       
  output$metaboliteUI = renderUI({
    if (input$pathway_boxplot == "Any") {
      selectInput("metabolite", "Select Metabolite", choices = unique(data$combined_dt$Name))
    } else {
      selectInput("metabolite", "Select Metabolite", choices = unique(data$combined_dt[Pathway == input$pathway_boxplot, Name]))
    }
  })
  
  # Plot heatmap
  output$heatmapPlot = renderPlot({
    req(input$pathway)
    if (input$pathway == "All") {
      ggplot(data$combined_dt_median_all, aes(x = annot_1, y = Name, fill = median_zscore)) +
        facet_wrap(~annot_2) + 
        geom_tile() +
        scale_fill_paletteer_c("pals::coolwarm", limits = c(min(data$combined_dt_median$median_zscore),
                                                            max(data$combined_dt_median$median_zscore))) +
        labs(fill = "Median Z-score", x = "", title = "All metabolites") +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    } else {
      ggplot(data$combined_dt_median[Pathway == input$pathway], aes(x = annot_1, y = Name, fill = median_zscore)) +
        facet_wrap(~annot_2) + 
        geom_tile() +
        scale_fill_paletteer_c("pals::coolwarm", limits = c(min(data$combined_dt_median$median_zscore),
                                                            max(data$combined_dt_median$median_zscore))) +
        labs(fill = "Median Z-score", x = "", title = input$pathway) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
  })

  
  # Plot boxplot
  output$boxplotPlot = renderPlot({
    req(input$pathway_boxplot, input$metabolite)
    # Make sure it's unique
    boxplot_dt = unique(data$combined_dt[Name == input$metabolite], by = "Spot index")
    
    ggplot(boxplot_dt, aes(x = annot_1, y = value)) +
      facet_wrap(~annot_2) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      labs(y = "Raw Intensity", x = "", title = input$metabolite) +
      scale_color_brewer(palette = "Set1") + 
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  })

  
  # Plot MSI
  output$msiImagePlot = renderPlot({
    req(input$pathway_boxplot, input$metabolite)
    msi_image = data$msi_images_dt[Name == input$metabolite, imagepath]
    if (length(msi_image) > 0) {
      plot.new()
      grid.raster(readPNG(msi_image))
    } else {
      plot.new()
      text(0.5, 0.5, "No MSI image available for this metabolite", cex = 1.5)
    }
  })

  
  # Download handler for heatmap
  output$downloadHeatmap = downloadHandler(
    filename = function() {
      paste("heatmap-", input$metabolite, ".", input$heatmapFileType, sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = input$heatmapFileType, width = input$heatmapWidth, height = input$heatmapHeight)
    }
  )
  
  # Download handler for boxplot
  output$downloadBoxplot = downloadHandler(
    filename = function() {
      paste("boxplot-", input$metabolite, ".", input$boxplotFileType, sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = input$boxplotFileType, width = input$boxplotWidth, height = input$boxplotHeight)
    }
  )
  # Download handler for MSI image
  output$downloadMSI = downloadHandler(
    filename = function() {
      paste("msi-", input$metabolite, ".png", sep = "")
    },
    content = function(file) {
      msi_image = data$msi_images_dt[Name == input$metabolite, imagepath]
      file.copy(msi_image, file)
    }
  )
  
  # Observe that outdir has been specified and enable save all button
  observeEvent(input$save_dir, {
    dir = parseDirPath(roots, input$save_dir)
    output$selected_path = renderText({ dir })
    shinyjs::enable("save_all")
    shinyjs::enable("split_msi")
  })
  
  
  # Run split MSI
  observeEvent(input$process_msi, {
    req(input$msi_images, input$metabolites_data, input$pathway_annotation)
    shinyjs::show("loading_spinner")
    output$loading_text = renderText("Processing and saving, please wait...")
    removeModal()

    # Load pathway annotation
    data$metabolites = fread(input$metabolites_data$datapath, fill = TRUE, skip = "m/z", header = TRUE, sep = ";")
    data$pathways = fread(input$pathway_annotation$datapath)
    
    meta_pathways = merge(data$metabolites, data$pathways, by = "Name")
    
    # Load MSI images
    data$msi_images_dt = lapply(1:nrow(input$msi_images), function(i) {
      mz_value = as.numeric(sub(" .*", "", input$msi_images$name[i]))
      data.table(name = input$msi_images$name[i], "m/z" = mz_value, imagepath = input$msi_images$datapath[i])
    }) |> rbindlist()
    
    # Match MSI image with name
    data$msi_images_dt = data$msi_images_dt[meta_pathways[, .(`m/z`, Name, Pathway)], on = "m/z", roll = "nearest"]

    # Create a temporary directory to store the files
    temp_dir = tempdir()
    msi_download_dir = file.path(temp_dir, "msi_download")
    dir.create(msi_download_dir, showWarnings = FALSE)
    
    for (pathway in unique(data$msi_images_dt$Pathway)) {
      sanitized_pathway = sanitize_dir(pathway)
      pathway_dir = file.path(msi_download_dir, sanitized_pathway)
      dir.create(paste0(pathway_dir, "/msi-images"), showWarnings = FALSE, recursive = TRUE)
      
      for (metabolite in unique(data$msi_images_dt[Pathway == pathway, Name])) {
        msi_image = data$msi_images_dt[Pathway == pathway & Name == metabolite, imagepath]
        
        if (length(msi_image) > 0) {
          file.copy(msi_image, paste0(pathway_dir, "/msi-images/", sanitize_dir(metabolite), ".png", sep = ""))
        }
      }
    }
    
    # Create a ZIP file containing all the processed files
    msi_zip_file = file.path(temp_dir, "msi_images.zip")
    zip::zip(
      zipfile = msi_zip_file,
      files = list.files(msi_download_dir, full.names = FALSE, recursive = TRUE),
      root = msi_download_dir)
    
    shinyjs::hide("loading_spinner")
    shinyjs::enable("downloadAllMSI")
    showNotification("All MSI images processed!", type = "message")
    
    # Provide the ZIP file for download
    output$downloadAllMSI = downloadHandler(
      filename = function() {
        "msi_images.zip"
      },
      content = function(file) {
        file.copy(msi_zip_file, file)
      },
      contentType = "application/zip"
    )
    
  })
  # Clean up temp directory when session ends
  session$onSessionEnded(function() {
    unlink(c(msi_download_dir, full_download_dir, msi_zip_file, zip_file), recursive = TRUE)
  })
}

shinyApp(ui = ui, server = server)