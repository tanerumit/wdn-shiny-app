

# Import R packages needed for the app here:
library(shiny)
library(DT)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(tidyr)

library(dplyr)
library(ggplot2)
library(ggmap)
library(ggnetwork)
library(ggnewscale)
library(scales)

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES = c('pip', 'pandas', 'numpy','casadi')

# Begin app server
shinyServer(function(input, output) {

  ################## DO NOT CHANGE THIS PART ###################################
  
  virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
  python_path = Sys.getenv('PYTHON_PATH')
  
  # Create virtual env and install dependencies --- DO NOT CHANGE
  reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
  reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=FALSE)
  reticulate::use_virtualenv(virtualenv_dir, required = T)
  
  ##############################################################################
  
  # Source R and PY functions  
  reticulate::source_python('calculate_network2.py')
  source("functions.R", local = TRUE)

  ## Stress test results upload (UI element)
  output$nodeDataUploadUI = renderUI({fileInput("nodeDataUpload",
    label = "Node table", multiple = F, accept = ".csv", width = '95%')
   })
  
  nodeData <- reactive({req(input$nodeDataUpload)
    read.csv(input$nodeDataUpload$datapath, header = T, sep = ",", stringsAsFactors = T, row.names = NULL)
  })
  
    ## Stress test results upload (UI element)
  output$pipeDataUploadUI = renderUI({fileInput("pipeDataUpload",
       label = "Connection table", multiple = F, accept = ".csv", width = '95%')
  })
  
  pipeData <- reactive({req(input$pipeDataUpload)
    read.csv(input$pipeDataUpload$datapath, header = T, sep = ",", stringsAsFactors = T, row.names = NULL)
  })

  mapData <- reactive({req(nodeData())
    
    # Create the map
    bbox <- c(left = min(nodeData()$lon)-0.20, right = max(nodeData()$lon)+0.20,
              bottom = min(nodeData()$lat)-0.20, top = max(nodeData()$lat)+0.20)
    
    get_stamenmap(bbox = bbox,
      maptype = "toner-lite", color = "bw", crop = TRUE, force = TRUE)
  })
  
  out_reactive <- eventReactive(input$execute,{

    req(nodeData(), pipeData())

    nodes_scn  <- nodeData()
    pipes_scn <- pipeData()
  
    dnodes <- nodes_scn$id[which(nodes_scn$type == "demand")]
    snodes <- nodes_scn$id[which(nodes_scn$type == "supply")] 
    
    nodes_scn$discharge[snodes] <- nodes_scn$discharge[snodes] * 1.05
    nodes_scn$disuse[which(nodes_scn$id %in% input$supplyFailure)] <- 1
    pipes_scn$disuse[which(pipes_scn$id %in% input$linkFailure)] <- 1
            
    simulateWDN(
        edev.change = 0.0005,
        tdev.change = 0.01,
        temp.change = 0.04,
        price.change = 0.01,
        pop.change  = 0.01,
        dnetwork.change = 0.01,
        wqual.change = 0.005,
        edev.elasticity = 1.0,
        price.elasticity = -0.2,
        temp.elasticity = 0.03,
        dom.peak.factor = 1,
        ind.peak.factor = 1,
        year.sim = 2021,
        year.ref = 2021,
        nodes.data = nodes_scn,
        pipes.data = pipes_scn,
        global.output = TRUE)
  })
  
  plotWDN <- reactive({

      req(out_reactive())

      visualizeWDN(
          nodes.data = out_reactive()$nodes,
          pipes.data = out_reactive()$pipes,
          node.fill.var = "rel",
          node.size.var = "discharge",
          edge.color.var = "usage",
          edge.size.var = "diameter",
          background.map = mapData())
  })

  output$plotWDNUI <- renderPlot(plotWDN(),  res = 90)
  
  # Test that the Python functions have been imported
  output$message <- DT::renderDataTable({out_reactive()$summary})
  output$nodeTable <- DT::renderDataTable({
    DT::datatable(nodeData(), options = list(pageLength = 10)) #%>%
     # formatRound(columns=1:ncol(nodeData()))
  })
  output$linkTable <- DT::renderDataTable({
      DT::datatable(pipeData(), options = list(pageLength = 10)) %>%
      formatRound(columns=1:ncol(pipeData()), digits=1)
  })

  # ~~ INTERACTIVE INPUT VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  output$supplyFailureUI <- renderUI({
      
      prod_sites <- nodeData() %>% filter(type == "supply")
      prod_sites_lst <- prod_sites$id
      names(prod_sites_lst) <- prod_sites$label
      
      pickerInput(
          inputId = "supplyFailure",
          label = "Production failures",
          choices = prod_sites_lst,
          selected = NULL,
          options = list(`actions-box` = TRUE),
          multiple = TRUE,
          width = '95%')
  })
  
  output$linkFailureUI <- renderUI({
      
      node_names <- nodeData()$label
      connection_lst <- pipeData()$id
      names(connection_lst) <- paste(node_names[pipeData()$start], node_names[pipeData()$end], sep ="_")
    
      pickerInput(
          inputId = "linkFailure",
          label = "Connection failures",
          choices = connection_lst,
          selected = NULL,
          options = list(`actions-box` = TRUE),
          multiple = TRUE,
          width = '95%')
  })

})