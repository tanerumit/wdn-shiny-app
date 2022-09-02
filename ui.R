# Import R packages needed for the UI
library(shiny)
library(shinycssloaders)
library(DT)
library(shinyWidgets)

# Begin UI for the R + reticulate example app
ui <- fluidPage(
  titlePanel("Friesland Water Distribution Network Analysis"),
  helpText("Simulation exercise to assess system resilience to water shortages"),
  br(),
  # ---------------- Sidebar panel with changeable inputs ----------------- #
  sidebarPanel(
      h1(strong("1) Upload datasets"), style = "font-size:18px;"),
      helpText("First upload the node and connection datasets in CSV format. 
                Uploaded datasets can be viewed in the main panel."),
      uiOutput('nodeDataUploadUI'),
      uiOutput('pipeDataUploadUI'),
      br(),
      h1(strong("2) Select system failures"), style = "font-size:18px;"),
      helpText("Select failures in production sites and the connections between the sites"),
      uiOutput('supplyFailureUI'),
      uiOutput('linkFailureUI'),
      br(),
      actionButton("execute", "Run"),
  ),
  # ---------------- Sidebar panel with changeable inputs ----------------- #
  mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Sites",
          helpText("This table displays main features of the production and demand sites"),
          br(),
          DT::dataTableOutput('nodeTable'),
        ),
        tabPanel("Connections",
          helpText("This table displays main features of the connections between the sites"),
          br(),
          DT::dataTableOutput('linkTable')
        ),
        tabPanel("Results plot",
          helpText("This visual shows the distribution of simulated water shortages accross the system"),
          br(),
          plotOutput("plotWDNUI", height = "800px", width = "900px") %>% 
          withSpinner()
        ),
        tabPanel("Results table",
          helpText("This table displays calculated metrics from the simulation"),
          br(),
          DT::dataTableOutput('message')
        )
      )
  ) #mainpanel close
)