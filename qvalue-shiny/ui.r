library(shiny)
library(qvalue, lib.loc= .libPaths()[1])

shinyUI(pageWithSidebar(
  headerPanel("q-value estimation for FDR control"),
   sidebarPanel(
    fileInput('file1', 'Choose file containing p-values',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE),
    
gsub("label class=\"radio\"", "label class=\"radio inline\"",radioButtons("sep", "Separator:",
             c("Comma" = ",",
               "Semicolon" = ";",
               "Tab" = "\t",
               "Space" = " ",
               "End-of-line" = "\n"), 'Comma')),
   # radioButtons('sep', 'Separator',
   #              c(Comma=',',
   #                Semicolon=';',
   #                Tab='\t',
   #                Space=' ',
   #                End_of_line='\n'),
   #              'Comma'),
    wellPanel(
      p( strong( HTML("&pi;<sub>0</sub>"), "estimate inputs")),
      selectInput("pi0method", p("Choose a", HTML("&pi;<sub>0</sub>"), "method:"), 
                   choices = c("smoother", "bootstrap")),
      sliderInput(inputId = "lambda", label = p(HTML("&lambda;"),"range"),
                  min = 0, max = 1, value = c(0, 0.95), step = 0.01),
      
      numericInput("step", p(HTML("&lambda;"),"step size:"), 0.05),
      numericInput("sdf", "smooth df:", 3.0),
      checkboxInput(inputId = "pi0log", label = p(HTML("&pi;<sub>0</sub>"), "smoother log"),     value = FALSE)
    ),
    wellPanel(
      p(strong("Local FDR inputs")),
      selectInput("transf", "Choose a transformation method:", 
                  choices = c("probit", "logit")),
      checkboxInput(inputId = "trunc", label = "truncate local FDR values",     value = TRUE),
      checkboxInput(inputId = "mono", label = "monotone",     value = TRUE),
      numericInput("adj", "adjust:", 1.5),
      numericInput("eps", "threshold:", 10^-8)
    ),
    wellPanel(
      p(strong("Output")),
      sliderInput("fdr", 
                  "FDR level:",
                  step = 0.01,
                  value = 0.05,
                  min = 0, 
                  max = 1),
      checkboxInput(inputId = "pfdr", label = "pfdr",     value = FALSE),
      downloadButton('downloadData', 'Download Output')
    )

  ),
  mainPanel(

    tabsetPanel(id="tabSelected",
      tabPanel("About", h4("Using the App"), uiOutput("about"), h4("References"), uiOutput("ref")),
     # tabPanel("Figures", h4("Plot"), plotOutput("qvaluePlot"), h4("Histogram"), plotOutput("qvalueHist"),   h4("Summary"),  verbatimTextOutput("summary") ),
      tabPanel("Output", uiOutput("subTabs")),
      tabPanel("Help", uiOutput("help")))

  )
))

