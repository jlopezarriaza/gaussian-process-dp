# ui.R
#ui.R
library('shinyjs')

shinyUI(fluidPage(
  theme = "yeti.css",
  titlePanel(h2("Non-parametric Multi-species Fisheries Management", 
                align = "center")),
  withMathJax(),
  fluidRow(
    column(width = 3,
           wellPanel(helpText(
             h5("This app simulates data from a multi-species model, fits
                       a GP to the data and performs dynamic programming to find 
                       the optimal harvest of each species",
                style = "font-size:13pt;")),
             
             h4("General Parameters:"),
             sliderInput("time.series.length",
                         label = "Time Series Length:",
                         min = 10, max = 50, value = 25),
             sliderInput("V",
                         label = "Process Stochasticity Variance:",
                         min = 0, max = 0.1, value = 0.01),
             selectInput("num", "Number of Species", 
                         choices = seq(from = 1, to = 2, by = 1)),
             shinyjs::useShinyjs(),
             actionButton(inputId = "submit_loc1", 
                          label = "Generate Data", 
                          icon("paper-plane"), 
                          width = '100%',
                          style = "color: #fff; background-color:
                          #337ab7; border-color: #2e6da4"),
             hidden(actionButton(inputId = "submit_loc2",
                                 label = "Run Dynamic Programming", 
                                 icon("ship"),
                                 width = '100%',
                                 style = "color: #fff; background-color:
                                 #337ab7; border-color: #2e6da4")),
             hidden(actionButton(inputId = "submit_loc3",
                                 label = "Do Forward Simulations", 
                                 icon("forward"),
                                 width = '100%',
                                 style = "color: #fff; background-color:
                                 #337ab7; border-color: #2e6da4")),
             hidden(actionButton(inputId = "submit_loc4",
                                 label = "Reset",
                                 icon("forward"),
                                 width = '100%',
                                 style = "color: #fff; background-color:
                                 #337ab7; border-color: #2e6da4")))),
    ###########################################################################
    column(width = 8, 
           conditionalPanel("input.num ==1",
                            helpText('$$x_{t+1}=x_t\\exp(r_1(1-\\frac{x_t}{k_1}))$$',
                                     style = "font-family: 'times'; font-size:16pt;")),
           conditionalPanel("input.num ==2",
                            helpText('$$x_{1,t+1}=x_{1,t}\\exp(r_1(1-\\frac{x_{1,t}}{k_1})+c_1 x_{2,t}  )$$',
                                     style = "font-family: 'times'; font-size:14pt;"),
                            helpText('$$x_{2,t+1}=x_{2,t}\\exp(r_2(1-\\frac{x_{2,t}}{k_2})+c_2 x_{1,t} )$$',
                                     style = "font-family: 'times'; font-size:14pt;"))),
    ###########################################################################
    column(width = 4,
           wellPanel(h5("Carrying Capacities:"),
                     uiOutput("input_ui"),
                     conditionalPanel("input.num > 1", 
                                      sliderInput("c1",
                                                  label = withMathJax("$$c_1$$"),
                                                  min = -1, 
                                                  max = 1, 
                                                  value = 0,
                                                  step = 0.1)))),
    column(width = 4,
           wellPanel(h5("Intrinsic Growth Rate:"),
                     uiOutput("input_ui2"),
                     conditionalPanel("input.num >1",
                                      sliderInput("c2",
                                                  label = withMathJax("$$c_2$$"),
                                                  min = -1,
                                                  max = 1, 
                                                  value = 0,
                                                  step = 0.1))))),
  ###########################################################################
  fluidRow(column(width = 8,
                  offset = 2,
                  plotOutput("graph1"))),
  fluidRow(column(width = 8,
                  offset = 2,
                  plotOutput("graph2"))),
  fluidRow(column(width = 8,
                  offset = 2,
                  plotOutput("graph3"))),
  fluidRow(conditionalPanel("input.num > 1",
                            column(width = 8,
                                   offset = 2,
                                   plotOutput("graph5")))),
  fluidRow(column(width = 8,
                  offset = 2,
                  plotOutput("graph4")))
))