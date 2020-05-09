#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    tags$head(tags$style(
      HTML('
         #sidebar {
          background-color: #d0d5f0;
         }
         ')
    )),
            
  # Application title
  titlePanel(h1("Classifying Graft Rejection Using Genomic Data",
                style = "font-weight: bold; color: #1685a0; font-family: sans-serif")),
  headerPanel(""),
  
  # Sidebar with a slider input for number of bins
  
  selectInput("classifier", "Please choose a model:",
              list("KNN", "Random Forests", "Naive Bayes", "Multilayer Perceptron", "Linear Support Vector Machine", "Decision Tree")
  ),
  
  
  sidebarLayout(
    sidebarPanel(id = "sidebar",
      sliderInput("nFolds", "Number of Folds for Cross Validation:", min = 2, max = 50, value = 5),
      sliderInput("repeats", "Number of Repeats:", min = 1, max = 50, value = 5),
      sliderInput("nComps", "Number of Principal Components:", min = 1, max = 88, value = 20),
      sliderInput("k2", "(KNN) Maximum k Value:", min = 1, max = 50, value = 5),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("accPlot")
    )
  ),
  actionButton(inputId = "modelConfirm", label = "Run Selected Model"),
  headerPanel(""),
  textOutput("result"),
  headerPanel(""),
  textOutput("accuracy")
))
