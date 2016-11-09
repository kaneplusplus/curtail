library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Critical Value Calculator for 2-Stage Design"),

  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    numericInput("n", "Total Number of Enrollees", value=20, step=1, min=2),
    sliderInput("p", "Response Rate for Stage 1", min=0.01, max=0.99, 
      value=c(0.2, 0.8))
  ),
      
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("minimax_size"),
    tableOutput("optimal_size")
  )
))
