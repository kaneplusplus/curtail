library(shiny)

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {

  get_inputs <- reactive({
    print("here")
    print(p)
    list(p1=input$p[1], p2=input$p[2], n=input$n)
  })

  output$minimax_size <- renderTable({
    print(input$p)
    inputs <- get_inputs()
    a <- inputs$n + 3
    iris[1:10,]
  })

  output$optimal_size <- renderTable({
    inputs <- get_inputs()
    a <- inputs$n + 3
    iris[1:10,]
  })

})

