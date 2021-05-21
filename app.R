#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinythemes)

source("analysis.R")

#load("data/musem_workspace.rda")
setup_workspace()

get_intro_text <- function(){
  div(h4("Welcome to Form Segmentation App"), 
         p("This app allows you visualize the data from a study on form perception in classical music",
           "that was carried out by the Max Planck Institute for empirical Aesthetics, Frankfurt/M., Germany"),
      p("Have fun!"),
      style = "width:50%:text-align:justify")
}
input_width <- 300
# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(tags$style(".butt{background-color:#add8e6;} .butt{color: #337ab7;}")),
  
   # Application title
   titlePanel("Form Perception"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        tags$style("#mode:{background-color: #72b573}"),
        selectInput("plot_type", "Plot Type:", c("Lines" = "lines", "Gaussification" = "gauss"), selected = "lines"), 

        width = 2
      ),      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(type ="tabs",
          tabPanel("Plot",
                   plotOutput("marker_plot", width = "900px", height = "600px")),
          tabPanel("Stats",
                   tableOutput("data_stats")),
          tabPanel("Info",
                   p(htmlOutput("introduction")),
          )
        )
      )
   )
)



# Define server logic required to draw a plot
server <- function(input, output, session) {
   message("*** STARTING APP***")
   output$introduction <- renderUI({
     get_intro_text()
   })
   output$data_stats <- renderTable({
     part2_meta %>% mutate(n = 1:n()) %>% select(n, everything())
   })
   output$marker_plot <- renderPlot({
     if(input$plot_type == "lines"){
        plot_marker() 
     }
     else{
       plot_gaussification()
     }
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

