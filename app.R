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
part1_result_dir <- "data/part1"
part2_result_dir <- "data/part2"
setup_workspace(part1_result_dir, part2_result_dir)
stimulus_choices <- sort(unique(both_parts$stimulus))
names(stimulus_choices) <- stimulus_choices %>% 
   str_replace(".wav", "") %>% 
   str_replace("part", "Part ") %>% 
   str_replace("_", ": ") 

get_intro_text <- function(){
  div(h3("Welcome to the Form Segmentation Analysis App"), 
         p("This app allows you visualize and inspect the data from a study on form and phrase perception in classical music",
           "that was carried out by the Max Planck Institute for empirical Aesthetics, Frankfurt/M., Germany"),
      p("Have fun!"),
      style = "width:50%;text-align:justify")
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
        selectInput("stimulus", "Stimulus:", stimulus_choices, selected = stimulus_choices[1]), 
        selectInput("plot_type", "Plot Type:", c("Lines" = "lines", "Gaussification" = "gauss"), selected = "lines"), 

        width = 2
      ),      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(type ="tabs",
               tabPanel("Info",
                        htmlOutput("introduction"),
                        h4("Summary"),
                        tableOutput("overall_stats")
                    ),
               tabPanel("Plot",
                        plotOutput("marker_plot", width = "900px", height = "600px")),
               tabPanel("Stats",
                        tableOutput("data_stats"))
        )
      )
   )
)



# Define server logic required to draw a plot
server <- function(input, output, session) {
   message("*** STARTING APP***")
   check_data1 <- reactiveFileReader(1000, session, part1_result_dir, read_part1_data)
   check_data2 <- reactiveFileReader(1000, session, part2_result_dir, read_part2_data)
   assign("both_parts", bind_rows(part1, part2), globalenv())
   output$introduction <- renderUI({
     get_intro_text()
   })
   output$overall_stats <- renderTable({
      check_data1()
      check_data2()
      assign("both_parts", bind_rows(part1, part2), globalenv())
      p_id_stats <- both_parts_meta %>% 
         group_by(part) %>% 
         distinct(p_id, gender, age, GMS.general) %>% 
         summarise(n_female = sum(gender == "female", na.rm = T), 
                   n_male = sum(gender == "male", na.rm = T), 
                   mean_age = mean(age, na.rm = T), 
                   mean_GMS = mean(GMS.general, na.rm = T), 
                   n_unique = n(), .groups = "drop")
      
      parts <- both_parts_meta %>% group_by(part) %>% 
         summarise(n = n(),  
                   mean_marker = mean(count),
                   mean_difficulty = mean(difficulty, na.rm = T),  
                   mean_liking = mean(liking, na.rm = T), .groups = "drop") %>% 
         rename(type = part) %>% 
         bind_cols(p_id_stats)
      stimuli <- both_parts_meta %>% 
         group_by(stimulus) %>% 
         summarise(n = n(),  
                   n_unique = n_distinct(p_id),
                   n_male = sum(gender == "male", na.rm = T),
                   n_female = sum(gender == "female", na.rm = T),
                   mean_marker = mean(count),
                   mean_age = mean(age, na.rm = T), 
                   mean_GMS = mean(GMS.general, na.rm = T), 
                   mean_difficulty = mean(difficulty, na.rm = T),  
                   mean_liking = mean(liking, na.rm = T), .groups = "drop") %>% 
         rename(type = stimulus)
      bind_rows(parts, stimuli) %>% 
         select(type, starts_with("n"), mean_age, mean_GMS, everything())
      })
   output$data_stats <- renderTable({
      check_data1()
      check_data2()
      assign("both_parts", bind_rows(part1, part2), globalenv())
         
      both_parts_meta %>% 
         filter(stimulus == input$stimulus) %>% 
         mutate(n = 1:n()) %>% 
         mutate(p_id = sprintf("%s...", substr(p_id, 1, 8))) %>% 
         select(n, everything())
   })
   
   output$marker_plot <- renderPlot({
      check_data1()
      check_data2()
      assign("both_parts", bind_rows(part1, part2), globalenv())
      data <- both_parts %>% filter(stimulus == input$stimulus)
      end <- 450
      if(str_detect(input$stimulus, "part1")){
         end <- 70
      }
      if(input$plot_type == "lines"){
         plot_marker(data = data) 
      } else{
       plot_gaussification(data = data, end = end)
     }
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

