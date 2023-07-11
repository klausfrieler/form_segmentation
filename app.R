#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinyjs)
library(shinythemes)
library(tidyverse)

#source("analysis.R")
source("workspace_part2.R")

setup_workspace()

stimulus_choices <- c("Louise Farrenc, Nonet, 1st mov" = "1", "John K. Paine, Symphony No. 1, 1st mov" = "2")
trial_choices <- c("All trials"  = "both",  "Trial 1" = "1" , "Trial 2" = "2")
plot_choices <- c("Histogram" = "histogram", "Densities" = "gauss", "DTW Alignment" = "dtw_alignment")

input_width <- 300
text_size <- 14

get_intro_text <- function(){
  div(h3("Welcome to the Form Segmentation Analysis App"), 
         p("This app allows you visualize and inspect the data from a study on form and phrase perception in classical music",
           "that was carried out by the Max Planck Institute for empirical Aesthetics, Frankfurt/M., Germany"),
      p("Have fun!"),
      style = "width:50%;text-align:justify")
}



# Define UI for application that draws a histogram
ui <- fluidPage(
   useShinyjs(),
   tags$head(tags$style(".butt{background-color:#add8e6;} .butt{color: #337ab7;}")),
   
   # Application title
   titlePanel("Form Perception"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         tags$style("#mode:{background-color: #72b573}"),
         selectInput("source", "Data Source:", choices = c("Lab" = "lab", "Online" = "online", "Lab + Online" = "both"), selected = "lab"), 
         selectInput("piece", "Piece:", choices = stimulus_choices, selected = stimulus_choices[1]), 
         selectInput("trial", "Trials", choices = trial_choices, selected = "both"),
         tags$hr(),
         selectInput("plot_type", "Plot Type:", plot_choices, selected = "gauss"), 
         sliderInput("bw", "Bandwidth: ", min = .5, max = 5, value = 1, step = .25),
         sliderInput("range", "Time Range: ", min = 0, max = 322, value = c(0, 322), step = 10),
         selectInput("threshold", "Threshold: ", choices = c("No" = "0", "Mean" = "mean", "Median" = "median"), selected = "0"),
         
         tags$hr(),
         checkboxInput("add_ground_truth", label = "Add Ground Truth", value = FALSE),
         checkboxInput("random_gt", label = "Random Ground Truth", value = FALSE),
         #actionButton("generate_random_ground_truth", label = "Go!", width = input_width/2),
         sliderInput("max_level", "Max Level: ", min = 1, max = 3, value = 3, step = 1),
         
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
                              plotOutput("marker_plot", width = "900px", height = "600px"),
                              tableOutput("sim_f1"))#,
                     #tabPanel("Stats",
                     #         tableOutput("data_stats"))
         )
      )
   )
)



# Define server logic required to draw a plot
server <- function(input, output, session) {
   message("*** STARTING APP***")
   g_random_ground_truth <- NULL
   random_ground_truth <- reactive({
      print("generated")
      input$random_gt
      simulate_ground_truth(ground_truth, 
                            piece = as.integer(input$piece),
                            max_level = as.integer(input$max_level),
                            theory = 1) 
      })
   
   output$introduction <- renderUI({
     get_intro_text()
   })
   
   observeEvent(input$generate_random_ground_truth, {
      print("generate_random_ground_truth")
      g_random_ground_truth <- random_ground_truth()
      
   })
   observeEvent(input$random_gt, {
      print("generate_random_ground_truth")
      if(input$random_gt){
         print("triggered")
         g_random_ground_truth <- random_ground_truth()
      }
      
   })
   observeEvent(input$plot_type, {
      if(input$plot_type == "dtw_alignment"){
         print("toggle")
         shinyjs::hide("add_ground_truth")
      }
      else{
         shinyjs::show("add_ground_truth")
         
      }
   })
   
   observeEvent(input$source, {
      if(input$source == "online"){
         updateSelectInput(session, inputId = "trial", selected = "1", choices = c("Trial 1" = "1"))
      }
      else{
         updateSelectInput(session, inputId = "trial", selected = "both", choices = trial_choices)
         
      }
   })
   observeEvent(input$piece, {
      curr_level <- as.integer(input$max_level)
      if(input$piece == "2"){
         updateSliderInput(session, inputId = "max_level", value = min(2, curr_level), min = 1, max = 2)
      }
      else{
         updateSliderInput(session, inputId = "max_level", value = min(3, curr_level), min = 1, max = 3)
         
      }
   })
   
   output$overall_stats <- renderTable({
      p_id_stats <- metadata %>% 
         group_by(source) %>% 
         distinct(p_id, gender, age, GMS.general) %>% 
         summarise(n_female = sum(gender == "female", na.rm = T), 
                   n_male = sum(gender == "male", na.rm = T), 
                   mean_age = mean(age, na.rm = T), 
                   mean_GMS = mean(GMS.general, na.rm = T), 
                   n_unique = n(), .groups = "drop")
      
      sources <- metadata %>% 
         group_by(source) %>% 
         summarise(n = n(),  
                   mean_marker = mean(count),
                   mean_difficulty = mean(difficulty, na.rm = T),  
                   mean_liking = mean(liking, na.rm = T), .groups = "drop") %>% 
         bind_cols(p_id_stats %>% select(-source))
      # stimuli <- metadata %>% 
      #    group_by(piece) %>% 
      #    summarise(n = n(),  
      #              n_unique = n_distinct(p_id),
      #              n_male = sum(gender == "male", na.rm = T),
      #              n_female = sum(gender == "female", na.rm = T),
      #              mean_marker = mean(count),
      #              mean_age = mean(age, na.rm = T), 
      #              mean_GMS = mean(GMS.general, na.rm = T), 
      #              mean_difficulty = mean(difficulty, na.rm = T),  
      #              mean_liking = mean(liking, na.rm = T), .groups = "drop") %>% 
      #    rename(type = stimulus)
      bind_rows(sources) %>% 
         select(source, starts_with("n"), mean_age, mean_GMS, everything())
   })
   
   output$data_stats <- renderTable({
      metadata %>% 
         filter(stimulus == input$stimulus) %>% 
         mutate(n = 1:n()) %>% 
         mutate(p_id = sprintf("%s...", substr(p_id, 1, 8))) %>% 
         select(n, everything())
   })
   
   current_data <- reactive({
      data <- all_boundaries %>% filter(piece == as.integer(input$piece))
      
      if(input$source != "both"){
         data <- data %>% filter(source == input$source)
         messagef("Filtered for source: %s", input$source)
      }
      
      if(input$trial != "both"){
         data <- data %>% filter(trial == as.integer(input$trial))
         messagef("Filtered for trial %s", input$trial)
      }
      data
   })
   
   current_ground_truth <- reactive({
      gt <- ground_truth
      if(input$add_ground_truth){
         if(input$random_gt){
            gt <- g_random_ground_truth
            if(is.null(gt)){
               gt <- random_ground_truth()
            }
         }
      }
      gt
   })
   
   current_threshold <- reactive({
      threshold <- input$threshold
      if(threshold == "0"){
         threshold <- 0
      }
      threshold
   })   
   
   output$marker_plot <- renderPlot({
      data <- current_data()
      end <- piece_durations[as.integer(input$piece)]
      external_markers <- NULL
      #input$generate_random_ground_truth
      gt <- current_ground_truth()
      if(input$add_ground_truth){
         external_markers <- gt %>% 
            filter(piece == as.integer(input$piece), 
                   level <= as.integer(input$max_level), 
                   theory == 1,
                   boundary_type != "ending")
      }
      #messagef("Max level (slider): %s, max level (data): %s", input$max_level, max(external_markers$level))
      if(input$plot_type == "gauss"){
         plot_gaussification(data = data, end = end, sigma = as.numeric(input$bw), 
                             threshold = current_threshold(),
                             with_marker = !is.null(external_markers), 
                             external_markers = external_markers) + 
            xlim(as.integer(input$range[1]), as.integer(input$range[2])) + 
            theme(axis.title = element_text(size = text_size),
                  axis.text.x =  element_text(size = text_size),
                  axis.text.y =  element_text(size = text_size))
      }
      else if(input$plot_type == "dtw_alignment"){
         threshold <- current_threshold()
         gt <- current_ground_truth()

         compare_segmentations(segs = data, 
                               gt = gt, 
                               piece = as.integer(input$piece), 
                               theory = 1, 
                               max_level = as.numeric(input$max_level),
                               start = as.numeric(input$range[1]),
                               end = as.numeric(input$range[2]),
                               sigma = as.numeric(input$bw), 
                               threshold = threshold,
                               only_plot = T) + 
            theme(axis.title = element_text(size = text_size),
                  axis.text.x =  element_text(size = text_size),
                  axis.text.y =  element_text(size = text_size))
      }
   })
   output$sim_f1 <- renderTable({
      #browser()
      compare_segmentations(segs = current_data(), 
                            gt = current_ground_truth(), 
                            piece = as.integer(input$piece), 
                            theory = 1, 
                            max_level = as.numeric(input$max_level),
                            start = as.numeric(input$range[1]),
                            end = as.numeric(input$range[2]),
                            sigma = as.numeric(input$bw), 
                            threshold = current_threshold(),
                            with_plot = F) %>% pluck("summary")
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

