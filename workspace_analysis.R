library(tidyverse)
library(janitor)

piece_durations <- c(322, 421)

read_new_ground_truth <- function(){
  gt  <- bind_rows(read_csv2("data/raw/boundaries_stim1_figureklaus.csv", col_types = "nnnnnn")  %>% 
                     mutate(time_in_s = time_in_s/1000, piece = 1), 
                   read_csv2("data/raw/boundaries_stim2_figureklaus.csv", col_types = "nnnnnn")  %>% 
                     mutate(time_in_s = time_in_s/1000, piece = 2)) %>% 
    mutate(boundary_type = factor(sprintf("%s%s%s", beginning, ending, both), 
                                  labels = c("both", "ending", "beginning")), 
           theory = 1) %>% 
    select(-c(beginning, ending, both))
  gt
}

segment_questions <- c("• Wie wichtig finden Sie diesen Abschnitt?" = "SEG.importance",
                       "• Wie stark ist Ihr Gefühl, dass hier etwas endet?" = "SEG.ending",
                       "• Wie stark ist Ihr Gefühl, dass hier etwas neues kommt?" = "SEG.beginning")
task_questions <- c("• Das Stück hat mir gefallen." = "liking",
                    "• Die Aufgabe war schwierig."  = "difficulty")

BMRI_questions <- 
  c("Wenn ich mit jemandem gemeinsam Musik höre, spüre ich eine besondere Verbindung zu dieser Person." = "BMR.q1",
    "In meiner Freizeit höre ich kaum Musik." = "BMR.q2",
    "Ich höre gerne Musik, die Emotionen enthält." = "BMR.q3",
    "Musik leistet mir Gesellschaft, wenn ich alleine bin." = "BMR.q4",
    "Ich tanze nicht gerne, auch nicht zu Musik, die ich sonst gerne höre." = "BMR.q5",                           
    "Durch Musik fühle ich mich anderen Menschen verbunden." = "BMR.q6",                
    "Ich informiere mich über Musik, die ich mag." = "BMR.q7",                               
    "Ich werde emotional, wenn ich bestimmte Stücke höre." = "BMR.q8", 
    "Musik beruhigt und entspannt mich." = "BMR.q9",                   
    "Musik verleitet mich oft zum Tanzen." = "BMR.q10",
    "Ich suche immer nach neuer Musik." = "BMR.q11",
    "Mir kommen die Tränen oder ich weine, wenn ich eine Melodie höre, die ich sehr mag." = "BMR.q12",
    "Ich singe oder musiziere gerne mit anderen Menschen." = "BMR.q13",
    "Musik hilft mir zu entspannen." = "BMR.q14",               
    "Ich kann nicht anders, als bei meiner Lieblingsmusik mitzusummen oder mitzusingen." = "BMR.q15",
    "Bei einem Konzert fühle ich mich den MusikerInnen und dem Publikum verbunden." = "BMR.q16",   
    "Ich gebe relativ viel Geld für Musik und verwandte Dinge aus." = "BMR.q17",        
    "Ich bekomme manchmal Gänsehaut, wenn ich eine Melodie höre, die mir gefällt." = "BMR.q18",
    "Musik tröstet mich." = "BMR.q19",
    "Wenn ich ein Stück höre, das ich sehr gerne mag, muss ich im Takt klopfen oder mich bewegen." = "BMR.q20")


FCQ_questions <-c(
  "FCQ01" = "FCQ.q1",
  "FCQ02" = "FCQ.q2",                                                                                            
  "FCQ03" = "FCQ.q3",                                                                                            
  "FCQ04" = "FCQ.q4",
  "FCQ05" = "FCQ.q5") 

read_lab_questions <- function(fname = "data/part2_lab_questions.csv"){
  browser()
  labq <- read_csv2(fname) %>% 
    janitor::clean_names() %>% 
    filter(filename != "P31-SegmentExp2.log") %>% 
    select(-c(filename, event, trial, frame_no)) %>% 
    rename(p_id = participant)
  
  q_names <- c(task_questions, segment_questions, BMRI_questions, FCQ_questions)
  labq <- labq %>% mutate(question = q_names[q_text]) %>% select(-q_text, -q_code)
  labq[is.na(labq$stimulus),]$stimulus <- "participant"
  #browser()
  p_meta <- labq %>% 
    filter(stimulus == "participant", question %in% BMRI_questions) %>% 
    pivot_wider(id_cols = p_id, names_from = question, values_from = rating)
  
  piece_meta <- labq %>% 
    filter(question %in% c("liking", "difficulty"), stimulus != "participant") %>% 
    pivot_wider(id_cols = c(p_id, stimulus), 
                names_from = question, 
                values_from = rating)
  
  labq <- labq %>% 
    filter(stimulus != "participant", 
           !(question %in% c("liking", "difficulty"))) %>% 
    pivot_wider(id_cols = c(p_id, stimulus, time_in_s, boundary), 
                names_from = question, 
                values_from = rating)
  labq %>% 
    left_join(p_meta, by = "p_id") %>% 
    left_join(piece_meta, by = c("p_id", "stimulus")) %>% 
    mutate(piece = stimulus %>% str_extract("0[0-9]+") %>% as.numeric(), 
           trial = 2, 
           source = "lab") %>% 
    select(-stimulus)
}

prepare_workspace <- function(){
  online_data <- list.files("data/part2/", pattern = "*.rds", full.names = T)
  map_dfr(online_data, function(x){
    tmp <- x %>% readRDS() 
    p_id <- tmp$session$p_id
    familiarity <- tmp$results$familiarity
    complete <- tmp$session$complete
    GMS.general <- tmp$GMS$General 
    DEG.gender <- c("female", "male", "rather not say", "other")[tmp$DEG$Gender]
    DEG.age <- tmp$DEG$Age/12 
    markers <- tmp %>% pluck("MSM") %>% pluck("q1") %>% pull(marker) %>% str_split(",") %>% pluck(1) %>% as.numeric()
    messagef("%s - %s", tmp$DEG$Gender, DEG.gender)
    if(length(tmp$DEG$Gender) == 0){
      return(NULL)
    }
    #print(mean(sign(diff(markers) < 0 )))
    tibble(time_in_s = markers/1000, p_id = p_id, GMS.general = GMS.general, DEG.gender = DEG.gender, DEG.age = DEG.age, familiarity = familiarity)
  })
}

fix_part2_data <- function(part2_data){
  part2_data <- part2_data %>% filter(count > 1) 
  #browser()
  #for unknown reasons some session had the same psychTestR p_id, which should not happen at all, so here we go
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 4.611111) < .01,]$p_id <- "OP05-1"
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 3.722222) < .01,]$p_id <- "OP05-2"
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 2.166667) < .01,]$p_id <- "OP05-3"

  part2_data[part2_data$p_id == "OP17" & abs(part2_data$GMS.general - 3.277778) < .01,]$p_id <- "OP17-1"
  part2_data[part2_data$p_id == "OP17" & abs(part2_data$GMS.general - 5.222222) < .01,]$p_id <- "OP17-2"
  
  part2_data <- part2_data %>% filter(p_id != "OP05" & p_id != "OP17")
  
  #also some marker seem to be IOI instead of absolute time points, we fixed those that can be fixed
  bad_ids <- part2_data %>% group_by(p_id) %>% summarise(s = mean(diff(marker) < 0, na.rm = T)) %>% filter(s > 0) %>% pull(p_id)
  map_dfr(unique(part2_data$p_id), function(pid){
    tmp <- part2_data %>% filter(p_id == pid)
    if(pid %in% bad_ids){
      fixed_marker <- cumsum(tmp$marker)  
      s <- mean(diff(fixed_marker) < 0)
      m <- max(fixed_marker)
      # if(m > 422){
      #   browser()
      # }
      tmp$marker <- fixed_marker
      if(m > 450){
        messagef("Removed %s with max = %.3f > 450)", pid, s, m)
        return(NULL)
      }
      messagef("Fixed marker for %s (s = %s, max = %.3f)", pid, s, m)
      
    }
    tmp
  }) %>% filter(marker < 422)
}

setup_workspace <- function(){
  all_online <- readr::read_csv("data/part2_all_online_clean.csv")  %>% 
    fix_part2_data()
  #browser()
  ground_truth <- readr::read_csv("data/part2_ground_truth.csv")

  boundaries_lab <- readr::read_csv("data/part2_boundaries_lab.csv") %>% 
    mutate(trial = ((trial - 1) %% 2) + 1)
  
  metadata <- bind_rows(
    readr::read_csv("data/part2_metadata_lab.csv") %>% 
      mutate(source = "lab", 
             GMS.general = GMS.general/18),
    readr::read_csv("data/part2_metadata_online.csv") %>% 
      mutate(source = "online", 
             piece = str_extract(stimulus, "01|02") %>% as.integer()),
  ) %>% select(-stimulus, -part)
  
  all_boundaries <- bind_rows(
    boundaries_lab %>% 
      mutate(source = "lab"), 
    all_online %>% 
      select(p_id, time_in_s = marker, piece) %>% 
      mutate(source = "online", trial = 1)    
  )  
  if(file.exists("data/seg_stats.rds")){
    seg_stats <- readRDS("data/seg_stats.rds")
  }
  else{
    seg_stats = get_segmentation_stats(ground_truth = ground_truth %>% 
                                         filter(level < 3), 
                                       band_widths = seq(0.5, 4, .125))
    saveRDS(seg_stats, "data/seg_stats.rds")
  }
  # labq <- read_lab_questions()
  # 
  # boundaries_lab_annotations <- all_boundaries %>% 
  #   filter(source == "lab", trial == 2) %>% 
  #   left_join(labq, by = c("time_in_s", "boundary", "p_id", "piece", "source", "trial"))
  boundaries_lab_annotations <- readRDS("data/boundaries_lab_annotations.rds")
  assign("seg_stats", seg_stats, globalenv())
  assign("all_online", all_online, globalenv())
  assign("all_boundaries", all_boundaries, globalenv())
  assign("ground_truth", ground_truth, globalenv())
  assign("boundaries_lab_annotations", boundaries_lab_annotations, globalenv())
  assign("boundaries_lab", boundaries_lab, globalenv())
  assign("metadata", metadata, globalenv())
  
} 