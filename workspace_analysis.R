library(tidyverse)

piece_durations <- c(322, 421)

setup_workspace <- function(){
  all_online <- readr::read_csv("data/part2_all_online_clean.csv")  %>% 
    filter(count != 0)
  
  ground_truth <- readr::read_csv("data/part2_ground_truth.csv")
  
  boundaries_lab <- readr::read_csv("data/part2_boundaries_lab.csv") %>% 
    mutate(trial = ((trial - 1) %% 2) + 1)
  
  metadata <- bind_rows(
    readr::read_csv("data/part2_metadata_lab.csv") %>% 
      mutate(source = "lab"),
    readr::read_csv("data/part2_metadata_online.csv") %>% 
      mutate(source = "online", piece = str_extract(stimulus, "01|02") %>% as.integer()),
  ) %>% select(-stimulus, -part)
  
  all_boundaries <- bind_rows(
    boundaries_lab %>% 
      mutate(source = "lab"), 
    all_online %>% 
      select(p_id, time_in_s = marker, piece) %>% 
      mutate(source = "online", trial = 1)    
  )
  assign("all_online", all_online, globalenv())
  assign("all_boundaries", all_boundaries, globalenv())
  assign("ground_truth", ground_truth, globalenv())
  
  assign("boundaries_lab", boundaries_lab, globalenv())
  assign("metadata", metadata, globalenv())
  
} 