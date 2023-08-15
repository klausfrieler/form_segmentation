library(tidyverse)

piece_durations <- c(322, 421)
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
  
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 4.611111) < .01,]$p_id <- "OP05-1"
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 3.722222) < .01,]$p_id <- "OP05-2"
  part2_data[part2_data$p_id == "OP05" & abs(part2_data$GMS.general - 2.166667) < .01,]$p_id <- "OP05-3"

  part2_data[part2_data$p_id == "OP17" & abs(part2_data$GMS.general - 3.277778) < .01,]$p_id <- "OP17-1"
  part2_data[part2_data$p_id == "OP17" & abs(part2_data$GMS.general - 5.222222) < .01,]$p_id <- "OP17-2"
  
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
             GMS.general = GMS.general/19),
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
  assign("seg_stats", seg_stats, globalenv())
  
  assign("all_online", all_online, globalenv())
  assign("all_boundaries", all_boundaries, globalenv())
  assign("ground_truth", ground_truth, globalenv())
  
  assign("boundaries_lab", boundaries_lab, globalenv())
  assign("metadata", metadata, globalenv())
  
} 