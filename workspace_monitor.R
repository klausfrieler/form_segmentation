library(tidyverse)
setup_workspace <- function(part1_results = "data/part1", part2_results = "data/part2"){
  read_part1_data(part1_results)
  read_part2_data(part2_results)
  assign("both_parts", bind_rows(part1, part2), globalenv())
  assign("both_parts_meta", bind_rows(part1_meta, part2_meta), globalenv())
}

read_part1_data <- function(result_dir = "data/part1"){
  messagef("Setting up part1 data from %s", result_dir)
  part1_all <- MSM::read_MSM_data(result_dir, expand_markers = T) 
  part1_all <- part1_all %>% 
    left_join(part1_all %>% 
                group_by(p_id) %>% mutate(d = c(diff(marker), NA)) %>% 
                summarise(d = mean(d, na.rm = T), m = max(marker), .groups ="drop"), 
              by = "p_id")
  part1_meta <- part1_all %>% 
    distinct(p_id, stimulus, difficulty, count, liking, age, gender, GMS.general) %>% 
    mutate(part = "PART1")
  part1 <- part1_all %>% 
    filter(!is.na(d))
  #combined_marker <- gaussification(part2$marker, sigma = 2, deltaT = .1)
  assign("part1", part1, globalenv())
  assign("part1_all", part1_all, globalenv())
  assign("part1_meta", part1_meta, globalenv())
  #assign("combined_marker", combined_marker, globalenv())
}

read_part2_data <- function(result_dir = "data/part2"){
  messagef("Setting up part2 data from %s", result_dir)
  part2_all <- MSM::read_MSM_data(result_dir, expand_markers = T) %>% 
    filter(lubridate::day(time_ended) > 20 | lubridate::month(time_ended) >= 6)
  part2_all <- part2_all %>% 
    left_join(part2_all %>% 
                group_by(p_id) %>% 
                mutate(d = c(diff(marker), NA)) %>% 
                summarise(d = mean(d, na.rm = T), 
                          s = mean(sign(d) < 0), 
                          m = max(marker), .groups ="drop"), 
              by = "p_id")
  browser()
  
  part2_meta <- part2_all %>% 
    distinct(p_id, stimulus, count, difficulty, liking, age, gender, GMS.general)%>% 
    mutate(part = "PART2")
  part2 <- part2_all %>% 
    filter(m > 200, !is.na(d), marker < 450)
  #combined_marker <- gaussification(part2$marker, sigma = 2, deltaT = .1)
  assign("part2", part2, globalenv())
  assign("part2_all", part2_all, globalenv())
  assign("part2_meta", part2_meta, globalenv())
  #assign("combined_marker", combined_marker, globalenv())
}

