library(tidyverse)
library(MSM)

messagef <- function(...) message(sprintf(...))
printf <- function(...) print(sprintf(...))

gaussification <- function(onsets, deltaT = 0.1, sigma = 2, weights = NULL, start = -1, end = 450){
  if(is.null(start)){
    start <- min(onsets) - 2 * sigma
  }
  if(is.null(end)){
    end <- max(onsets) + 2*sigma
  }
  messagef("start = %.2f, end = %.2f", start,  end)
  ret <- NULL
  if(is.null(weights)){
    weights <- rep(1, length(onsets))
  }
  for(t in seq(start, end, deltaT)){
    val <- sum(weights * exp( -.5 * (onsets - t) * (onsets - t)/sigma/sigma))
    ret <- rbind(ret, tibble(t = t, val = val))       
  }
  return(ret)
}

setup_workspace <- function(result_dir = "data/part2"){
  part2_all <- MSM::read_MSM_data(result_dir, expand_markers = T) %>% 
    filter(lubridate::day(time_ended) > 20)
  part2_all <- part2_all %>% 
    left_join(part2_all %>% 
                group_by(p_id) %>% mutate(d = c(diff(marker), NA)) %>% 
                summarise(d = mean(d, na.rm = T), m = max(marker), .groups ="drop"), 
              by = "p_id")
  part2_meta <- part2_all %>% distinct(p_id, difficulty, liking, age, gender, GMS.general)
  part2 <- part2_all %>% 
    filter(m > 200, !is.na(d), marker < 450)
  #combined_marker <- gaussification(part2$marker, sigma = 2, deltaT = .1)
  assign("part2", part2, globalenv())
  assign("part2_all", part2_all, globalenv())
  assign("part2_meta", part2_meta, globalenv())
  #assign("combined_marker", combined_marker, globalenv())
}
plot_marker <- function(data = part2){
  q <- data %>% ggplot(aes(x = marker, y = as.integer(factor(p_id)), color = p_id))
  q <- q + geom_vline(aes(xintercept = marker), alpha = .2)
  q <- q + geom_point(size = 2) 
  q <- q + theme_minimal()
  q <- q + theme(legend.position = "none")  
  q <- q + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
  q <- q + labs(x = "Time (s)", y ="Participant")
  
  q
}
plot_gaussification <- function(data = part2, sigma = 2, deltaT = .1){
  combined_marker <- gaussification(part2$marker, sigma = sigma, deltaT = deltaT)
  q <- combined_marker %>% ggplot(aes(x = t, y = val))
  q <- q + geom_line(color = "indianred") 
  q <- q + theme_minimal() 
  q <- q + scale_color_brewer(palette = "RdBu") 
  q <- q + labs(x = "Time (s)", y = "Combined segments")
  q  
}