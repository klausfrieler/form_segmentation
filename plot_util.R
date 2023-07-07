library(tidyverse)

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

plot_gaussification <- function(data = part2, sigma = 2, deltaT = .1, start = -1, end = 450, with_marker = F, only_markers = F){
  combined_marker <- gaussification(data$marker, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  print(get_gaussification_peaks(combined_marker))
  #peaks <- tibble(w = get_gaussification_peaks(combined_marker))
  q <- combined_marker %>% ggplot(aes(x = t, y = val))
  if(only_markers){
    with_marker <- T
  }
  if(!only_markers){
    q <- q + geom_line(color = "indianred") 
  }
  if(with_marker){
    peaks <- tibble(w = get_gaussification_peaks(combined_marker))
    q <- q + geom_vline(data = peaks, aes(xintercept = w), color = "lightblue4")
  }
  q <- q + theme_minimal() 
  q <- q + scale_color_brewer(palette = "RdBu") 
  q <- q + labs(x = "Time (s)", y = "Combined segments")
  q  
}

