library(tidyverse)
library(MSM)

messagef <- function(...) message(sprintf(...))
printf <- function(...) print(sprintf(...))

get_arg_max <- function(x){
  d_x <- sign(diff(x))
  dd_x <- diff(d_x)
  which(dd_x < 0) + 1
}

get_gaussification_peaks <- function(gauss_data, with_plot = F){
  if(is.list(gauss_data) & "sum" %in% names(gauss_data)){
    gauss_data <- gauss_data[["sum"]]
  }
  peak_pos <- get_arg_max(gauss_data$val)  
  t_max <- gauss_data[peak_pos,]$t
  if(with_plot){
    q <- ggplot(gauss_data, aes(x = t, y = val)) + geom_line() + geom_vline(data = tibble(x = t_max), aes(xintercept = x))
    print(q)
  }
  t_max
}

gaussification <- function(onsets, deltaT = 0.1, sigma = 2, weights = NULL, start = -1, end = 450, with_singles = F){
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
  time_points <- seq(start, end, deltaT)
  sig_fact <- 1/sigma/sigma 
  singles <- map(seq_along(onsets), function(i){
    weights[i] * exp( -.5 * (onsets[i] - time_points) * (onsets[i] - time_points)*sig_fact)
  })
  val <- reduce(singles, function(x,y) x + y)
  ret <- tibble(t = time_points, val = val)
  if(with_singles){
    dfs <- map_dfr(seq_along(singles), function(x) tibble(t = time_points, val = singles[[x]], onset = onsets[[x]]))
    ret <- list(sum = ret, singles = singles)
  }
  ret
}

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
    distinct(p_id, stimulus, difficulty, liking, age, gender, GMS.general) %>% 
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
    filter(lubridate::day(time_ended) > 20)
  part2_all <- part2_all %>% 
    left_join(part2_all %>% 
                group_by(p_id) %>% mutate(d = c(diff(marker), NA)) %>% 
                summarise(d = mean(d, na.rm = T), m = max(marker), .groups ="drop"), 
              by = "p_id")
  part2_meta <- part2_all %>% distinct(p_id, stimulus, difficulty, liking, age, gender, GMS.general)%>% 
    mutate(part = "PART2")
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

plot_gaussification <- function(data = part2, sigma = 2, deltaT = .1, start = -1, end = 450, with_marker = F){
  combined_marker <- gaussification(data$marker, sigma = sigma, deltaT = deltaT, end = end)
  peaks <- tibble(w = get_gaussification_peaks(combined_marker))
  q <- combined_marker %>% ggplot(aes(x = t, y = val))
  q <- q + geom_line(color = "indianred") 
  if(with_marker){
    peaks <- tibble(w = get_gaussification_peaks(combined_marker))
    q <- q + geom_vline(data = peaks, aes(xintercept = w), color = "lightblue4")
  }
  q <- q + theme_minimal() 
  q <- q + scale_color_brewer(palette = "RdBu") 
  q <- q + labs(x = "Time (s)", y = "Combined segments")
  q  
}

