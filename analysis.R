library(tidyverse)
library(mccr)
library(cvAUC)
library(MASS)
library(caret)
library(proxy)
library(ggridges)
library(psych)

#library(MSM)

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

rect_func <- function(t, t0, half_width = .5){
  ifelse(abs(t - t0) <= half_width, 1.0, 0.0)
}

gaussification <- function(onsets, deltaT = 0.1, sigma = 2, weights = NULL, start = -1, end = 450, use_rect_func = FALSE, with_singles = F){
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
  if(use_rect_func){
    singles <- map(seq_along(onsets), function(i){
      weights[i] * rect_func(time_points, onsets[i], half_width = sigma)
    })
    
  }
  else{
    sig_fact <- 1/sigma/sigma 
    singles <- map(seq_along(onsets), function(i){
      weights[i] * exp( -.5 * (onsets[i] - time_points) * (onsets[i] - time_points)*sig_fact)
    })
    
  }
  val <- reduce(singles, function(x,y) x + y)
  ret <- tibble(t = time_points, val = val)
  if(with_singles){
    dfs <- map_dfr(seq_along(singles), function(x) tibble(t = time_points, val = singles[[x]], onset = onsets[[x]]))
    ret <- list(sum = ret, singles = singles)
  }
  ret
}

get_sims <- function(time_points, ground_truth, start = 0, end = 330, bw = 1){
  #browser()
  bin_v <- hist(time_points, breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  bin_gt <- hist(ground_truth, breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  
  bin_v[bin_v > 1] <- 1
  bin_gt[bin_gt > 1] <- 1
  
  sim_f1 <- confusionMatrix(as.factor(bin_v), as.factor(bin_gt), mode = "everything", positive="1")$byClass[7] #function to return prediction vector based on bw
  if(is.na(sim_f1)){
    sim_f1 <- 0
  }  
  return(tibble(sim_f1 = sim_f1))
  sim_ppv <- confusionMatrix(as.factor(bin_v), as.factor(bin_gt), mode = "everything", positive="1")$byClass[3] #function to return prediction vector based on bw
  if(is.na(sim_ppv)){
    sim_ppv <- 0
  }  
  #browser()
  sim_AUC <- cvAUC::AUC(bin_v, bin_gt) 
  sim_mccr <- mccr::mccr(bin_v, bin_gt) 
  sim_euc_dist <- proxy::dist(rbind(bin_v, bin_gt), "euclidean")  %>% as.numeric()
  sim_cosine <- proxy::simil(rbind(bin_v, bin_gt), "cosine") %>% as.numeric() 
  sim_corr <- cor(bin_v, bin_gt) 
  return(tibble(sim_AUC = sim_AUC, 
                sim_mccr = sim_mccr, 
                sim_euc_dist = sim_euc_dist/sqrt(length(bin_v)), 
                sim_cosine = sim_cosine, 
                sim_corr = sim_corr, 
                sim_ppv = sim_ppv,
                sim_f1 = sim_f1, 
                n_bins = length(bin_v)))
  
}

#function taking test and GT vectors, testy type, and bin width - returning results df
#functional programming - everything is a function and you apply them to multiple vars
get_all_sims <- function(segment_data, vect_gt, bw_range = seq(.5, 1.5, .25)){
  p_ids <- unique(segment_data$p_id)
  
  map_dfr(p_ids, function(id){ #works like lapply; function is anonymous functions also called lamda
    
    map_dfr(bw_range, function(bw){ # like a nested loop 

      messagef("p_id = %s, bw = %.2f", id, bw)
      # a <-   get_sims(segment_data[segment_data$p_id ==id,]$time_in_s, vect_gt, bw) %>%  
      #   mutate(bw = bw, id = id)
      get_sims(segment_data[segment_data$p_id ==id,]$time_in_s, vect_gt, bw) %>%  
        mutate(bw = bw, p_id = id) # adds columns to data - also able to do operations on columns
    })
  })
} 

get_all_sims_by_trials <- function(segment_data, ground_truth, bw_range = seq(.75, .75, 0), max_level = 2){
  p_ids <- unique(segment_data$p_id)
  trials <- unique(segment_data$trial)
  pieces < unique(segment_data$piece)
  map_dfr(p_ids, function(id){ #works like lapply; function is anonymous functions also called lamda
    
    map_dfr(bw_range, function(bw){ # like a nested loop 
      
      messagef("p_id = %s, bw = %.2f", id, bw)
      # a <-   get_sims(segment_data[segment_data$p_id ==id,]$time_in_s, vect_gt, bw) %>%  
      #   mutate(bw = bw, id = id)
      get_sims(segment_data[segment_data$p_id ==id,]$time_in_s, vect_gt, bw) %>%  
        mutate(bw = bw, p_id = id) # adds columns to data - also able to do operations on columns
    })
  })
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
                group_by(p_id) %>% mutate(d = c(diff(marker), NA)) %>% 
                summarise(d = mean(d, na.rm = T), m = max(marker), .groups ="drop"), 
              by = "p_id")
  part2_meta <- part2_all %>% distinct(p_id, stimulus, count, difficulty, liking, age, gender, GMS.general)%>% 
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
  combined_marker <- gaussification(data$marker, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
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

