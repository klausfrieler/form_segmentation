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

plot_gaussification <- function(data = boundaries_lab, 
                                sigma = 2, 
                                deltaT = .1, 
                                start = -1, end = 450, 
                                threshold = NULL,
                                with_markers = F,
                                only_markers = F, 
                                external_markers = NULL){
  if(is.data.frame(data)){
    if(!("marker" %in% names(data))){
      data <- data %>% rename(marker = time_in_s)
    }
    
    combined_marker <- gaussification(data$marker, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  }
  else{
    combined_marker <- gaussification(data, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
    
  }
  #print(get_gaussification_peaks(combined_marker))
  #peaks <- tibble(w = get_gaussification_peaks(combined_marker))
  q <- combined_marker %>% ggplot(aes(x = t, y = val))
  if(only_markers || !is.null(external_markers)){
    with_markers <- T
  }
  if(!only_markers){
    q <- q + geom_line(color = "indianred") 
  }
  if(!is.null(external_markers)){
    if(is.data.frame(external_markers)){
      peaks <- tibble(w = external_markers$time_in_s)  
    }
    else{
      peaks <- tibble(w = external_markers)  
    }
  }
  else{
    peaks <- tibble(w = get_gaussification_peaks(combined_marker))
    
  }
  if(with_markers){
    q <- q + geom_vline(data = peaks %>% filter(w >= start, w <= end), aes(xintercept = w), color = "lightblue4")
  }
  if(!is.null(threshold)){
    #browser()
    peaks <- get_gaussification_peaks(combined_marker)
    peak_vals <- combined_marker %>% filter(t %in% peaks) %>% pull(val) 
    if(threshold == "median"){
      intercept <- median(peak_vals)
    }
    else if(threshold == "mean"){
      intercept <- mean(peak_vals)
    }
    else{
      intercept <- as.numeric(threshold)
      if(is.na(intercept)){
        intercept <- 0
      }
    }
    q <- q + geom_hline(yintercept = intercept)
  }
  q <- q + theme_minimal() 
  q <- q + scale_color_brewer(palette = "RdBu") 
  q <- q + labs(x = "Time (s)", y = "Combined segments")
  q <- q + xlim(start, end)
  q  
}

plot_marker_histogram <- function(data = boundaries_lab, 
                                  sigma = 2, 
                                  external_markers = NULL, 
                                  start = 0, 
                                  end = 430){
  q <- data %>% ggplot(aes(x = time_in_s, y = after_stat(count))) 
  q <- q + geom_histogram(binwidth = sigma, fill = "lightblue", color = "black")
  
  if(!is.null(external_markers)){
    if(is.data.frame(external_markers)){
      if("level" %in% names(external_markers)){
        peaks <- external_markers %>% mutate(level = factor(level))
      }
      else{
        peaks <- tibble(time_in_s = external_markers$time_in_s, level = "1")
      }
    }
    else{
      peaks <- tibble(time_in_s = external_markers, level = "1")  
    }
    q <- q + geom_vline(data = peaks %>% filter(time_in_s >= start,
                                                time_in_s <= end), 
                        aes(xintercept = time_in_s, linetype = level, color = level), 
                        linewidth = 1)
  }
  q <- q + theme_minimal() 
  q <- q + scale_color_manual(values = c("indianred", "darkgreen", "coral")) 
  #q <- q + scale_color_brewer(palette = "RdBu", direction = -1) 
  q <- q + labs(x = "Time (s)", y = "Count")
  q <- q + xlim(start, end)
  q  
}

plot_dtw_alignment <- function(x, y = NULL){
  if(class(x) == "dtw"){
    d <- x
    x <- d$query %>% as.vector()
    y <- d$reference %>% as.vector()
  }
  else{
    if(is.null(y)){
      stop("y must have value")
    }
    d <- dtw(x, y, keep.internals = T)
  }
  plot_df <- bind_rows(tibble(x = x, type = "query"), tibble(x = y, type = "reference"))
  plot_df2 <- tibble(x = x[d$index1], y = y[d$index2]) %>% mutate(d = x - y)
  
  #browser()
  q <- plot_df %>% ggplot(aes(x = x, y  = type, colour = type)) + geom_point(size = 5)
  q <- q + geom_segment(data = plot_df2, 
                        aes(x = x, y = "query", xend = y, yend = "reference"),
                        colour = "black", 
                        arrow = arrow(length = unit(0.30, "cm"), 
                                      ends = "last", 
                                      type = "closed"))
  q <- q + theme_minimal()
  q <- q + labs(x = "Time (s)", title = sprintf("DTW: dist = %.2f, norm = %.2f, d = %.2f, abs(d) = %.2f", 
                                                d$distance, d$normalizedDistance, 
                                                mean(plot_df2$d), mean(abs(plot_df2$d))))
  q <- q + theme(legend.title = element_blank())
  q
  
}


plot_bandwidth <- function(data = seg_stats, base_size  = 11, type = "F1"){
  tmp <- data %>% 
    mutate(max_level = factor(max_level), piece = sprintf("Piece %s", piece), psig = p_nd < .05, 
           sdf = (d_nd < 0), rel = 2 *(sdf) + psig, 
           relation = factor(rel, 
                             labels = c("smaller non-sig", "greater  (sig.)", "smaller (non-sig)", "greater (sig)")), 
           max_level = sprintf("Level %d", max_level), 
           sig_str = factor(psig, labels = c("Non-sig", "Sig"))) 
  
  if(type == "F1"){
    tmp <- tmp  %>%  
      select(sigma, max_level, piece, sim_f1, sim_f1_mean, relation, sig_str)
  }
  else{
    tmp <- tmp %>% 
      select(sigma, max_level, piece, sim_f1 = norm_dist, sim_f1_mean = norm_dist_mean, relation, sig_str) 
  }
    
  tmp <- tmp %>% 
    pivot_longer(-c(sigma, max_level, piece, relation, sig_str)) %>% 
    mutate(name = factor(name, labels = c("Theory 1", "Random (N = 25)")))
  
  q <- tmp %>%  
    ggplot(aes(x = sigma, y = value, color = name)) 
  q <- q + geom_point(aes(shape = sig_str)) 
  q <- q + geom_line() 
  q <- q + facet_grid(piece ~ factor(max_level)) 
  q <- q + theme_bw(base_size = base_size)
  q <- q + theme(legend.title = element_blank(), 
                 legend.position = "bottom", 
                 strip.background = element_rect(fill = "white")) 
  q <- q + labs(x = "Bandwidth (sec)", y = "Normalized Distance (sec)") 
  q <- q + scale_color_brewer(palette = "Set1")
  q <- q  + geom_smooth(method = "lm", formula = "y ~ poly(x,2)")
  q
}