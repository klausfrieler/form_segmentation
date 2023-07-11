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

plot_gaussification <- function(data = part2, 
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

