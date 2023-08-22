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
                                start = 0, end = 422,
                                threshold = NULL,
                                with_markers = F,
                                only_markers = F, 
                                external_markers = NULL,
                                alpha = NULL){
  #browser()
  if(is.data.frame(data)){
    if(!("marker" %in% names(data))){
      data <- data %>% rename(marker = time_in_s)
    }
    combined_markers <- gaussification(data$marker, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  }
  else{
    combined_markers <- gaussification(data, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
    
  }
  #print(get_gaussification_peaks(combined_marker))
  #peaks <- tibble(w = get_gaussification_peaks(combined_marker))
  q <- combined_markers %>% ggplot(aes(x = t, y = val))
  if(only_markers || !is.null(external_markers) || !is.null(threshold)){
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
    peaks <- tibble(w = get_gaussification_peaks(combined_markers))
    external_markers <- combined_markers
  }
  intercept <- 0
  #browser()
  if(!is.null(threshold)){
    #browser()
    peaks <- tibble(w = get_gaussification_peaks(combined_markers))
    peak_vals <- combined_markers %>% filter(t %in% peaks$w) %>% pull(val) 
    if(threshold == "median"){
      intercept <- median(peak_vals)
    }
    else if(threshold == "mean"){
      intercept <- mean(peak_vals)
    }
    else if(threshold == "q75"){
      intercept <- quantile(peak_vals)[4]
    }
    else{
      intercept <- as.numeric(threshold)
      if(is.na(intercept)){
        intercept <- 0
      }
    }
    q <- q + geom_hline(yintercept = intercept)
  }
  #browser()
  if(with_markers){
    if("val" %in% names(external_markers)){
      markers <- external_markers %>% 
        filter(val >= intercept, t %in% peaks$w, t >= start, t <= end)
    }
    else{
      markers <- external_markers %>% 
        filter(time_in_s  %in% peaks$w, time_in_s  >= start, time_in_s  <= end) %>% mutate(t = time_in_s)
      
    }
    if(is.null(alpha)){
      alpha <- pmax(.01, pmin(1, .25*(end - start)/nrow(markers)))
    }
    q <- q + geom_vline(data = markers, 
                        aes(xintercept = t), 
                        color = "lightblue4", 
                        alpha = alpha)
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
  labels <- NULL
  #browser()
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
                        aes(xintercept = time_in_s, 
                            linetype = level, 
                            color = level), 
                        linewidth = 1)
    if("label" %in% names(external_markers)){
      q <- q + ggrepel::geom_text_repel(data = external_markers %>% filter(time_in_s >= start, time_in_s <= end), 
                          aes(x = time_in_s + (end - start)*.00, 
                              y = 40, 
                              color = factor(level),
                              label = as.character(label)))
    }
    
  }
  q <- q + theme_minimal() 
  q <- q + scale_color_manual(values = c("indianred", "darkgreen", "coral")) 
  #q <- q + scale_color_brewer(palette = "RdBu", direction = -1) 
  q <- q + labs(x = "Time (s)", y = "Count")
  q <- q + xlim(start, end)
  q <- q + theme(legend.position = "none")
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
  q <- q + scale_color_brewer(palette = "Set1")
  q <- q + labs(x = "Time (s)", title = sprintf("DTW: dist = %.2f, norm = %.2f, d = %.2f, abs(d) = %.2f", 
                                                d$distance, d$normalizedDistance, 
                                                mean(plot_df2$d), mean(abs(plot_df2$d))))
  q <- q + theme(legend.title = element_blank())
  q
  
}


plot_bandwidth <- function(data = seg_stats, base_size  = 11, type = "F1"){
  tmp <- data %>% 
    mutate(max_level = factor(max_level), 
           piece = sprintf("Piece %s", piece), 
           psig = p_nd < .05, 
           sdf = (d_nd < 0), 
           rel = 2 *(sdf) + psig, 
           relation = factor(rel, levels = 0:3,
                             labels = c("smaller non-sig", 
                                        "greater  (sig.)", 
                                        "smaller (non-sig)", 
                                        "greater (sig)")), 
           max_level = sprintf("Level %d", max_level), 
           sig_str = factor(psig, labels = c("Non-sig", "Sig"))) 
  
  if(type == "F1"){
    tmp <- tmp  %>%  
      select(sigma, max_level, piece, sim_f1, sim_f1_mean, relation, sig_str)
    ylab <- "F1"
  }
  else{
    tmp <- tmp %>% 
      select(sigma, max_level, piece, sim_f1 = norm_dist, sim_f1_mean = norm_dist_mean, relation, sig_str) 
    ylab <- "Normalized Distance (sec)"
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
  q <- q + labs(x = "Bandwidth (sec)", y = ylab) 
  q <- q + scale_color_brewer(palette = "Set1")
  q <- q  + geom_smooth(method = "lm", formula = "y ~ poly(x,2)")
  q
}

plot_isi_dist <- function(data = all_boundaries, ground_truth = ground_truth ){
  bs <- get_boundary_stats() %>% 
    mutate(piece = piece_names[piece], 
           Trial = sprintf("Trial %s", trial))
        
  theory_means = bs %>% 
    filter(source == "theory", level == 1) %>% 
    group_by(piece) %>% 
    summarise(m_theory = mean(log_m), .groups = "drop") 
  
  q <- bs %>%
    filter(source != "theory") %>% 
    distinct() %>% 
    ggplot(aes(x = log_m, y = after_stat(count),  fill = Trial)) 
  
  q <- q + geom_histogram(color = "black")
  q <- q + facet_grid(source ~ piece) 
  q <- q + xlim(-3, 6)  
  q <- q + geom_vline(data = theory_means, aes(xintercept = m_theory), linewidth = 1, linetype = "dotted")
  q <- q + theme_bw() 
  q <- q + theme(legend.title = element_blank(), 
             strip.background = element_rect(fill = "white")) 
  q <- q + scale_fill_brewer(palette = "Set1") 
  q <- q + labs(x = "Log ISI")
  q
}

make_analytical_histograms <- function(sigma = 1, max_level = 1){
  # Fig1a: Stimulus 1 komplett
  tmp <- all_boundaries %>% filter(piece == 1)
  gt <- ground_truth %>% filter(piece == 1, level <= max_level) %>% mutate(label = 1:nrow(.))
  fig1a <- plot_marker_histogram(tmp, 
                                 sigma = sigma, 
                                 external_markers = gt,
                                 start = 0,
                                 end = piece_durations[1])
  browser()
  # Fig1b: Stimulus 1, 0-100 sec
  fig1b <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt, start = 0, end = 100)
  # Fig1c: Stimulus 1, 122-190 sec
  fig1c <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt, start = 122, end = 199)
  # Fig1d: Stimulus 1, 175-322 sec (end)
  fig1d <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt, start = 175, end = piece_durations[1])
  # 
  # Fig2a: Stimulus 2 komplett
  tmp <- all_boundaries %>% filter(piece == 2)
  gt <- ground_truth %>% filter(piece == 2, level <= max_level) %>% mutate(label = 1:nrow(.))
  fig2a <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt,  end = piece_durations[2])
  # Fig2b: Stimulus 2, 0-210 sec
  fig2b <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt,  start = 0, end = 210)
  # Fig2c: Stimulus 2, 200-380 sec (end)
  fig2c <- plot_marker_histogram(tmp, sigma = sigma, external_markers = gt,  start = 200, end = piece_durations[2])
  ggsave(plot = fig1a, filename = "figs/fig1a.png", dpi = 300)
  ggsave(plot = fig1b, filename = "figs/fig1b.png", dpi = 300)
  ggsave(plot = fig1c, filename = "figs/fig1c.png", dpi = 300)
  ggsave(plot = fig1d, filename = "figs/fig1d.png", dpi = 300)
  ggsave(plot = fig2a, filename = "figs/fig2a.png", dpi = 300)
  ggsave(plot = fig2b, filename = "figs/fig2b.png", dpi = 300)
  ggsave(plot = fig2c, filename = "figs/fig2c.png", dpi = 300)
}
