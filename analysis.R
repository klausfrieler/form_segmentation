library(tidyverse)
#library(mccr)
#library(cvAUC)
#library(MASS)
library(caret)
library(proxy)
library(psych)
library(dtw)

#library(MSM)

messagef <- function(...) message(sprintf(...))
printf <- function(...) print(sprintf(...))
select <- dplyr::select

remove_doublets<- function(str){
  str <- as.character(str)
  str[str == lag(str)] <- ""
  str
}

get_arg_max <- function(x){
  d_x <- sign(diff(x))
  dd_x <- diff(d_x)
  which(dd_x < 0) + 1
}

get_gaussification_peaks <- function(gauss_data, with_plot = F, min_value = 0, troughs = F, output = c("time", "values", "time_values")){
  #browser()
  output <- match.arg(output)
  if(is.list(gauss_data) & "sum" %in% names(gauss_data)){
    gauss_data <- gauss_data[["sum"]]
  }
  if(troughs){
    gauss_data$val <- -gauss_data$val
  }
  peak_pos <- get_arg_max(gauss_data$val)  
  values <- gauss_data[peak_pos,]$val
  threshold <- 0
  if(is.numeric(min_value)){
    threshold <- min_value
  }
  else{
    if(tolower(min_value) == "median"){
      threshold <- median(values)
    }
    else if(tolower(min_value) == "mean"){
      threshold <- mean(values)
    }
  }
  #browser()
  if(troughs){
    filtered_pos <- intersect(peak_pos, which(gauss_data$val <= threshold))
  }
  else{
    filtered_pos <- intersect(peak_pos, which(gauss_data$val >= threshold))
  }
  t_max <- gauss_data[filtered_pos,]$t
  
  if(troughs){
    gauss_data$val <- -gauss_data$val 
  }

  if(output == "time"){
    ret <- t_max    
  }
  else if(output == "values"){
    ret <- gauss_data[filtered_pos,]$val
  }
  else if(output == "time_values"){
    ret <- gauss_data[filtered_pos,]
  }
  if(with_plot){
    q <- ggplot(gauss_data, aes(x = t, y = val)) + geom_line() + geom_vline(data = tibble(x = t_max), aes(xintercept = x))
    print(q)
  }
  ret
}

rect_func <- function(t, t0, half_width = .5){
  ifelse(abs(t - t0) <= half_width, 1.0, 0.0)
}

resample_bin_vec <- function(bin_vec, window_size = 3, shift = 1){
  l <- length(bin_vec)
  pad_bin_vec <- c(bin_vec, rep(0, ceiling(l/window_size) * window_size - l))
  window_starts <- seq(1, length(pad_bin_vec), shift)
  window_starts <- window_starts[window_starts <= (length(pad_bin_vec) - window_size + 1)]
  #values <- rep(1, window_size) * any(pad_bin_vec[window_starts])
  ret <- rep(0, length(bin_vec))
  for(i in window_starts){
    value <- as.integer(any(pad_bin_vec[i:(i + window_size - 1)] == 1))
    ret[i:(i + window_size - 1)] <- value
  }
  ret[1:l]
}

gaussification <- function(onsets, deltaT = 0.1, sigma = 5, weights = NULL, start = 0, end = 450, use_rect_func = FALSE, with_singles = F){
  if(is.null(start)){
    start <- min(onsets) - 2 * sigma
  }
  if(is.null(end)){
    end <- max(onsets) + 2*sigma
  }
  #messagef("start = %.2f, end = %.2f", start,  end)
  ret <- NULL
  if(is.null(weights)){
    weights <- rep(1, length(onsets))
  }
  time_points <- seq(start, end, deltaT)
  if(use_rect_func){
    singles <- purrr::map(seq_along(onsets), function(i){
      weights[i] * rect_func(time_points, onsets[i], half_width = sigma)
    })
    
  }
  else{
    sig_fact <- 1/sigma/sigma 
    singles <- purrr::map(seq_along(onsets), function(i){
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

get_conditional_prob <- function(time_points, ground_truth, forward_win = 2, backward_win = forward_win){
  if(is.data.frame(time_points)){
    time_points <- time_points$time_in_s
  }
  if(is.data.frame(ground_truth)){
    ground_truth <- ground_truth$time_in_s
  }
  total_counts <- length(time_points)
  map_dfr(1:length(ground_truth), function(i){
    #browser()
    counts <- time_points[time_points >= (ground_truth[i] - forward_win) & time_points <= (ground_truth[i] + backward_win)]
    tibble(i = i, time_in_s = ground_truth[i], counts = length(counts), rel_freq = counts/total_counts)
  })
}

get_cond_prob_baseline <- function(boundary_data, ground_truth, piece = 1:2, max_level = 1:3, size = 100, win_sizes = c(.5, 1, 1.5, 2, 2.5)){
  
  map_dfr(piece, function(pi){
    max_dur <- piece_durations[pi]
    
    map_dfr(max_level, function(ml){
      gt <- ground_truth %>% 
        filter(piece == pi, 
               level <= max_level)
      num_segs <- nrow(gt)
      map_dfr(win_sizes, function(ws){
        messagef("Checking piece = %s, max_level = %s, window size = %.2f, max_dur = %.2f, num_seqs = %d", pi, ml, ws, max_dur, num_segs)
        mu <- get_conditional_prob(boundaries_lab %>% 
                                     filter(piece == pi), 
                                   gt, 
                                   forward_win = ws)  %>% 
          pull(counts) %>% 
          mean() 
        baseline <- map_dbl(1:size, function(i){
          get_conditional_prob(boundaries_lab %>% 
                                 filter(piece == pi), 
                               sample(seq(0, max_dur, 1), num_segs, replace = F), 
                               forward_win = ws)  %>% 
            pull(counts) %>% 
            mean()})
        
        t.test(baseline, mu = mu) %>% 
          broom::tidy() %>% 
          mutate(win_size = ws, piece = pi, max_level = ml, mu = mu)
      })
      
    })
    
  })
}

get_sims <- function(time_points, 
                     ground_truth, 
                     start = 0,  end = 330, 
                     bw = 1, 
                     resample_window = 0, 
                     resample_shift = 3, 
                     with_overlap = T){
  browser()
  bin_v <- hist(time_points[time_points>= start & time_points <= end], breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  bin_gt <- hist(ground_truth[ground_truth >= start & ground_truth <= end], breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  if(with_overlap){
    # print(ground_truth)
    # messagef("start = %.2f / %.2f end = %.2f / %.2f, maxgt = %.2f", start, start - bw/2, end, end + bw/2, max(ground_truth[ground_truth >= start & ground_truth <= end]))
    bin_v_2 <- hist(time_points[time_points>= start & time_points <= end], 
                    breaks = seq(from = start - bw/2, to = end + bw, by = bw), 
                    plot = FALSE)$counts
    bin_gt_2 <- hist(ground_truth[ground_truth >= start & ground_truth <= end], 
                     breaks = seq(from = start - bw/2, to = end + bw, by = bw), 
                     plot = FALSE)$counts
    #messagef("Length difference: %d last1 = %d, last2 = %d", length(bin_v) - length(bin_v_2), bin_v[length(bin_v)], bin_v_2[length(bin_v_2)])
    l <- min(length(bin_v), length(bin_v_2))
    bin_v <- bin_v[1:l]  + bin_v_2[1:l]
    bin_gt <- bin_gt[1:l]  + bin_gt_2[1:l]
    
  }
  bin_v[bin_v > 1] <- 1
  bin_gt[bin_gt > 1] <- 1
  if(resample_window != 0){
    if(resample_shift > resample_window){
      resample_shift <- resample_window
    } 
    bin_v <- resample_bin_vec(bin_v, window_size = resample_window, resample_shift)
    bin_gt <- resample_bin_vec(bin_gt, window_size = resample_window, resample_shift)
  }
  #print(sprintf("%s-%s", bin_v, bin_gt)[sprintf("%s-%s", bin_v, bin_gt) != "0-0"])
  #print(table(bin_v, bin_gt))
  #browser()
  if(length(unique(bin_v)) == 1 || length(unique(bin_gt)) == 1){
    return(tibble(sim_f1 = 1))
  }
  #print(confusionMatrix(as.factor(bin_v), as.factor(bin_gt), mode = "everything", positive="1"))
  sim_f1 <- confusionMatrix(as.factor(bin_v), as.factor(bin_gt), mode = "everything", positive="1")$byClass[7] #function to return prediction vector based on bw
  if(is.na(sim_f1)){
    #browser()
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

get_pairwise_sims <- function(segments, start = 0, end = NULL, bw = .75, window_size = 0, window_shift = 3){
  ids <- unique(segments$p_id)
  pieces <- unique(segments$piece) %>% na.omit()
  
  map_dfr(pieces, function(pi){
    trials <-   segments %>% filter(piece == pi) %>% pull(trial) %>% unique()
    map_dfr(trials, function(tri){
      map_dfr(1:(length(ids)-1), function(i){
        onsets1 <- segments %>% filter(p_id == ids[i], piece == pi, trial == tri) %>% pull(time_in_s)
        map_dfr((i + 1):length(ids), function(j){
          #browser()
          if(i %% 10 == 0 & j %% 10 == 0){
            messagef("Calculating: piece = %s, trial = %s, id1 = %s, id2 = %s", pi, tri, ids[i], ids[j])
          }
          onsets2 <- segments %>% filter(p_id == ids[j], piece == pi) %>% pull(time_in_s)
          f1 <- get_sims(onsets1, onsets2, start = 0, end = piece_durations[pi], bw = bw, window_size, window_shift)
          bind_rows(tibble(id1 = ids[i], id2 = ids[j], f1= f1$sim_f1),
                    tibble(id1 = ids[j], id2 = ids[i], f1= f1$sim_f1))
        })
      }) %>% mutate(bw = bw, window_size = window_size, windowm_shift = window_shift, piece = pi, trial = tri)
      
    })
  }) %>% arrange(piece, id1, id2)
}

get_sims_between_trials <- function(segments, start = 0, end = NULL, bw = .75, window_size = 0, window_shift = 3){
  ids <- unique(segments$p_id)
  pieces <- unique(segments$piece) %>% na.omit()
  
  map_dfr(pieces, function(pi){
    trials <-   segments %>% filter(piece == pi) %>% pull(trial) %>% unique()
    if(length(trials) == 1){
      return(NULL)
    }
    map_dfr(1:length(ids), function(i){
      onsets1 <- segments %>% filter(p_id == ids[i], piece == pi, trial == 1) %>% pull(time_in_s)
      #browser()
      onsets2 <- segments %>% filter(p_id == ids[i], piece == pi, trial == 2) %>% pull(time_in_s)
      f1 <- get_sims(onsets1, onsets2, start = 0, end = piece_durations[pi], bw = bw, window_size, window_shift)
      bind_rows(tibble(id = ids[i], f1= f1$sim_f1))
    }) %>% mutate(bw = bw, window_size = window_size, windowm_shift = window_shift, piece = pi)
    
  }) %>% arrange(piece, id)
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

get_sims_by_gaussification <- function(onsets1, onsets2, sigma = 1, deltaT = .1, start = 0, end = 300){
  g1 <- gaussification(onsets1, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  g2 <- gaussification(onsets2, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  peaks1 <- get_gaussification_peaks(g1)
  peaks2 <- get_gaussification_peaks(g2)
  get_sims(peaks1, peaks2, start = start, end = end, bw = sigma)
}

get_gaussification_sd <- function(onsets, sigma = 1, deltaT = .1, start = 0, end = 300){
  #browser()
  g <- gaussification(onsets, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  peaks <- get_gaussification_peaks(g)
  g %>% filter(t %in% peaks) %>% pull(val) %>% sd(na.rm = T) 
}

export_gaussification_peaks <- function(boundaries, piece = 1, source = NULL, 
                                        sigma = 2, deltaT = .1, start = 0, end = 300){
  onsets <- boundaries %>% filter(piece == !!piece)
  if(!is.null(source)){
    onsets <- onsets %>% filter(source %in% !!source)
  }
  g <- gaussification(onsets$time_in_s, 
                      sigma = sigma, 
                      deltaT = deltaT, end = piece_durations[piece], 
                      use_rect_func = FALSE)
  peaks <- get_gaussification_peaks(g, output = "time_values")
  peaks
}

export_all_gaussification_peaks <- function(data = all_boundaries, sigma = 2, fname = "gauss_peaks.xlsx"){
  N_1 <- data %>% distinct(piece, trial, source, p_id)  %>% filter(piece == 1) %>% nrow()
  N_2 <- data %>% distinct(piece, trial, source, p_id)  %>% filter(piece == 2) %>% nrow()
  gauss_peaks <- bind_rows(export_gaussification_peaks(all_boundaries) %>% 
                             mutate(piece = piece_names[1],
                                    N = N_1), 
                           export_gaussification_peaks(all_boundaries, piece = 2) %>% 
                             mutate(piece = piece_names[2],
                                    N = N_2)) 
  #browser()
  gauss_peaks <- gauss_peaks %>% 
    rename(peak_height = val) %>% 
    group_by(piece) %>% 
    mutate(rel_peak_height = peak_height/N,
           mean_peak_height = mean(peak_height), 
           median_peak_height = median(peak_height), 
           q75 = quantile(peak_height)[4],
           max = max(peak_height)) %>% 
    ungroup()
  if(!is.null(fname) && nchar(fname) > 0 ){
    writexl::write_xlsx(gauss_peaks, fname)
  }
  gauss_peaks
}

get_gaussification_peakiness <- function(onsets, n_rater = 1, sigma = 1, deltaT = .1, start = 0, end = 300){
  #browser()
  g <- gaussification(onsets, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  if(is.null(g) || is.na(g)|| length(g) == 0 || nrow(g) == 0 || n_rater < 1){
    return(NA)
  }
  peaks <- get_gaussification_peaks(g, output = "values")
  troughs <- get_gaussification_peaks(g, troughs = T, output = "values")
  tibble(m_peaks = mean(peaks), m_roughs = mean(troughs), full = (mean(peaks) - mean(troughs))/n_rater, simple = m_peaks/n_rater, n_rater = n_rater)
}

get_peakiness_from_segmentation_ratings <- function(seg_data, piece = 1, trial = 1, source = NULL, sigma = 1){
  seg_data <- seg_data %>% filter(piece == !!piece, trial == !!trial)
  if(!is.null(source)){
    seg_data <- seg_data %>% filter(source == !!source)
  }
  get_gaussification_peakiness(seg_data$time_in_s, n_rater = n_distinct(seg_data$p_id), end = piece_durations[piece])
}
#' get_best_alignment: Calc best DTW alignment between to timelines and optional features derived from this 
#'
#' @param query <dbl> List of onsets
#' @param target <dbl> List of onsets 
#' @param summary <lgl> Returned summary?
#'
#' @return Either a tibble object with a query -> target mapping or summary features
#' @export
#'
#' @examples
get_best_alignment <- function(query, target){
  #get the DTW alignment
  align <- dtw(query, target) 
  ret <- c()
  
  #find the best alignment by selecting only the closest time points in case of multiple mappings
  ba <- 
    map_dfr(1:length(query), function(i){
      target_idz <- align$index2[which(align$index1 == i)]
      candidates <- target[target_idz]
      candidate_pos <- paste(which(align$index1 == i), collapse = ";")
      
      abs_diff <- abs(candidates - query[i])
      best_idx <- which.min(abs_diff)
      best <- min(abs_diff)
      target_best <- target_idz[best_idx ]
      
      tibble(query_pos = i, 
             target_pos = target_best, 
             candidate_pos = candidate_pos, 
             query = query[i], 
             target = target[target_best], 
             d = best[1])
    })
  
  # MAE mean absolute error; mean of absolute dIntermediate statewm-ifferences between query time points and best target timepoints
  # MAS mean absolute error stddev; Sd of absolute differences between query time points and best target timepoints
  # norm_dist: normaliized DTW distance
  # d_n: difference in length between target and query, if positive, query contains *more* events than target
  summary <- tibble(MAE = mean(ba$d),
                    MAS = sd(ba$d), 
                    MAX = max(ba$d),
                    norm_dist = align$normalizedDistance, 
                    d_n = length(query) - length(target))
  list(dtw = align, raw = ba, summary = summary)  
}

simulate_segmentation <- function(mean_log_isi = 1.5, sd_log_isi = .5,  n_seg = 100, max_t = 420){
  #messagef("Simulate segmentation, mean = %.2f, sd = %.2f, max = %.2f", mean_log_isi, sd_log_isi, max_t)
  #browser()
  isi <- rnorm(n_seg, mean_log_isi, sd_log_isi) %>% exp()
  #print(isi)
  t <- cumsum(c(0, isi))
  t[t < max_t]
}

simulate_segmentation_from_groundtruth <- function(ground_truth, 
                                                   size = 50,
                                                   piece = 1, 
                                                   theory = 1, 
                                                   max_level = 3){
  gt_log_ioi <- ground_truth %>% 
    filter(level <= max_level, 
           piece == !!piece, 
           theory == !!theory) %>% 
    pull(time_in_s) %>% 
    diff() %>% 
    log()
  map_dfr(1:size, function(i){
    tibble(time_in_s = simulate_segmentation(mean(gt_log_ioi), 
                        sd(gt_log_ioi), 
                        round(1.1 * length(gt_log_ioi)), 
                        piece_durations[piece]),
           p_id = i)
    })
}

simulate_ground_truth <- function(ground_truth, 
                                  piece = 1, 
                                  theory = 1, 
                                  max_level = 3){
  gt_log_ioi <- ground_truth %>% 
    filter(level <= max_level, 
           piece == !!piece, 
           theory == !!theory) %>% 
    pull(time_in_s) %>% 
    diff() %>% 
    log()
  
  ret <-  tibble(time_in_s = simulate_segmentation(mean(gt_log_ioi), 
                                                   sd(gt_log_ioi), 
                                                   round(1.5 * length(gt_log_ioi)), 
                                                   piece_durations[piece]))
  messagef("d_l: %.1f", 100*(length(gt_log_ioi) - nrow(ret))/length(gt_log_ioi))
  ret %>% mutate(boundary = 1, 
                 beginning = 0, 
                 ending = 0, 
                 both = 1, 
                 level = max_level,
                 piece = piece, 
                 theory = theory, 
                 boundary_type = "both")
  
}

simulate_segmentation_from_data <- function(segments = boundaries_lab, 
                                            piece = 1,
                                            trial = 1,
                                            size = length(unique(segments$p_id)),
                                            sigma = 1, 
                                            dT = .1,
                                            min_value = 0)
{
  #browser()
  # log_part_isi <- segments %>% 
  #   filter(piece == !!piece, trial == !!trial) %>% 
  #   pull(time_in_s) %>% 
  #   gaussification(end = piece_durations[piece] + 1, sigma = sigma) %>% 
  #   get_gaussification_peaks(with_plot = F, min_value = min_value) %>% 
  #   diff() %>% 
  #   log()
  log_part_isi <- segments %>% 
    filter(piece == !!piece, trial == !!trial) %>%
    group_by(p_id) %>% 
    summarise(m = mean(log(diff(time_in_s)), na.rm = T)) %>% 
    pull(m) 
  #browser()
  messagef("Simulating data with %d peaks for bw = %.2f, mean log ISI = %.2f, mean ISI = %.2f, sd log ISI = %.2f", 
           length(log_part_isi), sigma, mean(log_part_isi), exp(mean(log_part_isi)), sd(log_part_isi))
  if(length(log_part_isi) == 0){
    browser()
  }
  map_dfr(1:size, function(i){
    ret <- 
      tibble(time_in_s = simulate_segmentation(mean(log_part_isi), 
                                             sd(log_part_isi), 
                                             round(1.1 * length(log_part_isi)), 
                                             piece_durations[piece]),
           p_id = sprintf("S%02d", i), 
           piece = piece, 
           trial = trial)
    ret
  })
}

generate_fake_data <- function(segments = boundaries_lab, size = length(unique(segments$p_id)), bw = .75){
  pieces <- unique(segments$piece) %>% na.omit()
  map_dfr(pieces, function(pi){
    trials <- unique(segments[segments$piece == pi,]$trial) %>% na.omit()
    map_dfr(trials, function(tri){
     ret <- simulate_segmentation_from_data(segments, piece = pi, trial = tri, sigma = bw) %>% 
       filter(time_in_s > 0)
     #browser()
     ret
    })
  })
}

compare_segmentations <- function(segs = boundaries_lab, 
                                  gt = ground_truth, 
                                  piece = 1,
                                  theory = 1,
                                  max_level = 3,
                                  start = 0, 
                                  end = NULL,
                                  sigma = 1,
                                  threshold = 0, 
                                  with_plot = F,
                                  only_plot = F){
  browser()
  if(is.null(end)){
    end <- piece_durations[piece] + 1
  }
  part_boundaries <- segs %>% 
    filter(piece == !!piece) %>% 
    pull(time_in_s) %>% 
    gaussification(start = start, end = end, sigma = sigma) %>% 
    get_gaussification_peaks(with_plot = F, min_value = threshold)
  
  gt_boundaries <- gt %>% 
    filter(level <= max_level, 
           piece == !!piece, 
           theory == !!theory, 
           time_in_s >= start, 
           time_in_s <= end) %>% 
    pull(time_in_s)
  save(part_boundaries, gt_boundaries, file = "tmp.rda")
  
  best <- get_best_alignment(part_boundaries, gt_boundaries)
  sims <- get_sims_by_gaussification(part_boundaries, gt_boundaries, sigma = sigma, start = start, end = end)
  if(with_plot && !only_plot) print(plot_dtw_alignment(part_boundaries, gt_boundaries) + xlim(start, end))
  #browser()
  if(only_plot){
    return(plot_dtw_alignment(part_boundaries, gt_boundaries) + xlim(start, end))
  }
  else{
    best$summary <- bind_cols(best$summary, sims) 
    best
    
  }
}

get_all_dtw_distances_by_particpant <- function(segs = all_boundaries, 
                                             gt = ground_truth, 
                                             theory = 1,
                                             max_level = 2,
                                             start = 0, 
                                             end = NULL){
  gt_boundaries <- gt_boundaries <- gt %>% 
    filter(level <= max_level, 
           theory == !!theory)
           
  map_dfr(1:2, function(piece){
    if(is.null(end)){
      end <- piece_durations[piece] + 1
    }
    gt_boundaries <- gt %>% 
      filter(level <= max_level, 
             piece == !!piece, 
             theory == !!theory, 
             time_in_s >= start, 
             time_in_s <= end) %>% 
      pull(time_in_s)
    ids <- unique(segs$p_id)
    
    segs <- segs %>% 
      filter(piece == !!piece) 
    
    map_dfr(ids, function(id){
      p_data  <- segs %>% filter(p_id == id) 
      trials <- p_data %>% pull(trial) %>% unique()
      map_dfr(trials, function(tr){
        part_boundaries <- p_data %>% 
          filter(trial == tr) %>% 
          pull(time_in_s) 
        best <- get_best_alignment(part_boundaries, gt_boundaries)
        tibble(p_id = id, trial = tr, norm_dist = best$summary$norm_dist[1])  
      })
    }) %>% mutate(piece = piece, level = max_level, theory = theory)
  })
}

get_baseline <- function(boundary_data, ground_truth, size, piece, max_level, start, end, sigma, threshold, summary = T){
  se <- function(x, na.rm = FALSE) {
    if (na.rm) x <- na.omit(x)
    sqrt(var(x)/length(x))
  }
  ci95 <- function(x, type = c("low", "up")){
    type <- match.arg(type)
    b <- ifelse(type == "low", -1.96, 1.96)
    mean(x, na.rm = T) + b * se(x, na.rm = T)
  }
  ci95_low <- function(x) ci95(x, "low")
  ci95_up  <- function(x) ci95(x, "up")
  ci95_str <- function(x){
    sprintf("[%.2f, %.2f]", ci95_low(x), ci95_up(x))
  }
  ret <- 
    map_dfr(1:size, function(i){
      gt <- simulate_ground_truth(ground_truth, 
                                  piece = as.integer(piece),
                                  max_level = as.integer(max_level),
                                  theory = 1) %>% mutate(type = "simulation")
      compare_segmentations(segs = boundary_data, 
                            gt = gt, 
                            piece = as.integer(piece), 
                            theory = 1, 
                            max_level = as.numeric(max_level),
                            start = as.numeric(start),
                            end = as.numeric(end),
                            sigma = as.numeric(sigma), 
                            threshold = threshold,
                            with_plot = F) %>% pluck("summary")
      
  })
  sim_f1 <-     compare_segmentations(segs = boundary_data, 
                                      gt = ground_truth, 
                                      piece = as.integer(piece), 
                                      theory = 1, 
                                      max_level = as.numeric(max_level),
                                      start = as.numeric(start),
                                      end = as.numeric(end),
                                      sigma = as.numeric(sigma), 
                                      threshold = threshold,
                                      with_plot = F) %>% 
    pluck("summary") 
  
  
  #if(max_level == 2)browser()
  if(summary){
    tt1 <- t.test(ret$sim_f1, mu = sim_f1$sim_f1) %>% 
      broom::tidy() %>% 
      select(d_f1 = estimate, t_f1 = "statistic",  df_f1 = parameter, p_f1 = p.value) %>% 
      mutate(d_f1 = d_f1 - sim_f1$sim_f1)
    
    tt2 <- t.test(ret$norm_dist, mu = sim_f1$norm_dist) %>% 
      broom::tidy() %>% 
      select(d_nd = estimate, t_nd = "statistic",  df_nd = parameter, p_nd = p.value) %>% 
      mutate(d_nd = d_nd - sim_f1$norm_dist)
    ret <- ret %>% 
      summarise(across(c("MAE","MAS","MAX", "d_n"), 
                       list("mean" = mean)
      ), 
      across(c("norm_dist", "sim_f1"),
             list("mean" = mean,  "ci95" = ci95_str)
      )) %>% bind_cols(tt1, tt2)  
    
  }
  ret
}

get_segmentation_stats <- function(boundary_data = all_boundaries, 
                                   ground_truth,
                                   band_widths = seq(0.75, 3.75, .5),
                                   baseline_size = 25){
  pieces <- unique(boundary_data$piece)
  start <- 0 
  #thresholds <- c("0", "Mean")
  thresholds <- c("mean")
  map_dfr(pieces, function(pi){
    end <- piece_durations[pi]
    gt <- ground_truth %>% filter(piece == pi)
    max_level <- max(gt$level)
    bd <- boundary_data %>% filter(piece == pi)
    #trials <- c("both", unique(bd$trial))
    trials <- "both"
    map_dfr(band_widths, function(bw){
      map_dfr(thresholds, function(thr){
        map_dfr(trials, function(tri){
          map_dfr(1:max_level, function(ml){
            #browser()
            if(tri != "both"){
              bd_loc <- bd %>% filter(trial == tri)
            }
            else{
              bd_loc <- bd
            }
            messagef("Checking bw = %.2f [size = %d], piece = %s, trial = %s, max_level = %d, threshold = %s", bw, baseline_size, pi, tri, ml, thr)
            bl <- get_baseline(bd_loc,
                               gt, 
                               baseline_size, 
                               pi, 
                               max_level = ml, 
                               start, end, 
                               sigma = bw, 
                               threshold = thr, 
                               summary = T)
            f1 <- compare_segmentations(bd_loc, 
                                        gt = gt, 
                                        piece = pi, 
                                        theory = 1, 
                                        max_level = ml,
                                        start = start,
                                        end = end,
                                        sigma = bw, 
                                        threshold = thr,
                                        with_plot = F) %>% pluck("summary")
            tibble(piece = pi, 
                   trial = tri,
                   threshold = thr,
                   theory = 1,
                   max_level = ml,
                   sigma = bw) %>% bind_cols(bl, f1)  
          })
        })
      })
    })
  })
  

}

get_boundary_stats <- function(){
  tmp <- 
    all_boundaries %>% 
    mutate(level = 1, 
           theory = 1,
           boundary_type = "both",
           trial_id = sprintf("%s_%s_%s_%s", p_id, piece, trial, level)) %>% 
    bind_rows(bind_rows(
      ground_truth %>% filter(level == 1),
      ground_truth %>% filter(level <= 2) %>% mutate(level = 2)
    )%>% 
      mutate(source = "theory", 
             trial = 1, 
             p_id = "T01", 
             trial_id = sprintf("%s_%s_%s_%s", p_id, piece, trial, level))) %>% 
    ungroup() %>% 
    group_by(trial_id) %>% 
    mutate(isi = c(diff(time_in_s), NA),
           m = mean(isi, na.rm = T), 
           s = sd(isi, na.rm = T),
           log_m = log(m),
           log_sd = log(s), .groups = "drop",
           label = 1:n())
  tmp
}

test_peakiness_values <- function(seg_data = all_boundaries, size = 100, sigma = 2, type = "simple"){
  pieces <- unique(seg_data$piece)
  map_dfr(pieces, function(p){
    trials <- unique(seg_data[seg_data$piece == p,]$trial)
    map_dfr(trials, function(tr){
      #browser()
      n_rater <- seg_data %>% filter(piece == p, trial == tr) %>% pull(p_id) %>% unique() %>% length()
      map_dfr(sigma, function(s){
        mu <- get_peakiness_from_segmentation_ratings(seg_data, piece = p, trial = tr, sigma = s) %>% pull(!!type)
        sim <- map_dbl(1:size, function(x){
          tmp <- simulate_segmentation_from_data(seg_data, piece = p, trial = tr)
          get_peakiness_from_segmentation_ratings(tmp, piece = p, trial = tr, sigma = s) %>% pull(!!type)
        })
        #browser()
        t.test(sim, mu = mu) %>% 
          broom::tidy() %>% 
          mutate(piece = p, trial = tr, sigma = s, size = size, mu = mu, n_rater = n_rater)
        
      })
   })
  })
}

find_annotated_peaks <- function(data = all_boundaries, annotations = boundaries_lab_annotations, window_size = 2, summary = T){
  all_peaks <- export_all_gaussification_peaks(data, fname = NULL)
    
  lab2 <- annotations %>% select(p_id, time_in_s, piece, starts_with("SEG"))
  ret <- 
    map_dfr(1:2, function(p){
    peak_positions <- all_peaks %>% filter(piece == piece_names[p]) 
    map2_dfr(peak_positions$t, peak_positions$rel_peak_height, function(t, val){
      #browser()
      lab2 %>% 
        filter(piece == p, abs(time_in_s - t) < window_size) %>% 
        mutate(peak_pos = t, 
               rel_peak_height = val)
    })
  }) %>% na.omit()
  
  if(summary){
    ret <- ret %>% 
      group_by(piece, peak_pos) %>% 
      summarise(beginning = mean(SEG.beginning, na.rm = T), 
                ending = mean(SEG.ending, na.rm = T), 
                importance = mean(SEG.importance), 
                rel_peak_height = mean(rel_peak_height), .groups = "drop") 
  }
  ret
}