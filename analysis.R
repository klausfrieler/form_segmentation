library(tidyverse)
library(mccr)
library(cvAUC)
library(MASS)
library(caret)
library(proxy)
library(psych)
library(dtw)

#library(MSM)

messagef <- function(...) message(sprintf(...))
printf <- function(...) print(sprintf(...))
select <- dplyr::select

get_arg_max <- function(x){
  d_x <- sign(diff(x))
  dd_x <- diff(d_x)
  which(dd_x < 0) + 1
}

get_gaussification_peaks <- function(gauss_data, with_plot = F, min_value = 0){
  if(is.list(gauss_data) & "sum" %in% names(gauss_data)){
    gauss_data <- gauss_data[["sum"]]
  }
  peak_pos <- get_arg_max(gauss_data$val)  
  values <- gauss_data[peak_pos,]$val
  threshold <- 0
  if(is.numeric(min_value)){
    threshold <- min_value
  }
  else{
    if(min_value == "median"){
      threshold <- median(values)
    }
    else if(min_value == "mean"){
      threshold <- mean(values)
    }
  }
  #browser()
  filtered_pos <- intersect(peak_pos, which(gauss_data$val >= threshold))
  t_max <- gauss_data[filtered_pos,]$t
  if(with_plot){
    q <- ggplot(gauss_data, aes(x = t, y = val)) + geom_line() + geom_vline(data = tibble(x = t_max), aes(xintercept = x))
    print(q)
  }
  t_max
}

rect_func <- function(t, t0, half_width = .5){
  ifelse(abs(t - t0) <= half_width, 1.0, 0.0)
}

resample_bin_vec <- function(bin_vec, window_size = 3, overlap = 1){
  l <- length(bin_vec)
  pad_bin_vec <- c(bin_vec, rep(0, ceiling(l/window_size) * window_size - l))
  window_starts <- seq(1, length(pad_bin_vec), overlap)
  window_starts <- window_starts[window_starts <= (length(pad_bin_vec) - window_size + 1)]
  #values <- rep(1, window_size) * any(pad_bin_vec[window_starts])
  ret <- rep(0, length(bin_vec))
  for(i in window_starts){
    value <- as.integer(any(pad_bin_vec[i:(i + window_size - 1)] == 1))
    ret[i:(i + window_size - 1)] <- value
  }
  ret[1:l]
}

gaussification <- function(onsets, deltaT = 0.1, sigma = 5, weights = NULL, start = -1, end = 450, use_rect_func = FALSE, with_singles = F){
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

get_sims <- function(time_points, ground_truth, start = 0, end = 330, bw = 1, resample_window = 3, resample_overlap = 3){
  bin_v <- hist(time_points, breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  bin_gt <- hist(ground_truth, breaks = seq(from = start, to = end + bw, by = bw), plot = FALSE)$counts
  
  bin_v[bin_v > 1] <- 1
  bin_gt[bin_gt > 1] <- 1
  if(resample_window != 0){
    bin_v <- resample_bin_vec(bin_v, window_size = resample_window, resample_overlap)
    bin_gt <- resample_bin_vec(bin_gt, window_size = resample_window, resample_overlap)
  }
  browser()
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

get_sims_by_gaussification <- function(onsets1, onsets2, sigma = 1, deltaT = .1, start = 0, end = 300){
  browser()
  g1 <- gaussification(onsets1, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  g2 <- gaussification(onsets2, sigma = sigma, deltaT = deltaT, end = end, use_rect_func = FALSE)
  peaks1 <- get_gaussification_peaks(g1)
  peaks2 <- get_gaussification_peaks(g2)
  get_sims(peaks1, peaks2, start = 0, end = end, bw = sigma)
  
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
  isi <- rnorm(n_seg, mean_log_isi, sd_log_isi) %>% exp()
  t <- cumsum(c(0, isi))
  t[t < max_t]
}

simulate_segmentation_from_groundtruth <- function(ground_truth = ground_truth, 
                                                   size = 50,
                                                   piece = 1, 
                                                   theory = 1, 
                                                   max_level = 3){
  gt_log_ioi <- ground_truth %>% 
    filter(level <= max_level, 
           piece == !!piece, 
           theory == !!theory, 
           boundary_type != "ending") %>% pull(time_in_s) %>% diff() %>% log()
  map_dfr(1:size, function(i){
    tibble(time_in_s = simulate_segmentation(mean(gt_log_ioi), 
                        sd(log(gt_log_ioi)), 
                        round(1.1 * length(gt_log_ioi)), 
                        piece_durations[piece]),
           p_id = i)
    })
}

simulate_segmentation_from_data <- function(segs = boundaries_lab, 
                                            piece = 1,
                                            size = 50,
                                            sigma = 1, 
                                            dT = .1,
                                            min_value = 0)
{
  log_part_isi <- segs %>% 
    filter(piece == !!piece) %>% 
    pull(time_in_s) %>% 
    gaussification(end = piece_durations[piece] + 1, sigma = sigma) %>% 
    get_gaussification_peaks(with_plot = F, min_value = min_value) %>% 
    diff() %>% 
    log()
  #browser()
  map_dfr(1:size, function(i){
    tibble(time_in_s = simulate_segmentation(mean(log_part_isi), 
                                             sd(log_part_isi), 
                                             round(1.1 * length(log_part_isi)), 
                                             piece_durations[piece]),
           p_id = i)
  })
}

compare_segmentations <- function(segs = boundaries_lab, 
                                  gt = ground_truth, 
                                  piece = 1,
                                  theory = 1,
                                  max_level = 3, 
                                  sigma = 1,
                                  threshold = 0, 
                                  with_plot = T){
  part_boundaries <- segs %>% 
    filter(piece == !!piece) %>% 
    pull(time_in_s) %>% 
    gaussification(end = piece_durations[piece] + 1, sigma = sigma) %>% 
    get_gaussification_peaks(with_plot = F, min_value = threshold)
  print(length(part_boundaries))
  gt_boundaries <- ground_truth %>% 
    filter(level <= max_level, 
           piece == !!piece, 
           theory == !!theory, 
           boundary_type != "ending") %>% 
    pull(time_in_s)
  best <- get_best_alignment(part_boundaries, gt_boundaries)
  sims <- get_sims_by_gaussification(part_boundaries, gt_boundaries, sigma = sigma)
  if(with_plot) print(plot_dtw_alignment(part_boundaries, gt_boundaries))
  browser()
  best$summary <- bind_cols(best$summary, sims) 
  best
}
