library(tidyverse)

#' rect_func, 
#' @param t: Vector of time points
#' @param t0: Center of rectangle
#' @param half_width: Half width of rectangle
rect_func <- function(t, t0, half_width = .5){
  ifelse(abs(t - t0) <= half_width, 1.0, 0.0)
}

#' triangle_func, 
#' @param t: Vector of time points
#' @param t0: Center of rectangle
#' @param half_width: Half width of rectangle
triangle_func <- function(t, t0, half_width = .5){
  tri <- function(t, t0, hw){
    ifelse(t <=  t0, (t-t0)/hw + 1, (t0-t)/hw + 1 )
    
  }
  ifelse(abs(t - t0) <= half_width, tri(t, t0, half_width), 0)
}

#' gaussification
#' 
#' Generates a "gaussifcation" from  onsets. All times should be measured in the same units (e.g., seconds, milliseconds)
#' @param onsets (vector numeric) Vector of onsets.
#' @param deltTa (scalar numeric) Time steps to use for the gaussification
#' @param sigma (scalar numeric) Standard deviation for the Gaussians
#' @param weights (vector numeric) Vector of the same length as onsets, can be used to weight onsets, default to NULL, which means equal weights.
#' @param start (scalar numeric) Start point of gaussification, if NULL set to a 2*sigma before first onsets
#' @param end (scalar numeric) End point of gaussification, if NULL set to a 2*sigma after last onsets
#' @param use_rect_func (scalar boolean or a function ) If not TRUE, a rectular window function with half_width equal to sigma is used instead of Gaussians. 
#' Can also be any other function f with signature f(time_points = seq(start, end, deltaT), onset = onset[i], half_width = sigma)
#' @param with_single (boolean) Should all single gaussificatio be returned?
#' 
#' @return Either a tibble with columns t and val, containing the sum of all Gaussians, or a list with the sum as well as every single gaussification. 

gaussification <- function(onsets, deltaT = 0.1, sigma = 5, weights = NULL, start = 0, end = 450, use_rect_func = FALSE, with_singles = F){
  if(is.null(start)){
    start <- min(onsets) - 2 * sigma
  }
  if(is.null(end)){
    end <- max(onsets) + 2*sigma
  }
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
      weights[i] * exp( -.5 * (onsets[i] - time_points) * (onsets[i] - time_points) * sig_fact)
    })
    
  }
  val <- purrr::reduce(singles, function(x,y) x + y)
  ret <- tibble(t = time_points, val = val)
  if(with_singles){
    dfs <- map_dfr(seq_along(singles), function(x) tibble(t = time_points, val = singles[[x]], onset = onsets[[x]]))
    ret <- list(sum = ret, singles = singles)
  }
  ret
}

#' get_arg_max: Get all maximum positions in a vector
#' @param x: numeric vector
#' 
#' @return indizes in vector x where local maxima are located
#' 
get_arg_max <- function(x){
  d_x <- sign(diff(x))
  dd_x <- diff(d_x)
  which(dd_x < 0) + 1
}

#' get_gaussification_peaks
#' 
#' @param gauss_data: A gaussification tibble as return by gaussification, if a list with singles only the sum is used.
#' @param with_plot  (boolean) Do you want a plot? 
#' @param min_value  (scalar numeric or  "mean"  "median") Defines a threshold for peaks to pick. The peak should be higher
#' as the numerical value of min_value or as mean/median of all peaks found.
#' @param troughs  (boolean) Return minima (troughs) instead of maxima? ( = maxima of the negative gaussfication)
#' @param output   (one of "time", "values", or "time_values"). Should we return either time points or peak values or both?
#' 
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

#' plot_gaussification: plot gaussification 
#' 
#' @param data: data frame containing time points in a column 'marker' or a column 'time_in_s'
#' @param  sigma (scalar numeric): Standard deviation for the gaussification
#' @param deltaT (scalar numeric): Time step  for the gaussification
#' @param start, end (scalar numeric): where to start and end the plot
#' @param threshold (scalar numeric or "median", "mean", "q75") Adds vertical lines for each peak above threshold
#' @param with_markers (boolean scalar) if TRUE adds marker for all peaks.
#' @param only_arkers (booelan scalar) if TRUE only markers  without gaussification will be plotted
#' @external_markers (data frame) If not NULL, it should be data frame with a column "time_in_s" containing marker positions to be added. Peak related markers take precedence if threshold is not NULL, i.e., external markers will be ignored.
#' @param alpha (scalar numeric between 0 and 1): Alpha value for markers, otherwise alpha will be set according to number of markers for better displays.
#' 
#' @return A ggplot  
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
  q <-  ggplot()
  if(only_markers || !is.null(external_markers) || !is.null(threshold)){
    with_markers <- T
  }
  if(!only_markers){
    q <- q + geom_line(data = combined_markers, aes(x = t, y = val), color = "indianred") 
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
  q <- q + labs(x = "Time (s)", y = "Gaussification Values")
  q <- q + xlim(start, end)
  #browser()
  q  
}
