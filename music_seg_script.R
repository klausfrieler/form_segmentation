
#packages

library(tidyverse)
library(mccr)
install.packages("cvAUC")
library(cvAUC)
install.packages('caret')
library(caret)
install.packages('lsa')
library(TSdist)
library(lsa)


# read data


df_all <- read.csv("dataLONG.csv")

################### Clean data and create DF for different experimental conditions #################

# trying to do all the data cleaning with pipes

# piece 1 trial 1
df_p1_t1 <- subset(df_all, Stimulus=="part2_01.avi") %>% 
  filter(Event=="BOUNDARY_REGISTERED") %>% 
  subset(Trial==1) %>% 
  select(ID, FrameNo, Time_in_s)

# piece 1 trial 2
df_p1_t2 <- subset(df_all, Stimulus=="part2_01.avi") %>% 
  filter(Event=="BOUNDARY_REGISTERED") %>% 
  subset(Trial==2) %>% 
  select(ID, FrameNo, Time_in_s)

# piece 2 trial 1
df_p2_t1 <- subset(df_all, Stimulus=="part2_02.avi") %>% 
  filter(Event=="BOUNDARY_REGISTERED") %>% 
  subset(Trial==3) %>% 
  select(ID, FrameNo, Time_in_s)

# piece 2 trial 2
df_p2_t2 <- subset(df_all, Stimulus=="part2_02.avi") %>% 
  filter(Event=="BOUNDARY_REGISTERED") %>% 
  subset(Trial==4) %>% 
  select(ID, FrameNo, Time_in_s)


#read in file with theory predictions and levels 
df_qual <- read.csv("boundaries_complete.csv", sep=",")
names(df_qual)[1] <- "Time"


##################### Latest code - functions that creae results table ##############################

#some data prep that I can intergrate into functions

#make data into list

df_p1_t1_list = split(df_p1_t1, f = df_p1_t1$ID)

#change df into list of vectors of time
obs_list <-    lapply(df_p1_t1_list, "[", , 3)

vect_gt <- df_qual$Time


####### Master function#####

sum_func <- function(vect1, vect_gt, sim, bw){
  if (sim == 1) { #AUC
    #binning vectors
    vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    
    match <- AUC(vect1, vect_gt) #function to return prediction vector based on bw
    
  } else if (sim == 2) { # Matthew correlation
    #binning vectors
    vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    
    match <- mccr(vect1, vect_gt) #function to return prediction vector based on bw
    
    
    #  } else if (sim == 3) { #F1 score
    #binning vectors
    #  vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    # vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    #creating factors out of biined vectors
    # vect1 <- as.factor(vect1)
    #vect_gt <-as.factor(vect_gt)
    
    
    #match <- confusionMatrix(vect1, vect_gt, mode = "everything", positive="1")$byClass[7] #function to return prediction vector based on bw
    #getting just the value not working. Ask Klaus
    
  } else if (sim==4) {# Euclidean dist
    
    vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    
    match <- EuclideanDistance(vect1, vect_gt) #function to return prediction vector based on bw
    
  } else if (sim==5) { #cosine similarity
    #binning vectors
    vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    
    match <- cosine(vect1, vect_gt) #function to return prediction vector based on bw
    
  } else if (sim==6) { # correlation
    #binning vectors
    vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
    
    match <- cor(vect1, vect_gt) #function to return prediction vector based on bw
    
  }
  
  return(match)
}


# wrapper function 

#function taking test and GT vectors, testy type, and bin width - returning results df

lapp_df <- function (test, list, vect_gt, bw){
  res_df <- lapply(list, sum_func, vect_gt=vect_gt, sim=test, bw=bw)
  res_df <- do.call(rbind.data.frame, res_df)
  
  res_df[1:43,2] <- test
  res_df[1:43,3] <- bw
  res_df[1:43,4] <- "main"
  res_df[1:43,5] <- 1:43
  
  
  
  colnames(res_df) <- c("result", "test", "bw", "gt", "ID")
  
  return(res_df)
}



#WORKING!
test_lapp <- lapp_df(1, obs_list, vect_gt, 0.5)


# testing multiple tests function
sappl_df <- function (list, vect_gt, bw) {
  res_all_tests_df <- lapp_df(1, list, vect_gt, bw) %>% 
    rbind(lapp_df(2, list, vect_gt, bw)) %>% 
    rbind(lapp_df(4, list, vect_gt, bw)) %>% 
    rbind(lapp_df(5, list, vect_gt, bw)) %>% 
    rbind(lapp_df(6, list, vect_gt, bw)) %>% 
    
    
    return(res_all_tests_df)
}

#test - working!!! 
test_lapp <- sappl_df(obs_list, vect_gt, 0.5)

# multiple bw function 

bw_function <- function (list, vect_gt) {
  res_all_tests_bw_df <- sappl_df(list, vect_gt, 0.1) %>%
    rbind(sappl_df(list, vect_gt, 0.3)) %>%
    rbind(sappl_df(list, vect_gt, 0.5)) %>%
    rbind(sappl_df(list, vect_gt, 0.7))
  
  return(res_all_tests_bw_df)
}


#test - working! 
test_bw <- bw_function(obs_list, vect_gt)



################################################################################




########################################################################################################


#################### Old redundant code ###########################
#select only piece 1 and create data frame
df_piece1 <- subset(df_all, Stimulus=="part2_01.avi")
#select only piece 2 and create data frame
df_piece2 <- subset(df_all, Stimulus=="part2_02.avi")


#filter out irrelevant rows 

df_piece1 <- filter(df_piece1, Event=="BOUNDARY_REGISTERED")
df_piece2 <- filter(df_piece2, Event=="BOUNDARY_REGISTERED")


# make Trial variable a factor which you need when you are plotting both trials on one plot
df_piece1$Trial <- factor(df_piece1$Trial)
df_piece1$Trial <- factor(df_piece1$Trial)

# data per trial 
df_p1_t1 <- subset(df_piece1, Trial==1)
df_p1_t2 <- subset(df_piece1, Trial==2)



#keeping only relevant vairables
df_p1_t1 <- select(df_p1_t1, ID, FrameNo, Time_in_s)

#####################################################################################


# create prediction timeline - for now just the broadest theory with all the boundaries

pred_th1 <- list(hist(df_qual$Time, breaks = seq(from=0, to=330, by=5), plot= FALSE)$counts)

# function that takes in a data set and a binning value and returns a prediction list



#########################OLD CODE BEYOND HERE#######################################



####################### Functions ##################################################

# Function to bin theory time series 

# function that takes in a data frame and a binning value and returns a prediction vector

# NOTE: i will update this function to return all different theories vectors in a list
predX <- function(df, x){
  predVect <- hist(df$Time, breaks = seq(from=0, to=330, by=x), plot= FALSE)$counts
  return(predVect)
}

test <- predX(df_qual, 10)



#function that returns list of  histogram based vectors from a df and binning value

# QUESTION FOR KLAUS: WHAT ABOUT INSTANCES WHEN A BIN HAS MORE THAN ONE REPORT?

histX <- function(df, x){
  vectors <- list()
  for(i in 1:length(unique(df$ID))) {
    
    vectors[i] <- list(hist(subset(df, ID==i)$Time_in_s, breaks = seq(from=0, to=330, by=x), plot = FALSE)$counts)
  }
  return(vectors)
}

#example usage  

temp1 <- histX (df_p1_t1, 10)



#function that returns list of density vectors based on df and bin width

#QUESTION FOR KLAUS: HISTOGRAM OUTPUT HAS DENSITY ESTIMATES, 
#BUT AM GUESSING WE WANT TO WRITE A SEPARATE FUNCTION USING DENSITY FUNCTION. iS THAT RIGHT?
#i had a go with some of the output from the density function, but i don't think i am doing it right

densX <- function(df, x){
  vectors <- list()
  for(i in 1:length(unique(df$ID))) {
    vectors[i] <- list(density(subset(df, ID==i)$Time_in_s, bw=x)$x)
  }
  return(vectors)
}


temp1 <- density(df_p1_t1$Time_in_s)


#example usage  

temp2 <- densX(df_p1_t1, 10)








#matching function that takes a list of vectors and a theoretical vector and returns a DATA FRAME of matching scores

match_Func_DF <- function(vecs, x){
  vect_df <- data.frame()
  for(i in 1:length(vecs)) {
    vect_df[i,1] <- i
    vect_df[i,2] <- "math"
    vect_df[i,3] <- 330/unlist(length(unlist(vecs[1])))
    vect_df[i,4] <- mccr(x, unlist(vecs[i]))  
  }
  names(vect_df) <- c('ID', "match_func",'bin_width',  "match_score")
  return(vect_df)
}

test5 <- matchFuncDF(temp1, test)
mccr(test, vect1)

#matching function that takes a list of vectors and a theoretical vector and returns a list of matching scores
# THE PREVIOUS FUNCTION MAKES THIS ONE REDUNDANT, BUT KEEPING IT HERE IN CASE IT IS USED LATER

matchFunc <- function(vecs, x){
  matchScr <- list()
  for(i in 1:length(vecs)) {
    matchScr[i] <- mccr(x, unlist(vecs[i]))  
  }
  return(matchScr)
}

test2 <- matchFunc(temp1, test)

######## variables created for testing

vectObs <- subset(df_p1_t1, ID==1) %>% 
  select(Time_in_s)

vectObs <-  vectObs$Time_in_s

vect_gt <- df_qual$Time

vect2 <- hist(vectObs, breaks = seq(from=0, to=330, by=10), plot = FALSE)$counts
vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=10), plot = FALSE)$counts

mccr(vect2, vect_gt)
test <- as.factor(test)
vect1 <- as.factor(vect1)
res <- confusionMatrix(actual, pred, mode = "everything", positive="1")
res2 <- res$byClass:F1

res2 <- res$byClass[7]

actual <- factor(rep(c(1, 0), times=c(160, 240)))
pred <- factor(rep(c(1, 0, 1, 0), times=c(120, 40, 70, 170)))

####################################



sumFunc(vectObs, vectGT, 6, 2)
sumFunc(vect1, vectGT, 1, 1)

sumFunc(vect1, vectGT, 2, 10)

