
#packages

library(tidyverse)
library(mccr)



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

#example usage  

temp2 <- densX(df_p1_t1, 10)






#creating matching function that takes a list of vectors and a theoretical vector and returns a DATA FRAME of matching scores

matchFuncDF <- function(vecs, x){
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
