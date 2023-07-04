
#packages

library(tidyverse)
library(mccr)
library(cvAUC)
library(MASS)
library(caret)
library(TSdist)
library(lsa)
library(ggridges)
library(psych)
library(gtools)
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
df_qual <- read.csv("boundaries_complete2.csv", sep=",")
names(df_qual)[1] <- "Time"

# recoding starts and ends into one variable

df_qual$start_end[df_qual$beginning==1]=1
df_qual$start_end[df_qual$ending==1]=2
df_qual$start_end[df_qual$both==1]=3

df_qual <- select(df_qual, Time, boundary_all, boundary_1, start_end, level)

#selecting only level 1 and 2 boundaries
df_qual <- subset(df_qual, level < 3)


##################### Latest code - functions that creae results table ##############################

#some data prep that I can integrate into functions

#make data into list

df_p1_t1_list = split(df_p1_t1, f = df_p1_t1$ID)

#change df into list of vectors of time
obs_list <-    lapply(df_p1_t1_list, "[", , 3)

vect_gt <- df_qual$Time


####### Master function#####

# 
#   if (sim[1] == "AUC") { #AUC
#     #binning vectors
#     vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
# 
#     match <- AUC(vect1, vect_gt) #function to return prediction vector based on bw
# 
#   } else if (sim == 2) { # Matthew correlation
#     #binning vectors
#     vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
# 
#     match <- mccr(vect1, vect_gt) #function to return prediction vector based on bw
# 
# 
#     #  } else if (sim == 3) { #F1 score
#     #binning vectors
#     #  vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     # vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     #creating factors out of biined vectors
#     # vect1 <- as.factor(vect1)
#     #vect_gt <-as.factor(vect_gt)
# 
# 
#     #match <- confusionMatrix(vect1, vect_gt, mode = "everything", positive="1")$byClass[7] #function to return prediction vector based on bw
#     #getting just the value not working. Ask Klaus
# 
#   } else if (sim==4) {# Euclidean dist
# 
#     vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
# 
#     match <- EuclideanDistance(vect1, vect_gt) #function to return prediction vector based on bw
# 
#   } else if (sim==5) { #cosine similarity
#     #binning vectors
#     vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
# 
#     match <- cosine(vect1, vect_gt) #function to return prediction vector based on bw
# 
#   } else if (sim==6) { # correlation
#     #binning vectors
#     vect1 <- hist(vect1, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
#     vect_gt <- hist(vect_gt, breaks = seq(from=0, to=330, by=bw), plot = FALSE)$counts
# 
#     match <- cor(vect1, vect_gt) #function to return prediction vector based on bw
# 
#   }
# 
#   return(match)




get_sims <- function(vect1, vect_gt, bw){
  match = NA
  vect1 <- hist(vect1, breaks = seq(from = 0, to = 330, by = bw), plot = FALSE)$counts
  vect_gt <- hist(vect_gt, breaks = seq(from = 0, to = 330, by = bw), plot = FALSE)$counts

  vect1[vect1 > 1] <- 1
  
  vect_gt[vect_gt > 1] <- 1
  

  sim_f1 <- confusionMatrix(as.factor(vect1), as.factor(vect_gt), mode = "everything", positive="1")$byClass[7] #function to return prediction vector based on bw
  if(is.na(sim_f1)){
    sim_f1 <- 0
  }  

  sim_AUC <- auc(vect1, vect_gt) #function to return prediction vector based on bw
  sim_mccr <- mccr(vect1, vect_gt) #function to return prediction vector based on bw
  sim_euc_dist <- EuclideanDistance(vect1, vect_gt) #function to return prediction vector based on bw
  sim_cosine <- cosine(vect1, vect_gt) %>% as.numeric() #function to return prediction vector based on bw
  sim_corr <- cor(vect1, vect_gt) #function to return prediction vector based on bw
  
  return(tibble(sim_AUC=sim_AUC, sim_mccr= sim_mccr, sim_euc_dist = sim_euc_dist, sim_cosine=sim_cosine, sim_corr = sim_corr, sim_f1=sim_f1))
  
}





#function taking test and GT vectors, testy type, and bin width - returning results df
# functional programming - everything is a function and you apply them to multiple vars
get_all_sims <- function(segment_data, vect_gt, bw_range = seq(.5, 1.5, .25)){
  p_ids <- unique(segment_data)$p_id
  map_dfr(p_ids, function(id){ #works like lapply; function is anonymous functions also called lamda
    map_dfr(bw_range, function(bw){ # like a nested loop 
      get_sims(segment_data[segment_data$p_id ==id,], vect_gt, bw) %>%  # 
        mutate(bw = bw, id = id) # adds columns to data - also able to do operations on columns
    })
  })
} 


# testing f1 errors

vect1 <-  hist(obs_list[[2]][[1]], breaks = seq(from=0, to=330, by=.5), plot = FALSE)$counts
vect2 <-hist(vect_gt, breaks = seq(from=0, to=330, by=.5), plot = FALSE)$counts

vect1[vect1 > 1] <- 1

vect2[vect2 > 1] <- 1

sim_f1 <- confusionMatrix(as.factor(vect1), as.factor(vect2), mode = "everything", positive="1")$byClass[7]

test_res <- confusionMatrix(as.factor(vect1), as.factor(vect2), mode = "everything", positive="1")


################# HACK JOB ######################################

# after a lot of attempts to do anything neat, here comes the hack job

# creating GT vectors of interest
#broad theory
vect_gt <- df_qual$Time

# levels 
vect_gt_lv1 <- subset(df_qual, level==1)$Time
vect_gt_lv2 <- subset(df_qual, level==2)$Time

#beginnings, endings, and both
vect_gt_b <- subset(df_qual, beginning==1)$Time
vect_gt_e <- subset(df_qual, ending==1)$Time
vect_gt_bth <- subset(df_qual, both==1)$Time

#narrow theory
vect_gt_2 <- subset(df_qual, boundary_1==1)$Time


# levels 
vect_gt__2_lv1 <- subset(df_qual, boundary_1==1 & level==1)$Time
vect_gt__2_lv2 <- subset(df_qual, boundary_1==1 & level==2)$Time


#beginnings, endings, and both
vect_gt_2_b <- subset(df_qual, boundary_1==1 & beginning==1)$Time #same as theory 1 and therefore not something that can differentiate theories
vect_gt_2_e <- subset(df_qual, boundary_1==1 & ending==1)$Time
vect_gt_2_bth <- subset(df_qual, boundary_1==1 & both==1)$Time


vect_gt_2 <- subset(df_qual, boundary_1==1)$Time
vect_gt_b <- subset(df_qual, beginning==1)$Time
vect_gt_e <- subset(df_qual, ending==1)$Time
vect_gt_bth <- subset(df_qual, both==1)$Time

vect_gt_2_b <- subset(df_qual, boundary_1==1 & beginning==1)$Time #same as theory 1 and therefore not something that can differentiate theories
vect_gt_2_e <- subset(df_qual, boundary_1==1 & ending==1)$Time
vect_gt_2_bth <- subset(df_qual, boundary_1==1 & both==1)$Time


#broad theory results
th1_res <- get_all_sims(obs_list, vect_gt)
th1_res[1:215,9] <- "broad"
names(th1_res)[9] <- "theory"
th1_res[1:215,10] <- 1
names(th1_res)[10] <- "level1"
th1_res[1:215,11] <- 1
names(th1_res)[11] <- "level2"
th1_res[1:215,12] <- 1
names(th1_res)[12] <- "beg"
th1_res[1:215,13] <- 1
names(th1_res)[13] <- "end"
th1_res[1:215,14] <- 1
names(th1_res)[14] <- "both"

#broad theory level 1 only 
th1_lv1 <- get_all_sims(obs_list, subset(df_qual, level==1)$Time)
th1_lv1[1:215,9] <- "broad"
names(th1_lv1)[9] <- "theory"
th1_lv1[1:215,10] <- 1
names(th1_lv1)[10] <- "level1"
th1_lv1[1:215,11] <- 0
names(th1_lv1)[11] <- "level2"
th1_lv1[1:215,12] <- 1
names(th1_lv1)[12] <- "beg"
th1_lv1[1:215,13] <- 1
names(th1_lv1)[13] <- "end"
th1_lv1[1:215,14] <- 1
names(th1_lv1)[14] <- "both"

#broad theory level 2 only 
th1_lv2 <- get_all_sims(obs_list, subset(df_qual, level==2)$Time)
th1_lv2[1:215,9] <- "broad"
names(th1_lv2)[9] <- "theory"
th1_lv2[1:215,10] <- 0
names(th1_lv2)[10] <- "level1"
th1_lv2[1:215,11] <- 1
names(th1_lv2)[11] <- "level2"
th1_lv2[1:215,12] <- 1
names(th1_lv2)[12] <- "beg"
th1_lv2[1:215,13] <- 1
names(th1_lv2)[13] <- "end"
th1_lv2[1:215,14] <- 1
names(th1_lv2)[14] <- "both"

#broad theory beginnings only 
th1_beg <- get_all_sims(obs_list, subset(df_qual, beginning==1)$Time)
th1_beg[1:215,9] <- "broad"
names(th1_beg)[9] <- "theory"
th1_beg[1:215,10] <- 1
names(th1_beg)[10] <- "level1"
th1_beg[1:215,11] <- 1
names(th1_beg)[11] <- "level2"
th1_beg[1:215,12] <- 1
names(th1_beg)[12] <- "beg"
th1_beg[1:215,13] <- 0
names(th1_beg)[13] <- "end"
th1_beg[1:215,14] <- 0
names(th1_beg)[14] <- "both"


#broad theory endings only 
th1_end <- get_all_sims(obs_list, subset(df_qual, ending==1)$Time)
th1_end[1:215,9] <- "broad"
names(th1_end)[9] <- "theory"
th1_end[1:215,10] <- 1
names(th1_end)[10] <- "level1"
th1_end[1:215,11] <- 1
names(th1_end)[11] <- "level2"
th1_end[1:215,12] <- 0
names(th1_end)[12] <- "beg"
th1_end[1:215,13] <- 1
names(th1_end)[13] <- "end"
th1_end[1:215,14] <- 0
names(th1_end)[14] <- "both"

#broad theory both end and beginning 
th1_both <- get_all_sims(obs_list, subset(df_qual, both==1)$Time)
th1_both[1:215,9] <- "broad"
names(th1_both)[9] <- "theory"
th1_both[1:215,10] <- 1
names(th1_both)[10] <- "level1"
th1_both[1:215,11] <- 1
names(th1_both)[11] <- "level2"
th1_both[1:215,12] <- 0
names(th1_both)[12] <- "beg"
th1_both[1:215,13] <- 0
names(th1_both)[13] <- "end"
th1_both[1:215,14] <- 1
names(th1_both)[14] <- "both"


#narrow theory results
th2_res <- get_all_sims(obs_list, vect_gt_2)
th2_res[1:215,9] <- "narrow"
names(th2_res)[9] <- "theory"
th2_res[1:215,10] <- 1
names(th2_res)[10] <- "level1"
th2_res[1:215,11] <- 1
names(th2_res)[11] <- "level2"
th2_res[1:215,12] <- 1
names(th2_res)[12] <- "beg"
th2_res[1:215,13] <- 1
names(th2_res)[13] <- "end"
th2_res[1:215,14] <- 1
names(th2_res)[14] <- "both"


# narrow theory level 1 only 
th2_lv1 <- get_all_sims(obs_list, subset(df_qual, boundary_1==1 & level==1)$Time)
th2_lv1[1:215,9] <- "narrow"
names(th2_lv1)[9] <- "theory"
th2_lv1[1:215,10] <- 1
names(th2_lv1)[10] <- "level1"
th2_lv1[1:215,11] <- 0
names(th2_lv1)[11] <- "level2"
th2_lv1[1:215,12] <- 1
names(th2_lv1)[12] <- "beg"
th2_lv1[1:215,13] <- 1
names(th2_lv1)[13] <- "end"
th2_lv1[1:215,14] <- 1
names(th2_lv1)[14] <- "both"


# narrow theory level 2 only 
th2_lv2 <- get_all_sims(obs_list, subset(df_qual, boundary_1==1 & level==2)$Time)
th2_lv2[1:215,9] <- "narrow"
names(th2_lv2)[9] <- "theory"
th2_lv2[1:215,10] <- 0
names(th2_lv2)[10] <- "level1"
th2_lv2[1:215,11] <- 1
names(th2_lv2)[11] <- "level2"
th2_lv2[1:215,12] <- 1
names(th2_lv2)[12] <- "beg"
th2_lv2[1:215,13] <- 1
names(th2_lv2)[13] <- "end"
th2_lv2[1:215,14] <- 1
names(th2_lv2)[14] <- "both"

#narrow theory beginnings only 
th2_beg <- get_all_sims(obs_list, subset(df_qual,  boundary_1==1 & beginning==1)$Time)
th2_beg[1:215,9] <- "broad"
names(th2_beg)[9] <- "theory"
th2_beg[1:215,10] <- 1
names(th2_beg)[10] <- "level1"
th2_beg[1:215,11] <- 1
names(th2_beg)[11] <- "level2"
th2_beg[1:215,12] <- 1
names(th2_beg)[12] <- "beg"
th2_beg[1:215,13] <- 0
names(th2_beg)[13] <- "end"
th2_beg[1:215,14] <- 0
names(th2_beg)[14] <- "both"


#broad theory endings only 
th2_end <- get_all_sims(obs_list, subset(df_qual,  boundary_1==1 & ending==1)$Time)
th2_end[1:215,9] <- "broad"
names(th2_end)[9] <- "theory"
th2_end[1:215,10] <- 1
names(th2_end)[10] <- "level1"
th2_end[1:215,11] <- 1
names(th2_end)[11] <- "level2"
th2_end[1:215,12] <- 0
names(th2_end)[12] <- "beg"
th2_end[1:215,13] <- 1
names(th2_end)[13] <- "end"
th2_end[1:215,14] <- 0
names(th2_end)[14] <- "both"

#broad theory both end and beginning 
th2_both <- get_all_sims(obs_list, subset(df_qual,  boundary_1==1 & both==1)$Time)
th2_both[1:215,9] <- "broad"
names(th2_both)[9] <- "theory"
th2_both[1:215,10] <- 1
names(th2_both)[10] <- "level1"
th2_both[1:215,11] <- 1
names(th2_both)[11] <- "level2"
th2_both[1:215,12] <- 0
names(th2_both)[12] <- "beg"
th2_both[1:215,13] <- 0
names(th2_both)[13] <- "end"
th2_both[1:215,14] <- 1
names(th2_both)[14] <- "both"





#ask Klaus - many NaN's in F1 test; AUC function stopped working see tests above

res <- rbind(th1_res, th1_lv1) %>% 
  rbind(th2_res) %>% 
  rbind(th2_lv1) %>% 
  rbind(th1_lv2) %>% 
  rbind(th2_lv2) %>% 
  rbind(th1_beg) %>% 
  rbind(th1_end) %>% 
  rbind(th1_both) %>% 
  rbind(th2_beg) %>% 
  rbind(th2_end) %>% 
  rbind(th2_both)
  
  
#################################################





#############trying all sorts of things to reduce the mess above ##########################



get_all_gt <- function(truth_var = c("boundary_all","boundary_1", "start_end", "level"), gt_df) {
  temp_list2 =NA
  n = 0
  #all combinations of gt indicators 
  for(i in 1:length(truth_var)){
    print(paste("i is",i))
    comb_per_r_temp = combinations(length(truth_var), i, v = truth_var, set = TRUE, repeats.allowed = FALSE)
    for(j in 1:dim(comb_per_r_temp)[1]){
      print(paste("j is",j))
      comb_temp = comb_per_r_temp[j,]
      temp_ind2 = paste0(comb_temp[1],":",gt_df[[comb_temp[1]]])
      if(!length(comb_temp) == 1){
        for(k in 2:length(comb_temp)){
          print(paste("k is",k))
          temp_ind2 = paste0(temp_ind2[[1]],"-",paste0(comb_temp[k],":",gt_df[[comb_temp[k]]]))  
        }
        temp_list2 = c(temp_list2, temp_ind2)
      }
    }
  }
  return(temp_list2[-1])
}
    #temp_ind2 <-paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$start_end, df_qual$level)  
  
  
  temp_iter2 = unique(temp_ind2)
  for (jj in 1:length(temp_iter2)){
    print(paste("jj is",jj))
    n = n+1
    #subset df based on whether row matches boolean list
    temp_list2[n] <- subset(df_qual, temp_ind2==temp_iter2[jj])["Time"]
    #set name based on group membership e.g. 10 vs 01  
    names(temp_list2)[n] <- temp_iter2[jj] 
  #naming


grand_truth_list = get_all_gt(gt_df = df_qual)
# removing boundary_1:0
grand_truth_list = grand_truth_list[-grep("boundary_1:0",names(grand_truth_list))]



# create character vector of all the possible combinations of the variables specified in brackets
temp_ind <-paste0(df_qual$boundary_all, df_qual$boundary_1)
#get all the unique values from the character vector
temp_iter <- unique(paste0(df_qual$boundary_all, df_qual$boundary_1))

#create a list to be used in the loop
temp_list <- list()

#temp_ind==temp_iter[i] - creates Boolean vector matching given char to list of chars - run it to understand

for (i in 1:length(temp_iter)){
  #subset df based on whether row matches boolean list
  temp_list[i] <- subset(df_qual, temp_ind==temp_iter[i])["Time"]
  #set name based on group membership e.g. 10 vs 01  
  names(temp_list)[i] <- temp_iter[i] 
}


names(temp_list) <- c("broad", "narrow")

# create character vector of all the possible combinations of the variables specified in brackets
temp_ind2 <-paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$start_end, df_qual$level)
#get all the unique values from the character vector
temp_iter2 <- unique(paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$start_end, df_qual$level))

temp_list2 <- list()


for (i in 1:length(temp_iter2)){
  #subset df based on whether row matches boolean list
  temp_list2[i] <- subset(df_qual, temp_ind2==temp_iter2[i])["Time"]
  #set name based on group membership e.g. 10 vs 01  
  names(temp_list2)[i] <- temp_iter2[i] 
}


# create character vector of all the possible combinations of the variables specified in brackets
temp_ind3 <-paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$start_end)
#get all the unique values from the character vector
temp_iter3 <- unique(paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$start_end))

# create character vector of all the possible combinations of the variables specified in brackets
temp_ind4 <-paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$level)
#get all the unique values from the character vector
temp_iter4 <- unique(paste0(df_qual$boundary_all, df_qual$boundary_1, df_qual$level))





#creating a list of test vectors

#Theory 1 vectors


#creating list of ground truth vectors

gt_list <- list()
gt_list_name <- list(c("main", "alt", "main_b", "main_e", "main_bth",  ))

gt_list[[1]] <-  df_qual$Time
gt_list[[2]] <- subset(df_qual, boundary_1==1)$Time
gt_list[[3]] <- subset(df_qual, beginning==1)$Time
gt_list[[4]] <- subset(df_qual, ending==1)$Time
gt_list[[5]] <- subset(df_qual, both==1)$Time
gt_list[[6]] <- subset(df_qual, boundary_1==1 & beginning==1)$Time
gt_list[[7]] <- subset(df_qual, boundary_1==1 & ending==1)$Time
gt_list[[8]] <- subset(df_qual, boundary_1==1 & both==1)$Time
gt_list[[9]] <-
gt_list[[10]] <-
gt_list[[11]] <-
  
#function that returns conditional vectors
  # theory: 1, 2 or both
  #end_beg 1, 2, 3, 
  #level 
  #the above are inputs - still thinking about it
  
gt_vects <- function(theory, )

  
  
  
##################### OLD stuff#######################    

lapp_df <- function (gt_df){
  res_df <- lapply(list, sum_func, vect_gt=vect_gt, sim=test, bw=bw)
  res_df <- do.call(rbind.data.frame, res_df)
  
  res_df[1:43,2] <- test
  res_df[1:43,3] <- bw
  res_df[1:43,4] <- "main"
  res_df[1:43,5] <- 1:43
  
  
  
  colnames(res_df) <- c("result", "test", "bw", "gt", "ID")
  
  return(res_df)
}
#all
vect_gt <- df_qual$Time
#beginnings, endings, and both
vect_gt_b <- subset(df_qual, beginning==1)$Time
vect_gt_e <- subset(df_qual, ending==1)$Time
vect_gt_bth <- subset(df_qual, both==1)$Time


vect_gt_2 <- subset(df_qual, boundary_1==1)$Time
vect_gt_2_b <- subset(df_qual, boundary_1==1 & beginning==1)$Time #same as theory 1 and therefore not something that can differentiate theories
vect_gt_2_e <- subset(df_qual, boundary_1==1 & ending==1)$Time
vect_gt_2_bth <- subset(df_qual, boundary_1==1 & both==1)$Time


vect_gt_2 <- subset(df_qual, boundary_1==1)$Time
vect_gt_b <- subset(df_qual, beginning==1)$Time
vect_gt_e <- subset(df_qual, ending==1)$Time
vect_gt_bth <- subset(df_qual, both==1)$Time

vect_gt_2_b <- subset(df_qual, boundary_1==1 & beginning==1)$Time #same as theory 1 and therefore not something that can differentiate theories
vect_gt_2_e <- subset(df_qual, boundary_1==1 & ending==1)$Time
vect_gt_2_bth <- subset(df_qual, boundary_1==1 & both==1)$Time


# plots
test %>% ggplot(aes(x = result, y = factor(bw), fill=bw)) + geom_density_ridges(alpha= 2) + facet_wrap(~test)

test %>% filter(bw==1) %>% select(-c(id, bw)) %>% psych::pairs.panels()


test %>% pivot_longer(-c(id, bw)) %>%  ggplot(aes(x = value, y = factor(bw), fill=bw)) + geom_density_ridges(alpha= 2) + facet_wrap(~name, scale="free_x")

# wide to long 

test %>% pivot_longer(-c(id, bw))






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




# learning 

a <- 1 
> a
[1] 1
> a[1]
[1] 1
> a[2] <2
[1] NA
> a[2] <- 2
> a
[1] 1 2
> names(a) <- c("col 1")
> a
col 1  <NA> 
  1     2 
> a[Col 1]
Error: unexpected numeric constant in "a[Col 1"
> a["Col 1"]
<NA> 
  NA 
> a["col 1"]
col 1 
1 
> 

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

