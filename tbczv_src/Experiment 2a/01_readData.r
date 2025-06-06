# ********** Start of header **************
# Title: Pizza DCE SAT 
#
# Author: Guy Hawkins
# Date: 3rd June 2018
#
# *********** End of header ****************

# read data files from discrete choice task #
rm(list=ls())

require(tidyverse)
#Pizza SAT directory
setwd("/Users/jpc/Documents/PhD/Experiments/Exp1b-Pizza/01_instructions/analysis/data/")

# insert name that has all participant data combined from JATOS
results_file_name <- "exp1b-sat.txt"
raw_data <- read.delim(results_file_name, sep=",")

names(raw_data)[1]<-"block"
tmp<-as.character(raw_data$block)
new_subj <- c(1,grep("&fn",tmp)) # &fn is the sona ID

raw_data <- raw_data %>% 
  mutate(Condition = factor(Condition, levels = c(" 1", " 2"), labels = c("1", "2")))

read_data <- function(x) {
  # tidy data file # Change the data in the following columns to numeric values - making sure the data in those columns AREN'T factors
  for(i in c("block","round","MaxTrialTime","fixationTime","rt")) x[[i]] <- as.numeric(as.character(x[[i]]))

  # remove semicolon from selected response
  x$selected <- sapply(strsplit(as.character(x$selected),";"),function(x) as.numeric(x))
  # reference blocks and trials to 1 (instead of 0)
  x$block <- x$block+1
  x$round <- x$round+1

  # generate wide format data frame
  n_blocks <- max(unique(x$block))
  trials_per_block <- max(unique(x$round))
  # maximum number of options, to set up wide format
  max_options <- max(table(x$block,x$round))
  n_trials <- nrow(table(x$block,x$round))*ncol(table(x$block,x$round))
  # order to store each attribute_ get text that uniquely identifies each attribute level
  att <- c("$","sides","minutes")
  att_names <- c("cost","sides","deliveryTime")
  n_att <- length(att_names)
  for(i in paste0("feature_",0:(n_att-1))) x[[i]] <- factor(as.character(x[[i]]))
  # identify which feature position (0, 1, 2) contains the sides attribute
  use_feat <- NULL
  for(i in 0:(n_att-1)) if(length(grep("sides",levels(x[[paste0("feature_",i)]])))) use_feat <- i
  # rename and reorder levels of sides attribute
  tmp <- levels(x[[paste0("feature_",use_feat)]])
  tmp[tmp==" No sides"] <- "0 sides"
  tmp[tmp==" Drink OR Garlic Bread"] <- "1 sides"
  tmp[tmp==" Drink AND Garlic Bread"] <- "2 sides"
  levels(x[[paste0("feature_",use_feat)]]) <- tmp
  x[[paste0("feature_",use_feat)]] <- factor(as.character(x[[paste0("feature_",use_feat)]]))

  ### variables to code experimental conditions

  #   attributeOrder -        fixed / random
  #   noptions -              number of options in trial
  #   SA -                    1=speed emphasis, 2=accuracy emphasis
  #   maxTrialTime -          time permitted for trial; codes countdown timer value, or max trial time when countdown timer absent
  #   recommendationPresent - TRUE / FALSE
  #   recommendedOption -     digit from 1:n_options if recommendationPresent==TRUE, else NA
  ### code row position of each attribute (1:n_att) where 1 is top and n_att is bottom
  ### generate variables that code for which option had the best level for:
  #   each attribute (cost, camera, memory, battery; irrespective of attribute position), AND
  #   each attribute position (ie, row - 1, 2, 3, 4; irrespective of attribute)

  wide_data_names <- c("block","trial","attributeOrder","noptions","SA","maxTrialTime","recommendationPresent","recommendedOption","response","rt","fixationTime","trialError",paste0("pos",att_names),paste0("pos",1:n_att),paste0(att_names,"BestOption"),paste0("pos",1:n_att,"BestOption"))
  for(i in 1:max_options) {
    for(j in 1:n_att) {
      wide_data_names <- c(wide_data_names,paste0("opt",i,att_names[j]))
    }
  }
  # data frame to hold participant's data
  wide_data <- data.frame(matrix(nrow=n_trials,ncol=length(wide_data_names),
                              dimnames=list(NULL,wide_data_names),data=NA))

  # setup to work for fixed number of options in each trial
  row_id <- 0
  for(blk in 1:n_blocks) {
    for(trl in 1:trials_per_block) {
      row_id <- row_id+1
      y <- subset(x,x$block==blk & x$round==trl)
      n_options <- nrow(y)
      wide_data$block[row_id] <- blk
      wide_data$trial[row_id] <- trl
      wide_data$noptions[row_id] <- n_options
      wide_data$rt[row_id] <- y$rt[1]/1000  # code in seconds
      wide_data$fixationTime[row_id] <- y$fixationTime[1]/1000  # code in seconds
      wide_data$SA[row_id] <- y$Condition[1]  # speed vs accuracy instructions
      wide_data$maxTrialTime[row_id] <- y$MaxTrialTime[1]/1000  # duration of countdown timer
      # check if there was a trial error
      if(all(as.character(y$Error)==" null")) {
        wide_data$trialError[row_id] <- FALSE
      } else {
        wide_data$trialError[row_id] <- as.character(y$Error[1])
      }
      # code recommendation present/absent; if present, code which option was recommended
      if(any(y$recommended==1)) {
        wide_data$recommendationPresent[row_id] <- TRUE
        wide_data$recommendedOption[row_id] <- which(y$recommended==1)
      } else {
        wide_data$recommendationPresent[row_id] <- FALSE
        wide_data$recommendedOption[row_id] <- NA
      }
      # include NAs in some columns if number of options on current trial is less than the maximum shown to
      #   if number of options is the maximum shown to this participant, don't include any NAs
      for(i in 1:nrow(y)) {   # loop over options
        for(j in 1:n_att) {   # loop over feature_0 - feature_[n_att-1]
          # pull out level corresponding to current option and feature number
          feat_ij <- as.character(y[[paste0("feature_",j-1)]][i])
          # if looping over the first option in a trial, then check and code for attribute order
          for(k in 1:n_att) {
            if((length(grep(paste0("\\",att[k]),feat_ij)) + length(grep(paste0(" ",att[k]),feat_ij)))>0) {
              # save attribute level for option i in trial row_id
              wide_data[[paste0("opt",i,att_names[k])]][row_id] <- feat_ij
              # save position of attribute j
              wide_data[[paste0("pos",att_names[k])]][row_id] <- j
              # save attribute of position k
              wide_data[[paste0("pos",j)]][row_id] <- att_names[k]
            }
          }
        }
        if(y$selected[i]==1) wide_data$response[row_id] <- i
      }
      # code for best option for each attribute AND attribute position
      for(i in 1:n_att) {
        current_attribute <- wide_data[row_id,paste0("pos",i)]
        att_levels <- wide_data[row_id,paste0("opt",1:n_options,current_attribute)]
        tmp <- as.numeric(gsub("[^[:digit:] ]","",unlist(att_levels)))
        # code best option on current_attribute and the position of current_attribute
        if(current_attribute=="sides") {
          best_option <- list(which(tmp==max(tmp)))
        } else {
          best_option <- list(which(tmp==min(tmp)))
        }
        wide_data[[paste0(current_attribute,"BestOption")]][row_id]  <- 
        wide_data[[paste0("pos",wide_data[[paste0("pos",current_attribute)]][row_id],"BestOption")]][row_id]  <-  best_option
      }
    }
  }
  # code whether attribute order was fixed or randomised across trials
  if(length(unique(wide_data$poscost))>1) {
    wide_data$attributeOrder <- "random"
  } else {
    wide_data$attributeOrder <- "fixed"
  }
  wide_data
}


# data frame to hold data from all participants
# store first participant, then loop over participants 2:N, accounting for JATOS concatenation of files
data <- cbind(subj=strsplit(as.character(raw_data$block[new_subj[2]]),"=")[[1]][2],
                        read_data(raw_data[new_subj[1]:(new_subj[2]-1),]))
for(i in 3:length(new_subj)) {
  data <- rbind(data,cbind(subj=strsplit(as.character(raw_data$block[new_subj[i]]),"=")[[1]][2],
                        read_data(raw_data[(new_subj[i-1]+2):(new_subj[i]-1),])))
}


# Save the data file for use in R Markdown  
# save(data, file = "/Users/jpc/Documents/PhD/Experiments/Exp1b-Pizza/01_instructions/analysis/dataObjects/exp1bSAT.RData")

# Things to do:
# 1. recode RT for trials with slow (missed) response as NA
data$rt[data$trialError == "Slow"] <- NA
# 2. As #1 but for fast trials
data$rt[data$trialError == "Fast"] <- NA

# SA column is incorrectly coded. Contains 2 = speed, 3 = accuracy. 
# Something going on with readData script because raw text file contains 1s and 2s
# recoded with lines below

# data$SA[data$SA == "2"] <- "1"
# # 2. As #1 but for fast trials
# data$SA[data$SA == "3"] <- "2"

# Recode 
data$SA <- factor(data$SA)
levels(data$SA) <- c("Speed", "Accuracy")
