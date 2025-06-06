# ********** Start of header **************
# Title: Mobile Phone DCE Deadline 
#
# Author: Guy Hawkins
# Date: 
#
# *********** End of header ****************

# read data files from discrete choice task #
rm(list=ls())

read.data=function(x) {
  # tidy data file
  # remove semicolon from selected response
  x$selected=sapply(strsplit(as.character(x$selected),";"),function(x) as.numeric(x))
  # reference blocks and trials to 1 (instead of 0)
  x$block=x$block+1
  x$round=x$round+1

  # generate wide format data frame
  n.blocks=max(unique(x$block))
  trials.per.block=max(unique(x$round))
  # maximum number of options, to set up wide format
  max.options=max(table(x$block,x$round))
  n.trials=nrow(table(x$block,x$round))*ncol(table(x$block,x$round))
  # order to store each attribute. get text that uniquely identifies each attribute level
  att=c("$","MP","GB","hours")
  att.names=c("cost","camera","memory","battery")
  n.att=length(att.names)

  ### variables to code experimental conditions

  #   attributeOrder -        fixed / random
  #   noptions -              number of options in trial
  #   SA -                    1=speed emphasis, 2=accuracy emphasis
  #   maxTrialTime -          time permitted for trial; codes countdown timer value, or max trial time when countdown timer absent
  #   recommendationPresent - TRUE / FALSE
  #   recommendedOption -     digit from 1:n.options if recommendationPresent==TRUE, else NA
  ### code row position of each attribute (1:n.att) where 1 is top and n.att is bottom
  ### generate variables that code for which option had the best level for:
  #   each attribute (cost, camera, memory, battery; irrespective of attribute position), AND
  #   each attribute position (ie, row - 1, 2, 3, 4; irrespective of attribute)

  wide.data.names=c("block","trial","attributeOrder","noptions","SA","maxTrialTime","recommendationPresent","recommendedOption","response","rt","fixationTime","trialError",paste0("pos",att.names),paste0("pos",1:n.att),paste0(att.names,"BestOption"),paste0("pos",1:n.att,"BestOption"))
  for(i in 1:max.options) {
    for(j in 1:n.att) {
      wide.data.names=c(wide.data.names,paste0("opt",i,att.names[j]))
    }
  }
  # data frame to hold participant's data
  wide.data=data.frame(matrix(nrow=n.trials,ncol=length(wide.data.names),
                              dimnames=list(NULL,wide.data.names),data=NA))


  # setup to work for fixed number of options in each trial
  row.id=0
  for(blk in 1:n.blocks) {
    for(trl in 1:trials.per.block) {
      row.id=row.id+1
      y=subset(x,x$block==blk & x$round==trl)
      n.options=nrow(y)
      wide.data$block[row.id]=blk
      wide.data$trial[row.id]=trl
      wide.data$noptions[row.id]=n.options
      wide.data$rt[row.id]=y$rt[1]/1000  # code in seconds
      wide.data$fixationTime[row.id]=y$fixationTime[1]/1000  # code in seconds
      wide.data$SA[row.id]=y$Condition[1]  # speed vs accuracy instructinos
      wide.data$maxTrialTime[row.id]=y$MaxTrialTime[1]/1000  # duration of countdown timer
      # check if their was a trial error
      if(all(as.character(y$Error)==" null")) {
        wide.data$trialError[row.id]=FALSE
      } else {
        wide.data$trialError[row.id]=as.character(y$Error[1])
      }
      # code recommendation present/absent; if present, code which option was recommended
      if(any(y$recommended==1)) {
        wide.data$recommendationPresent[row.id]=TRUE
        wide.data$recommendedOption[row.id]=which(y$recommended==1)
      } else {
        wide.data$recommendationPresent[row.id]=FALSE
        wide.data$recommendedOption[row.id]=NA
      }
      # include NAs in some columns if number of options on current trial is less than the maximum shown to  
      #   if number of options is the maximum shown to this participant, don't include any NAs
        for(i in 1:nrow(y)) {   # loop over options
          for(j in 1:n.att) {   # loop over feature_0 - feature_[n.att-1]
            # pull out level corresponding to current option and feature number
            feat.ij=as.character(y[[paste0("feature_",j-1)]][i])
            # if looping over the first option in a trial, then check and code for attribute order
            for(k in 1:n.att) {
              if(length(grep(paste0("\\",att[k]),feat.ij))==1) {
                # save attribute level for option i in trial row.id
                wide.data[[paste0("opt",i,att.names[k])]][row.id]=feat.ij
                # save position of attribute j
                wide.data[[paste0("pos",att.names[k])]][row.id]=j
                # save attribute of position k
                wide.data[[paste0("pos",j)]][row.id]=att.names[k]
              } 
            }
          }
          if(y$selected[i]==1) wide.data$response[row.id]=i
        }
      # code for best option for each attribute AND attribute position
      for(i in 1:n.att) {
        current.attribute=wide.data[row.id,paste0("pos",i)]
        att.levels=wide.data[row.id,paste0("opt",1:n.options,current.attribute)]
        tmp=as.numeric(gsub("[^[:digit:] ]","",unlist(att.levels)))
        # code best option on current.attribute and the position of current.attribute
        if(current.attribute=="cost") {
          best.option=list(which(tmp==min(tmp)))
        } else {
          best.option=list(which(tmp==max(tmp)))
        }
        wide.data[[paste0(current.attribute,"BestOption")]][row.id] =
        wide.data[[paste0("pos",wide.data[[paste0("pos",current.attribute)]][row.id],"BestOption")]][row.id] = best.option
      }
    }
  }
  # code whether attribute order was fixed or randomised across trials
  if(length(unique(wide.data$poscost))>1) {
    wide.data$attributeOrder="random"
  } else {
    wide.data$attributeOrder="fixed"
  }
  wide.data
}


# move to directory data.set
setwd("data/")
# get experimental conditions
conditions=dir()

# data frame to hold data from all participants
data=NULL

# loop over data files in each condition
for(cn in conditions) {
  cat("\n\n",cn,"\n")
  setwd(cn)
  file.names=dir(pattern="[0-9][0-9][0-9][0-9][0-9]_")
  for(f in file.names) {
    cat("\n",f)
    data=rbind(data,cbind(subj=strsplit(f,"_")[[1]][1],read.data(read.delim(f,sep=","))))
  }
  setwd("../")
}
setwd("../")

# 1. recode RT for trials with slow (missed) response as NA
data$rt[data$trialError == "Slow"] <- NA
# 2. As #1 but for fast trials
data$rt[data$trialError == "Fast"] <- NA
# 3. recode 'id' column to 'subj
# data <- rename(.data = data, subj = id )
# 4. recode SA condition
data$SA <- as.factor(data$SA)
levels(data$SA) <- c("Speed", "Accuracy")

# Save raw data as R data object (data frame)
save(data, 
     file = "exp2bSAT.RData")

