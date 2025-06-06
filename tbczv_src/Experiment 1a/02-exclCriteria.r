# ********** Start of header **************
# Title: Exp1b Deadline Pizza: Exclusion criteria
#
# Author: J-P Cavallaro (jon-paul.cavallaro@uon.edu.au)
# Date: 2nd June 2020
#
# *********** End of header ****************


# Load and install required packages
pkgs = c("tidyverse", "plotly", "ggforce", "reshape2")
for (pkg in pkgs){
  if (!require(pkg, character.only = T)){
    install.packages(pkg)
    library(pkg)
  }
}

# Load data data object - i.e. output from read-data.r script
load(file = "exp1aDead.RData")

##--------------------REMOVE PARTICIPANTS WHO COMPLETED FEWER THAN 10 BLOCKS-----------------------##
# Find subjs with fewer than 10 blocks
complete.subj <- with(data, tapply(block, subj, max))
# count number of subjs
dim(complete.subj)
# store subjs with 12 blocks from complete.subj 
keepers <- names(complete.subj)[complete.subj == 10]
# Keep subjs with 12 blocks 
new.data <- data[as.character(data$subj) %in% keepers,]
# Check all subjs have completed 12 blocks
table(new.data$subj)
# Drop levels of subjs with 0 (Subjs removed w/ < 10 blocks)
new.data <- droplevels(new.data)
# Check
table(new.data$subj)
# Count no. remaining participants
dim(table(new.data$subj))
# Assign to new variable
data <- new.data


# Porportion of trials missed within each block, for each subj
miss.blk <- with(data, tapply(trialError, list(subj, data$block),
                                       function(x) sum(x!="FALSE") / length(x)))
# ggplot - requires long shaped data
miss.blk <- melt(miss.blk)
# Convert block (Var2) and subj (Var1) to factor to display 
miss.blk$Var1 <- as.factor(miss.blk$Var1)
miss.blk$Var2 <- as.factor(miss.blk$Var2)
# Check structure of df
str(miss.blk)
#Change column names
colnames(miss.blk)[colnames(miss.blk)=="Var1"] <- "Participant"
colnames(miss.blk)[colnames(miss.blk)=="Var2"] <- "Block"
# Plot proportion of missed trials, across blocks for each participant
miss.trials.plot <- ggplot(miss.blk, aes(x = Block, y = value, group = Participant, colour = Participant)) +
  geom_line(stat="identity") + geom_point() + ylim(0,1) + ylab("Proportion of missed trials\n") +
  theme(axis.title.y = element_text(family = "Helvetica",size = 18, colour = "#333333"), 
        axis.text.y = element_text(family = "Helvetica",size = 14, colour = "#333333"), 
        axis.title.x = element_text(family = "Helvetica",size = 18, colour = "#333333"), 
        axis.text.x = element_text(family = "Helvetica",size = 14, colour = "#333333")) 
        # legend.text = element_text(family = "Helvetica", size= 12, colour = "#333333",  face = "italic"),  
        # legend.key.size = unit(0.5, 'cm'), 
        # legend.title = element_text(family = "Helvetica", size= 14, colour = "#333333",  face = "italic"),
        # legend.title.align = 0.5)
# Make plot interactive with 
ggplotly(miss.trials.plot)


# Remove participants w/ >50% missed trials in ANY block 
# find subj ids with >50% missed trials / block
excl.part <- miss.blk$Participant[miss.blk$value >= 0.50]
# drop unused levels of subj factor - i.e string of subj ids with >50% missed trials
excl.part <- unique(droplevels.factor(excl.part))
# Create 'cleaned' data frame with subjs <50% missed trials/block
clean.data <- subset(data,
                     !(data$subj 
                       %in% excl.part))
# drop unused levels of subj factor in clean.data df
clean.data$subj <- droplevels(clean.data$subj)
# N = 45
unique(clean.data$subj)

#---------------------------------Response Bias---------------------------------#
# Check for response bias (Bias to right or left key)
r.b.l <- with(subset(clean.data,!is.na(response)),
                         tapply(response, subj, 
                                function(x) 
                                  sum(x==1) / length(x)))
r.b.m <- with(subset(clean.data,!is.na(response)),
              tapply(response, subj, 
                     function(x) 
                       sum(x==2) / length(x)))

r.b.r <- with(subset(clean.data,!is.na(response)),
                           tapply(response, subj, 
                                  function(x)
                                    sum(x==3) / length(x)))

# Column bind left and right bias arrays
bias <- cbind(r.b.l, r.b.m, r.b.r)
# convert to data frame
bias <- as.data.frame(bias)
# melt for ggplot
bias <- melt(t(bias))
# convert subj to factor
bias$Var2 <- as.factor(bias$Var2)
# check structure - need factors for X and Y
str(bias)
# Change column names
names(bias)[1] <- "Bias"
names(bias)[2] <- "Participant"
# Change legend text
bias$Bias <- str_replace_all(bias$Bias, 'r.b.l', 'Opt1')
bias$Bias <- str_replace_all(bias$Bias, 'r.b.m', 'Opt2')
bias$Bias <- str_replace_all(bias$Bias, 'r.b.r', 'Opt3')

# Plot response bias x subj
optBiasPlt <- ggplot(bias, aes(fill = Bias, y = value, x = Participant)) + 
  geom_bar(stat = "identity") + 
  ylab("Proportion of trials") +
  theme(axis.text.x = element_text(angle = 90))

# Make plot interactive with 
ggplotly(optBiasPlt)


# # Look at response/option bias for each subj within blocks
rb.blk.l <- with(subset(clean.data,!is.na(response)),
              tapply(response, list(subj, block),
                     function(x)
                       sum(x==1) / length(x)))
rb.blk.m <- with(subset(clean.data,!is.na(response)),
              tapply(response, list(subj, block),
                     function(x)
                       sum(x==2) / length(x)))

rb.blk.r <- with(subset(clean.data,!is.na(response)),
              tapply(response, list(subj, block),
                     function(x)
                       sum(x==3) / length(x)))
# Add opt column name for grouping
rb.blk.l <- as.data.frame(rb.blk.l) %>% add_column(opt = "opt1")
rb.blk.m <- as.data.frame(rb.blk.m) %>% add_column(opt = "opt2")
rb.blk.r <- as.data.frame(rb.blk.r) %>% add_column(opt = "opt3")
# get subjids
subjids <- rownames(rb.blk.l)
# Row bind left, middle and right bias lists
blk.bias <- rbind(rb.blk.l, rb.blk.m, rb.blk.r) %>%
  pivot_longer(cols = 1:10, names_to = "block", values_to = "Bias") %>%
  add_column(subject = rep(subjids,each = 10, 3))
# Block as factor
blk.bias$block <- as.factor(blk.bias$block) %>%
  factor(levels = 1:10)
# Plot bias for each particpant w/in each block
part.blk.bias.plt <- ggplot(blk.bias, aes(fill = opt, y = Bias, x = block)) +
    geom_col() +
  ylab("Proportion of trials") +
  facet_wrap(~subject)

ggplotly(part.blk.bias.plt)

# Check for participants with >2/3 response in one or more blocks
duds <- subset(blk.bias, Bias > 0.66)

data <- clean.data
# save the data set 
save(data, file = "exp1aDeadReady.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Trial errors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prepare dataset for modelling
# Too fast trials - RTs < 1 second
# remove these trials - not a 'true' decision
data <- data[!data$trialError == " fast",] 
# RDM dataset - all trials remain in the data set and their likelihood is according to the tail probability
# Export data frame as data object with all too slow trials, but code as trialError = TRUE
data$trialError[data$trialError == " Slow"] <- TRUE
data$trialError <- as.logical(data$trialError)
save(data, file = "exp1aDeadRDM.RData")
# Too slow trials - fail to respond before deadline
# MNL dataset - remove trials where participants didn’t respond before the deadline and model all other responses
data <- data[!data$trialError == TRUE,] 
save(data, file = "exp1aDeadMNL.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
