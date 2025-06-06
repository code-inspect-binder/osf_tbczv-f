# Produce posterior predictive data
require(tidyverse)

# convert subject ids to numbers
data$subject <- as.numeric(factor(data$subject))
# Function to generate posterior predictive data
gen_post <- function(sampled, n) {
  n_post <- n # Number of parameter samples from posterior distribution.
  pp_data <- list()
  S <- sampled$n_subjects
  sampled_stage <- length(sampled$samples$stage[sampled$samples$stage == "sample"])
  for (s in 1:S) {
    cat(s," ")
    iterations <- round(seq(from = (sampled$samples$idx - sampled_stage),
                            to = sampled$samples$idx,
                            length.out = n_post))
    for (i in 1:length(iterations)) {
      x <- sampled$samples$alpha[,s,iterations[i]]
      tmp <- sampled$ll_func(x = x, data = data[data$subject == s, ], sample = TRUE)
      if (i == 1) {
        pp_data[[s]] <- cbind(i, tmp)
      } else {
        pp_data[[s]] <- rbind(pp_data[[s]], cbind(i, tmp))
      }
    }
  }
  return(pp_data)
}

# Generate posterior predictive data and bind rows (100 samples)
ppData <- gen_post(sampled, 100)
# gen_post produces posterior data as a list
temp <- do.call(rbind, ppData)
# Combine original data and ppdata
data <- bind_cols(source = rep("data", nrow(data)), data)
# source to be changed - RDM or MNL
ppData <- bind_cols(source = rep("RDM", nrow(temp)), temp)
allData <- bind_rows(data, ppData)

# Save posterior predictive data and original dataset - set location
save(allData, 
     file = "")
