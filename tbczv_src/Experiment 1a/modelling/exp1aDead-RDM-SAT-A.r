# hierarchical RDM fit with PMwG to exp1aSAT
rm(list=ls())

# load required packages
require(tidyverse)
require(mvtnorm)
require(pmwg)

### Identifiers for experiment, model, parameterisation
expt <- "exp1aDead"
model <- "RDM"
model.par <- "SAT-A"
jobid <- Sys.getenv()["PBS_JOBID"]

# load data and likelihood
load(paste0("data/", expt, model, ".RData"))
fnam <- paste0(expt, "-", model, "-", model.par, "-", jobid, ".RData")
data <- rename(data, subject = subj)

## Name the parameters:
par_names <- c("u_cost_5", "u_sides_5", "u_deliveryTime_5", 
               "u_cost_10", "u_sides_10", "u_deliveryTime_10", 
               "u_cost_20", "u_sides_20", "u_deliveryTime_20", 
               "alpha_E_5", "alpha_E_10", "alpha_E_20",    # evidence threshold
               "tau_E")                                    # non-decision time


# source DIC function for pmwg sampler object
source("pmwgDIC.r")

# source TRDM function for pmwg sampler object
source("TRDM-functions.r")


# convert attributes from character to numeric
for(i in 1:3) {
  data[[paste0("opt",i,"cost")]] <- as.numeric(substring(data[[paste0("opt",i,"cost")]], 3, 7))
  data[[paste0("opt",i,"sides")]] <- as.numeric(substring(data[[paste0("opt",i,"sides")]], 1, 2))
  data[[paste0("opt",i,"deliveryTime")]] <- as.numeric(substring(data[[paste0("opt",i,"deliveryTime")]], 2, 3))
}

# Mean centre attributes
for(i in 1:3) {
  data[[paste0("opt",i,"mcCost")]] <- data[[paste0("opt",i,"cost")]] - mean(data[[paste0("opt",i,"cost")]])
  data[[paste0("opt",i,"mcSides")]] <- data[[paste0("opt",i,"sides")]] - mean(data[[paste0("opt",i,"sides")]]) 
  data[[paste0("opt",i,"mcDeliveryTime")]] <- data[[paste0("opt",i,"deliveryTime")]] - mean(data[[paste0("opt",i,"deliveryTime")]])
}

# log density of the observed choices under subject-level parameters (RDM, SAT-A)
ll <- function(x, data, sample=FALSE) {
  names(x) <- par_names
  x[-(1:9)] <- exp(x[-(1:9)])
  n.options <- 3
  n.trials <- nrow(data)
  x["alpha_T"] <- 1  # timing threshold not estimated
  x["tau_T"] <- 60  # fixed offset for timer from Hawkins & Heathcote (in press)
  x["gamma_T"] <- .001 # small value to make sure it does not terminate
  
  # assume additive utility representation - like the MNL
  # multiply each attribute level by the corresponding utility and sum across attributes
  deadlines <- sort(unique(data$maxTrialTime))
  u <- drifts <- array(dim=c(nrow(data), n.options), data=NA)
  for (d in deadlines) {
    u[data$maxTrialTime==d,] <- as.matrix(
      (data[paste0("opt", 1:n.options, "mcCost")][data$maxTrialTime==d,] * x[paste0("u_cost_",d)]) +
        (data[paste0("opt", 1:n.options, "mcSides")][data$maxTrialTime==d,] * x[paste0("u_sides_",d)]) +
        (data[paste0("opt", 1:n.options, "mcDeliveryTime")][data$maxTrialTime==d,] * x[paste0("u_deliveryTime_",d)]))  
  }
  
  # drift rates are exponent of utilities (Hawkins et al. Integrating cog processes, 2014)
  for (d in deadlines) {
    drifts[data$maxTrialTime==d,] <- exp(u[data$maxTrialTime==d,])
  }
  
  # convenience vectors: easier access to parameters that vary across trials
  # unname all alpha_E parameters for each trial. alpha_E (threshold parameter) is determined by deadline
  # value for the corresponding trial. Executing the code below unnames the alpha_E for each trial 
  alpha_E <- unname(x[paste0("alpha_E_", data$maxTrialTime)])
  # contaminant mixture process
  p.cont <- .05
  max.RT <- 20
  
  if(sample == FALSE) {
    like <- rep(NA, n.trials)
    for(i in 1:n.trials) {
      if(data$trialError[i] == TRUE) {
        # GH: likelihood is the product of the survivor functions for all accumulators
        # note: don't need to re-order drift rates because there is no response 
        like[i] <- prod(fS_E(t=data$maxTrialTime[i]-x["tau_E"], 
                             drift=unlist(drifts[i,]),
                             threshold=rep(alpha_E[i], n.options),
                             n.acc=n.options)$S)
      } else {
        like[i] <- dTRDM(t = data$rt[i],
                         x = c(x["alpha_T"], x["tau_T"], x["tau_E"], x["gamma_T"]),
                         # Re-order drift rates
                         drift = unlist(drifts[i, c(data$response[i], (1:3)[-data$response[i]])]),
                         threshold = rep(alpha_E[i], n.options))
      }
    }
    # contaminant mixture process operates differently for deadline data sets - different terms for trials before and after the deadline
    like[!data$trialError] <- (1-p.cont)*like[!data$trialError] + p.cont*(dunif(data$rt[!data$trialError],0,max.RT)/n.options)
    like[data$trialError] <- (1-p.cont)*like[data$trialError] + p.cont*(((max.RT-data$maxTrialTime[data$trialError])/max.RT)/n.options)
    # Return sum log likelihood. Include protection against log(0) problems
    return(sum(log(pmax(like, 1e-10))))
    
  } else {
    # data is subjects data - modified with model predictions and returned
    #    remove responses and also code for response generated from the
    #    timing process (TRUE) or evidence response (FALSE)
    data$response <- data$rt <- data$timer.response <- data$trialError <- NA
    for(i in 1:n.trials) {
      tmp <- rTRDM(n=1,
                   x=c(x["alpha_T"], x["tau_T"], x["tau_E"], x["gamma_T"]),
                   drift=unlist(drifts[i,]),
                   threshold=rep(alpha_E[i], n.options))
      data$response[i] <- tmp$R
      data$rt[i] <- tmp$RT
      data$timer.response[i] <- FALSE
    }
    # contaminant mixture process for all trials
    replace.trials <- runif(n.trials) < p.cont
    data$rt[replace.trials] <- runif(sum(replace.trials), 0, max.RT)
    data$response[replace.trials] <- sample(1:n.options, sum(replace.trials), replace=TRUE)  	
    
    # now replace trials that missed the trial deadline w/ TRUE
    data$trialError <- data$rt > data$maxTrialTime
    # GH: I've guessed that missed responses and RTs in the observed data are coded as NA. 
    #   If that's not right, edit the next two lines so they are consistent with the formatting in the observed data.  
    data$response[data$trialError] <- NA 
    data$rt[data$trialError] <- NA
    
    return(data)
  }
}


# estimation settings
n.cores <- 32
n.posterior <- 25
burn.sets <- c(n=1500, particles=100)
adapt.sets <- c(n=5000, particles=100)   # set very high as it can terminate early
sample.sets <- c(n=10000, particles=25)
epsilon <- .3

## Some good starts for a chain.
# raw scale for utilities, log scale for TRDM parameters
start_points <- list(
  mu = c(-.02, .025, .01,
         -.02, .025, .01,  
         -.02, .025, .01,  
         log(c(5, 5, 7, .2))),
  sig2 = diag(rep(.01, length(par_names)))
)

## Some priors.
priors <- list(
  theta_mu_mean = c(rep(0, 9), 2, 2, 2,0),
  theta_mu_var = diag(rep(1, length(par_names)))
)


## Initialise things.
sampler <- pmwgs(data = data, pars = par_names, prior = priors, ll_func = ll)

## Initialse more stuff, with starts above
sampled <- init(sampler, start_mu=start_points$mu, start_sig=start_points$sig2)

# Run some burn-in samples
sampled <- run_stage(sampled, stage="burn", iter=burn.sets[1], particles=burn.sets[2],
                     epsilon=epsilon, n_cores=n.cores)
save.image(paste0("modelFits/",fnam))

# adaptation phase
sampled <- run_stage(sampled, stage="adapt", iter=adapt.sets[1], particles=adapt.sets[2],
                     epsilon=epsilon, n_cores=n.cores)
save.image(paste0("modelFits/",fnam))

# sampling phase
sampled <- run_stage(sampled, stage="sample", iter=sample.sets[1], particles=sample.sets[2],
                     epsilon=epsilon, n_cores=n.cores)
save.image(paste0("modelFits/",fnam))

# calculate DIC and then perform final save
DIC <- pmwg_DIC(sampled)
save.image(paste0("modelFits/",fnam))