# hierarchical MNL fit with PMwG to exp1bDead-SAT-A

# load required packages
require(tidyverse)
require(mvtnorm)
require(pmwg)

# source DIC function for pmwg sampler object
source("pmwgDIC.r")

### Identifiers for experiment, model, parameterisation
expt <- "exp1bDead"
model <- "MNL"
model.par <- "SAT-A"
jobid <- Sys.getenv()["PBS_JOBID"]

## Name the parameters:
par_names <- c("u_cost_5", "u_camera_5", "u_memory_5", "u_battery_5", 
               "u_cost_10", "u_camera_10", "u_memory_10", "u_battery_10",
               "u_cost_20","u_camera_20", "u_memory_20", "u_battery_20"
               )

# load data
load(paste0("data/", expt, model, ".RData"))
fnam <- paste0(expt, "-", model, "-", model.par, "-", jobid, ".RData")
data <- rename(data, subject = subj)

# convert attributes from character to numeric for MNL 
for(i in 1:3) {
  data[[paste0("opt",i,"cost")]] <- as.numeric(substring(data[[paste0("opt",i,"cost")]], 3, 4))
  data[[paste0("opt",i,"camera")]] <- as.numeric(gsub(pattern = "MP","",x = data[[paste0("opt", i,"camera")]]))
  data[[paste0("opt",i,"memory")]] <- as.numeric(substring(data[[paste0("opt",i,"memory")]], 1, 3))
  data[[paste0("opt",i,"battery")]] <- as.numeric(substring(data[[paste0("opt",i,"battery")]], 1, 3))
}

# Mean centre attributes
for(i in 1:3) {
  data[[paste0("opt",i,"mcCost")]] <- data[[paste0("opt",i,"cost")]] - mean(data[[paste0("opt",i,"cost")]])
  data[[paste0("opt",i,"mcCamera")]] <- data[[paste0("opt",i,"camera")]] - mean(data[[paste0("opt",i,"camera")]]) 
  data[[paste0("opt",i,"mcMemory")]] <- data[[paste0("opt",i,"memory")]] - mean(data[[paste0("opt",i,"memory")]])
  data[[paste0("opt",i,"mcBattery")]] <- data[[paste0("opt",i,"battery")]] - mean(data[[paste0("opt",i,"battery")]])
}


# log density of the observed choices under subject-level parameters (MNL, SAT-A)
ll <- function(x, data, sample = FALSE) {
  names(x) <- par_names
  n.options <- 3
  # MNL assumes additive utility representation.
  # multiply each attribute level by the corresponding utility and sum across attributes
  # Separate the deadline trials for independent utility parameters
  deadlines <- sort(unique(data$maxTrialTime))
  u <- array(dim=c(nrow(data), n.options), data=NA)
  for (d in deadlines) {
    u[data$maxTrialTime==d,] <- as.matrix(
      (data[paste0("opt", 1:n.options, "mcCost")][data$maxTrialTime==d,] * x[paste0("u_cost_",d)]) +
        (data[paste0("opt", 1:n.options, "mcMemory")][data$maxTrialTime==d,] * x[paste0("u_memory_",d)]) +
        (data[paste0("opt", 1:n.options, "mcCamera")][data$maxTrialTime==d,] * x[paste0("u_camera_",d)]) +
        (data[paste0("opt", 1:n.options, "mcBattery")][data$maxTrialTime==d,] * x[paste0("u_battery_",d)]))  
  }
  
  # convert utilities to choice probabilities following Luce's choice rule
  p <- apply(u, 1, function(x) exp(x)/sum(exp(x)))
  
  if(sample == FALSE) {
    # likelihood is the probability of the observed response
    like <- diag(p[data$response, ])
    # Return sum log likelihood. Include protection against log(0) problems
    return(sum(log(pmax(like, 1e-10))))
    
  } else {
    # data is subjects data - modified with model predictions and returned
    #    remove responses
    data$response <- NA
    for(i in 1:nrow(data)) {
      data$response[i] <- sample(x = 1:n.options, size = 1, prob = p[, i])
    }
    return(data)
  }
}



# estimation settings
n.cores <- 32
n.posterior <- 25
burn.sets <- c(n=1000, particles=100)
adapt.sets <- c(n=5000, particles=100)   # set very high as it can terminate early
sample.sets <- c(n=10000, particles=25)
epsilon <- .3

## Some good starts for a chain.
start_points <- list(
  mu = rep(0, length(par_names)),
  sig2 = diag(rep(.01, length(par_names)))
)

## Some priors
priors <- list(
  theta_mu_mean = rep(0, length(par_names)),
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