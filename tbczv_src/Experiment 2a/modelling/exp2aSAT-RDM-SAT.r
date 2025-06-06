# hierarchical RDM fit with PMwG to exp2aSAT
rm(list=ls())

# load required packages
require(tidyverse)
require(mvtnorm)
require(pmwg)

### Identifiers for experiment, model, parameterisation
expt <- "exp2aSAT"
model <- "RDM"
model.par <- "SAT"
jobid <- Sys.getenv()["PBS_JOBID"]

# source DIC function for pmwg sampler object
source("pmwgDIC.r")

# source TRDM functions
source("TRDM-functions.r")

## Name model parameters:
# by default, we assume effects of SAT on evidence threshold, timing rate, and non-decision time
# threshold in speed-emphasis condition is the scaling parameter (alpha_E_speed = 1)
par_names <- c("u_cost_speed", "u_sides_speed", "u_deliveryTime_speed",
               "u_cost_accuracy", "u_sides_accuracy", "u_deliveryTime_accuracy",
               "alpha_E_speed", "alpha_E_accuracy", # evidence threshold
               "tau_E")        # non-decision time


# load data
load(paste0("data/", expt, model, ".RData"))
fnam <- paste0(expt, "-", model, "-", model.par, "-", jobid, ".RData")
data <- rename(data, subject=subj)
data$SA <- factor(data$SA)
data$SA <- recode(data$SA, "2"="speed", "3"="accuracy")

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


# log density of the observed choices under subject-level parameters (TRDM, SAT)
ll <- function(x, data, sample=FALSE) {
  names(x) <- par_names
  x[-(1:6)] <- exp(x[-(1:6)])
  n.options <- 3
  n.trials <- nrow(data)
  x["alpha_T"] <- 1  # timing threshold not estimated
  x["tau_T"] <- 60  # fixed offset for timer from Hawkins & Heathcote (in press)
  x["gamma_T"] <- .001 # small value to make sure it does not terminate
  
  # assume additive utility representation - like the MNL
  # multiply each attribute level by the corresponding utility and sum across attributes
  is.speed <- data$SA == "speed"
  u <- drifts <- array(dim=c(nrow(data), n.options), data=NA)
  u[is.speed,] <- as.matrix(
    (data[paste0("opt", 1:n.options, "mcCost")][is.speed,] * x["u_cost_speed"]) +
      (data[paste0("opt", 1:n.options, "mcSides")][is.speed,] * x["u_sides_speed"]) +
      (data[paste0("opt", 1:n.options, "mcDeliveryTime")][is.speed,] * x["u_deliveryTime_speed"]))
  
  u[!is.speed,] <- as.matrix(
    (data[paste0("opt", 1:n.options, "mcCost")][!is.speed,] * x["u_cost_accuracy"]) +
      (data[paste0("opt", 1:n.options, "mcSides")][!is.speed,] * x["u_sides_accuracy"]) +
      (data[paste0("opt", 1:n.options, "mcDeliveryTime")][!is.speed,] * x["u_deliveryTime_accuracy"]))
  
  
  # drift rates given as exponent of utility
  drifts[is.speed,] <- exp(u[is.speed,])
  drifts[!is.speed,] <- exp(u[!is.speed,])
  # convenience vectors: easier access to parameters that vary across trials
  alpha_E <- unname(x[paste0("alpha_E_", data$SA)])
  #tau_E <- unname(x[paste0("tau_E_", data$SA)])
  
  # contaminant mixture process
  p.cont <- .05
  max.RT <- 40  
  
  if(sample == FALSE) {
    like <- rep(NA, n.trials)
    for(i in 1:n.trials) {
      like[i] <- dTRDM(t=data$rt[i],
                       x=c(x["gamma_T"], x["alpha_T"], x["tau_T"], x["tau_E"]),
                       drift=unlist(drifts[i, c(data$response[i], (1:3)[-data$response[i]])]),
                       threshold=rep(alpha_E[i], n.options))
    }
    like <- (1-p.cont)*like + p.cont*(dunif(data$rt,0,max.RT)/n.options)
    # Return sum log likelihood. Include protection against log(0) problems
    return(sum(log(pmax(like, 1e-10))))
    
  } else {
    # data is subjects data - modified with model predictions and returned
    #    remove responses and also code for response generated from the
    #    timing process (TRUE) or evidence response (FALSE)
    data$response <- data$rt <- data$timer.response <- NA
    for(i in 1:n.trials) {
      tmp <- rTRDM(n=1,
                   x=c(x["gamma_T"], x["alpha_T"], x["tau_T"], x["tau_E"]),
                   drift=unlist(drifts[i,]),
                   threshold=rep(alpha_E[i], n.options))
      data$response[i] <- tmp$R
      data$rt[i] <- tmp$RT
      data$timer.response[i] <- tmp$timer.response
    } # Contaminant mixture process
    replace.trials <- runif(n.trials) < p.cont
    data$rt[replace.trials] <- runif(sum(replace.trials), 0, max.RT)
    data$response[replace.trials] <- sample(1:n.options, sum(replace.trials), replace=TRUE)
    data$timer.response[replace.trials] <- FALSE
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
# raw scale for utilities, log scale for TRDM parameters
start_points <- list(
  mu = c(-.1, .4, -.05, -.1, .4, -.05, 
         log(c(5, 7, .2))),
  sig2 = diag(rep(.01, length(par_names)))
)

## Some priors.
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