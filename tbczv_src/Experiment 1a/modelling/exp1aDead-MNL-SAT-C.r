# hierarchical MNL fit with PMwG to exp1bDead
rm(list=ls())

# load required packages
require(tidyverse)
require(mvtnorm)
require(pmwg)

### Identifiers for experiment, model, parameterisation
expt <- "exp1aDead"
model <- "MNL"
model.par <- "SAT-C"
jobid <- Sys.getenv()["PBS_JOBID"]


## Name the parameters:
par_names <- c("u_cost_5", "u_sides_5", "u_deliveryTime_5",
               "u_cost_10_20", "u_sides_10_20", "u_deliveryTime_10_20")

# source DIC function for pmwg sampler object
source("pmwgDIC.r")

# load data and likelihood
load(paste0("data/", expt, model, ".RData"))
fnam <- paste0(expt, "-", model, "-", model.par, "-", jobid, ".RData")
data <- rename(data, subject = subj)

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