pmwg_DIC <- function(sampled, pD = FALSE){
  # Identify number of subjects
  nsubj <- length(unique(sampled$data$subject))

  # Mean log-likelihood of the overall (sampled-stage) model, for each subject
  mean_ll <- apply(sampled$samples$subj_ll[, sampled$samples$stage == "sample"],
                   1,
                   mean)

  # Mean of each parameter across iterations.
  # Keep dimensions for parameters and subjects
  mean_pars <- t(apply(sampled$samples$alpha[,, sampled$samples$stage == "sample"], 1:2, mean))

  # Name 'mean_pars' so it can be used by the log_like function
  colnames(mean_pars) <- sampled$par_names

  # log-likelihood for each subject using their mean parameter vector
  mean_pars_ll <- numeric(ncol(mean_pars))

  data <- transform(sampled$data,
                    subject = match(subject, unique(subject)))

  for (j in 1:nsubj) {
    mean_pars_ll[j] <- sampled$ll_func(mean_pars[j, ],
                                       data = data[data$subject == j,],
                                       sample = FALSE)
  }

  # Effective number of parameters
  pD <- sum(-2 * mean_ll + 2 * mean_pars_ll)

  # Deviance Information Criterion
  DIC <- sum(-4 * mean_ll + 2 * mean_pars_ll)

  if (pD){
    return(c("DIC " = DIC, " Effective parameters" = pD))
  }else{
    return(DIC)
  }

}

