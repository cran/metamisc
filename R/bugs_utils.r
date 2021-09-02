# Random effects meta-analysis of regression coefficients
run_Bayesian_REMA <- function(x, pars, n.chains, verbose, 
                              FUN_generate_bugs, # Function to generate BUGS code
                              ...) {
  # Identify number of studies
  numstudies = length(x$theta)
  
  # Perform a Bayesian meta-analysis
  model <- FUN_generate_bugs(pars = pars, ...)
  
  # Construct the hyperparameters
  model.pars <- generateHyperparametersMA(pars)
  
  # Generate initial values from the relevant distributions
  inits <- generateMCMCinits(n.chains = n.chains, model.pars = model.pars)
  
  
  jags.model <- runjags::run.jags(model = model$model.text, 
                                  monitor = c(model$model.pars, "PED"), 
                                  data = x, 
                                  confidence =  pars$level , # Which credibility intervals do we need?
                                  n.chains = n.chains,
                                  silent.jags = !verbose,
                                  inits = inits,
                                  ...)
  
  # Check if model converged
  psrf.ul <-  jags.model$psrf$psrf[,2]
  psrf.target <- jags.model$psrf$psrf.target
  
  if(sum(psrf.ul > psrf.target)>1) {
    warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                  psrf.target, "for the parameters", 
                  paste(rownames(jags.model$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                        round(jags.model$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
                  ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'.", sep=""))
  }
  
  fit <- jags.model$summaries
  
  #Extract PED
  fit.dev <- runjags::extract(jags.model,"PED")
  txtLevel <- (pars$level*100)
  
  out <- list(numstudies = numstudies, 
              fit = jags.model, 
              PED = sum(fit.dev$deviance) + sum(fit.dev$penalty),
              est = fit[model$model.pars["mu_t"], "Median"],
              se = fit[model$model.pars["mu_t"], "SD"],
              tau2 = fit[model$model.pars["tau2"], "Median"],
              se.tau2 = fit[model$model.pars["tau2"], "SD"],
              ci.lb  = fit[model$model.pars["mu_t"], paste("Lower", txtLevel, sep = "")],
              ci.ub  = fit[model$model.pars["mu_t"], paste("Upper", txtLevel, sep = "")],
              pi.lb  = fit[model$model.pars["theta_new_t"],  paste("Lower", txtLevel, sep = "")],
              pi.ub  = fit[model$model.pars["theta_new_t"],  paste("Upper", txtLevel, sep = "")])
  out
}
