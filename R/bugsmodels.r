generateBugsOE <- function(extrapolate=F,
                           pars) {
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 

  if (pars$model.oe =="normal/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
    out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
    out <- paste(out, " }\n")
    out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
    if (pars$hp.tau.dist=="dunif") {
      out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else if (pars$hp.tau.dist=="dhalft") {
      out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else {
      stop("Specified prior not implemented")
    }
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
    out <- paste(out, "  mu.obs <- exp(mu.tobs)\n", sep="")
    out <- paste(out, "  pred.obs <- exp(pred.tobs)\n", sep="")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep="")

  } else if (pars$model.oe=="poisson/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    obs[i] ~ dpois(mu[i])\n")
    out <- paste(out, "    mu[i] <- exc[i] * theta[i]\n")
    #out <- paste(out, "    exc[i] ~ dbinom (q[i], N[i])\n")
    #out <- paste(out, "    logit(q[i]) <- alphaQ[i]\n")
    #out <- paste(out, "    alphaQ[i] ~ dnorm(0.0,1.0E-6)\n")
    out <- paste(out, "    theta[i] <- exp(alpha[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
    out <- paste(out, " }\n")
    out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
    if (pars$hp.tau.dist=="dunif") {
      out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else if (pars$hp.tau.dist=="dhalft") {
      out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else {
      stop("Specified prior not implemented")
    }
    
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
    out <- paste(out, "  mu.obs <- exp(mu.tobs)\n", sep="")
    out <- paste(out, "  pred.obs <- exp(pred.tobs)\n", sep="")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep="")
  }
  
  out <- paste(out, "}", sep="")
  return(out)
}