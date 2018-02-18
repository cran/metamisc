generateBUGS.OE.discrete <- function(N.type1, N.type2, N.type3, pars, ...) {
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 

  # Likelihood of studies providing O, E and N
  if (N.type1 > 0) {
    out <- paste(out, "for (j in 1:", N.type1, ")\n  {\n", sep="")
    out <- paste(out, "    O[s1[j]] ~ dbinom(pobs[j], N[s1[j]])\n")
    out <- paste(out, "    OE[j] <- exp(theta[s1[j]])\n")
    
    # Make sure that 'pobs' is always between 0 and 1
    out <- paste(out, "    pobs[j] <- min(OE[j], N[s1[j]]/(E[s1[j]]+1)) * E[s1[j]]/N[s1[j]] \n")
    out <- paste(out, " }\n")
  }
  
  # Likelihood of studies providing O and E
  if (N.type2 > 0) {
    out <- paste(out, "for (j in 1:", N.type2, ")\n  {\n", sep="")
    out <- paste(out, "    O[s2[j]] ~ dpois(lambda[j]) \n")
    out <- paste(out, "    lambda[j] <- exp(theta[s2[j]])*E[s2[j]]\n")
    out <- paste(out, " }\n")
  }
  
  # Likelihood for studies providing log OE ratio
  if (N.type3>0) {
    out <- paste(out, "for (j in 1:", N.type3, ")\n  {\n", sep="")
    out <- paste(out, "    logOE[s3[j]] ~ dnorm(theta[s3[j]], wsprec[j])\n")
    out <- paste(out, "    wsprec[j] <- 1/(logOE.se[s3[j]]*logOE.se[s3[j]])\n")
    out <- paste(out, " }\n")
  }
  
  # Between-study distribution of the logOR
  out <- paste(out, "for (j in 1:", (N.type1+N.type2+N.type3), ")\n  {\n", sep="")
  out <- paste(out, "    theta[j] ~ dnorm(mu.logoe, bsprec.logoe)\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec.logoe <- 1/(bsTau*bsTau)\n")
  if (pars$hp.tau.dist=="dunif") {
    out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
  } else if (pars$hp.tau.dist=="dhalft") {
    out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
  } else {
    stop("Specified prior not implemented")
  }
  out <- paste(out, "  mu.logoe ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
  out <- paste(out, "  mu.oe <- exp(mu.logoe)\n", sep="")
  out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep="")
  out <- paste(out, "  pred.logoe ~ dnorm(mu.logoe, bsprec.logoe)\n", sep="")
  out <- paste(out, "}", sep="")
  return(out)
}


generateBugsOE <- function(extrapolate=F,
                           pars,...) {
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 

  if (pars$model.oe =="normal/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.logoe, bsprec)\n")
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
    out <- paste(out, "  mu.logoe ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
    out <- paste(out, "  mu.oe <- exp(mu.logoe)\n", sep="")
    out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep="")
    out <- paste(out, "  pred.logoe ~ dnorm(mu.logoe, bsprec)\n", sep="")

  } else if (pars$model.oe=="poisson/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    obs[i] ~ dpois(mu[i])\n")
    out <- paste(out, "    mu[i] <- exc[i] * theta[i]\n")
    #out <- paste(out, "    exc[i] ~ dbinom (q[i], N[i])\n")
    #out <- paste(out, "    logit(q[i]) <- alphaQ[i]\n")
    #out <- paste(out, "    alphaQ[i] ~ dnorm(0.0,1.0E-6)\n")
    out <- paste(out, "    theta[i] <- exp(alpha[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.logoe, bsprec)\n")
    out <- paste(out, " }\n")
    out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
    if (pars$hp.tau.dist=="dunif") {
      out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else if (pars$hp.tau.dist=="dhalft") {
      out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
    } else {
      stop("Specified prior not implemented")
    }
    
    out <- paste(out, "  mu.logoe ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
    out <- paste(out, "  mu.oe <- exp(mu.logoe)\n", sep="")
    out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep="")
    out <- paste(out, "  pred.logoe ~ dnorm(mu.logoe, bsprec)\n", sep="")
  }
  
  out <- paste(out, "}", sep="")
  return(out)
}