valmeta <- function(measure="cstat", cstat, cstat.se, cstat.95CI, OE, OE.se, OE.95CI, citl, citl.se,
                    N, O, E, Po, Po.se, Pe, t.val, t.ma, t.extrapolate=FALSE, method="REML", test="knha", 
                    verbose=FALSE, slab, n.chains = 4, pars, ...) {
  pars.default <- list(hp.mu.mean = 0, 
                       hp.mu.var = 1E6,
                       hp.tau.min = 0,
                       hp.tau.max = 2,
                       hp.tau.sigma = 0.5,
                       hp.tau.dist = "dunif", 
                       hp.tau.df = 3, 
                       method.restore.c.se="Newcombe.4",
                       model.cstat = "normal/logit", #Alternative: "normal/identity"
                       model.oe = "normal/log") #Alternative: "poisson/log" or "normal/identity"
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  if (!is.element(measure, c("cstat","OE")))
    stop("Unknown 'measure' specified.")
  
  #######################################################################################
  # Check if we need to load runjags
  #######################################################################################
  if (method=="BAYES") {
    if (!requireNamespace("runjags", quietly = TRUE)) {
      stop("The package 'runjags' is currently not installed!")
    } 
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The package 'rjags' is currently not installed!")
    } 
    if (n.chains<1 | n.chains%%1!=0) {
      stop("Invalid number of chains specified for the Gibbs sampler!")
    }
  }
    
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  if (measure=="cstat") {
    if (!missing(cstat)) {
      k <- length(cstat)
    } else if (!missing(cstat.se)) {
      k <- length(cstat.se)
    } else if (!missing(cstat.95CI)) {
      k <- dim(cstat.95CI)[2]
    }
  } else if (measure=="OE") {
    if (!missing(OE)) {
      k <- length(OE)
    } else if (!missing(E)) {
      k <- length(E)
    } else if (!missing(Pe)) {
      k <- length(Pe)
    } else if (!missing(citl)) {
      k <- length(citl)
    }
  }
  
  #######################################################################################
  # Prepare data
  #######################################################################################
  if (missing(O)) {
    O <- rep(NA, length=k)
  }
  if (missing(Po)) {
    Po <- rep(NA, length=k)
  }
  if (missing(N)) {
    N <- rep(NA, length=k)
  }
  
  theta.var.source <- rep("-", k)
  
  #######################################################################################
  # Prepare object
  #######################################################################################
  out <- list()
  out$call <- match.call()
  out$measure <- measure
  out$method <- method
  class(out) <- "valmeta"
  
  
  #######################################################################################
  # Assign study labels
  #######################################################################################
  if(missing(slab)) {
    out$slab <- paste("Study",seq(1, k))
  } else {
    out$slab <- as.character(slab)
  }
  

  #######################################################################################
  # Meta-analysis of the c-statistic
  #######################################################################################
  if (measure=="cstat") {
    if (missing(cstat.95CI)) {
      cstat.95CI <- array(NA, dim=c(k,2))
    }
    if (is.null(dim(cstat.95CI))) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(k,2))
    }
    if (dim(cstat.95CI)[2] != 2 | dim(cstat.95CI)[1] != k) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(k,2))
    }
    if (missing(cstat.se)) {
      cstat.se <- array(NA, dim=k)
    }

    out$model <- pars.default$model.cstat
    
    if (verbose) message("Extracting/computing estimates of the concordance statistic ...")
    
    
    
    # Apply necessary data transformations
    if (pars.default$model.cstat == "normal/identity") {
      theta <- cstat
      theta.var <- (cstat.se)**2
      theta.var.source[!is.na(theta.var)] <- "Standard Error"
      theta.cil <- cstat.95CI[,1]
      theta.ciu <- cstat.95CI[,2]
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var.source[is.na(theta.var) & !is.na(theta.var.CI)] <- "Confidence Interval"
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var) #Prioritize reported SE
    } else if (pars.default$model.cstat == "normal/logit") {
      theta <- log(cstat/(1-cstat))
      theta.var <- (cstat.se/(cstat*(1-cstat)))**2
      theta.var.source[!is.na(theta.var)] <- "Standard Error"
      theta.cil <- logit(cstat.95CI[,1])
      theta.ciu <- logit(cstat.95CI[,2])
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2
      theta.var.source[is.na(theta.var) & !is.na(theta.var.CI)] <- "Confidence Interval"
      theta.var <- ifelse(is.na(theta.var.CI), theta.var, theta.var.CI) #Prioritize variance from 95% CI
    } else {
      stop(paste("No appropriate meta-analysis model defined: '", pars.default$model.cstat, "'", sep=""))
    }
    
    # Restore missing standard errors
    if (NA %in% theta.var) {
      
      # Calculate O and N from other information if possible
      O <- ifelse(is.na(O), Po*N, O)
      N <- ifelse(is.na(N), O/Po, N)

      # Restore missing estimates of the standard error of the c-statistic using information on c, N and O
      theta.var.hat <- restore.c.var(cstat=cstat, N.subjects=N, N.events=O, 
                                     restore.method=pars.default$method.restore.c.se, 
                                     model=pars.default$model.cstat)
      theta.var.source[is.na(theta.var) & !is.na(theta.var.hat)] <- pars.default$method.restore.c.se
      theta.var <- ifelse(is.na(theta.var), theta.var.hat, theta.var)

      # Replace remaining missing values in theta.var by very large values
      # Omitted because it now gives error in metafor
      # Ratio of largest to smallest sampling variance extremely large. Cannot obtain stable results.
      #theta.var <- ifelse(is.na(theta.var), 10e6, theta.var)
    }
    
    #Only calculate 95% CI for which no original values were available
    theta.cil[is.na(theta.cil)] <- (theta+qnorm(0.025)*sqrt(theta.var))[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm(0.975)*sqrt(theta.var))[is.na(theta.ciu)]
    
    
    ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu, NA)
    colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu", "theta.blup")
    
    
    if (method != "BAYES") { # Use of rma
      
      # Identify which studies can be used for meta-analysis
      selstudies <- which(!is.na(ds[,"theta"]) & !is.na(theta.var))
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, test=test, slab=out$cstat$slab, ...) 
      preds <- predict(fit)
      
      ds[selstudies, "theta.blup"] <- blup(fit)$pred
      
      results <- as.data.frame(array(NA, dim=c(1,5)))
      if (pars.default$model.cstat == "normal/logit") {
        cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
        cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
        results <- c(inv.logit(coefficients(fit)), inv.logit(c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub)))
      } else if (pars.default$model.cstat == "normal/identity") {
        cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
        cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub))
      } else {
        stop ("Meta-analysis model not implemented yet")
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$rma <- fit
      out$numstudies <- fit$k
      out$results <- results
    } else {
      # All data are used!
      out$numstudies <- dim(ds)[1]
      
      # Perform a Bayesian meta-analysis
      model <- .generateBugsCstat(pars=pars.default, ...)
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      model.pars[[2]] <- list(param="bsTau", param.f=runif, param.args=list(n=1, min=pars.default$hp.tau.min, max=pars.default$hp.tau.max))
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      mvmeta_dat <- list(theta = theta,
                         theta.var = theta.var,
                         Nstudies = length(theta))
      jags.model <- runjags::run.jags(model=model, 
                                      monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "PED"), 
                                      data = mvmeta_dat, 
                                      n.chains = n.chains,
                                      silent.jags = !verbose,
                                      inits=inits,
                                      ...)
      fit <- jags.model$summaries
      
      
      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$runjags <- jags.model
      out$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$results <- results
    }
    
    out$data <- ds
    out$se.source <- theta.var.source
    return(out)
  }
  #######################################################################################
  # Meta-analysis of the total OE ratio
  #######################################################################################
  if (measure=="OE") {
    t.ma <- ifelse(missing(t.ma), NA, t.ma)
    
    if(missing(t.val)) {
      t.val <- rep(NA, length=k)
    }
    if (missing(E)) {
      E <- rep(NA, length=k)
    }
    if (missing(Po.se)) {
      Po.se <- rep(NA, length=k)
    }
    if (missing(Pe)) {
      Pe <- rep(NA, length=k)
    }
    if (missing(OE)) {
      OE <- rep(NA, length=k)
    }
    if (missing(OE.se)) {
      OE.se <- rep(NA, length=k)
    }
    if (missing(citl)) {
      citl <- rep(NA, length=k)
    }
    if (missing(citl.se)) {
      citl.se <- rep(NA, length=k)
    }
    if (missing(OE.95CI)) {
      OE.95CI <- array(NA, dim=c(k,2))
    }
    if (is.null(dim(OE.95CI))) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(k,2))
    }
    if (dim(OE.95CI)[2] != 2 | dim(OE.95CI)[1] != k) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(k,2))
    }
    
    # Check if the length of all relevant arguments is consistent
    if (length(unique(c(length(N), length(O), length(E), length(Po), length(Po.se), 
                        length(Pe), length(OE), length(OE.se), length(citl),
                        length(citl.se), dim(OE.95CI)[1]))) > 1) {
      stop("Dimension mismatch")
    }
    
    out$model <- pars.default$model.oe

    if (verbose) message("Extracting/computing estimates of the total O:E ratio ...")
    
    ds <- generateOEdata(O=O, E=E, Po=Po, Po.se=Po.se, Pe=Pe, OE=OE, OE.se=OE.se, OE.95CI=OE.95CI, 
                         citl=citl, citl.se=citl.se, N=N, t.ma=t.ma, t.val=t.val, t.extrapolate=t.extrapolate,
                         pars=pars.default, verbose=verbose)
    
    out$numstudies <- length(which(rowMeans(!is.na(ds))==1))
      
    if (method != "BAYES") { # Use of rma
      
      if (pars.default$model.oe=="normal/identity") {
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, test=test, slab=out$slab, ...) 
        preds <- predict(fit)
        cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
        cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub))
        out$rma <- fit
        out$numstudies <- fit$k
      } else if (pars.default$model.oe=="normal/log") {
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, test=test, slab=out$slab, ...) 
        preds <- predict(fit)
        cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
        cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
        results <- c(exp(coefficients(fit)), exp(c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub)))
        out$rma <- fit
        out$numstudies <- fit$k
      } else if (pars.default$model.oe=="poisson/log" && method!="FE") {
        if (method!="ML") warning("The poisson/log model was fitted using ML.")
        if (test=="knha") warning("The Sidik-Jonkman-Hartung-Knapp correction cannot be applied")
        
        fit <- glmer(O~1|Study, offset=log(E), family=poisson(link="log"), data=ds)
        preds.ci <- confint(fit, quiet=!verbose, ...)
        preds.cr <- lme4::fixef(fit) + qt(c(0.025, 0.975), df=(lme4::ngrps(fit)-2))*sqrt(vcov(fit)[1,1]+(as.data.frame(lme4::VarCorr(fit))["vcov"])[1,1])
        results <- c(exp(lme4::fixef(fit)), exp(c(preds.ci["(Intercept)",], preds.cr)))
        out$lme4 <- fit
        out$numstudies <- nobs(fit)
      } else if (pars.default$model.oe=="poisson/log" && method=="FE") {
        fit <- glm(O~1, offset=log(E), family=poisson(link="log"), data=ds)
        preds.ci <- confint(fit, quiet=!verbose, ...)
        preds.cr <- c(NA, NA)
        results <- c(exp(coefficients(fit)), exp(c(preds.ci, preds.cr)))
        out$glm <- fit
        out$numstudies <- nobs(fit)
      } else {
        stop("Model not implemented yet!")
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$results <- results
    } else {
      if (pars.default$model.oe=="normal/log") {
        i.select <- which(!is.na(ds$theta.se)) #omit non-informative studies
        
        mvmeta_dat <- list(theta=ds$theta[i.select],
                           theta.var=(ds$theta.se[i.select])**2,
                           Nstudies = length(i.select))
      } else if (pars.default$model.oe =="poisson/log") {
        
        # Truncate hyper parameter variance
        pars.default$hp.mu.var = min(pars.default$hp.mu.var, 100)
        
        i.select <- which(!is.na(ds$E)) #omit non-informative studies
        
        mvmeta_dat <- list(obs=round(ds$O[i.select]),
                           exc=ds$E[i.select],
                           Nstudies = length(i.select))
      } else {
        stop("Model not implemented yet!")
      }
      out$numstudies <- mvmeta_dat$Nstudies
      
      # Perform a Bayesian meta-analysis
      model <- generateBugsOE(pars=pars.default, ...)
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      model.pars[[2]] <- list(param="bsTau", param.f=runif, param.args=list(n=1, min=pars.default$hp.tau.min, max=pars.default$hp.tau.max))
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      
      jags.model <- runjags::run.jags(model=model, 
                                      monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "PED"), 
                                      data = mvmeta_dat, 
                                      n.chains = n.chains,
                                      silent.jags = !verbose,
                                      inits=inits,
                                      ...)
      fit <- jags.model$summaries
      

      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      
      results <- c(fit["mu.obs","Mean"], fit["mu.obs", c("Lower95", "Upper95")], fit["pred.obs", c("Lower95", "Upper95")])
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$runjags <- jags.model
      out$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$results <- results
    }
    
    out$data <- ds
    out$se.source <- theta.var.source
    
    return(out)
  }

}

.generateBugsCstat <- function(pars, 
                                ...) # standard deviation for student T prior
{

  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
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
  
  if (pars$model.cstat  == "normal/logit") {
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep="")
    out <- paste(out, "  mu.obs <- 1/(1+exp(-mu.tobs))\n", sep="")
    out <- paste(out, "  pred.obs <- 1/(1+exp(-pred.tobs))\n", sep="")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep="")
  } else {
    stop("Specified link function not implemented")
  }
  out <- paste(out, "}", sep="")
  
  return(out)
}

print.valmeta <- function(x, ...) {
  if (x$measure=="cstat") {
    cat("Model results for the c-statistic:\n\n")
  } else if (x$measure=="OE") {
    cat("Model results for the total O:E ratio:\n\n")
  }
  
  print(x$results)
  if (!is.null(x$runjags)) {
    #Print penalized expected deviance
    cat(paste("\nPenalized expected deviance: ", round(x$PED,2), "\n"))
    
    # Check if model converged
    psrf.ul <-  x$runjags$psrf$psrf[,2]
    psrf.target <- x$runjags$psrf$psrf.target
    
    if(sum(psrf.ul > psrf.target)>0) {
      warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                    psrf.target, "for the parameters", 
                    paste(rownames(x$runjags$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                          round(x$runjags$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
                    ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'."  ))
    }
  }
  
  cat("\n")
  cat(paste("Number of studies included: ", x$numstudies))
  if (x$measure=="cstat") {
    se.sources <- c("Hanley","Newcombe.2","Newcombe.4") 
    num.estimated.var.c <- sum(x$se.source %in% se.sources)
    if (num.estimated.var.c > 0) {
      restore.method <- (se.sources[se.sources %in% x$se.source])[1]
      cat(paste("\nNote: For ", num.estimated.var.c, " validation(s), the standard error of the concordance statistic was estimated using method '", restore.method, "'.\n", sep=""))
    }
  }
}



plot.valmeta <- function(x, ...) {
  if (x$measure=="cstat") {
    plotForest(x, xlab="c-statistic", refline=NULL, ...)
  } else if (x$measure=="OE") {
    plotForest(x, xlab="Total O:E ratio", refline=1, ...)
  }
}

