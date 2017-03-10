#TODO: A package listed in "Suggests" or "Enhances" should be used conditionally 
# in examples or tests if it cannot straightforwardly be installed on the major R platforms.

valmeta <- function(cstat, cstat.se, cstat.95CI, OE, OE.95CI,
                    N, O, E, method="REML", knha=TRUE, verbose=FALSE, 
                    method.restore.c.se="Newcombe.4", scale.c = "logit", 
                    scale.oe = "log", n.chains = 4,
                    ...) {
  out <- list()
  out$call <- match.call()
  out$method <- method
  class(out) <- "valmeta"
  
  N.studies.OE <- 0
  
  if (!missing(OE)) {
    N.studies.OE <- length(OE)
  } else if (!missing(E)) {
    N.studies.OE <- length(E)
  }
  
  

  inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }
  logit <- function(x) { log(x/(1-x)) }
  
  # Check if we need to load runjags
  if (method=="BAYES") {
    if (!requireNamespace("runjags", quietly = TRUE)) {
      stop("The package 'runjags' is currently not installed!")
    } 
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The package 'rjags' is currently not installed!")
    } 
  }
  
  
  #Update SE(c.index) using Method 4 of Newcombe
  restore.c.var<- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", scale=scale.c) {
    n <- N.events #Number of events
    m <- N.subjects-N.events #Number of non-events
    
    if (missing(restore.method)) {
      restore.method <- "Newcombe.4"
    }
    
    if (restore.method=="Hanley" | restore.method=="Newcombe.2") {
      mstar <- m-1
      nstar <- n-1
    } else if (restore.method=="Newcombe.4") {
      mstar <- nstar <- N.subjects/2-1
    } else {
      stop ("Method not implemented yet!")
    }
    
    if (scale=="logit") {
      out <- (((1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n*cstat*(1-cstat)))
    } else {
      out <- ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
    }
    
    return(out)
  }
  
  
  if(!missing(cstat)) {
    if (missing(cstat.se) & missing(cstat.95CI)) {
      stop("No sampling error was provided for the c-statistic!")
    }
    if (missing(cstat.95CI)) {
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (is.null(dim(cstat.95CI))) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (dim(cstat.95CI)[2] != 2 | dim(cstat.95CI)[1] != length(cstat)) {
      warning("Invalid dimension for 'cstat.95CI', argument ignored.")
      cstat.95CI <- array(NA, dim=c(length(cstat),2))
    }
    if (missing(cstat.se)) {
      cstat.se <- array(NA, dim=length(cstat))
    }
    
    out$cstat <- list()
    out$cstat$method.restore.se <- method.restore.c.se 
    out$cstat$scale <- scale.c
    class(out$cstat) <- "vmasum"
    
    # Apply necessary data transformations
    if (scale.c == "identity") {
      theta <- cstat
      theta.var <- (cstat.se)**2
      theta.cil <- cstat.95CI[,1]
      theta.ciu <- cstat.95CI[,2]
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var) #Prioritize reported SE
    } else if (scale.c == "logit") {
      theta <- log(cstat/(1-cstat))
      theta.var <- (cstat.se/(cstat*(1-cstat)))**2
      theta.cil <- logit(cstat.95CI[,1])
      theta.ciu <- logit(cstat.95CI[,2])
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2
      theta.var <- ifelse(is.na(theta.var.CI), theta.var, theta.var.CI) #Prioritize variance from 95% CI
    } else {
      stop(paste("No appropriate transformation defined: '", scale.c, "'", sep=""))
    }
    
    num.estimated.var.c <- 0
    
    # Restore missing standard errors
    if (NA %in% theta.var) {
      # Restore missing estimates of the standard error of the c-statistic using information on c, N and O
      if (!missing(O) & !missing(N)) {
        if (verbose) cat("Attempting to restore missing information on the standard error of the c-statistic\n")
        theta.var.hat <- restore.c.var(cstat=cstat, N.subjects=N, N.events=O, restore.method=method.restore.c.se, scale=scale.c)
        num.estimated.var.c <- length(which(is.na(theta.var) & !is.na(theta.var.hat)))
        theta.var <- ifelse(is.na(theta.var), theta.var.hat, theta.var)
      }
      # Replace remaining missing values in theta.var by very large values
      theta.var <- ifelse(is.na(theta.var), 10e6, theta.var)
    }
    
    #Only calculate 95% CI for which no original values were available
    theta.cil[is.na(theta.cil)] <- (theta+qnorm(0.025)*sqrt(theta.var))[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm(0.975)*sqrt(theta.var))[is.na(theta.ciu)]
    

    ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu)
    colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu")
    out$cstat$data <- ds
    out$cstat$slab <- paste("Study",seq(1, length(theta)))
    out$cstat$num.estimated.var.c <- num.estimated.var.c
    
    if (method != "BAYES") { # Use of rma
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, knha=knha, ...) 
      preds <- predict(fit)
      
      results <- as.data.frame(array(NA, dim=c(1,5)))
      if (scale.c == "logit") {
        results <- c(inv.logit(coefficients(fit)), inv.logit(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
      } else {
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$cstat$rma <- fit
      out$cstat$results <- results
    } else {
      # Perform a Bayesian meta-analysis
      model <- .generateBugsCstat(link=scale.c, ...)
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, param.args=list(n=1, mean=0, sd=sqrt(1E6)))
      model.pars[[2]] <- list(param="bsTau", param.f=runif, param.args=list(n=1, min=0, max=2))
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      mvmeta_dat <- list(theta = theta,
                         theta.var = theta.var,
                         Nstudies = length(theta))
      jags.model <- runjags::run.jags(model=model, 
                             monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "priorTau", "PED"), 
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
      
      out$cstat$runjags <- jags.model
      out$cstat$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$cstat$results <- results
    }
  }
  
  ##### Prepare data for OE ratio
  if(N.studies.OE > 0) {
    if (missing(N)) {
      N <- array(NA, dim=N.studies.OE)
    }
    if (missing(OE)) {
      OE <- array(NA, dim=N.studies.OE)
    }
    if (missing(OE.95CI)) {
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    if (is.null(dim(OE.95CI))) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    if (dim(OE.95CI)[2] != 2 | dim(OE.95CI)[1] != N.studies.OE) {
      warning("Invalid dimension for 'OE.95CI', argument ignored.")
      OE.95CI <- array(NA, dim=c(N.studies.OE,2))
    }
    
    out$oe$scale <- scale.oe
    class(out$oe) <- "vmasum"
    
    #TODO: allow confidence intervals of OE ratio
    #TODO: allow E/O ratio
    # Apply necessary data transformations
    if (scale.oe == "identity") {
      theta <- OE
      theta <- ifelse(is.na(theta), O/E, theta)
      theta.var <- O*(1-(O/N))/(E**2) #BMJ eq 20 (binomial var)
      theta.var <- ifelse(is.na(theta.var), (O/(E**2)), theta.var) #BMJ eq 30 (Poisson var)
      theta.cil <- OE.95CI[,1]
      theta.ciu <- OE.95CI[,2]
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2 #Derive from 95% CI
      theta.var <- ifelse(is.na(theta.var), theta.var.CI, theta.var) #Prioritize reported SE
    } else if (scale.oe == "log") {
      theta <- log(OE)
      theta <- ifelse(is.na(theta), log(O/E), theta)
      theta.var <- (1-(O/N))/O #BMJ eq 27 (binomial var)
      theta.var <- ifelse(is.na(theta.var), (1/O), theta.var) #BMJ eq 36 (Poisson var)
      theta.cil <- log(OE.95CI[,1])
      theta.ciu <- log(OE.95CI[,2])
      theta.var.CI <- ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2
      theta.var <- ifelse(is.na(theta.var.CI), theta.var, theta.var.CI) #Prioritize variance from 95% CI
    } else {
      stop(paste("No appropriate transformation defined: '", scale.oe, "'", sep=""))
    }
    
    #Only calculate 95% CI for which no original values were available
    theta.cil[is.na(theta.cil)] <- (theta+qnorm(0.025)*sqrt(theta.var))[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm(0.975)*sqrt(theta.var))[is.na(theta.ciu)]
    
    ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu)
    colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu")
    out$oe$data <- ds
    out$oe$slab <- paste("Study",seq(1, length(theta)))

    if (method != "BAYES") { # Use of rma
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, vi=theta.var, data=ds, method=method, knha=knha, ...) 
      preds <- predict(fit)
      
      results <- as.data.frame(array(NA, dim=c(1,5)))
      if (scale.oe == "log") {
        results <- c(exp(coefficients(fit)), exp(c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub)))
      } else {
        results <- c(coefficients(fit), c(preds$ci.lb, preds$ci.ub, preds$cr.lb, preds$cr.ub))
      }
      names(results) <- c("estimate", "95CIl", "95CIu", "95PIl", "95PIu")
      
      out$oe$rma <- fit
      out$oe$results <- results
    } else {
      stop("Bayesian method not implemented yet!")
    }
  }
  
  
  
  return(out)
}

.generateBugsCstat <- function(link="logit", #Choose between 'log', 'logit' and 'binom'
                               prior="dunif", #Choose between dunif (uniform) or dhalft (half student T)
                               prior.bound=c(0,2), #boundaries for uniform prior
                               prior.sigma=0.5, ...) # standard deviation for student T prior
  {

  prior.prec <- 1/(prior.sigma*prior.sigma)
  out <- "model {\n " 
  out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
  out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
  out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
  out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
  
  if (prior=="dunif") {
    out <- paste(out, "  priorTau ~ dunif(", prior.bound[1], ",", prior.bound[2], ")\n", sep="") 
    out <- paste(out, "  bsTau ~ dunif(", prior.bound[1], ",", prior.bound[2], ")\n", sep="") 
  } else if (prior=="dhalft") {
    out <- paste(out, "  priorTau ~ dt(0,", prior.prec, ",3)T(0,10)\n", sep="") 
    out <- paste(out, "  bsTau ~ dt(0,", prior.prec, ",3)T(0,10)\n", sep="") 
    
  } else {
    stop("Specified prior not implemented")
  }
  
  if (link == "logit") {
    out <- paste(out, "  mu.tobs ~ dnorm(0.0,1.0E-6)\n", sep="")
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
  if (!is.null(x$cstat$results)) {
    cat("Model results for the c-statistic:\n\n")
    print(x$cstat)
    if (x$cstat$num.estimated.var.c > 0)
      cat(paste("\nNote: For ", x$cstat$num.estimated.var.c, " validation(s), the standard error was estimated using method '", x$cstat$method.restore.se, "'.\n", sep=""))
    
    if (!is.null(x$oe$results)) {
      cat("\n\n")
    }
    }
  if (!is.null(x$oe$results)) {
    cat("Model results for the total O:E ratio:\n\n")
    print(x$oe)
  }
}

print.vmasum <- function(x, ...) {
  print(x$results)
  if (!is.null(x$runjags)) {
    #Print penalized expected deviance
    cat(paste("\nPenalized expected deviance: ", round(x$PED,2), "\n"))
    
    # Check if model converged
    psrf.ul <-  x$runjags$psrf$psrf[,"Upper C.I."]
    psrf.target <- x$runjags$psrf$psrf.target
    
    if(sum(psrf.ul > psrf.target)>1) {
      warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                    psrf.target, "for the parameters", 
              paste(rownames(x$runjags$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                    round(x$runjags$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
              ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'."  ))
    }
  }
    
}

plot.valmeta <- function(x, ...) {
  inv.logit <- function(x) {1/(1+exp(-x)) }
  
  if (!is.null(x$cstat)) {
    if (!is.null(x$cstat$rma)) {
      # Forest plot for the c-statistic
      if (x$cstat$scale=="logit") {
        forest(x$cstat$rma, transf=inv.logit, xlab="c-statistic", addcred=T, ...)
      } else {
        forest(x$cstat$rma, transf=NULL, xlab="c-statistic", addcred=T, ...)
      }
    } else {
      col <- c("black", "gray50")
      border <- "black"
      lty <- c("solid", "dotted", "solid")
      cex <- 0.8
      efac <- 1
      efac <- rep(efac, 2)
      xlim <- c(-0.5, 1.5)
      
      par.usr <- par("usr")
      height <- par.usr[4] - par.usr[3]
      
      k <- dim(x$cstat$data)[1]
      slab <- c(x$cstat$slab, "RE Model")
      yi <- x$cstat$data[,"theta"]
      ci.lb <- x$cstat$data[,"theta.95CIl"]
      ci.ub <- x$cstat$data[,"theta.95CIu"]
      
      if (x$cstat$scale=="logit") {
        yi <- sapply(yi, inv.logit)
        ci.lb <- sapply(ci.lb, inv.logit)
        ci.ub <- sapply(ci.ub, inv.logit)
      }
      
      #Add the meta-analysis summary to the results
      #Note that no transormations are needed here, as summaries are always presented on original scale
      b <- x$cstat$results["estimate"]
      yi <- c(yi, b)
      b.ci.lb <- x$cstat$results["95CIl"]
      b.ci.ub <- x$cstat$results["95CIu"]
      b.cr.lb <- x$cstat$results["95PIl"]
      b.cr.ub <- x$cstat$results["95PIu"]
      ci.lb <- c(ci.lb, b.ci.lb)
      ci.ub <- c(ci.ub, b.ci.ub)
      
      
      rows <- c(seq(k,1),-1)
      
      annotext <- round(cbind(yi, ci.lb, ci.ub), 2)
      annotext <- matrix(apply(annotext, 2, format, nsmall = 2), ncol = 3)
      annotext <- paste(annotext[,1], "[", annotext[,2], ",", annotext[,3], "]")
      
      
      par.mar <- par("mar")
      par.mar.adj <- par.mar - c(0, 3, 1, 1)
      par.mar.adj[par.mar.adj < 0] <- 0
      par(mar = par.mar.adj)
      on.exit(par(mar = par.mar))
      
      par.usr <- par("usr")
      height <- par.usr[4] - par.usr[3]
      lheight <- strheight("O")
      cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
      cex <- par("cex") * cex.adj
      
      plot(NA, NA, xlim=xlim, ylim=c(-2,k), ylab="", xlab="c-statistic",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
      for (i in 1:k) {
        points(yi[i], rows[i], pch = 15, ...)
        
        segments(ci.lb[i], rows[i], ci.ub[i], rows[i], ...)
        
        segments(ci.lb[i], rows[i] - (height/150) * cex * 
                   efac[1], ci.lb[i], rows[i] + (height/150) * cex * 
                   efac[1], ...)
        
        segments(ci.ub[i], rows[i] - (height/150) * cex * 
                   efac[1], ci.ub[i], rows[i] + (height/150) * cex * 
                   efac[1], ...)
      }
      text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
      text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, ...)
      
      # Add prediction interval
      segments(b.cr.lb, -1, b.cr.ub, -1, lty = lty[2], col = col[2], ...)
      segments(b.cr.lb, -1 - (height/150) * cex * efac[1], 
               b.cr.lb, -1 + (height/150) * cex * efac[1], 
               col = col[2], ...)
      segments(b.cr.ub, -1 - (height/150) * cex * efac[1], 
               b.cr.ub, -1 + (height/150) * cex * efac[1], 
               col = col[2], ...)
      
      # Add diamond for summary estimate
      polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
                                                     (height/100) * cex * efac[2], -1, -1 - (height/100) * 
                                                     cex * efac[2]), col = col[1], border = border, ...)
      

      
      # Add separation line between forest plot and meta-analysis results
      abline(h = 0, lty = 1, ...)

      axis(side = 1, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 1, ...)
    }
    if (!is.null(x$oe)) {
      readline(prompt="Press [enter] to continue")
    }
  }
  if (!is.null(x$oe)) {
    if (!is.null(x$oe$rma)) {
      # Forest plot for the c-statistic
      if (x$oe$scale=="log") {
        forest(x$oe$rma, transf=exp, xlab="Total O:E ratio", addcred=T, refline=1, ...)
      } else {
        forest(x$oe$rma, transf=NULL, xlab="Total O:E ratio", addcred=T, refline=1, ...)
      }
    } else {
      warning("Plot function not implemented yet for Bayesian MA")
    }
  }
}

