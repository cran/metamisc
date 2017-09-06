#Multivariate meta-analyse: http://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix/ (directly take log-determinant)

#TODO: allow data transformations
uvmeta <- function(r, r.se, method="REML", test="knha", labels, na.action,
                   n.chains=4, pars, verbose=FALSE, ...) 
  UseMethod("uvmeta")

uvmeta.default <- function(r, r.se, method="REML", test="knha", labels, na.action, 
                           n.chains, pars, verbose=FALSE, ...)
{

  
  pars.default <- list(level = 0.95,
                       hp.mu.mean = 0, 
                       hp.mu.var = 1000) 
  
  # Check if we need to load runjags
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
  
  
  if (length(r)!=length(r.se)) {
    stop("The vectors 'r' and 'r.se' have different lengths!")
  }
  
  ds <- as.data.frame(cbind(as.vector(r),as.vector(r.se)))
  colnames(ds) <- c("theta","theta.se")
  
  if (!missing(labels)) {
    if (length(labels) != length(r))
      stop("The vectors 'labels' and 'r' have different lengths!")
  } else {
    labels <- paste("Study",seq(1, length(r)))
  }
  rownames(ds) = labels
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  quantiles <- c((1-pars.default$level)/2, 0.50, (1-((1-pars.default$level)/2)))
  
  out <- list()
  out$call <- match.call()
  out$method <- method
  class(out) <- "uvmeta"
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies <- dim(ds)[1]
  dfr <- numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  
  if (method != "BAYES") { 
    fit <- rma(yi=r, sei=r.se, method=method, test=test, slab=labels, ...) 
    preds <- predict(fit)
    cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
    cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
    results <- c(coefficients(fit), sqrt(vcov(fit)), fit$tau2, fit$se.tau2, c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub))
    names(results) <- c("estimate", "SE", "tau2", "se.tau2", "95CIl", "95CIu", "95PIl", "95PIu")
    
    out$rma <- fit
    out$numstudies <- fit$k
    out$results <- results
    
  } else if (method == "BAYES") { 
    results.overview <- as.data.frame(array(NA,dim=c(2, length(quantiles)+2)))
    colnames(results.overview) <- c("Estimate","Var",paste(quantiles*100,"%",sep=""))
    rownames(results.overview) <- c("mu", "tausq")
    
    modelfile <- system.file(package="metamisc", "model", "uvmeta_ranef.bug")
    uvmeta_dat <- list('r' = ds$theta,
                       'vars' = ds$theta.se**2,
                       'k' = numstudies,
                       'hp.mu.mean' = pars.default$hp.mu.mean,
                       'hp.mu.prec' = 1/pars.default$hp.mu.var)
    
    model.pars <- list()
    model.pars[[1]] <- list(param="mu", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
    model.pars[[2]] <- list(param="tau", param.f=runif, param.args=list(n=1, min=0, max=100))
    inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
    
    jags.model <- runjags::run.jags(model=modelfile, 
                                    monitor = c("mu", "tausq", "theta.new", "PED"), 
                                    data = uvmeta_dat, 
                                    silent.jags = !verbose,
                                    inits=inits,
                                    ...)
    fit <- jags.model$summaries
    
    #Extract PED
    fit.dev <- runjags::extract(jags.model,"PED")
    
    results <- c(fit["mu",c("Mean","SD")], fit["tausq",c("Mean","SD")], fit["mu", c("Lower95", "Upper95")], fit["theta.new", c("Lower95", "Upper95")])
    names(results) <- c("estimate", "SE", "tau2", "se.tau2", "95CIl", "95CIu", "95PIl", "95PIu")

    out$runjags <- jags.model
    out$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
    out$results <- results
  }
  #attr(out$results,"level") <- pars.default$level
  out$data <- ds
  out$numstudies <- dim(ds)[1]
  out$na.action <- na.action
  return(out)
}

plot.uvmeta <- function(x, ...) {
  level <- 0.95 #attr(x$results,"level")
  quantiles <- c((1-level)/2, (1-((1-level)/2)))

  ci <- x$data[,"theta"]+t(qnorm(quantiles)*matrix(rep(x$data[,"theta.se"],length(quantiles)),nrow=(length(quantiles)), ncol=dim(x$data)[1],byrow=T))
  
  xlim <- c(min(ci),max(ci))
  ylim <- c(2,(x$numstudies+5))

  par(mfrow=c(1,1), mar=( c(5, 12, 4, 4) + 0.1))
  
  loc = c((x$numstudies+4):3)
    
  lcol = "#EBEBEB"
  plot(-500,-500,xlim=xlim,ylim=ylim,xlab="", ylab="",yaxt="n",...)
  axis(2,at=c((x$numstudies+4):5,3),labels=c(rownames(x$data),"Pooled estimate"),las=1)
  
  abline(v=0.00,lty=2,col=lcol)
  for (i in 1:x$numstudies) {
    yloc = loc[i]
    points(x$data[i,"theta"],yloc,pch=15)
    lines(c(min(ci[i,]),max(ci[i,])),c(yloc,yloc))
    for (j in 1:dim(ci)[2])
      lines(c(ci[i,j],ci[i,j]),c((yloc-0.1),(yloc+0.1)),pch=3)
  }
  
  ci.bounds <- x$results[c("95CIl", "95CIu")]
  lines(rep(x$results["estimate"],2), c(2.9,3.1),pch=3)
  lines(c(min(ci.bounds),max(ci.bounds)),c(3,3))
  points(x$results["estimate"],3,pch=23,bg="white")
  
  #Add prediction interval
  if (x$method != "FE") {
    pi.bounds <- x$results[c("95PIl", "95PIu")]
    lines(c(min(pi.bounds),min(ci.bounds)),c(3,3), lty=2)
    lines(c(max(ci.bounds),max(pi.bounds)),c(3,3), lty=2)
    lines(c(min(ci.bounds),min(ci.bounds)),c((3-0.2),(3+0.2)),pch=3)
    lines(c(max(ci.bounds),max(ci.bounds)),c((3-0.2),(3+0.2)),pch=3)
    lines(c(min(pi.bounds),min(pi.bounds)),c((3-0.1),(3+0.1)),pch=3)
    lines(c(max(pi.bounds),max(pi.bounds)),c((3-0.1),(3+0.1)),pch=3)
  }
  
  box()
}



print.uvmeta <- function(x, ...)
{
  out <- (x$results)[c("estimate", "95CIl", "95CIu")]
  text.model <- if (x$method=="FE") "Fixed" else "Random"
  text.method <- if(x$method=="BAYES") "credibility" else "confidence"
  cat(paste(text.model,"effects estimates with corresponding", text.method, "intervals:\n\n"))
	print(out)
  if (x$method!="FE") {
    cat(paste("\n\nPrediction interval for mu:\n\n"))
    print((x$results)[c("estimate", "95PIl", "95PIu")])
  }
  if (x$method=="BAYES") {
    cat(paste("\nPenalized expected deviance: ", round(x$PED,3),"\n"))
    
    # Check if model converged
    psrf.ul <-  x$runjags$psrf$psrf[,"Upper C.I."]
    psrf.target <- x$runjags$psrf$psrf.target
    
    if(sum(psrf.ul > psrf.target)>1) {
      warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                    psrf.target, "for the parameters", 
                    paste(rownames(x$runjags$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                          round(x$runjags$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
                    ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'.", sep=""))
    }
  }
  
	out
}


summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$method=="FE")  cat(paste("\nFixed effects summary:\t",round(object$results["estimate"],5))," (SE: ",round(object$results["SE"],5), ")",sep="")
    if (object$method!="FE") {
        cat(paste("\nRandom effects summary:\t",round(object$results["estimate"],5))," (SE: ",round(object$results["SE"],5), ")",sep="")
        cat(paste("\nTau squared: \t\t",round(object$results["tau2"],5))," (SE: ",round(object$results["se.tau2"],5), ")",sep="")
    }
}



