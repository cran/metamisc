#Multivariate meta-analyse: http://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix/ (directly take log-determinant)
#TODO: allow data transformations

#' Univariate meta-analysis
#' 
#' This function summarizes multiple estimates for a single parameter by assuming a fixed (i.e. common) 
#' effect or random effects across studies. The summary estimate is obtained by calculating a weighted mean that accounts for
#' sample size and (in case random effects are assumed) for between-study heterogeneity. 
#' 
#' @param r Vector of numerics containing the effect size of each study
#' @param r.se Vector of numerics containing the standard error of the effect sizes
#' @param method Character string specifying whether a fixed-effect or a random-effects model should be fitted. 
#' A fixed-effect model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal 
#' to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"EB"}, 
#' \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.
#' @param test Optional character string when \code{method!="BAYES"} to specify how test statistics and confidence intervals 
#' for the fixed effects should be computed. By default (\code{test="knha"}), the method by Knapp and Hartung (2003) is used for adjusting test statistics 
#' and confidence intervals.  Type '\code{?rma}' for more details.
#' @param labels Vector of characters containing the labels for the studies
#' @param na.action A function which indicates what should happen when the data contain NAs. 
#' Defaults to \code{"na.fail"}, other options are \code{"na.omit"}, \code{"na.exclude"} or \code{"na.pass"}.
#' @param n.chains Optional numeric specifying the number of chains to use in the Gibbs sampler (\code{method="BAYES"}). 
#' More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly. 
#' The default number of chains is 4.
#' @param pars Optional list with additional arguments. The width of confidence, credibility and prediction intervals is 
#' defined by \code{level} (defaults to 0.95). 
#' The following parameters configure the MCMC sampling procedure:  
#' \code{hp.mu.mean} (Hyperparameter: mean of the prior distribution of the fixed/random effects model, defaults to zero), 
#' \code{hp.mu.var} (Hyperparameter: variance of the prior distribution of the fixed/random effects model, defaults to 1000),
#' \code{hp.tau.min} (Hyperparameter: mininum value for the between-study standard deviation, defaults to 0),
#' \code{hp.tau.max} (Hyperparameter: maximum value for the between-study standard deviation, defaults to 100).
#' @param verbose If TRUE then messages generated during the fitting process will be displayed.
#' @param \dots Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).
#' 
#' @details Unless specified otherwise, all meta-analysis models assume random effects and are fitted  using restricted 
#' maximum likelihood estimation with the \pkg{metafor} package (Viechtbauer 2010).  Further, confidence intervals for 
#' the average performance are based on the Hartung-Knapp-Sidik-Jonkman method, to better account for the uncertainty 
#' in the estimated between-study heterogeneity (Debray 2016). A Bayesian meta-analysis can be performed by specifying 
#' \code{method="BAYES"}. In that case, the R packages \pkg{runjags} and \pkg{rjags} must be installed.]\cr
#' \cr
#' For random-effects models, a prediction interval for the pooled effect size is displayed. This interval predicts in what 
#' range future effect sizes will fall given what has already been observed (Higgins 2009, Riley 2011).  
#' 
#' \subsection{Bayesian meta-analysis models}{
#' For Bayesian meta-analysis models that involve the Gibbs sampler (\code{method="BAYES"}), the R packages \code{runjags} 
#' and \code{rjags} must be installed. The Bayesian approach uses an uninformative Normal prior for the mean and a 
#' uniform prior for the between-study variance of the pooled effect size (Higgins 2009). By default, the Normal prior 
#' has a mean of 0 and a variance of 1000. These hyperparameters can, however, be altered through the 
#' variables \code{hp.mu.mean} and \code{hp.mu.var} in the argument \code{pars}. The prior distribution of the between-study 
#' standard deviation is given by a uniform distribution, by default bounded between 0 and 100. 
#' }
#' 
#' @return An object of the class \code{uvmeta} for which many standard methods are available. If \code{method="BAYES"}, 
#' the results contain an object of class \link[runjags]{runjags}. Otherwise, the results contain an object of 
#' class \link[metafor]{rma}
#' 
#' @references \itemize{
#' \item Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model 
#' in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.
#' \item Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects 
#' models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.
#' \item DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.
#' \item Graham PL, Moran JL. Robust meta-analytic conclusions mandate the provision of prediction intervals in 
#' meta-analysis summaries. \emph{Journal of Clinical Epidemiology} 2012; \bold{65}: 503--510.
#' \item Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; 
#' \bold{21}: 1539--1558.
#' \item Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. 
#' \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} 2009, \bold{172}: 137--159.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. 
#' \emph{British Medical Journal} 2011; \bold{342}: d549.
#' \item Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 
#' 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}
#' }
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @examples 
#' data(Roberts)
#' 
#' # Frequentist random-effects meta-analysis
#' fit1 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts)))
#' summary(fit1)
#' plot(fit1) #show a forest plot
#' fit1
#' 
#' \dontrun{
#' # Bayesian random effects meta-analysis 
#' fit2 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts), method="BAYES"))
#' plot(fit2)
#' }
#' 
#' @keywords univariate fixed-effect random-effects meta-analysis heterogeneity
#' 
#' @import metafor
#' @export

uvmeta <- function(r, r.se, method="REML", test="knha", labels, na.action, 
                   n.chains=4, pars, verbose=FALSE, ...) 
  UseMethod("uvmeta")

#' @export
uvmeta.default <- function(r, r.se, method="REML", test="knha", labels, na.action, 
                           n.chains=4, pars, verbose=FALSE, ...)
{
  out <- list()
  out$call <- match.call()
  out$method <- method
  out$test <- test
  class(out) <- "uvmeta"

  
  pars.default <- list(level = 0.95,
                       hp.mu.mean = 0, 
                       hp.mu.var = 1000,
                       hp.tau.min = 0,
                       hp.tau.max = 100) 
  
  # Load default parameters
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  
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
    if (pars.default$hp.tau.min < 0) {
      stop("Invalid value for hyperparameter 'hp.tau.min")
    }
    if (pars.default$hp.tau.max < pars.default$hp.tau.min) {
      stop("Invalid value for hyperparameter 'hp.tau.max")
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
    
    out$slab <- make.unique(as.character(labels))
  } else {
    out$slab <- paste("Study",seq(1, length(r)))
  }
  rownames(ds) = out$slab
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  
  
  # Define quantiles for calculating intervals
  out$level <- pars.default$level
  quantiles <- c((1-pars.default$level)/2, 0.50, (1-((1-pars.default$level)/2)))
  
  if (out$level < 0 | out$level > 1) {
    stop ("Invalid value for 'level'!")
  } 
  
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies <- dim(ds)[1]
  dfr <- numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  
  if (method != "BAYES") { 
    fit <- rma(yi=r, sei=r.se, method=method, test=test, slab=out$slab, ...) 
    preds <- predict(fit, level=pars.default$level)
    cr.lb <- ifelse(method=="FE", NA, preds$cr.lb)
    cr.ub <- ifelse(method=="FE", NA, preds$cr.ub)
    results <- c(coefficients(fit), sqrt(vcov(fit)), fit$tau2, fit$se.tau2, c(preds$ci.lb, preds$ci.ub, cr.lb, cr.ub))
    names(results) <- c("estimate", "SE", "tau2", "se.tau2", "CIl", "CIu", "PIl", "PIu")
    
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
                       'hp.mu.prec' = 1/pars.default$hp.mu.var,
                       'hp.tau.min' = pars.default$hp.tau.min,
                       'hp.tau.max' = pars.default$hp.tau.max)
    
    model.pars <- list()
    model.pars[[1]] <- list(param="mu", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
    model.pars[[2]] <- list(param="tau", param.f=runif, param.args=list(n=1, min=pars.default$hp.tau.min, max=pars.default$hp.tau.max))
    inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
    
    jags.model <- runjags::run.jags(model=modelfile, 
                                    monitor = c("mu", "tausq", "theta.new", "PED"), 
                                    data = uvmeta_dat, 
                                    silent.jags = !verbose,
                                    inits=inits,
                                    confidence=out$level,
                                    ...)
    fit <- jags.model$summaries
    
    #Extract PED
    fit.dev <- runjags::extract(jags.model,"PED")
    
    txtLevel <- (out$level*100)
    
    results <- c(fit["mu",c("Mean","SD")], fit["tausq",c("Mean","SD")], fit["mu", paste(c("Lower", "Upper"), txtLevel, sep="")], fit["theta.new",  paste(c("Lower", "Upper"), txtLevel, sep="")])
    names(results) <- c("estimate", "SE", "tau2", "se.tau2", "CIl", "CIu", "PIl", "PIu")

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

#' Forest Plots
#' 
#' Function to create forest plots for objects of class \code{"uvmeta"}.
#' 
#' @param x An object of class \code{"uvmeta"}
#' @param sort By default, studies are ordered by ascending effect size (\code{sort="asc"}). For study ordering by descending
#' effect size, choose \code{sort="desc"}. For any other value, study ordering is ignored.
#' @param \dots Additional arguments which are passed to \link{forest}.
#' 
#' @details The forest plot shows the performance estimates of each validation with corresponding confidence 
#' intervals. A polygon is added to the bottom of the forest plot, showing the summary estimate based on the model. 
#' A 95\% prediction interval is added by default for random-effects models,  the dotted line indicates its (approximate) bounds
#' 
#' @note Full lines indicate confidence intervals or credibility intervals (in case of a Bayesian meta-analysis). Dashed
#' lines indicate prediction intervals. The width of all intervals is defined by the significance level chosen during 
#' meta-analysis. 
#' 
#' @references \itemize{
#' \item Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{BMJ}. 2011 342:d549--d549.
#' }
#' 
#' @examples 
#' data(Roberts)
#' 
#' # Frequentist random-effects meta-analysis
#' fit <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts)))
#' plot(fit) 
#' 
#' @keywords meta-analysis forest
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @importFrom graphics par plot axis polygon points lines abline
#' @import ggplot2
#' 
#' @method plot uvmeta
#' @export
plot.uvmeta <- function(x, sort="asc", ...) {
  level <- x$level
  quantiles <- c((1-level)/2, (1-((1-level)/2)))
  
  # Reconstruct confidence intervals for study effects with the level used during calculation of summary effect
  yi.ci <- x$data[,"theta"]+t(qnorm(quantiles)*matrix(rep(x$data[,"theta.se"],length(quantiles)),nrow=(length(quantiles)), ncol=dim(x$data)[1],byrow=T))
  
  #Extract data
  yi <- c(x$data[,"theta"])
  yi.slab <- c(as.character(x$slab))
  
  forest(theta=yi, theta.ci=yi.ci, theta.slab=yi.slab, 
         theta.summary=x$results["estimate"], 
         theta.summary.ci=x$results[c("CIl","CIu")], 
         theta.summary.pi=x$results[c("PIl","PIu")], 
         sort=sort, ...)
}


#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print uvmeta
#' @export
print.uvmeta <- function(x, ...)
{
  out <- (x$results)[c("estimate", "CIl", "CIu")]
  
  text.model <- if (x$method=="FE") "Fixed" else "Random"
  text.ci <- if(x$method=="BAYES") "credibility" else "confidence"
  text.pi <- if(x$method=="BAYES") "" else "(approximate)"
  
  if (x$method!="FE") {
    cat(paste("Summary estimate with ", x$level*100, "% ", text.ci, " and ", text.pi, " ", 
              x$level*100, "% prediction interval:\n\n", sep=""))
    print((x$results)[c("estimate", "CIl", "CIu", "PIl", "PIu")])
  } else {
    cat(paste("Summary estimate with ", x$level*100, "% ", text.ci, ":\n\n", sep=""))
    print((x$results)[c("estimate", "CIl", "CIu")])
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

#' Summarizing Univariate Meta-Analysis Models
#' 
#' This function provides summary estimates of a fitted univariate meta-analysis model.
#' 
#' @method summary uvmeta
#' 
#' @param  object  An object of class \code{"uvmeta"}
#' @param \dots Optional arguments to be passed on to other functions
#' 
#' @references  \itemize{
#' \item Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.
#' \item DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.
#' }
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @seealso \code{\link{uvmeta}}
#' 
#' @export
#' @keywords DerSimonian  Laird univariate random-effects meta-analysis 
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



