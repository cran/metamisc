#' Meta-analysis of prediction model performance
#'
#' This function provides summary estimates for the concordance statistic, the total observed-expected ratio 
#' or the calibration slope. Where appropriate, data transformations are applied and missing information 
#' is derived from available quantities. Unless specified otherwise, all meta-analysis models assume random effects 
#' and are fitted using restricted maximum likelihood estimation with the \pkg{metafor} package (Viechtbauer 2010).  
#' Further, confidence intervals for the average performance are based on the Hartung-Knapp-Sidik-Jonkman method. 
#' When conducting a Bayesian meta-analysis, the R packages \pkg{runjags} and \pkg{rjags} must be installed.
#' 
#' @param  measure A character string indicating which summary performance measure should be calculated. Options are
#' \code{"cstat"} (meta-analysis of the concordance statistic) and \code{"OE"} 
#' (meta-analysis of the total observed-expected ratio). See `Details' for more information.
#' @param cstat Optional vector with the estimated c-statistic for each valiation
#' @param cstat.se Optional vector with the standard error of the estimated c-statistics
#' @param cstat.95CI Optional 2-dimensional array with the lower (first column) and upper (second column) boundary 
#' of the 95\% confidence interval of the estimated c-statistics
#' @param sd.LP Optional vector with the standard deviation of the linear predictor (prognostic index)
#' @param OE Optional vector with the estimated ratio of total observed versus total expected events
#' @param OE.se Optional vector with the standard errors of the estimated O:E ratios
#' @param OE.95CI Optional 2-dimensional array with the lower (first column) and upper (second column) boundary 
#' of the 95\% confidence interval of the total O:E ratios
#' @param citl Optional vector with the estimated calibration-in-the-large for each valiation
#' @param citl.se Optional vector with the standard error of the estimated calibration-in-the-large statistics
#' @param N Optional vector with the total number of participants for each valiation
#' @param O Optional vector with the total number of observed events for each valiation
#' (if specified, during time \code{t.val})
#' @param E Optional vector with the total number of expected events for each valiation 
#' (if specified, during time \code{t.val})
#' @param Po Optional vector with the (cumulative) observed event probability for each valiation
#' (if specified, during time \code{t.val})
#' @param Po.se Optional vector with the standard errors of \code{Po}.
#' @param Pe Optional vector with the (cumulative) expected event probability for each validation
#' (if specified, during time \code{t.val})
#' @param t.val Optional vector specifying the time period for which \code{cstat}, \code{O}, \code{E}, \code{Po} or
#' \code{Pe} are applicable. Also specifies the time point at which \code{OE} and \code{CITL} have been calculated.
#' @param t.ma Optional numeric value, specifying the target time period (of time point) of the meta-analysis
#' @param t.extrapolate Optional logical indicating whether calibration performance of the prognostic model 
#' should be extrapolated to time \code{t.ma}. Otherwise, studies where \code{t.val!=t.ma} will be omitted
#' when meta-analysing the total O:E ratio. 
#' @param method Character string specifying whether a fixed- or a random-effects model should be fitted. 
#' A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method 
#' equal to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, 
#' \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.
#' @param test Optional character string specifying how test statistics and confidence intervals for the fixed effects 
#' should be computed. By default (\code{test="knha"}), the method by Knapp and Hartung (2003) is used for 
#' adjusting test statistics and confidence intervals. Type '\code{?rma}' for more details.
#' @param verbose If TRUE then messages generated during the fitting process will be displayed.
#' @param slab Optional vector specifying the label for each study
#' @param n.chains Optional numeric specifying the number of chains to use in the Gibbs sampler 
#' (if \code{method="BAYES"}). More chains will improve the sensitivity of the convergence diagnostic, but will 
#' cause the simulation to run more slowly. The default number of chains is 4.
#' @param pars A list with additional arguments.  See 'Details' for more information. The following parameters configure the MCMC sampling procedure:  
#' \code{hp.mu.mean} (mean of the prior distribution of the random effects model, defaults to 0), 
#' \code{hp.mu.var} (variance of the prior distribution of the random effects model, defaults to 1E6), 
#' \code{hp.tau.min} (minimum value for the between-study standard deviation, defaults to 0), 
#' \code{hp.tau.max} (maximum value for the between-study standard deviation, defaults to 2), 
#' \code{hp.tau.sigma} (standard deviation of the prior distribution for the between-study standard-deviation), 
#' \code{hp.tau.dist} (prior distribution for the between-study standard-deviation. Defaults to \code{"dunif"}), 
#' \code{hp.tau.df} (degrees of freedom for the prior distribution for the between-study standard-deviation. 
#' Defaults to 3). Other arguments are \code{method.restore.c.se} (method for restoring missing estimates for the standard error 
#' of the c-statistic. See \code{\link{ccalc}} for more information), \code{model.cstat} (The likelihood/link for modeling 
#' the c-statistic; see "Details"), \code{model.oe} (The likelihood/link for modeling the O:E ratio; see "Details")
#' @param \ldots Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).
#' 
#' @details 
#' \subsection{Meta-analysis of the concordance statistic}{
#' A summary estimate for the concorcance (c-) statistic can be obtained by specifying \code{measure="cstat"}.
#' The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to 
#' distinguish between patients developing and not developing the outcome. The c-statistic typically ranges 
#' from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). When missing, the c-statistic 
#' and/or its standard error are derived from other reported information. 
#' See \code{\link{ccalc}} for more information.
#' 
#' By default, it is assumed that the logit of the c-statistic is Normally distributed within and across studies
#' (\code{pars$model.cstat = "normal/logit"}). Alternatively, it is possible to assume that the raw c-statistic 
#' is Normally distributed across studies \code{pars$model.cstat = "normal/identity"}. 
#' } 
#' 
#' \subsection{Meta-analysis of the total observed versus expected ratio}{
#' A summary estimate for the total observed versus expected (O:E) ratio can be obtained by specifying
#' \code{measure="OE"}. The total O:E ratio provides a rough indication of the overall model calibration (across the 
#' entire range of predicted risks). 
#' 
#' For frequentist meta-analysis, within-study variation can either be modeled using a Normal (\code{model.oe = "normal/log"} 
#' or \code{model.oe = "normal/identity"}) or a Poisson distribution (\code{model.oe = "normal/log"}). 
#' 
#' When performing a Bayesian meta-analysis, all data are modeled using a one-stage random effects (hierarchical related regression) model.
#' In particular, a binomial distribution (if \code{O}, \code{E} and \code{N} is known), a poisson distribution 
#' (if only \code{O} and \code{E} are known) or a Normal distribution (if \code{OE} and \code{OE.se} or \code{OE.95CI} are known) is selected separately for each study.
#' 
#' For meta-analysis of prognostic models, it is recommended to provide information on the time period 
#' (\code{t.val}) during which calibration was assessed in the validation study. When the time period of 
#' the validation study does not correspond to the time period of interest (\code{t.ma}), corresponding 
#' studies are omitted from the meta-analysis. It is possible to extrapolate observed
#' event rates by setting \code{t.extrapolate=TRUE}. This approach is currently only supported for
#' frequentist meta-analysis models with \code{model.oe = "normal/log"} or \code{model.oe = "normal/identity"}. 
#' }
#' 
#' \subsection{Bayesian meta-analysis}{
#' All Bayesian meta-analysis models assume random effects by default. Results are based on the posterior median. 
#' Credibility and prediction intervals are directly obtained from the corresponding posterior quantiles.
#' 
#' The prior distribution for the (transformed) summary estimate is always modeled using a Normal distribution, 
#' with mean \code{hp.mu.mean} (defaults to 0) and variance 
#' \code{hp.mu.var} (defaults to 1E6). For meta-analysis of the total O:E ratio, the maximum value for \code{hp.mu.var} is 100.
#' 
#' By default, the prior distribution for the between-study standard deviation is modeled using a uniform distribution 
#' (\code{hp.tau.dist="dunif"}), with boundaries \code{hp.tau.min} and \code{hp.tau.max}. Alternatively, it is possible
#' to specify a truncated Student-t distribution (\code{hp.tau.dist="dhalft"}) with a mean of \code{hp.tau.mean}, 
#' a standard deviation of \code{hp.tau.sigma} and \code{hp.tau.df} degrees of freedom. This distribution is again 
#' restricted to the range \code{hp.tau.min} to \code{hp.tau.max}.
#' }
#' 
#' @note The width of calculated confidence, credibility and prediction intervals can be specified 
#' using \code{level} in the \code{pars} argument (defaults to 0.95).
#' 
#' @return An object of class \code{valmeta} with the following elements:
#' \describe{
##'  \item{"data"}{array with (transformed) data used for meta-analysis, and method(s) used for restoring missing information. }
##'  \item{"lme4"}{a fitted object of class \code{glmerMod} (if \pkg{lme4} was used for meta-analysis).}
##'  \item{"measure"}{character string specifying the performance measure that has been meta-analysed.}
##'  \item{"method"}{character string specifying the meta-analysis method.}
##'  \item{"model"}{character string specifying the meta-analysis model (link function).}
##'  \item{"results"}{numeric vector containing the meta-analysis results}
##'  \item{"rma"}{a fitted object of class \link[metafor]{rma} (if \pkg{metafor} was used for meta-analysis).}
##'  \item{"runjags"}{a fitted object of class \code{runjags} (if \pkg{runjags} was used for meta-analysis).}
##'  \item{"slab"}{vector specifying the label of each study.}
##' }
#' @references 
#' \itemize{
#' \item Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review 
#' and meta-analysis of prediction model performance. \emph{BMJ}. 2017; 356:i6460.
#' \item Stijnen T, Hamza TH, Ozdemir P. Random effects meta-analysis of event outcome in the framework of 
#' the generalized linear mixed model with applications in sparse data. \emph{Stat Med}. 2010; 29(29):3046--67.
#' \item Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 
#' 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}
#' }
#'   
#' @seealso \code{\link{ccalc}} to calculate concordance statistics and corresponding standard errors, 
#' \code{\link{plot.valmeta}} to generate forest plots
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Meta-analysis of the c-statistic (random effects)
#' fit <- with(EuroSCORE, valmeta(cstat=c.index, cstat.se=se.c.index, 
#'                                cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
#'                                N=n, O=n.events, slab=Study))
#' plot(fit)
#' 
#' # Nearly identical results when we need to estimate the SE
#' with(EuroSCORE, valmeta(cstat=c.index,  N=n, O=n.events, slab=Study))
#' 
#' # Two-stage meta-analysis of the total O:E ratio (random effects)
#' with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events, N=n))    
#' with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events))        
#' with(EuroSCORE, valmeta(measure="OE", Po=Po, Pe=Pe, N=n))
#' 
#' # One-stage meta-analysis of the total O:E ratio (random effects)
#' with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events, method="ML",
#'                         pars=list(model.oe="poisson/log")))
#' 
#' \dontrun{
#' # Bayesian random effects meta-analysis of the c-statistic
#' fit2 <- with(EuroSCORE, valmeta(cstat=c.index, cstat.se=se.c.index, 
#'                                 cstat.95CI=cbind(c.index.95CIl,c.index.95CIu),
#'                                 N=n, O=n.events, method="BAYES", slab=Study))
#' 
#' # Bayesian one-stage random effects meta-analysis of the total O:E ratio
#' # Consider that some (but not all) studies do not provide information on N
#' # A Poisson distribution will be used for studies 1, 2, 5, 10 and 20
#' # A Binomial distribution will be used for the remaining studies
#' EuroSCORE.new <- EuroSCORE
#' EuroSCORE.new$n[c(1, 2, 5, 10, 20)] <-  NA
#' pars <- list(hp.tau.dist="dhalft",   # Prior for the between-study standard deviation
#'              hp.tau.sigma=1.5,       # Standard deviation for 'hp.tau.dist'
#'              hp.tau.df=3,            # Degrees of freedom for 'hp.tau.dist'
#'              hp.tau.max=10)          # Maximum value for the between-study standard deviation
#' fit3 <- with(EuroSCORE.new, valmeta(measure="OE", O=n.events, E=e.events, N=n, 
#'         method="BAYES", slab=Study, pars=pars))
#' plot(fit3)
#' print(fit3$runjags$model) # Inspect the JAGS model
#' print(fit3$runjags$data)  # Inspect the JAGS data
#' } 
#' 
#' ######### Validation of prediction models with a time-to-event outcome #########
#' data(Framingham)
#' 
#' # Meta-analysis of total O:E ratio after 10 years of follow-up
#' with(Framingham, valmeta(measure="OE", Po=Po, Pe=Pe, N=n, t.val=t.val, t.ma=10))
#' with(Framingham, valmeta(measure="OE", Po=Po, Pe=Pe, N=n, t.val=t.val, t.ma=10, t.extrapolate=TRUE))
#' 
#' @keywords meta-analysis discrimination  calibration
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' @import metafor
#' @import mvtnorm
#' @importFrom lme4 glmer
#' @importFrom stats coef coefficients dnorm glm nobs optim pchisq qnorm qt pt rnorm runif confint poisson
#' predict vcov as.formula formula model.frame model.frame.default update.formula family

valmeta <- function(measure="cstat", cstat, cstat.se, cstat.95CI, sd.LP, OE, OE.se, OE.95CI, citl, citl.se,
                    N, O, E, Po, Po.se, Pe, t.val, t.ma, t.extrapolate=FALSE, method="REML", test="knha", 
                    verbose=FALSE, slab, n.chains = 4, pars, ...) {
  pars.default <- list(level = 0.95,
                       hp.mu.mean = 0, 
                       hp.mu.var = 1E6,
                       hp.tau.min = 0,
                       hp.tau.max = 2,
                       hp.tau.mean = 0,
                       hp.tau.sigma = 0.5,
                       hp.tau.dist = "dunif", 
                       hp.tau.df = 3, 
                       correction = 0.5,
                       method.restore.c.se=4,
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
  
  if (pars.default$level < 0 | pars.default$level > 1) {
    stop ("Invalid value for 'level'!")
  } 
  
  if (measure=="OE" & pars.default$model.oe=="poisson/log" & t.extrapolate) {
    t.extrapolate <- FALSE
    warning("Extrapolation not implemented yet for poisson/log models!")
  }else if (measure=="OE" & method=="BAYES" & t.extrapolate) {
    t.extrapolate <- FALSE
    warning("Extrapolation not implemented yet for Bayesian models!")
  }
  
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
    } else if (!missing(sd.LP)) {
      k <- length(sd.LP)
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
  if (missing(cstat.se)) {
    cstat.se <- rep(NA, length=k)
  }
  if (missing(O)) {
    O <- rep(NA, length=k)
  }
  if (missing(Po)) {
    Po <- rep(NA, length=k)
  }
  if (missing(N)) {
    N <- rep(NA, length=k)
  }
  if (missing(sd.LP)) {
    sd.LP <- rep(NA, length=k)
  }
  
  #######################################################################################
  # Prepare results
  #######################################################################################
  results <- data.frame(estimate=NA, CIl=NA, CIu=NA, PIl=NA, PIu=NA)
  
  #######################################################################################
  # Prepare object
  #######################################################################################
  out <- list()
  out$call <- match.call()
  out$measure <- measure
  out$method <- method
  out$level <- pars.default$level 
  class(out) <- "valmeta"
  
  
  #######################################################################################
  # Assign study labels
  #######################################################################################
  if(missing(slab)) {
    out$slab <- paste("Study", seq(1, k))
  } else {
    out$slab <- make.unique(as.character(slab))
  }

  #######################################################################################
  # Meta-analysis of the c-statistic
  #######################################################################################
  if (measure=="cstat") {
    out$model <- pars.default$model.cstat
    
    if (out$model == "normal/identity") {
      g <- NULL
    } else if (out$model=="normal/logit") {
      g <- "log(cstat/(1-cstat))"
    } else {
      stop (paste("Meta-analysis model currently not supported: '", out$model, '"', sep=""))
    }
      
    ds <- ccalc(cstat=cstat, cstat.se=cstat.se, cstat.CI.lower=cstat.95CI[,1],
                cstat.CI.upper=cstat.95CI[,2], sd.LP=sd.LP, 
                N=N, O=O, Po=Po, slab=out$slab, g=g, level=pars.default$level, 
                approx.se.method=pars.default$method.restore.c.se) 
    
    
    if (method != "BAYES") { # Use of rma
      
      # Identify which studies can be used for meta-analysis
      selstudies <- which(!is.na(ds$theta) & !is.na(ds$theta.se))
      
      # Apply the meta-analysis
      fit <- rma(yi=theta, sei=theta.se, data=ds, method=method, test=test, slab=out$slab, ...) 
      preds <- predict(fit, level=pars.default$level)
      
      ds[selstudies, "theta.blup"] <- blup(fit)$pred
      
      if (pars.default$model.cstat == "normal/logit") {
        results$estimate <- inv.logit(coefficients(fit))
        results$CIl <- inv.logit(preds$ci.lb)
        results$CIu <- inv.logit(preds$ci.ub)
        results$PIl <- ifelse(method=="FE", NA, inv.logit(preds$cr.lb))
        results$PIu <- ifelse(method=="FE", NA, inv.logit(preds$cr.ub))
      } else if (pars.default$model.cstat == "normal/identity") {
        results$estimate <- coefficients(fit)
        results$CIl <- preds$ci.lb
        results$CIu <- preds$ci.ub
        results$PIl <- ifelse(method=="FE", NA, preds$cr.lb)
        results$PIu <- ifelse(method=="FE", NA, preds$cr.ub)
      } else {
        stop ("There is no implementation for the specified meta-analysis model!")
      }
      
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
      model.pars[[1]] <- list(param="mu.tobs", param.f=rnorm, 
                              param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      
      if (pars.default$hp.tau.dist=="dunif") {
        model.pars[[2]] <- list(param="bsTau", param.f=runif, 
                                param.args=list(n=1, min=pars.default$hp.tau.min, 
                                                max=pars.default$hp.tau.max))
      } else if (pars.default$hp.tau.dist=="dhalft") {
        model.pars[[2]] <- list(param="bsTau", param.f=rstudentt, 
                                param.args=list(n=1, mean=pars.default$hp.tau.mean, 
                                                sigma=pars.default$hp.tau.sigma, 
                                                df=pars.default$hp.tau.df,
                                                lower=pars.default$hp.tau.min, 
                                                upper=pars.default$hp.tau.max))
      } else {
        stop("Invalid distribution for 'hp.tau.dist'!")
      }

      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      mvmeta_dat <- list(theta = ds$theta,
                         theta.var = ds$theta.se**2,
                         Nstudies = length(ds$theta))
      jags.model <- runjags::run.jags(model=model, 
                                      monitor = c("mu.tobs", "mu.obs", "pred.obs", "bsTau", "PED"), 
                                      data = mvmeta_dat, 
                                      confidence = out$level, # Which credibility intervals do we need?
                                      n.chains = n.chains,
                                      silent.jags = !verbose,
                                      inits=inits,
                                      ...)
      fit <- jags.model$summaries
      
      
      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      txtLevel <- (out$level*100)
      
      results$estimate  <- fit["mu.obs", "Median"]
      results$CIl       <- fit["mu.obs", paste("Lower", txtLevel, sep="")]
      results$CIu       <- fit["mu.obs", paste("Upper", txtLevel, sep="")]
      results$PIl       <- fit["pred.obs", paste("Lower", txtLevel, sep="")]
      results$PIu       <- fit["pred.obs", paste("Upper", txtLevel, sep="")]
      
      out$runjags <- jags.model
      out$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$results <- results
    }
    
    out$data <- ds
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
    
    
    #####################################################################################
    # Calculate the OE ratio when model!=poisson/log
    # Omit studies where t.val != t.ma (unless extrapolation is required)
    #####################################################################################
    t.O.E.N   <- restore.oe.O.E.N(O=O, E=E, N=N, correction = pars.default$correction, 
                               t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.O.Pe.N  <- restore.oe.OPeN(O=O, Pe=Pe, N=N, correction = pars.default$correction, 
                                t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe)
    t.E.Po.N  <- restore.oe.EPoN(E=E, Po=Po, N=N, correction = pars.default$correction, 
                                t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe)
    T.Po.Pe.N <- restore.oe.PoPeN(Po=Po, Pe=Pe, N=N, correction = pars.default$correction, 
                                 t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.O.Po.E  <- restore.oe.OPoE(O=O, Po=Po, E=E, correction = pars.default$correction, 
                                t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.O.Pe.E  <- restore.oe.OPeE(O=O, Pe=Pe, E=E, correction = pars.default$correction, 
                                t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.OE.SE   <- restore.oe.OE(OE=OE, OE.se=OE.se, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
                              model=pars.default$model.oe)
    t.OE.CI   <- restore.oe.OE.95CI(OE=OE, OE.95CI=OE.95CI, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
                                   model=pars.default$model.oe)
    t.O.E     <- restore.oe.O.E(O=O, E=E, correction = pars.default$correction, 
                              t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.pope   <- restore.oe.PoPe(Po=Po, Pe=Pe, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=pars.default$model.oe) 
    t.citl   <- restore.oe.citl(citl=citl, citl.se=citl.se, O=O, Po=Po, N=N, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
                                model=pars.default$model.oe) 
    
    # Select appropriate estimate for 'theta' and record its source
    t.method <- c("O, E and N",    "O, Pe and N", "E, Po and N", "Po, Pe and N", "O, E and Po", "O, E and Pe", "OE and SE(OE)", "OE and CI(OE)", "O and E", "Po and Pe", "CITL")
    dat.est  <- cbind(t.O.E.N[,1], t.O.Pe.N[,1], t.E.Po.N[,1], T.Po.Pe.N[,1], t.O.Po.E[,1], t.O.Pe.E[,1], t.OE.SE[,1], t.OE.CI[,1], t.O.E[,1], t.pope[,1],      t.citl[,1]) 
    dat.se   <- cbind(t.O.E.N[,2], t.O.Pe.N[,2], t.E.Po.N[,2], T.Po.Pe.N[,2], t.O.Po.E[,2], t.O.Pe.E[,2], t.OE.SE[,2], t.OE.CI[,2], t.O.E[,2], t.pope[,2],      t.citl[,2]) 
    myfun = function(dat) { which.min(is.na(dat)) }
    sel.theta <- apply(dat.est, 1, myfun)
    theta <- dat.est[cbind(seq_along(sel.theta), sel.theta)] 
    theta.se <- dat.se[cbind(seq_along(sel.theta), sel.theta)] 
    theta.source <-  t.method[sel.theta]
    
    # Define confidence intervals
    theta.cil <- theta.ciu <- rep(NA, k)
    
    if (pars.default$model.oe == "normal/identity") {
      if (pars.default$level==0.95) {
        theta.cil <- OE.95CI[,1]
        theta.ciu <- OE.95CI[,2]
      }
    } else if (pars.default$model.oe == "normal/log" | pars.default$model.oe == "poisson/log") {
      if (pars.default$level==0.95) {
        theta.cil <- log(OE.95CI[,1])
        theta.ciu <- log(OE.95CI[,2])
      }
    } else {
      stop("Undefined link function!")
    }
    theta.cil[is.na(theta.cil)] <- (theta+qnorm((1-pars.default$level)/2)*theta.se)[is.na(theta.cil)]
    theta.ciu[is.na(theta.ciu)] <- (theta+qnorm((1+pars.default$level)/2)*theta.se)[is.na(theta.ciu)]
    
    # Store results, and method for calculating SE
    ds <- data.frame(theta=theta, theta.se=theta.se, theta.CIl=theta.cil, theta.CIu=theta.ciu, theta.source=theta.source)
    
    #####################################################################################
    # Calculate O, E and N for meta-analyses using GLMER models
    # Omit studies where t.val != t.ma 
    #####################################################################################
    if (pars.default$model.oe == "poisson/log" | method == "BAYES") {
      O.new <- O
      E.new <- E
      N.new <- N
      N.new <- ifelse(is.na(N.new), O.new/Po, N.new)
      N.new <- ifelse(is.na(N.new), E.new/Pe, N.new)
      O.new <- ifelse(is.na(O.new), OE*E.new, O.new)
      O.new <- ifelse(is.na(O.new), Po*N.new, O.new)
      E.new <- ifelse(is.na(E.new), O.new/OE, E.new)
      E.new <- ifelse(is.na(E.new), Pe*N.new, E.new)
      
      O.new <- round(O.new) #round O to improve convergence of the models

      Study.new <- c(1:dim(ds)[1])
      ds <- data.frame(Study=Study.new, theta=theta, theta.se=theta.se, theta.CIl=theta.cil, theta.CIu=theta.ciu, 
                       O=O.new, E=E.new, N=N.new, theta.source=theta.source)
      
      # Omit values for O, E and N where t.val!=t.ma
      if (!is.na(t.ma) & class(t.val)=="numeric") {
        ds[which(t.val!=t.ma),c("O","E","N")] <- NA
      } 
    }
    
    out$numstudies <- length(which(rowMeans(!is.na(ds))==1))
      
    if (method != "BAYES") { # Use of rma
      
      if (pars.default$model.oe=="normal/identity") {
        if(verbose) print("Performing two-stage meta-analysis...")
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, test=test, slab=out$slab, ...) 
        preds <- predict(fit, level=pars.default$level)
        
        results$estimate <- coefficients(fit)
        results$CIl      <- preds$ci.lb
        results$CIu      <- preds$ci.ub
        results$PIl      <- ifelse(method=="FE", NA, preds$cr.lb)
        results$PIu      <- ifelse(method=="FE", NA, preds$cr.ub)
        
        out$rma <- fit
        out$numstudies <- fit$k
      } else if (pars.default$model.oe=="normal/log") {
        if(verbose) print("Performing two-stage meta-analysis...")
        fit <- rma(yi=ds$theta, sei=ds$theta.se, data=ds, method=method, test=test, slab=out$slab, ...) 
        preds <- predict(fit, level=pars.default$level)
        
        results$estimate <- exp(coefficients(fit))
        results$CIl      <- exp(preds$ci.lb)
        results$CIu      <- exp(preds$ci.ub)
        results$PIl      <- ifelse(method=="FE", NA, exp(preds$cr.lb))
        results$PIu      <- ifelse(method=="FE", NA, exp(preds$cr.ub))
        
        out$rma <- fit
        out$numstudies <- fit$k
      } else if (pars.default$model.oe=="poisson/log") { 
        if(verbose) print("Performing one-stage meta-analysis...")
        if (method=="ML") { 
          if (test=="knha") warning("The Sidik-Jonkman-Hartung-Knapp correction cannot be applied")
          fit <- glmer(O~1|Study, offset=log(E), family=poisson(link="log"), data=ds)
          preds.ci <- confint(fit, level=pars.default$level, quiet=!verbose, ...)
          
          results$estimate <- exp(lme4::fixef(fit))
          results$CIl      <- exp(preds.ci["(Intercept)",1])
          results$CIu      <- exp(preds.ci["(Intercept)",2])
          results$PIl      <- exp(lme4::fixef(fit) + qt((1-pars.default$level)/2, df=(lme4::ngrps(fit)-2))*sqrt(vcov(fit)[1,1]+(as.data.frame(lme4::VarCorr(fit))["vcov"])[1,1]))
          results$PIu      <- exp(lme4::fixef(fit) + qt((1+pars.default$level)/2, df=(lme4::ngrps(fit)-2))*sqrt(vcov(fit)[1,1]+(as.data.frame(lme4::VarCorr(fit))["vcov"])[1,1]))
          
          out$lme4 <- fit
          out$numstudies <- nobs(fit)
        } else if (method=="FE") { #one-stage fixed-effects meta-analysis
          fit <- glm(O~1, offset=log(E), family=poisson(link="log"), data=ds)
          preds.ci <- confint(fit, level=pars.default$level, quiet=!verbose, ...)
          
          results$estimate <- exp(coefficients(fit))
          results$CIl      <- exp(preds.ci[1])
          results$CIu      <- exp(preds.ci[2])
          results$PIl      <- results$PIu <- NA
          
          out$glm <- fit
          out$numstudies <- nobs(fit)
        } else {
          stop(paste("No implementation found for ", method, " estimation (model '", pars.default$model.oe, "')!", sep=""))
        }
      } else {
        stop("Model not implemented yet!")
      }
      names(results) <- c("estimate", "CIl", "CIu", "PIl", "PIu")
      
      out$results <- results
    } else {
      if(verbose) print("Performing Bayesian one-stage meta-analysis...")
      
      out$model <- "hierarchical related regression"
        
      # Truncate hyper parameter variance
      pars.default$hp.mu.var = min(pars.default$hp.mu.var, 100)
      
      # Select studies where we have info on O, E and N
      i.select1 <- which(!is.na(ds$O) & !is.na(ds$E) & !is.na(ds$N))
      
      # Select studies where we only have info on O and E
      i.select2 <- which(!is.na(ds$O) & !is.na(ds$E) & is.na(ds$N))
      
      # Select studies where we have (estimated) information on log(OE) and its standard error
      i.select3 <- which(!is.na(ds$theta) & !is.na(ds$theta.se) & is.na(ds$O) & is.na(ds$E))
      
      mvmeta_dat <- list(O=ds$O, E=ds$E)
      
      if (length(i.select1)>0) {
        mvmeta_dat$s1 <- i.select1
        mvmeta_dat$N <- ds$N
      }
      if (length(i.select2)>0)
        mvmeta_dat$s2 <- i.select2
      if (length(i.select3)>0) {
        mvmeta_dat$s3 <- i.select3
        mvmeta_dat$logOE <- ds$theta
        mvmeta_dat$logOE.se <- ds$theta.se
      }
      
      # Generate model
      model <- generateBUGS.OE.discrete(N.type1=length(i.select1), 
                                        N.type2=length(i.select2),
                                        N.type3=length(i.select3),
                                        pars=pars.default, ...)
      
      out$numstudies <- length(c(i.select1, i.select2, i.select3))
     
      
      
      # Generate initial values from the relevant distributions
      model.pars <- list()
      model.pars[[1]] <- list(param="mu.logoe", param.f=rnorm, param.args=list(n=1, mean=pars.default$hp.mu.mean, sd=sqrt(pars.default$hp.mu.var)))
      
      if (pars.default$hp.tau.dist=="dunif") {
        model.pars[[2]] <- list(param="bsTau", param.f=runif, 
                                param.args=list(n=1, min=pars.default$hp.tau.min, 
                                                max=pars.default$hp.tau.max))
      } else if (pars.default$hp.tau.dist=="dhalft") {
        model.pars[[2]] <- list(param="bsTau", param.f=rstudentt, 
                                param.args=list(n=1, mean=pars.default$hp.tau.mean, 
                                                sigma=pars.default$hp.tau.sigma, 
                                                df=pars.default$hp.tau.df,
                                                lower=pars.default$hp.tau.min, 
                                                upper=pars.default$hp.tau.max))
      } else {
        stop("Invalid distribution for 'hp.tau.dist'!")
      }
      
      inits <- generateMCMCinits(n.chains=n.chains, model.pars=model.pars)
      
      jags.model <- runjags::run.jags(model=model, 
                                      monitor = c("mu.logoe", "mu.oe", "pred.oe", "bsTau", "PED"), 
                                      data = mvmeta_dat, 
                                      n.chains = n.chains,
                                      confidence = out$level, # Which credibility intervals do we need?
                                      silent.jags = !verbose,
                                      inits=inits,
                                      ...)
        
      fit <- jags.model$summaries
      

      #Extract PED
      fit.dev <- runjags::extract(jags.model,"PED")
      txtLevel <- (out$level*100)
      
      results$estimate <- fit["mu.oe", "Median"]
      results$CIl      <- fit["mu.oe", paste("Lower", txtLevel, sep="")]
      results$CIu      <- fit["mu.oe", paste("Upper", txtLevel, sep="")]
      results$PIl      <- fit["pred.oe", paste("Lower", txtLevel, sep="")]
      results$PIu      <- fit["pred.oe", paste("Upper", txtLevel, sep="")]
      
      out$runjags <- jags.model
      out$PED <- sum(fit.dev$deviance)+sum(fit.dev$penalty)
      out$results <- results
    }
    
    if ("Study" %in% colnames(ds)) {
      ds <- ds[,-which(colnames(ds)=="Study")]
    }
    
    out$data <- ds
    
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
    out <- paste(out, "  bsTau ~ dt(", pars$hp.tau.mean," ,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep="") 
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

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print valmeta
#' @export
print.valmeta <- function(x, ...) {
  if (x$measure=="cstat") {
    text.stat <- "c-statistic"
  } else if (x$measure=="OE") {
    text.stat <- "O:E ratio"
  }
  
  text.model <- if (x$method=="FE") "Fixed" else "Random"
  text.ci <- if(x$method=="BAYES") "credibility" else "confidence"
  text.pi <- if(x$method=="BAYES") "" else "(approximate)"
  
  
  if (x$method!="FE") {
    cat(paste("Summary ", text.stat, " with ", x$level*100, "% ", text.ci, " and ", text.pi, " ", x$level*100, "% prediction interval:\n\n", sep=""))
    print(x$results)
  } else {
    cat(paste("Summary ", text.stat, " with ", x$level*100, "% ", text.ci, " interval:\n\n", sep=""))
    print((x$results)[c("estimate", "CIl", "CIu")])
  }
  
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
}

#' Forest Plots
#' 
#' Function to create forest plots for objects of class \code{"valmeta"}.
#' 
#' @param x An object of class \code{"valmeta"}
#' @param sort By default, studies are ordered by ascending effect size (\code{sort="asc"}). For study ordering by descending
#' effect size, choose \code{sort="desc"}. For any other value, study ordering is ignored.
#' @param \ldots Additional arguments which are passed to \link{forest}.
#' 
#' @details The forest plot shows the performance estimates of each validation with corresponding confidence 
#' intervals. A polygon is added to the bottom of the forest plot, showing the summary estimate based on the model. 
#' A 95\% prediction interval is added by default for random-effects models,  the dotted line indicates its (approximate) bounds.
#' 
#' @references \itemize{
#' \item Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review 
#' and meta-analysis of prediction model performance. \emph{BMJ}. 2017;356:i6460.
#' \item Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{BMJ}. 2011 342:d549--d549.
#' }
#' 
#' @examples 
#' data(EuroSCORE)
#' fit <- with(EuroSCORE, valmeta(cstat=c.index, cstat.se=se.c.index, 
#'             cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), N=n, O=n.events))
#' plot(fit)
#' 
#' library(ggplot2)
#' plot(fit, theme=theme_grey())
#' 
#' @keywords meta-analysis discrimination calibration forest
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @import metafor
#' @import ellipse
#' @import ggplot2
#' @importFrom stats reorder
#' @return An object of class \code{ggplot}
#' 
#' @method plot valmeta
#' @export
plot.valmeta <- function(x, sort="asc", ...) {
  k <- dim(x$data)[1]
  yi.slab <- c(as.character(x$slab))
  yi <- c(x$data[,"theta"])
  ci.lb <- c(x$data[,"theta.CIl"])
  ci.ub <- c(x$data[,"theta.CIu"])
  
  if (x$model=="normal/logit") {
    yi <- sapply(yi, inv.logit)
    ci.lb <- sapply(ci.lb, inv.logit)
    ci.ub <- sapply(ci.ub, inv.logit)
  } else if (x$model == "normal/log" | x$model == "poisson/log") {
    yi <- sapply(yi, exp)
    ci.lb <- sapply(ci.lb, exp)
    ci.ub <- sapply(ci.ub, exp)
  }
  
  yi.ci <- cbind(ci.lb, ci.ub)
  
  if (x$measure=="cstat") {
    metamisc::forest(theta=yi, theta.ci=yi.ci, theta.slab=yi.slab, 
           theta.summary=as.numeric(x$results["estimate"]), 
           theta.summary.ci=as.numeric(x$results[c("CIl","CIu")]), 
           theta.summary.pi=as.numeric(x$results[c("PIl","PIu")]), 
           xlim=c(0,1),
           refline=0.5, xlab="c-statistic", sort=sort, ...)
  } else if (x$measure=="OE") {
    metamisc::forest(theta=yi, theta.ci=yi.ci, theta.slab=yi.slab, 
           theta.summary=as.numeric(x$results["estimate"]), 
           theta.summary.ci=as.numeric(x$results[c("CIl","CIu")]), 
           theta.summary.pi=as.numeric(x$results[c("PIl","PIu")]), 
           xlim=c(0,NA),
           refline=1, xlab="Total O:E ratio", sort=sort, ...)
  }
}

