\name{valmeta}
\alias{valmeta}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-analysis of prediction model performance
}
\description{
This function allows to meta-analyze the performance of a prediction model, and to obtain summary estimates of the c-statistic, the total observed-expected ratio and the calibration slope. Where appropriate, data transformations are applied and missing information is derived from available quantities.
}
\usage{
valmeta(cstat, cstat.se, cstat.95CI, OE, OE.se, OE.95CI, citl, citl.se,
        N, O, E, Po, Po.se, Pe, t.val, t.ma, method="REML", knha=TRUE, verbose=FALSE, 
        slab, n.chains = 4, pars, \ldots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cstat}{vector with the estimated c-statistic for each valiation}
  \item{cstat.se}{vector with the standard error of the estimated c-statistics}
  \item{cstat.95CI}{2-dimensional array with the lower (first column) and upper (second column) boundary of the 95\% confidence interval of the estimated c-statistics}
  \item{OE}{vector with the estimated ratio of total observed versus total expected events}
  \item{OE.se}{vector with the standard error of the estimated O:E ratios}
  \item{OE.95CI}{2-dimensional array with the lower (first column) and upper (second column) boundary of the 95\% confidence interval of the total O:E ratio}
  \item{citl}{vector with the estimated calibration-in-the-large for each valiation}
  \item{citl.se}{vector with the standard error of the estimated calibration-in-the-large statistics}
  \item{N}{vector with the total number of participants for each valiation.}
  \item{O}{vector with the total number of observed events for each valiation. }
  \item{E}{vector with the total number of expected events for each valiation}
  \item{Po}{vector with the (cumulative) observed event probability for each valiation}
  \item{Po.se}{vector with the standard errors of \code{Po}. Note that when \code{Po} is derived from Kaplan-Meier estimates, \code{Po.se} is equal to the standard error of the (cumulative) observed survival probability. }
  \item{Pe}{vector with the (cumulative) expected event probability for each validation}
  \item{t.val}{optional vector for prognostic models, containing the time period for which performance was assessed}
  \item{t.ma}{optional numeric value for prognostic models, containing the time period of primary interest for meta-analysis}
  \item{method}{character string specifying whether a fixed- or a random-effects model should be fitted. A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.}
  \item{knha}{logical specifying whether the method by Knapp and Hartung (2003) should be used for adjusting test statistics and confidence intervals (default is TRUE). Type '\code{?rma}' for more details.}
  \item{verbose}{if TRUE then messages generated during the fitting process will be displayed.}
  \item{slab}{optional vector with labels for the \emph{k} studies.}
  \item{n.chains}{the number of chains to use in the Gibbs sampler (\code{method="BAYES"}). More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly. The default number of chains is 4.}
  \item{pars}{A list with additional arguments.  The following parameters configure the MCMC sampling procedure:  \code{hp.mu.mean} (mean of the prior distribution of the random effects model, defaults to 0), \code{hp.mu.var} (variance of the prior distribution of the random effects model, defaults to 1E6), \code{hp.tau.min} (minimum value for the between-study standard deviation, defaults to 0), \code{hp.tau.max} (maximum value for the between-study standard deviation, defaults to 2), \code{hp.tau.sigma} (standard deviation of the prior distribution for the between-study standard-deviation), \code{hp.tau.dist} (prior distribution for the between-study standard-deviation. Defaults to \code{"dunif"}), \code{hp.tau.df} (degrees of freedom for the prior distribution for the between-study standard-deviation. Defaults to 3), \code{method.restore.c.se} (method for restoring missing estimates for the standard error of the c-statistic. So far, only \code{"Newcombe.2"} and \code{"Newcombe.4"} are supported. These methods have been described by Newcombe in 2006.), \code{model.cstat} (The likelihood/link for modeling the c-statistic; see "Details"), \code{model.oe} (The likelihood/link for modeling the O:E ratio; see "Details") }
  \item{\ldots}{Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).}
}
\details{
Unless specified otherwise, all meta-analysis models assume random effects and are fitted  using restricted maximum likelihood estimation with the \pkg{metafor} package (Viechtbauer 2010).  Further, confidence intervals for the average performance are based on the Hartung-Knapp-Sidik-Jonkman method, to better account for the uncertainty in the estimated between-study heterogeneity (Debray 2016). A Bayesian meta-analysis can be performed by specifying \code{method="BAYES"}. In that case, the R packages \pkg{runjags} and \pkg{rjags} must be installed. 


\bold{Meta-analysis of the concordance statistic (c-statistic)}

The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to distinguish between patients developing and not developing the outcome. The c-statistic typically ranges from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). A meta-analysis for the c-statistic will be performed if the c-statistics (\code{cstat}) and their respective standard errors (\code{cstat.se}) are defined. For studies where the standard error is unknown, it can be derived from the 95\% confidence interval, or from \code{cstat}, \code{O} and \code{N} (Newcombe 2006). By default, the meta-analysis model assumes Normality for the logit of the c-statistic (\code{model.cstat = "normal/logit"}). Alternatively, it is possible to summarize raw estimates of the c-statistic by setting \code{model.cstat = "normal/identity"}. 

\bold{Meta-analysis of the total observed versus expected ratio (O:E ratio)}

The total O:E ratio provides a rough indication of the overall model calibration (across the entire range of predicted risks). Currently, three methods have been implemented to obtain a summary estimate of the total O:E ratio. By default, the meta-analysis model assumes   Normality for the (natural) logarithm of the O:E ratios (\code{model.oe = "normal/log"}). Continuity corrections are applied when necessary by adding 0.5 to \code{O}, \code{E} and \code{N}. Alternatively, it is possible to model the total number of observed and expected events using a Poisson likelihood (Stijnen 2010). The resulting model does not require continuity corrections for \code{O} and can be implemented by setting \code{model.oe = "poisson/log"} (note that \code{hp.mu.var} is truncated to a maximum value of 100 for \code{method="BAYES"}). Finally, it is possible to summarize raw estimates of the O:E ratio by setting \code{model.oe = "normal/identity"}. 

When unkown, the standard error of the O:E ratio will be approximated in the following order from  (1) the 95\% confidence interval, (2) the standard error of \code{Po}, (3)  the error variance of the binomial distribution, (4) the error variance of the Poisson distribution, or from (5)  the calibration-in-the-large statistic.    

For meta-analysis of prognostic models, it is recommended to provide information on the time period (\code{t.val}) during which calibration was assessed in the validation study. When the time period of the validation study does not correspond to the time period of interest (\code{t.ma}), observed and expected survival probabilities will be extrapolated using Poisson distributions. Currently, extrapolation of event rates is only supported for \code{model.oe = "normal/log"} and \code{model.oe = "normal/identity"}. Note that  values for \code{O} and \code{N} should take the presence of drop-out into account. This implies that \code{O} is ideally based on Kaplan-Meier estimates, or that \code{N} should represent the total number of participants with complete follow-up.

\bold{Bayesian meta-analysis}

The prior distribution for the between-study standard deviation can be specified by \code{hp.tau.dist}, and is always truncated by \code{hp.tau.min} and \code{hp.tau.max}. Initial values for the between-study standard deviation are sampled from a uniform distribution with aformentioned boundaries. The following distributions are supported for modeling the prior of the between-study standard deviation: Uniform distribution (\code{hp.tau.dist="dunif"}; default), truncated Student-t distribution  (\code{hp.tau.dist="dhalft"}).

}
\value{
An object of class \code{valmeta} with the following elements:
\item{cstat}{object of class \code{vmasum} containing the meta-analysis results of the c-statistic.}
\item{oe}{object of class \code{vmasum} containing the meta-analysis results of the total O:E ratio.}
\item{method}{character string specifying whether a fixed- or a random-effects model was fitted.}
%\item{results}{vector with the summary c-statistic, the corresponding standard error, 95\% confidence interval and 95\% prediction interval.}
%\item{scale}{the scale used for meta-analysis.}
%\item{data}{2-dimensional array with the effect estimates (\code{theta}) and (restored) error variances (\code{theta.var}) of each validation. For \code{method="BAYES"}, the estimates represent draws from the posterior distribution.}
%\item{rma}{if applicable, fitted \code{rma} object of the meta-analysis.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017; 356:i6460. 

Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. \emph{Radiology}. 1982; 143(1):29--36. 

Newcombe RG. Confidence intervals for an effect size measure based on the Mann-Whitney statistic. Part 2: asymptotic methods and evaluation. \emph{Stat Med}. 2006; 25(4):559--73. 

Stijnen T, Hamza TH, Ozdemir P. Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Stat Med}. 2010; 29(29):3046--67. 

Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}
}
\author{Thomas Debray <thomas.debray@gmail.com>}

\seealso{
\code{\link{plot.valmeta}}
}

\examples{
######### Validation of prediction models with a binary outcome #########
data(EuroSCORE)
attach(EuroSCORE)

# Meta-analysis of the c-statistic (random effects)
fit <- valmeta(cstat=c.index, cstat.se=se.c.index, 
               cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
               N=n, O=n.events, slab=Study)
plot(fit)

# Meta-analysis of the total O:E ratio (random effects)
valmeta(O=n.events, E=e.events, N=n)     # Option 1
valmeta(O=n.events, E=e.events)          # Option 2 
valmeta(Po=Po, Pe=Pe, N=n)               # Option 3
valmeta(O=n.events, E=e.events, pars=list(model.oe="poisson/log"))  # Option 4


# Meta-analysis of the c-statistic and the total O:E ratio (random effects)
valmeta(cstat=c.index, cstat.se=se.c.index, 
        cstat.95CI=cbind(c.index.95CIl,c.index.95CIu), 
        N=n, O=n.events, E=e.events, slab=Study)
        
\dontrun{
# Bayesian meta-analysis of the c-statistic (random effects)
fit2 <- valmeta(cstat=c.index, cstat.se=se.c.index, 
                cstat.95CI=cbind(c.index.95CIl,c.index.95CIu),
                N=n, O=n.events, method="BAYES", slab=Study)
plot(fit2)

# Bayesian meta-analysis of the O:E ratio
pars <- list(model.oe="poisson/log", # Use a Poisson-Normal model
             hp.tau.dist="dhalft",   # Prior for the between-study standard deviation
             hp.tau.sigma=1.5,       # Standard deviation for 'hp.tau.dist'
             hp.tau.df=3,            # Degrees of freedom for 'hp.tau.dist'
             hp.tau.max=10)          # Maximum value for the between-study standard deviation
valmeta(O=n.events, E=e.events, N=n, method="BAYES", slab=Study, pars=pars)
} 

# Detach the dataset
detach(EuroSCORE)

######### Validation of prediction models with a time-to-event outcome #########
data(Framingham)
attach(Framingham)

# Meta-analysis of total O:E ratio after 10 years of follow-up
valmeta(Po=Po, Pe=Pe, N=n, t.val=t.val, t.ma=10)

# Detach the dataset
detach(Framingham)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
\keyword{ discrimination }
\keyword{ calibration }
