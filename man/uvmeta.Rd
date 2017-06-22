\name{uvmeta}
\alias{uvmeta}
\alias{uvmeta.default}
\alias{plot.uvmeta}
\alias{print.uvmeta}

\title{Univariate meta-analysis.}

\description{This function performs a univariate meta-analysis by assuming fixed or random effects.  Whereas the fixed effects model assumes that all studies in the analysis share a common effect size, the random-effects model allows different study-specific effect sizes. Concretely, if we move from fixed-effect weights to random-effects weights, large studies lose influence and small studies gain influence (Borenstein 2010). }

\usage{
uvmeta(r, r.se, method="REML", test="knha", labels, na.action, pars, 
       verbose=FALSE, ...)}

\arguments{
  \item{r}{vector of numerics containing the effect sizes}
  \item{r.se}{vector of numerics containing the standard error of the effect sizes}
  \item{method}{character string specifying whether a fixed- or a random-effects model should be fitted. A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.}
  \item{test}{character string specifying how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="knha"}), the method by Knapp and Hartung (2003) is used for adjusting test statistics and confidence intervals.  Type '\code{?rma}' for more details.}
  \item{labels}{vector of characters containing the labels for the studies}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. Defaults to "na.fail", other options are "na.omit", "na.exclude" or "na.pass".}
  \item{pars}{A list with additional arguments. The following parameters configure the MCMC sampling procedure:  \code{hp.mu.mean} (Hyperparameter: mean of the prior distribution of the fixed/random effects model, defaults to zero), \code{hp.mu.var} (Hyperparameter: variance of the prior distribution of the fixed/random effects model, defaults to 1000), \code{n.chains} (specifies the number of parallel chains)}
  \item{verbose}{if TRUE then messages generated during the fitting process will be displayed.}
  \item{\dots}{Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).}
}

\details{Unless specified otherwise, all meta-analysis models assume random effects and are fitted  using restricted maximum likelihood estimation with the \pkg{metafor} package (Viechtbauer 2010).  Further, confidence intervals for the average performance are based on the Hartung-Knapp-Sidik-Jonkman method, to better account for the uncertainty in the estimated between-study heterogeneity (Debray 2016). A Bayesian meta-analysis can be performed by specifying \code{method="BAYES"}. In that case, the R packages \pkg{runjags} and \pkg{rjags} must be installed. 

For random effects models, a prediction interval for the pooled effect size is displayed. This interval predicts in what range future effect sizes will fall given what has already been observed (Higgins 2009, Riley 2011).  

\bold{Bayesian meta-analysis models}


For Bayesian meta-analysis models that involve the Gibbs sampler (\code{method="BAYES"}), the R packages \code{runjags} and \code{rjags} must be installed. The Bayesian approach uses an uninformative Normal prior for the mean and a uniform prior for the between-study variance of the pooled effect size (Higgins 2009). By default, the Normal prior has a mean of 0 and a variance of 1000. These hyperparameters can, however, be altered through the variables \code{hp.mu.mean} and \code{hp.mu.var} in the argument \code{pars}. The prior distribution of the between-study standard deviation is given by a uniform distribution, bounded between 0 and 100. Future versions of \code{metamisc} will allow to alter these boundaries.
}

\value{ An object of the class \code{uvmeta} for which many standard methods are available. If \code{method="BAYES"}, the results contain an object of the class \code{runjags}. }

\references{
Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.

DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Graham PL, Moran JL. Robust meta-analytic conclusions mandate the provision of prediction intervals in meta-analysis summaries. \emph{Journal of Clinical Epidemiology} 2012; \bold{65}: 503--510.

Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; \bold{21}: 1539--1558.

Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} 2009, \bold{172}: 137--159.

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.

Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\seealso{
\code{uvmeta-class}
}

\examples{
data(Roberts)

# Frequentist random-effects meta-analysis
fit1 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts)))
summary(fit1)
plot(fit1, main="Forest plot") #show a forest plot
fit1

\dontrun{
# Bayesian random effects meta-analysis 
fit2 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts), method="BAYES"))
plot(fit2)
}
}


\keyword{ DerSimonian }
\keyword{ Laird }
\keyword{ univariate }
\keyword{ random-effects }
\keyword{ meta-analysis }

