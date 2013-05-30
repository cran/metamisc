\name{uvmeta}
\alias{uvmeta}
\alias{uvmeta.default}
\alias{plot.uvmeta}
\alias{print.uvmeta}

\title{Univariate meta-analysis.}

\description{This function performs a univariate meta-analysis by assuming fixed or random effects.  Whereas the fixed effects model assumes that all studies in the analysis share a common effect size, the random-effects model allows different study-specific effect sizes. Concretely, if we move from fixed-effect weights to random-effects weights, large studies lose influence and small studies gain influence (Borenstein 2010). }

\usage{
uvmeta(r, vars, model="random", method="MOM", labels, na.action, pars, 
       verbose=FALSE, ...)}

\arguments{
  \item{r}{vector of numerics containing the effect sizes}
  \item{vars}{vector of numerics containing the error variance of the effect sizes}
  \item{model}{Assume "random" or "fixed" effects.}
  \item{method}{Estimation method: use "MOM" to implement the non-parametric method-of-moment estimator from DerSimonian and Laird, "ml" to implement the maximum-likelihood estimator, "pl" to use the profile-likelihood estimator or "bayes" to implement a Bayesian meta-analysis assuming normality of the random effects (Higgins 2009).}
  \item{labels}{vector of characters containing the labels for the studies}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. Defaults to "na.fail", other options are "na.omit", "na.exclude" or "na.pass".}
  \item{pars}{A list with additional arguments. Use "level" to specify the level of confidence or credibility intervals. The following parameters configure the MCMC sampling procedure and are ignored if \code{method="MOM"}:  \code{hp.mu.mean} (Hyperparameter: mean of the prior distribution of the fixed/random effects model, defaults to zero), \code{hp.mu.var} (Hyperparameter: variance of the prior distribution of the fixed/random effects model, defaults to 1000), \code{n.chains} (specifies the number of parallel chains), \code{n.adapt} (specifies the number of iterations for adaptation), \code{n.init} (number of iterations to run for initializing the Markov chain) and \code{n.iter} (number of iterations to monitor).}
  \item{verbose}{if TRUE then messages generated during the fitting process will be displayed.}
  \item{\dots}{arguments to be passed on to other functions}
}

\details{The Bayesian approach uses an uninformative Normal prior for the mean and an uninformative uniform prior for the variance of the pooled effect size (Higgins 2009).  For random effects models, a prediction interval for the pooled effect size is displayed. This interval predicts in what range future effect sizes will fall given what has already been observed (Higgins 2009, Riley 2011).  
}

\value{ An object of the class \code{uvmeta} for which many standard methods are available. }

\references{
Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.

DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Gasparrini A, Armstrong B, Kenward, MG. Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine} 2012; \bold{31}: 3821--3839.

Graham PL, Moran JL. Robust meta-analytic conclusions mandate the provision of prediction intervals in meta-analysis summaries. \emph{Journal of Clinical Epidemiology} 2012; \bold{65}: 503--510.

Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; \bold{21}: 1539--1558.

Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} 2009, \bold{172}: 137--159.

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.

Rucker G, Schwarzer G, Carpenter JR, Schumacher M . Undue reliance on I2 in assessing heterogeneity may mislead. \emph{BMC Medical Research Methodology} 2008; \bold{8}: 79.

}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\seealso{
\code{uvmeta-class}
}

\examples{
data(Roberts)

#Extract effect size and error variance
r <- Roberts$SDM
vars <- Roberts$SE**2

#Frequentist random-effects meta-analysis
fit1 <- uvmeta(r,vars,labels=rownames(Roberts))
plot(fit1, main="Forest plot") #show a forest plot
fit1


#Bayesian random-effects meta-analysis
#fit2 <- uvmeta(r,vars,method="bayes")
#fit2
}

\keyword{ DerSimonian }
\keyword{ Laird }
\keyword{ univariate }
\keyword{ random-effects }
\keyword{ meta-analysis }

