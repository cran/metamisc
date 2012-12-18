\name{uvmeta}
\alias{uvmeta}
\alias{uvmeta.default}
\alias{uvmetaMOM}
\alias{print.uvmeta}

\title{Perform a univariate meta-analysis.}

\description{This function performs a univariate meta-analysis by assuming fixed or random effects. }

\usage{ uvmeta(r, vars, model="random", method="MOM", pars=list(quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), n.chains=4, n.adapt=5000, n.init=1000, n.iter=10000), ...)}

\arguments{
  \item{r}{vector of numerics containing the effect sizes}
  \item{vars}{vector of numerics containing the error variance of the effect sizes}
  \item{model}{Specifies whether \code{random} (default) or \code{fixed} effects should be assumed}
  \item{method}{estimation method, where \code{MOM} implements the method-of-moment estimator proposed by DerSimonian and Laird. Altnernatively, \code{bayes} implements a Bayesian approach that assumes normality of the random effects.}
  \item{pars}{A list with additional arguments. Use \code{quantiles} to specify which quantiles should be calculated for parameter estimates. The following parameters configure the MCMC sampling procedure and are ignored if method is \code{MOM}:  \code{n.chains} (specifies the number of parallel chains), \code{n.adapt} (specifies the number of iterations for adaptation), \code{n.init} (number of iterations to run for initializing the Markov chain) and \code{n.iter} (number of iterations to monitor).}
  \item{\dots}{arguments to be passed on to other functions}
}

\details{ The defining feature of the fixed-effect model is that all studies in the analysis share a common effect size. Conversely, the random-effects model assumes a distribution of true effect sizes and estimates the mean of this distribution. Concretely, if we move from fixed-effect weights to random-effects weights, large studies lose influence and small studies gain influence (Borenstein 2010). Note that JAGS need to be installed in order to run a Bayesian meta-analysis, see http://mcmc-jags.sourceforge.net/. }

\value{ An object of the class \code{uvmeta} for which many standard methods are available.}

\references{
DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.

Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; \bold{21}: 1539--1558.

Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} 2009, \bold{172}: 137--159.

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\seealso{
\code{\link{summary.uvmeta}} \code{\link{predict.uvmeta}}
}

\examples{
data(Roberts)

#Extract effect size and error variance
r <- Roberts[,1]
vars <- Roberts[,2]**2

#Frequentist random-effects meta-analysis
fit1 <- uvmeta(r,vars)
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

