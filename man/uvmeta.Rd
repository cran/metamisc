\name{uvmeta}
\alias{uvmeta}
\alias{uvmeta.default}
\alias{uvmetaMOM}
\alias{print.uvmeta}

\title{Perform a univariate meta-analysis.}

\description{This function performs a univariate meta-analysis by assuming fixed or random effects. }

\usage{ uvmeta(r, vars, method="MOM", ...)}

\arguments{
  \item{r}{vector of numerics containing the effect sizes}
  \item{vars}{vector of numerics containing the error variance of the effect sizes}
  \item{method}{estimation method, where \code{MOM} implements the method-of-moment estimator proposed by DerSimonian and Laird.}
  \item{\dots}{arguments to be passed on to other functions}
}

\details{ The defining feature of the fixed-effect model is that all studies in the analysis share a common effect size. Conversely, the random-effects model assumes a distribution of true effect sizes and estimates the mean of this distribution. Concretely, if we move from fixed-effect weights to random-effects weights, large studies lose influence and small studies gain influence (Borenstein 2010). }

\value{ An object of the class \code{uvmeta} for which many standard methods are available. }

\references{
DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.

Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; \bold{21}: 1539--1558.

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\seealso{
\code{\link{summary.uvmeta}} \code{\link{predict.uvmeta}}
}

\keyword{ DerSimonian }
\keyword{ Laird }
\keyword{ univariate }
\keyword{ random-effects }
\keyword{ meta-analysis }

