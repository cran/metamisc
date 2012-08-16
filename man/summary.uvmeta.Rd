\name{summary.uvmeta}
\alias{summary.uvmeta}

\title{Parameter summaries}

\description{This function provides summary estimates of univariate random-effects models with their corresponding confidence intervals. The model parameters are given as mu (overall treatment effect) and tau squared (between-study variance).}

\usage{ \method{summary}{uvmeta}(object, level = 0.95, \ldots)}

\arguments{
  \item{object}{ a \code{uvmeta} object.}
  \item{level}{ numeric, the level for calculations of confidence intervals }
  \item{\dots}{ arguments to be passed on to other functions }
}

\value{ Array containing confidence regions for the overall treatment effect and heterogeneity measures.}

\references{
DerSimonian, R., & Laird, N. (1986). 
\dQuote{Meta-analysis in clinical trials.} \emph{Controlled Clinical Trials}, \bold{7}, 177--188.

Biggerstaff, BJ., & Tweedie, RL. (1997). 
\dQuote{Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis.} \emph{Statistics in Medicine}, \bold{16},  753--768.

Higgins, JPT., & Thompson, SG. (2002). 
\dQuote{Quantifying heterogeneity in a meta-analysis.} \emph{Statistics in Medicine}, \bold{21}, 1539--1558.

Borenstein, M., Hedges, LV., Higgins, JPT., & Rothstein HR. (2010). 
\dQuote{A basic introduction to fixed-effect and random-effects models for meta-analysis.} \emph{Research Synthesis Methods}, \bold{1}, 97--111.

Riley, RD., Higgins, JPT., & Deeks, JJ. (2011). 
\dQuote{Interpretation of random effects meta-analyses.} \emph{British Medical Journal}, \bold{342}, d549.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\note{ There are no confidence intervals for the tau-parameter, as it is considered fixed. }

\seealso{
\code{\link{uvmeta}}
}


\keyword{ DerSimonian }
\keyword{ Laird }
\keyword{ univariate }
\keyword{ random-effects }
\keyword{ meta-analysis }

