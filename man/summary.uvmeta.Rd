\name{summary.uvmeta}
\alias{summary.uvmeta}

\title{Parameter summaries}

\description{This function provides summary estimates of a univariate meta-analysis model.}

\usage{ \method{summary}{uvmeta}(object, \ldots)}

\arguments{
  \item{object}{ a \code{uvmeta} object.}
  \item{\dots}{ arguments to be passed on to other functions }
}

\value{ The model parameters are given as \code{mu} (overall treatment effect),\code{tausq} (between-study variance if random effects were assumed), \code{Q} (Cochran's Q statistic) and \code{Isq} (I-square index).}

\references{
DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}:  753--768.

Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; \bold{21}: 1539--1558.

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.
}

\author{ Thomas Debray <thomas.debray@gmail.com> }

\note{ There are no confidence intervals for \code{tausq} when estimated with a frequentistic approach, as it is considered fixed. }

\seealso{
\code{\link{uvmeta}}
}


\keyword{ DerSimonian }
\keyword{ Laird }
\keyword{ univariate }
\keyword{ random-effects }
\keyword{ meta-analysis }

