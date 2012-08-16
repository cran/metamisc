\name{summary.riley}
\alias{summary.riley}
\alias{print.summary.riley}

\title{Parameter summaries}

\description{Provides the summary estimtes of the alternative model for bivariate random-effects meta-analysis by Riley et al. (2008) with their corresponding confidence intervals. The model parameters are given as mu1 (logit sensitivity), mu2 (logit FPR), psi1, psi2 and rho; the summary sensitivity and false positive rate are added for completeness. Confidence intervals are derived from the inverse Hessian. }

\usage{\method{summary}{riley}(object, level = 0.95, \ldots)}

\arguments{
  \item{object}{a \code{riley} object}
  \item{level}{numeric, the level for calculations of confidence intervals}
  \item{\dots}{arguments to be passed on to other functions}
}

\value{ array with confidence intervals for the estimated model parameters, as well for the resulting mean sensitivity (Sens) and false positive rate (FPR).}

\references{
Riley, RD., Thompson, JR., & Abrams, KR. (2008).
\dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}


\seealso{
\code{\link{riley}}
\code{\link{plot.riley}}
}

\keyword{riley}
\keyword{confidence}
\keyword{bivariate}
