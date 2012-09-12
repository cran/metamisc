\name{summary.riley}
\alias{summary.riley}
\alias{print.summary.riley}

\title{Parameter summaries}

\description{Provides the summary estimtes of the alternative model for bivariate random-effects meta-analysis by Riley et al. (2008) with their corresponding confidence intervals. The model parameters are given as \code{beta1}, \code{beta2}, \code{psi1}, \code{psi2} and \code{rho}.  Confidence intervals are derived from the inverse Hessian.}

\usage{\method{summary}{riley}(object, level = 0.95, \ldots)}

\arguments{
  \item{object}{a \code{riley} object}
  \item{level}{numeric, the level for calculations of confidence intervals}
  \item{\dots}{arguments to be passed on to other functions}
}

\details{For diagnostic test accuracy data, \code{beta1} equals the logit sensitivity (Sens) and \code{beta2} equals the logit false positive rate (FPR). The summary sensitivity and FPR are added for completeness.}


\value{ array with confidence intervals for the estimated model parameters. For diagnostic test accuracy data, the resulting summary sensitivity and false positive rate are included.}

\references{
Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}


\seealso{
\code{\link{riley}}
\code{\link{plot.riley}}
}

\keyword{riley}
\keyword{confidence}
\keyword{bivariate}
