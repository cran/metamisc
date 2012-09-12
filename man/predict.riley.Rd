\name{predict.riley}
\alias{predict.riley}
\concept{prediction interval}

\title{Prediction Interval}

\description{Calculates a prediction interval for the summary parameters of Riley's alternative model for bivariate random-effects meta-analysis. This interval predicts in what range future observations will fall given what has already been observed.}

\usage{\method{predict}{riley}(object, level = 0.95, \ldots)}

\arguments{
  \item{object}{ a \code{riley} object.}
  \item{level}{ numeric, the level for calculations of confidence intervals }
  \item{\dots}{ arguments to be passed on to other functions }
}

\details{ Prediction intervals are based on Student's t-distribution with (numstudies - 5) degrees of freedom.}

\value{Array containing prediction intervals for the summary estimates \code{beta1} and \code{beta2} (for effect size data), or for the mean sensitivity and false positive rate (diagnostic test accuracy data).}

\author{Thomas Debray <thomas.debray@gmail.com> }


\seealso{\code{\link{riley}}}

\keyword{ riley }
\keyword{ prediction }
\keyword{ bivariate }
